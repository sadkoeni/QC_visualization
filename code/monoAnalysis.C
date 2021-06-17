/*
 gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");

  .L $NOTES/JIRA/ATO-432/code/monoAnalysis.C+g

  AliMonopolRecParam recoParam;
  recoParam.fMinQSeedInit=200;
  doReconstruction(recoParam,-1);

  initTree();
  initClusters();
loadClusters(0,100)
  // calgrind fast options
AliMonopolRecParam recoParam;
  recoParam.fMinQSeedInit=500;
  recoParam.fMinQSeed=300;

*/
///
/// construct array
/// fill array in tree loop
/// per event
/// Sort array - create index arrays (sorting by z )
/// TMath::BinarySearch

#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TVectorF.h>
#include <TVector.h>
#include <AliXRDPROOFtoolkit.h>
#include "AliTPCclusterMI.h"
#include "AliTPCROC.h"
#include "TMath.h"
#include "TObject.h"
#include "AliDrawStyle.h"
#include "TStyle.h"
#include "ROOT/RVec.hxx"
#include "TTreeStream.h"
#include "AliRieman.h"
#include "AliExternalTrackParam.h"
#include "TDatabasePDG.h"

using namespace ROOT::VecOps;
#pragma link C++ class RVec<AliTPCclusterMI>+;
#pragma link C++ class RVec<AliExternalTrackParam>+;
#pragma link C++ class RVec<float>+;

TChain * tree = 0;
TArrayD gidArray(0);
TArrayI eventIndeces(0);
TTreeSRedirector *pcstream = 0;
TTree * treeSeeds0=0;
TTree * treeSeeds1=0;

const AliTPCROC*roc = AliTPCROC::Instance();
const Int_t kMaxCl=100000;
//const Int_t nRows= roc->GetNRows(0)+roc->GetNRows(36);
const Int_t nRows=159;
const Float_t Bz=5;   // to be read from database for given run
RVec<AliTPCclusterMI> clusterArray[nRows] ;
UInt_t   nClusters[nRows];
TVectorF * arrayGX[nRows], *arrayGY[nRows], *arrayGZ[nRows], *arrayPhi[nRows];
AliTPCclusterMI *pcluster=0;

class AliMonopolRecParam: public TObject{
  public:
  AliMonopolRecParam();
  // seeding parameters
  Int_t fMinQSeed;               // minimal threshold for cluster charge  to connect O(100);
  Int_t fMinQSeedInit;           // minimal threshold for cluster charge  for 3 initial clusters to seed O(200)
  Int_t fMinClustersIter1;       // minimal ammount of clusters in 1 iteration 3 clusters +-delta road O(5)
  const Int_t fDRows;            // number of delta rows for seed Iter1  - O(2)
  const float fDeltaZseed0;      // deltaZ  in initial 3 cluster seeding O(5 cm)
  const float fDeltaYseed0;      // deltaY  in initial 3 cluster seeding O(5 cm)
  const float fDeltaZseed1;      // deltaZ  in  seeding1 phase (linear interpolation  O(2 cm)
  const float fDeltaYseed1;      // deltaY  in  seeding1 phase (linear interpolation  O(2 cm)
  // kalman parameters
  const int fKalmanStep;         // number of rows to update kalman
  const int fKalmanStepMin;      // minimal number of clusters to update  if not present tracking stops
  // ToT parameters
  Float_t   fDeltaRPhiTot;       //  search window for TOT cluster    default O(2 cm)
  Float_t   fZStepTot;           //  search step for Tot clucter      default O(0.5 cm) ~ 2 time bins
  Float_t   fZCutTot;            //  minimal length without cluster   default O(3 cm) ~ 10 time bins for RUN1,2 with 100 ns
  Float_t   fChi2Cut;            //  chi2 cut for track/cluster association -err y, err z is calcuated as weighted mean of the local and global rms
  // Gas parameters as obtained from GeoManager ad Atima for ArCO2
  Float_t   fX0;                 //  radiation length =7.8350968e-05 for CorrectForMaterial
  Float_t   fXRho;               //  density =0.0016265266
  Float_t    fdEdxMIP;            //  1.523 10^-6 GeV/(g/cm2))
   ClassDef(AliMonopolRecParam, 1); // monopole seed
};

AliMonopolRecParam::AliMonopolRecParam():
  TObject(),
  fMinQSeed(100),
  fMinQSeedInit(200),
  fMinClustersIter1(5),
  fDRows(2),
  fDeltaZseed0(5),
  fDeltaYseed0(5),
  fDeltaZseed1(2),
  fDeltaYseed1(2),
  fKalmanStep(4),
  fKalmanStepMin(2),
  fDeltaRPhiTot(2),
  fZStepTot(0.5),
  fZCutTot(3),
  fChi2Cut(16),
  fX0(7.8350968e-05),    //  radiation length =7.8350968e-05 of ArCO2 per cm for CorrectForMaterial
  fXRho(0.0016265266),   //  density  ArCO2 =0.0016265266
  fdEdxMIP(1.523e-6)
{
}


class AliTPCseedMonopol : public TObject {
public:
  AliTPCseedMonopol();
  Int_t  fnSeedIter1;                      //
  Float_t fAlpha;                          // rotation coordinate system
  Int_t   fRowSeed;                        // middle row for seeding
  Float_t fQMeanSeed1;                     // reconstructed mean charge in seed1 phase
  Float_t fQMedianSeed1;
  RVec<AliTPCclusterMI>  fSeedIter1;       // cluster array assingin seed 1 stage
  RVec<float>  fSeedIter1Tot;
  RVec<float>  fSeedIter1QTot;
  Int_t fNclusters;                        // number of clusters
  RVec<AliTPCclusterMI>  fClustersTrack;   // clusters assigned to full track (index is pad row)
  RVec<AliTPCclusterMI>  fClustersTot;     // clusters i saturation  (clusters not sorted)
  AliRieman  fRieman;                      // rieman filter for seed 1 stage
  AliExternalTrackParam fParamSeed1;       // parameter from seed1 at radius of middle cluster
  //
  AliExternalTrackParam fParamInner;         // parameters at lower R
  AliExternalTrackParam fParamOuter;         // parameter at highest R
  RVec<AliExternalTrackParam>  fParamArrayIn;   // parameters along trajectory inner direction
  RVec<AliExternalTrackParam>  fParamArrayOut;  // parameters along trajectory outer direction

  Float_t               fQTotMedian;         // median Qtot dEdx estimator
  Float_t               fQMaxMedian;         // median Qtot dEdx estimator
  Float_t               fRmsYMedian0;        // median RMS in iteration 0  - local RMS
  Float_t               fRmsZMedian0;        // median RMS in iteration 0  - local RMS
  Float_t               fChi2Inner;          // chi2 of track, cluster association
  Float_t               fChi2Outer;          // chi2 of track, cluster association
  Float_t               fNClustersInner;     // number of clusters assigned in inner refit
  Float_t               fNClustersOuter;     // number of clusters assigned in outer refit
  ClassDef(AliTPCseedMonopol, 1); // monopole seed
};

AliTPCseedMonopol::AliTPCseedMonopol():
  fnSeedIter1(0),                       //
  fAlpha(0),                            // rotation coordinate system
  fRowSeed(0),                          // middle row for seeding
  fQMeanSeed1(0),                       // reconstructed mean charge in seed1 phase
  fQMedianSeed1(0),
  fSeedIter1(),                         // cluster array assingin seed 1 stage
  fSeedIter1Tot(),
  fSeedIter1QTot(),
  fNclusters(0),                        // number of clusters
  fClustersTrack(),                     // clusters assigned to full track (index is pad row)
  fClustersTot(),                       // clusters i saturation  (clusters not sorted)
  fRieman(160),                         // rieman filter for seed 1 stage
  fParamSeed1(),                        // parameter from seed1 at radius of middle cluster
  //
  fParamInner(),                        // parameters at lower R
  fParamOuter(),                        // parameter at highest R
  fParamArrayIn(),                      // parameters along trajectory inner direction
  fParamArrayOut(),                     // parameters along trajectory outer direction
  fQTotMedian(0),                       // median Qtot dEdx estimator
  fQMaxMedian(0),                       // median Qtot dEdx estimator
  fRmsYMedian0(0),                      // median RMS in iteration 0  - local RMS
  fRmsZMedian0(0),                      // median RMS in iteration 0  - local RMS
  fChi2Inner(0),                        // chi2 of track, cluster association
  fChi2Outer(0),                        // chi2 of track, cluster association
  fNClustersInner(0),                   // number of clusters assigned in inner refit
  fNClustersOuter(0)                   // number of clusters assigned in outer refit
{

}
///
// sort using a custom function object
struct {
        bool operator()(AliTPCclusterMI a, AliTPCclusterMI b) const{ return a.GetZ() < b.GetZ();}
} customLessCluster;

int  clusterLowerBound(RVec<AliTPCclusterMI> clArray, float z){
  AliTPCclusterMI clusterSort; clusterSort.SetZ(z);
  auto lower = std::lower_bound(clArray.begin(), clArray.end(), clusterSort, customLessCluster);
  return lower - clArray.begin();
}

/// find nearest cluster (L2 metric) close to the y0, z0 with tolerance dy0, dz0
/// WE ARE USING AliTPCClusterMI structure with function GetX,GetY,GetZ but internaly R,Phi,Z is stored
/// Distance
/// \param clArray
/// \param y0     - position y0 (in our case we store in y0 phi angle) - to get distance delta = dphi *R
/// \param z0
/// \param dy0    - tolerance dy
/// \param dz0    - tolerance dz
/// \return index of the cluster
int findNearest(RVec<AliTPCclusterMI> clArray, float y0, float z0, float dy0, float dz0){
  uint index0= clusterLowerBound(clArray, z0-dz0);
  int nCl=0;
  float dist2Min=dy0*dy0+dz0*dz0;
  int iClosest=-1;
  for (uint iCl=index0; iCl<clArray.size(); iCl++){
    AliTPCclusterMI &cl = clArray[iCl];
    float dy=abs(cl.GetY()-y0)*cl.GetX();
    if (dy>dy0) continue;
    float dz=abs(cl.GetZ()-z0);
    if (dz>dz0) break;
    nCl++;
    float  dist2 = dy*dy+dz*dz;
    if (dist2<dist2Min){
      iClosest=iCl;
      dist2Min=dist2;
    }
  }
  return iClosest;
}


void initClusters(Int_t capacity=kMaxCl) {
  for (Int_t i = 0; i < nRows; i++) {
    //arrayClusters[i] = new TClonesArray("AliTPCclusterMI", capacity);
    //arrayClusters[i]->ExpandCreateFast(capacity);
    arrayGX[i] = new TVectorF(capacity);
    arrayGY[i] = new TVectorF(capacity);
    arrayGZ[i] = new TVectorF(capacity);
    arrayPhi[i] = new TVectorF(capacity);
    //indecesSorted[i]=new TArrayI(capacity);
  }
}


/// load cluster above recoParam.fMinQseed threshold  for eventNumber  to RVec array
/// * AliTPCCluster.fRow  -  changed and is indexed from 0-159 (not as in ALICE tracking)
/// * AliTPCCluster.fY    - chaged and is used for phi coordinate - to be able to search in global space
/// \param eventNumber
/// \param recoParam
void loadClusters(Int_t eventNumber, const AliMonopolRecParam&recoParam){
  const float reserveFactor=2.;
  for (Int_t i=0; i<nRows; i++) nClusters[i]=0;
  //for (Int_t i=0; i<nRows; i++) clusterArray[i].resize(0);
  Float_t gx,gy,gz;
  tree->SetBranchAddress("Cluster.",&pcluster);
  tree->SetBranchAddress("gx",&gx);
  tree->SetBranchAddress("gy",&gy);
  tree->SetBranchAddress("gz",&gz);
  for (Int_t i=0; i<nRows; i++) {
    clusterArray[i].resize(kMaxCl);
  }
  for (Int_t i=eventIndeces[eventNumber]; i<eventIndeces[eventNumber+1]; i++){
    tree->GetEntry(i);
    //pcluster->Dump();
    if (pcluster->GetMax()<recoParam.fMinQSeed) continue;
    UInt_t   row = pcluster->GetRow()+((pcluster->GetDetector()>=36)?roc->GetNRows(0):0);
    if (row>=nRows) {
      ::Error("loadClusters","Invalid row");
      continue;
    }else{
      //::Info("loadClusters","%d",i);
    }
    RVec<AliTPCclusterMI> &clusterArrayRow = (clusterArray[row]);
    if (nClusters[row]>=clusterArrayRow.capacity()) clusterArrayRow.reserve((nClusters[row]+1)*reserveFactor);
    if (nClusters[row]>=clusterArrayRow.size())     clusterArrayRow.resize(nClusters[row]);

    AliTPCclusterMI       & cluster= clusterArrayRow[nClusters[row]];
    cluster= *pcluster;
    Double_t radius    = TMath::Sqrt(gx*gx+gy*gy);
    Double_t phi      = TMath::ATan2(gy,gx);
    cluster.SetZ(gz);
    cluster.SetX(radius);
    cluster.SetY(phi);
    cluster.SetRow(row);
    if (nClusters[row]>=kMaxCl) continue;   ///TODO expand arrays instead of continue
    //(*arrayGX[row])[nClusters[row]]=gx;
    //(*arrayGY[row])[nClusters[row]]=gy;
    //(*arrayGZ[row])[nClusters[row]]=gz;
    //(*arrayPhi[row])[nClusters[row]]=TMath::ATan2(gy,gx);
    nClusters[row]++;
  }
  for (Int_t i=0; i<nRows; i++) {
    clusterArray[i].resize(nClusters[i]);
    std::sort(clusterArray[i].begin(), clusterArray[i].end(), customLessCluster);
  }
}


void initTree(){
  for (int iRow=0; iRow<nRows; iRow++) {
    clusterArray[iRow].reserve(kMaxCl);
    clusterArray[iRow].resize(kMaxCl);
  }
  tree  = AliXRDPROOFtoolkit::MakeChainRandom("cluster.list","clusterDumpFull",0,1000);
  tree->SetCacheSize(100000000);
  tree->SetEstimate(tree->GetEntries());
  Int_t entries = tree->Draw("gid>>his(100,0,1)","1","goff");
  gidArray.Set(entries);
  eventIndeces.Set(entries);
  Double_t gidOld=-1;
  Double_t *gids = tree->GetV1();
  Int_t nEvents=0;
  for (Int_t i=0; i<entries; i++){
    //if (i%100==0) printf("%d\n",i);
    if (gidOld!=gids[i]){
      gidArray[nEvents]=gids[i];
      gidOld=gids[i];
      eventIndeces[nEvents]=i;
      nEvents++;
    }
  }
  eventIndeces[nEvents]=entries;
  nEvents++;
  eventIndeces.Set(nEvents);
  gidArray.Set(nEvents);
  tree->SetMarkerStyle(21);
  tree->SetMarkerSize(0.5);
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
}

/// find clusters in range cl0-2 -> cl2+2
///      linear interpolation between tracks used
/// \param cl0        - seeding cluster lowest
/// \param cl1        - seeding cluster lowest
/// \param cl2        - seeding cluster lowest
/// \param recoParam  - monopole reconstruction parameters
/// \param seed
/// \return number of cluster associated + sets clusters in seed
int makeSeed1(AliTPCclusterMI&  cl0, AliTPCclusterMI&  cl1, AliTPCclusterMI&  cl2, const AliMonopolRecParam& recoParam, AliTPCseedMonopol &seed){
  RVec<AliTPCclusterMI> &seedIter1=seed.fSeedIter1;
  int row0 = cl0.GetRow();
  int row1 = cl1.GetRow();
  int row2 = cl2.GetRow();
  int nClusters =0;
  //
  for (int row=row0-recoParam.fDRows; row<row2+recoParam.fDRows; row++){
    int rowIndex=row-(row0-2);
   seedIter1[rowIndex].SetDetector(-1);
    if (row<0) continue;
    if (row>=nRows) continue;
    float dX,dphi,dz=0;
    float x0,y0,z0,dx0,dy0,dz0;
    // linear interpolation of the x,y,z position at row
    AliTPCclusterMI &rCl0=(row>row1)?cl1:cl0;
    AliTPCclusterMI &rCl1=(row>row1)?cl2:cl1;
    int dRow=row-rCl0.GetRow();
    float dRow0=rCl1.GetRow()-rCl0.GetRow();
    dx0   =  rCl1.GetX()-rCl0.GetX();
    dy0   =  rCl1.GetY()-rCl0.GetY();
    dz0   =  rCl1.GetZ()-rCl0.GetZ();
    float xRow  =  rCl0.GetX()+dRow*dx0/dRow0;
    float yRow  =  rCl0.GetY()+dy0*(xRow-rCl0.GetX())/dx0;
    float zRow  =  rCl0.GetZ()+dz0*(xRow-rCl0.GetX())/dx0;
    // find nearest =0;
    RVec<AliTPCclusterMI>& clusterRaw = clusterArray[row];
    int nearestIndex = findNearest(clusterRaw, yRow,zRow,recoParam.fDeltaYseed1,recoParam.fDeltaZseed1);
    if (nearestIndex<0) continue;
    nClusters++;
    seedIter1[rowIndex]=clusterRaw[nearestIndex];
    // define tot
  }
  seed.fnSeedIter1=nClusters;
  return nClusters;
}

/// Open rphi road around seed/cluster  interpolation
//   * loop in z direction in increasing time direction with  some step
//   * find nearest cluster
//   ( stop in case no cluster for some time interval ...
//integrate charge and time
/// \param cluster    - TPC luster in r,phi z coordinate
/// \param phiTrack
/// \param zTrack
/// \param deltaPhi   - range to look for new cluster
/// \param zStep      - step in which we look to nearest cluster
/// \param zCut       - maximal length with cluster in z
float accumulateTOT(AliTPCclusterMI& cluster, Float_t phiTrack, Float_t zTrack, Float_t deltaPhi, Float_t zStep=0.5, Float_t zCut=3){
  ///
  float z0= cluster.GetZ();
  float sign=(cluster.GetDetector()%36<18) ? -1.:1.;  // on A side V drift has other direction tha in other side
  int clusterId=1;
  UInt_t   row = cluster.GetRow();
  RVec<AliTPCclusterMI>& clusterRaw = clusterArray[row];
  float zLast=z0;
  Float_t length=0;
  for (float z=z0 +zStep*sign; abs(z)<250; z+=sign*zStep){
    int nearestIndex = findNearest(clusterRaw, cluster.GetY(),z,deltaPhi,zStep);
    if (nearestIndex<0) {
      if (TMath::Abs(z-zLast) > zCut) break;
      continue;
    }
    float zClusterT=clusterArray[row][nearestIndex].GetZ();
    length=abs(zClusterT-z0);
  }
  return length;
}

/// Refit seed in local rotated coordinate system - with rotation -media of phi cluster position
/// \param seed
void refitTrackRiemanSeed1(AliTPCseedMonopol & seed, const AliMonopolRecParam & recoParam){
  // Debug symbols:
  //      xOrig.fData[0]@seed.fRieman.fN
  //      yOrig.fData[0]@seed.fRieman.fN
  //      zOrig.fData[0]@seed.fRieman.fN
  //seed.fRieman= AliRieman(seed.fSeedIter1.size());
  seed.fRieman.Reset();
  const Float_t kMinRadius=80;
  RVec<float> alpha(seed.fSeedIter1.size());
  RVec<float> xOrig(seed.fSeedIter1.size());
  RVec<float> yOrig(seed.fSeedIter1.size());
  RVec<float> zOrig(seed.fSeedIter1.size());
  RVec<float> qOrig(seed.fSeedIter1.size());
  Int_t nClUsed=0;
  for (uint i=0; i<seed.fSeedIter1.size(); i++) {
    if (seed.fSeedIter1[i].GetX()<kMinRadius) continue;
    alpha[nClUsed]= seed.fSeedIter1[i].GetY();
    xOrig[nClUsed]= seed.fSeedIter1[i].GetX();
    yOrig[nClUsed]= seed.fSeedIter1[i].GetY();
    zOrig[nClUsed]= seed.fSeedIter1[i].GetZ();
    qOrig[nClUsed]= seed.fSeedIter1[i].GetMax();
    nClUsed++;
  }
  float alphaMedian = TMath::Median(nClUsed,alpha.data());
  seed.fQMeanSeed1=TMath::Mean(nClUsed, qOrig.data());
  seed.fQMedianSeed1=TMath::Median(nClUsed, qOrig.data());

  for (uint i=0; i<seed.fSeedIter1.size(); i++) {
    float dAlpha = seed.fSeedIter1[i].GetY()-alphaMedian;
    float x,y,z;
    float R= seed.fSeedIter1[i].GetX();
    if (R<kMinRadius) continue;
    float dPhi= seed.fSeedIter1[i].GetY()-alphaMedian;
    x= cos(dPhi)*R;
    y= sin(dPhi)*R;
    z= seed.fSeedIter1[i].GetZ();
    seed.fRieman.AddPoint(x,y,z,1,1);
  }
  seed.fAlpha=alphaMedian;
  seed.fRieman.Update();
  AliRieman * delta = seed.fRieman.MakeResiduals();
  Double_t xref= TMath::Median(nClUsed,xOrig.data()), params[5], covar[15];
  seed.fRieman.GetExternalParameters(xref, params,covar);   // curvatue from 1/cm  to 1/GeV
  params[4]/=Bz*kB2C;
  covar[14]/=(Bz*kB2C)*(Bz*kB2C);
  seed.fParamSeed1=AliExternalTrackParam(xref, alphaMedian, params, covar);
  // accumulate ToT in the region above track
  for (uint i=0; i<seed.fSeedIter1.size(); i++) {
    float dAlpha = seed.fSeedIter1[i].GetY()-alphaMedian;   // not userd info for debugger
    float x,y,z;
    float R= seed.fSeedIter1[i].GetX();
    if (R<kMinRadius) continue;
    Double_t xyz[3];
    seed.fParamSeed1.GetXYZatR(R,Bz, xyz);
    //
    float dAlphaTrack=TMath::ATan2(xyz[1],xyz[0]);
    float phiTrack=dAlphaTrack+alphaMedian;
    AliTPCclusterMI &cluster=seed.fSeedIter1[i];
    seed.fSeedIter1Tot[i]=accumulateTOT(cluster, phiTrack,xyz[2], recoParam.fDeltaRPhiTot, recoParam.fZStepTot, recoParam.fZCutTot);
  }
  //seed.fRieman.Dump();
  delete delta;
}

/// in the refit track we will use kalman filter to take energy loss and multiple scattering into account
///    * Kalman parameters are local
///    * track momentum along trajectory will change
///    * not clear in which direction particle is flying In (add energy) or Out (subtract energy)
///  Algorithm:
///    * propagate track to next fKalmanStep (5?)  rows and assing clusters
///    * calculate robust estimator - e.g. median to characterize deltas in rphi, z direction
///    * update track in "median"position with median delta
///    * store track parameters in array with given step fKalman
///    * difference between in and out along trajectory to define "kinks" - split tracks
/// \param seed          - seed to propagate
/// \param recoParam     - reconstruction parameters
void refitTrackKalman(AliTPCseedMonopol & seed, const AliMonopolRecParam & recoParam){
  const float delta0=5;
  const float resol02=0.1*0.1;
  const float massProton= TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  const float kMIP=20;   // very rough aerage dEdx normalization to MIP - normaly angular and diffusion correction needed
  // create track parameters from Rieman step
  const int nTrackRows= nRows/recoParam.fKalmanStep;
  AliExternalTrackParam paramArray[nTrackRows];

  //
  float deltaRPhiStep[recoParam.fKalmanStep];
  float deltaZStep[recoParam.fKalmanStep];
  Int_t rowM=TMath::Nint(seed.fRowSeed/float(recoParam.fKalmanStep))*recoParam.fKalmanStep;
  float rmsyArray[nRows]={0};
  float rmszArray[nRows]={0};
  //
  int nClustersAll=0;
  int rowMin=nRows;    // first and last pad row where we assingend clusters
  int rowMax=0;
  AliExternalTrackParam paramInner0, paramOuter0;
  for (Int_t dir=-1; dir<=1; dir+=2) {
    int rowEnd = (dir < 0) ? 0 : nRows;
    AliExternalTrackParam param(seed.fParamSeed1);
    int nCustersDir=0;
    for (int iRow = rowM; iRow<nRows ;iRow += dir * recoParam.fKalmanStep) {
      if (iRow < 0) iRow = 0;
      if (iRow > nRows) iRow = nRows - 1;
      // find clusters in kalman steps
      Double_t yat, zat = 0;
      int counter = 0;
      double x[recoParam.fKalmanStep], dy[recoParam.fKalmanStep], dz[recoParam.fKalmanStep], row[recoParam.fKalmanStep], indexCl[recoParam.fKalmanStep];
      for (int dRow = 0; dRow < recoParam.fKalmanStep; dRow++) {
        // 1.) propagate to radius
        // 2.) correct for  material - energy loss
        // 3.) find nearest cluster and store in array
        int jRow = iRow + dir * dRow;
        if (jRow < 0 || jRow >= nRows) continue;
        float phi0 = param.GetY() / param.GetX() + param.GetAlpha();
        int i0 = findNearest(clusterArray[jRow], phi0, param.GetZ(), recoParam.fDeltaYseed1 + delta0, recoParam.fDeltaZseed1 + delta0);
        if (i0 < 0) continue;
        AliTPCclusterMI &cl0 = clusterArray[jRow][i0];
        param.GetYAt(cl0.GetX(), Bz, yat);
        param.GetZAt(cl0.GetX(), Bz, zat);
        float phi1 = yat / cl0.GetX() + param.GetAlpha();
        int i1 = findNearest(clusterArray[jRow], phi1, zat, recoParam.fDeltaYseed1, recoParam.fDeltaZseed1);
        if (i1 < 0) continue;
        AliTPCclusterMI &cl1 = clusterArray[jRow][i1];
        x[counter] = cl1.GetX();
        dy[counter] = (cl1.GetY() - param.GetAlpha()) * cl1.GetX() - yat;
        dz[counter] = cl1.GetZ() - zat;
        row[counter] = jRow;
        indexCl[counter] = i1;
        seed.fClustersTrack[jRow] = cl1;
        counter++;
        if (jRow > rowMax) rowMax = jRow;
        if (jRow < rowMin) rowMin = jRow;
        nClustersAll++;
        nCustersDir++;
      }
      if (counter < recoParam.fKalmanStepMin) {
        if (dir>0) paramOuter0=param;
        if (dir<0) paramInner0=param;
        break;
      }
      // calculate robust mean
      float xmean = TMath::Mean(counter, x);
      float dymean = TMath::Median(counter, dy);
      float dzmean = TMath::Median(counter, dz);
      float dyrms = TMath::RMS(counter, dy);
      float dzrms = TMath::RMS(counter, dz);
      param.PropagateTo(xmean, Bz);
      double pos[2] = {param.GetY() + dymean, param.GetZ() + dzmean};
      double cov[3] = {dyrms * dyrms / counter + resol02, 0, dzrms * dzrms / counter + resol02};
      param.Update(pos, cov);
      for (int dRow = 0; dRow < recoParam.fKalmanStep; dRow++){
        rmsyArray[iRow+dRow] = dyrms;
        rmszArray[iRow+dRow] = dzrms;
      }
      if (dir>0) paramOuter0=param;
      if (dir<0) paramInner0=param;
      if (iRow == 0 || iRow == nRows - 1) break;
    }
  }
  seed.fNclusters=nClustersAll;
  // calculate dEdx estimator
  float qMaxArray[nRows], qTotArray[nRows];
  int counterCl=0;
  for (int i=0; i<nRows;i++){
    if (seed.fClustersTrack[i].GetMax()==0) continue;
    qMaxArray[counterCl]=  seed.fClustersTrack[i].GetMax();
    qTotArray[counterCl]=  seed.fClustersTrack[i].GetQ();
    counterCl++;
  }
  seed.fQTotMedian=TMath::Median(counterCl,qTotArray);
  seed.fQMaxMedian=TMath::Median(counterCl,qMaxArray);
  seed.fRmsYMedian0=TMath::Median(rowMax-rowMin, & rmsyArray[rowMin]);
  seed.fRmsZMedian0=TMath::Median(rowMax-rowMin, & rmszArray[rowMin]);
  //
  // refit in and refit out - find kinks
  //
  Int_t debugSwitch=1;  //debug switch only for gdb short tracks debugging
  for (Int_t dir=-1; dir<=1; dir+=2) {
    int rowStart = (dir < 0) ? rowMax : rowMin;
    int rowEnd = (dir < 0) ? 0 : nRows;
    AliExternalTrackParam param = (dir < 0) ? paramOuter0 : paramInner0;
    param.ResetCovariance(5);
    Int_t rowStep=(dir<0)? -1:1;
    for (int iRow = rowStart; (iRow < nRows) && (iRow >= 0); iRow += rowStep) {
      AliTPCclusterMI &cl = seed.fClustersTrack[iRow];
      if (cl.GetQ() <= 0) continue;            // could be not assigned
      Double_t xOld[3], x[3];
      param.GetXYZ(xOld);
      param.PropagateTo(cl.GetX(), Bz);
      param.GetXYZ(x);
      Double_t stepLength = TMath::Sqrt((x[0] - xOld[0]) * (x[0] - xOld[0]) + (x[1] - xOld[1]) * (x[1] - xOld[1]) + (x[2] - xOld[2]) * (x[2] - xOld[2]));
      param.CorrectForMeanMaterialdEdx(stepLength * recoParam.fX0, -dir * stepLength * recoParam.fXRho, massProton, (seed.fQMaxMedian / kMIP) * recoParam.fdEdxMIP);    // assuming
      if (iRow % recoParam.fKalmanStep == 0) {    //
        if (dir < 0) seed.fParamArrayIn[iRow / recoParam.fKalmanStep] = param;
        if (dir > 0) seed.fParamArrayOut[iRow / recoParam.fKalmanStep] = param;
      }
      float errY2 = 0.5 * resol02 + 0.5 * seed.fRmsYMedian0 * seed.fRmsYMedian0 + 0.5 * rmsyArray[iRow] * rmsyArray[iRow];
      float errZ2 = 0.5 * resol02 + 0.5 * seed.fRmsZMedian0 * seed.fRmsZMedian0 + 0.5 * rmszArray[iRow] * rmszArray[iRow];
      float y = (cl.GetY() - param.GetAlpha()) * cl.GetX();
      double pos[2] = {y, cl.GetZ()};
      double cov[3] = {errY2, 0, errZ2};
      double chi2 = param.GetPredictedChi2(pos, cov);
      if (chi2 > recoParam.fChi2Cut) continue;
      param.Update(pos, cov);
      if (dir < 0) {
        seed.fChi2Inner += chi2;
        seed.fNClustersInner++;
        seed.fParamInner = param;
      } else {
        seed.fChi2Outer += chi2;
        seed.fNClustersOuter++;
        seed.fParamOuter = param;
      }
    }
    if (dir > 0 && seed.fNClustersInner > 20 && seed.fNClustersOuter < seed.fNClustersInner * 0.5 && debugSwitch == 1) {
      ::Error("", "Outer fit failed %d\t%d\t%d", nClustersAll, seed.fNClustersInner, seed.fNClustersOuter);
      debugSwitch = 0;
      //dir = 1;
    }
  }


}

//void findTracks(Int_t eventID,int qMaxSeedThr, int minClusters=4){

void findTracks(Int_t eventID, const AliMonopolRecParam &recoParam){
  /// loop over rows
  ///     WE ARE USING -  GetY() to store phi coordinate, GetX() to store radius  - to be be considered in calculations
  //
  for (Int_t iRow0=nRows-recoParam.fDRows-1; iRow0>=0; iRow0-=recoParam.fDRows){ //seeding row
      UInt_t iRowP=iRow0+recoParam.fDRows;
      UInt_t iRowM=iRow0-recoParam.fDRows;
      if(iRowP>=nRows) continue;
      if(iRowM>=nRows) continue;
      RVec<AliTPCclusterMI>& clusterRaw0 = clusterArray[iRow0];
      RVec<AliTPCclusterMI>& clusterRawP = clusterArray[iRowP];
      RVec<AliTPCclusterMI>& clusterRawM = clusterArray[iRowM];
      Int_t nSeedClusters=0;
      Int_t nSeeds0=0;
      Int_t nSeeds1=0;
      /// seeding at Row0
      for (UInt_t iCl0=0; iCl0<clusterRaw0.size(); iCl0++){
        if (clusterRaw0[iCl0].GetMax()<recoParam.fMinQSeedInit) continue;
        AliTPCclusterMI& cl0=clusterRaw0[iCl0];
        nSeedClusters++;
        Int_t indexP0 = clusterLowerBound(clusterRawP,cl0.GetZ()-recoParam.fDeltaZseed0);
        for (UInt_t iClP=indexP0; iClP<clusterRawP.size(); iClP++) {
          AliTPCclusterMI& clP=clusterRawP[iClP];
          if (clP.GetMax()<recoParam.fMinQSeedInit) continue;
          if (clP.GetZ()-cl0.GetZ()>recoParam.fDeltaZseed0) break;
          if (TMath::Abs(cl0.GetY()-clP.GetY())*cl0.GetX()>recoParam.fDeltaYseed0) continue;
          for (UInt_t  iClM=0; iClM<clusterRawM.size(); iClM++) {
            AliTPCclusterMI& clM=clusterRawM[iClM];
            if (TMath::Abs(cl0.GetZ()-clM.GetZ())>recoParam.fDeltaZseed1) continue;
            if (TMath::Abs(cl0.GetY()-clM.GetY())*cl0.GetX()>recoParam.fDeltaYseed1) continue;
            //
            if (TMath::Abs(0.5*(clP.GetZ()+clM.GetZ())-cl0.GetZ())>recoParam.fDeltaZseed1) continue;
            if (TMath::Abs(0.5*(clP.GetY()+clM.GetY())-cl0.GetY())*cl0.GetX()>recoParam.fDeltaYseed1) continue;
            nSeeds0++;
            //::Info("seedZ","%d\t%f\t%f\t%f\n",nSeeds, clM.GetZ(), cl0.GetZ(),clP.GetZ());
            //::Info("seedY","%d\t%f\t%f\t%f\n",nSeeds, clM.GetY(), cl0.GetY(),clP.GetY());
            //
            AliTPCseedMonopol seed;
            seed.fRowSeed=clM.GetRow();
            seed.fSeedIter1.resize(1+(recoParam.fDRows+2)*2);
            seed.fSeedIter1Tot.resize(1+(recoParam.fDRows+2)*2);
            seed.fClustersTrack.resize(nRows);
            seed.fClustersTot.resize(nRows);
            seed.fParamArrayIn.resize(nRows/recoParam.fKalmanStep+1);
            seed.fParamArrayOut.resize(nRows/recoParam.fKalmanStep+1);
            int nClustersSeed1 = makeSeed1(clM,cl0,clP,recoParam,seed);
            (*pcstream)<<"seed0"<<
                       "eventID="<<eventID<<                   // event ID
                       "seedID="<<nSeeds0<<                     // see ID
                       "nClusterSeed1="<<nClustersSeed1<<      // number of clusters after search
                       "cl0.="<<&cl0<<                         // center cluster for seed
                       "clP.="<<&clP<<                         // upper radius seed cluster
                       "clM.="<<&clM<<                         // lower radius  seed cluster
                       "\n";

            if (nClustersSeed1<recoParam.fMinClustersIter1) continue;
            ::Info("seedZ","%d\t%f\t%f\t%f\t%d\n",nSeeds1, clM.GetZ(), cl0.GetZ(),clP.GetZ(),nClustersSeed1);
            ::Info("seedY","%d\t%f\t%f\t%f\t%d\n",nSeeds1, clM.GetY(), cl0.GetY(),clP.GetY(),nClustersSeed1);
            seed.fnSeedIter1=nClustersSeed1;
            refitTrackRiemanSeed1(seed,recoParam);
            refitTrackKalman(seed, recoParam);
            nSeeds1++;
            AliRieman * rieman1Delta = seed.fRieman.MakeResiduals();
            (*pcstream)<<"seed1"<<
              "eventID="<<eventID<<                   // event ID
              "seedID0="<<nSeeds0<<                     // see ID 0
              "seedID1="<<nSeeds1<<                     // see ID 1
              "nClustersSeed1="<<nClustersSeed1<<        // number of cluster assoctated to seed
              "cl0.="<<&cl0<<                         // center cluster for seed
              "clP.="<<&clP<<                         // upper radius seed cluster
              "clM.="<<&clM<<                         // lower radius  seed cluster
              "seed.="<<&seed<<                       //
              "rieman1Delta.="<<rieman1Delta<<        // rieman delta
              "\n";
              delete rieman1Delta;
          }
        }
      }
      ::Info("findTracks ", "Event\t%d\tnSeeds0\t %d\t nSeeds\t%d",eventID, nSeeds0,nSeeds1);
  }
}

///
/// \param maxQ           - threshold for cluster charge
/// \param qMaxSeedThr    - trheshold to start seeding
void doReconstruction(const AliMonopolRecParam &recoParam, Int_t nEvents=-1, Int_t firstEvent=0){   /// 500 ADC seeding
  const Int_t kMaxClustersRow=10000;
  initTree();
  initClusters();
  pcstream=new TTreeSRedirector("trackingDebug.root","recreate");
  if (nEvents<0) nEvents=eventIndeces.GetSize()-1;
  if (firstEvent<0 | firstEvent<nEvents) firstEvent=0;
  //loadClusters(0, maxQ);
  for (int iEvent=firstEvent; iEvent<nEvents; iEvent++) {
    ::Info("Processing event","%d",iEvent);
    loadClusters(iEvent, recoParam);
    findTracks(iEvent, recoParam);
  }
  delete pcstream;
}


void loadSeeds() {
  TFile *f = TFile::Open("trackingDebug.root");
  //
  treeSeeds0 = (TTree *) f->Get("seed0");
  treeSeeds0->SetAlias("gy","cos(cl0.fY)*cl0.fX");
  treeSeeds0->SetAlias("gx","sin(cl0.fY)*cl0.fX");
  treeSeeds0->SetAlias("gySeed","cos(seed.fSeedIter1.fData.fY)*seed.fSeedIter1.fData.fX");
  treeSeeds0->SetAlias("gxSeed","sin(seed.fSeedIter1.fData.fY)*seed.fSeedIter1.fData.fX");
  treeSeeds0->SetAlias("gzSeed","seed.fSeedIter1.fData.fZ");
  treeSeeds0->SetMarkerStyle(21);
  treeSeeds0->SetMarkerSize(0.5);
  //
  treeSeeds1 = (TTree *) f->Get("seed1");
  treeSeeds1->SetAlias("gySeed","cos(seed.fSeedIter1.fData.fY)*seed.fSeedIter1.fData.fX");
  treeSeeds1->SetAlias("gxSeed","sin(seed.fSeedIter1.fData.fY)*seed.fSeedIter1.fData.fX");
  treeSeeds1->SetAlias("gzSeed","seed.fSeedIter1.fData.fZ");
  treeSeeds1->SetAlias("rSeed","seed.fSeedIter1.fData.fX");
  treeSeeds1->SetAlias("qSeed","seed.fSeedIter1.fData.fMax");
  treeSeeds1->SetAlias("chi2N","sqrt(seed.fRieman.fChi2/seed.fRieman.fN)");
  treeSeeds1->SetAlias("fQMeanSeed1","seed.fQMeanSeed1");
  treeSeeds1->SetAlias("fQMedianSeed1","seed.fQMedianSeed1");
  treeSeeds1->SetAlias("chi2YN","sqrt(seed.fRieman.fChi2Y/seed.fRieman.fN)");
  treeSeeds1->SetAlias("chi2ZN","sqrt(seed.fRieman.fChi2Z/seed.fRieman.fN)");
  //
  treeSeeds1->SetAlias("seed1Tot","seed.fSeedIter1Tot.fData");
  // full track aliases

  // kink aliases
  treeSeeds1->SetAlias("nclMin","min(seed.fNClustersInner,seed.fNClustersOuter)");
  treeSeeds1->SetAlias("trackletOK","(seed.fParamArrayOut.fData[].fX*seed.fParamArrayOut.fData[].fX)>0");
  treeSeeds1->SetAlias("deltaP2","(seed.fParamArrayOut.fData[].fP[2]-seed.fParamArrayIn.fData[].fP[2])");   //
  treeSeeds1->SetAlias("deltaP3","(seed.fParamArrayOut.fData[].fP[3]-seed.fParamArrayIn.fData[].fP[3])");   //
  treeSeeds1->SetAlias("pullP2","(seed.fParamArrayOut.fData[].fP[2]-seed.fParamArrayIn.fData[].fP[2])/sqrt(seed.fParamArrayOut[].fData.fC[5]+seed.fParamArrayIn.fData[].fC[5])");   //
  treeSeeds1->SetAlias("pullP3","(seed.fParamArrayOut.fData[].fP[3]-seed.fParamArrayIn.fData[].fP[3])/sqrt(seed.fParamArrayOut.fData[].fC[9]+seed.fParamArrayIn.fData[].fC[9])");   //

  treeSeeds1->SetMarkerStyle(21);
  treeSeeds1->SetMarkerSize(0.5);

  //AliTPCseedMonopol seed,*ppseed=&seed;
  //treeSeeds1->SetBranchAddress("seed.",&ppseed);

}

void monoAnalysis(Int_t minQseed=200){
  AliMonopolRecParam recoParam;
  recoParam.fMinQSeedInit=200;
  doReconstruction(recoParam,-1);
}
