/// \brief  Analysis of the triggerred ESD data
///
/*
 gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
 .L $NOTES/JIRA/ATO-432/code/drawHighdEdx.C+
 // LoadChain("filtered.list",10000);
 // MakeHistograms()
 // MakeMaps();
 // MakeNDLocalRegressions()
 //
 */

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "AliXRDPROOFtoolkit.h"

#include "AliTreePlayer.h"
#include "TTreeStream.h"
#include "TROOT.h"
#include "THn.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TStatToolkit.h"
#include "TPRegexp.h"
#include "AliNDLocalRegression.h"
#include "AliESDtrack.h"
#include "AliAnalysisTaskFilteredTree.h"
#include "AliMathBase.h"
#include "AliTMinuitToolkit.h"
#include "TStopwatch.h"


TChain * tree= NULL;
TChain * treeEv=NULL;
TChain * treeV0=NULL;
TChain * treeV01=NULL;
TChain * treePt=NULL;
TChain * treeCosmic=NULL;
TChain * treeCosmic1=NULL;
TStopwatch timer;
//
TTree * treeMap= NULL;
TTree * treeMapID= NULL;
AliNDLocalRegression *regRMS[12]={0};
AliNDLocalRegression *regMeanG[12]={0};

AliESDtrack *track=NULL;

const char *regNames[12]={"hisRatioTotMax0Dist","hisRatioTotMax1Dist","hisRatioTotMax2Dist","hisRatioTotMax2Dist",\
        "hisRatioTot01Dist", "hisRatioTot12Dist","hisRatioTot02Dist","",\
        "hisRatioMax01Dist", "hisRatioMax12Dist","hisRatioMax02Dist",""};


void LoadChain(const char *filteredList="filtered.list",Int_t nChunks=1);
void MakeHistograms(Int_t nEvents=100000);
void MakeMaps();
void LoadSummaryTrees();
TTree* LoadSummaryTree(TString sregexp="hisRatio*", Bool_t addName=kFALSE);
void LoadNDLocal();
void makeNDLocalFit(TString varName, TString customCut, TString errorVar);
void MakeNDLocalRegressions();
TF1 *likeGausCachy = new TF1("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike,-10,10,2);

///  \brief Set aliases
/// \param tree
void SetAliases(TChain * treeA, Int_t v0Index=-1, Int_t cosmicIndex=-1){
  treeA->SetAlias("ratioTotMax0", "fTPCdEdxInfo.GetSignalTot(0)/fTPCdEdxInfo.GetSignalMax(0)");
  treeA->SetAlias("ratioTotMax1", "fTPCdEdxInfo.GetSignalTot(1)/fTPCdEdxInfo.GetSignalMax(1)");
  treeA->SetAlias("ratioTotMax2", "fTPCdEdxInfo.GetSignalTot(2)/fTPCdEdxInfo.GetSignalMax(2)");
  treeA->SetAlias("ratioTotMax3", "fTPCdEdxInfo.GetSignalTot(3)/fTPCdEdxInfo.GetSignalMax(3)");
  //
  treeA->SetAlias("logRatioTot03", "log(fTPCdEdxInfo.GetSignalTot(0)/fTPCdEdxInfo.GetSignalTot(3))");
  treeA->SetAlias("logRatioTot13", "log(fTPCdEdxInfo.GetSignalTot(1)/fTPCdEdxInfo.GetSignalTot(3))");
  treeA->SetAlias("logRatioTot23", "log(fTPCdEdxInfo.GetSignalTot(2)/fTPCdEdxInfo.GetSignalTot(3))");
  treeA->SetAlias("logRatioMax03", "log(fTPCdEdxInfo.GetSignalMax(0)/fTPCdEdxInfo.GetSignalMax(3))");  //
  treeA->SetAlias("logRatioMax13", "log(fTPCdEdxInfo.GetSignalMax(1)/fTPCdEdxInfo.GetSignalMax(3))");  //
  treeA->SetAlias("logRatioMax23", "log(fTPCdEdxInfo.GetSignalMax(2)/fTPCdEdxInfo.GetSignalMax(3))");  //
  //
  treeA->SetAlias("logRatioTot01", "log(fTPCdEdxInfo.GetSignalTot(0)/fTPCdEdxInfo.GetSignalTot(1))");
  treeA->SetAlias("logRatioTot12", "log(fTPCdEdxInfo.GetSignalTot(1)/fTPCdEdxInfo.GetSignalTot(2))");
  treeA->SetAlias("logRatioTot02", "log(fTPCdEdxInfo.GetSignalTot(0)/fTPCdEdxInfo.GetSignalTot(2))");
  treeA->SetAlias("logRatioMax01", "log(fTPCdEdxInfo.GetSignalMax(0)/fTPCdEdxInfo.GetSignalMax(1))");  //
  treeA->SetAlias("logRatioMax12", "log(fTPCdEdxInfo.GetSignalMax(1)/fTPCdEdxInfo.GetSignalMax(2))");  //
  treeA->SetAlias("logRatioMax02", "log(fTPCdEdxInfo.GetSignalMax(0)/fTPCdEdxInfo.GetSignalMax(2))");  //
  //
  treeA->SetAlias("nCross0","esdTrack.fTPCdEdxInfo.GetNumberOfCrossedRows(0)");
  treeA->SetAlias("nCross1","esdTrack.fTPCdEdxInfo.GetNumberOfCrossedRows(1)");
  treeA->SetAlias("nCross2","esdTrack.fTPCdEdxInfo.GetNumberOfCrossedRows(2)");
  treeA->SetAlias("nFraction0","esdTrack.GetTPCClusterInfo(1,0,0,62)");
  treeA->SetAlias("nFraction1","esdTrack.GetTPCClusterInfo(1,0,63,127)");
  treeA->SetAlias("nFraction2","esdTrack.GetTPCClusterInfo(1,0,127,159)");
  treeA->SetAlias("nFraction3","esdTrack.GetTPCClusterInfo(1,0,0,159)");
  treeA->SetAlias("n3Fraction0","esdTrack.GetTPCClusterInfo(3,0,0,62)");
  treeA->SetAlias("n3Fraction1","esdTrack.GetTPCClusterInfo(3,0,63,127)");
  treeA->SetAlias("n3Fraction2","esdTrack.GetTPCClusterInfo(3,0,127,159)");
  treeA->SetAlias("n3Fraction3","esdTrack.GetTPCClusterInfo(3,0,0,159)");

  //
  treeA->SetAlias("pileUpZ","esdTrack.fzTPC*sign(esdTrack.fP[3])"); // to be used for pileup correction  (pp) - better definition to be made A side/C side by delta Z
  treeA->SetAlias("nclCut", "esdTrack.fTPCdEdxInfo.GetNumberOfCrossedRows(0)>15&&esdTrack.fTPCdEdxInfo.GetNumberOfCrossedRows(1)>15&&esdTrack.fTPCdEdxInfo.GetNumberOfCrossedRows(2)>15");
  treeA->SetAlias("nclCutGold", "esdTrack.GetTPCClusterInfo(2,1,0,63)>30&&esdTrack.GetTPCClusterInfo(2,1,128,159)>20");
  //
  treeA->SetAlias("dcaCut", "abs(esdTrack.fdTPC)<6&&abs(esdTrack.fzTPC)<30");      //DCAxy and DCAz cut
  treeA->SetAlias("resolCut","sqrt(esdTrack.fIp.fC[14])/abs(esdTrack.fIp.fP[4])<0.25");
  treeA->SetAlias("chi2Cut","esdTrack.fTPCchi2/esdTrack.fTPCncls<5");               // high dEdx has lower chi2 ~1/dEdx
  treeA->SetAlias("chi2TPC","esdTrack.fTPCchi2/esdTrack.fTPCncls");
  treeA->SetAlias("ratioMax3","50/fTPCdEdxInfo.GetSignalMax(3)");
  treeA->SetAlias("qP","sign(esdTrack.fIp.fP[4])/(esdTrack.fIp.P())");
  treeA->SetAlias("qPt","sign(esdTrack.fIp.fP[4])/(esdTrack.fIp.Pt())");
  treeA->SetAlias("tgl","esdTrack.fIp.fP[3]");
  treeA->SetAlias("atgl","abs(esdTrack.fIp.fP[3]+0)");
  treeA->SetAlias("logdEdxMax","log(fTPCdEdxInfo.GetSignalMax(3)+0)");
  treeA->SetAlias("logdEdxTot","log(fTPCdEdxInfo.GetSignalTot(3)+0)");
  treeA->SetAlias("dEdxMax","(fTPCdEdxInfo.GetSignalMax(3)+0)");
  treeA->SetAlias("sdEdxMax","sqrt(fTPCdEdxInfo.GetSignalMax(3)+0)");
  treeA->SetAlias("sdEdxMax0","sqrt(fTPCdEdxInfo.GetSignalMax(0)+0)");
  treeA->SetAlias("sdEdxMax1","sqrt(fTPCdEdxInfo.GetSignalMax(1)+0)");
  treeA->SetAlias("sdEdxMax2","sqrt(fTPCdEdxInfo.GetSignalMax(2)+0)");
  treeA->SetAlias("dEdxMax0","(fTPCdEdxInfo.GetSignalMax(0)+0)");
  treeA->SetAlias("dEdxMax1","(fTPCdEdxInfo.GetSignalMax(1)+0)");
  treeA->SetAlias("dEdxMax2","(fTPCdEdxInfo.GetSignalMax(2)+0)");

  treeA->SetAlias("pileUpCut","(esdTrack.fTRDncls/15.+esdTrack.fITSncls)>2");
  // First rough cuts - ratios to be parameterized - n sigma cut to be used
  treeA->SetAlias("logRatioCutTot","abs(logRatioTot03)<0.5&&abs(logRatioTot13)<0.5&&abs(logRatioTot23)<0.5");
  treeA->SetAlias("logRatioCutMax","abs(logRatioMax03)<0.5&&abs(logRatioMax13)<0.5&&abs(logRatioMax23)<0.5");
  treeA->SetAlias("multA","ntracks");   // by mistake ntracks was not writen diretly in some trees - event information used
  if (v0Index<0) {
    treeA->SetAlias("ratioMaxCut", "(ratioTotMax0-1)<0.3&&(ratioTotMax1-1)<0.3&&abs(ratioTotMax2-1)<0.3");              // checked
  }else{
    treeA->SetAlias("ratioMaxCut", "(ratioTotMax0-1)<0.35&&(ratioTotMax1-1)<0.35&&(ratioTotMax2-1)<0.35");              // checked
  }
  treeA->SetAlias("dEdxCutRough","logRatioCutTot&&logRatioCutMax&&ratioMaxCut");
  // pulls
  //
  treeA->SetAlias("pullTotMax0","(ratioTotMax0-hisRatioTotMax0Dist.meanGFit)/hisRatioTotMax0Dist.rmsGFit");
  treeA->SetAlias("pullTotMax1","(ratioTotMax1-hisRatioTotMax1Dist.meanGFit)/hisRatioTotMax1Dist.rmsGFit");
  treeA->SetAlias("pullTotMax2","(ratioTotMax2-hisRatioTotMax2Dist.meanGFit)/hisRatioTotMax2Dist.rmsGFit");
  treeA->SetAlias("pullTotMax","sqrt(pullTotMax0**2+pullTotMax1**2+pullTotMax2**2)");
  treeA->SetAlias("pullTot01","(logRatioTot01-hisRatioTot01Dist.meanGFit)/hisRatioTot01Dist.rmsGFit");
  treeA->SetAlias("pullTot12","(logRatioTot12-hisRatioTot12Dist.meanGFit)/hisRatioTot12Dist.rmsGFit");
  //
  // nsgima cut to come ....
  treeA->SetAlias("isSelected0","nclCut&&dcaCut&&resolCut&&chi2Cut");
  treeA->SetAlias("isSelected","nclCut&&dcaCut&&resolCut&&chi2Cut&&dEdxCutRough");
  treeA->SetAlias("downscale","(1+0)");
  if (v0Index>-1){

    TList * aliasList=treeA->GetListOfAliases();
    for (Int_t i=0;i<aliasList->GetEntries();i++) {
      TString aliasName=aliasList->At(i)->GetName();
      TString aliasTitle=aliasList->At(i)->GetTitle();
      if (v0Index==0) {
        aliasTitle.ReplaceAll("fTPCdEdxInfo.", "track0.fTPCdEdxInfo.");
        aliasTitle.ReplaceAll("esdTrack.", "track0.");
        treeA->SetAlias("downscale","sqrt(track0.fD^2+track0.fZ^2)>(20*rndm/(1+5*((type==1)*abs(track1.fTPCsignal-78)<10)))");

      }else{
        aliasTitle.ReplaceAll("fTPCdEdxInfo.", "track1.fTPCdEdxInfo.");
        aliasTitle.ReplaceAll("esdTrack.", "track1.");
        treeA->SetAlias("downscale","sqrt(track1.fD^2+track1.fZ^2)>(5*rndm/(1+3*((type==1)*abs(track0.fTPCsignal-78)<10)))");
      }
      treeA->SetAlias(aliasName.Data(),aliasTitle.Data());
    }
  }
  if (cosmicIndex>-1){
    TList * aliasList=treeA->GetListOfAliases();
    for (Int_t i=0;i<aliasList->GetEntries();i++) {
      TString aliasName=aliasList->At(i)->GetName();
      TString aliasTitle=aliasList->At(i)->GetTitle();
      if (cosmicIndex==0) {
        aliasTitle.ReplaceAll("fTPCdEdxInfo.", "t0.fTPCdEdxInfo.");
        aliasTitle.ReplaceAll("esdTrack.", "t0.");
      }else{
        aliasTitle.ReplaceAll("fTPCdEdxInfo.", "t1.fTPCdEdxInfo.");
        aliasTitle.ReplaceAll("esdTrack.", "t1.");
      }
      treeA->SetAlias(aliasName.Data(),aliasTitle.Data());
    }
  }


}

/*

*/
///'brief Load chain and define aliases
//// \param list
/// \param nChunks
void LoadChain(const char *filteredList,Int_t nChunks){
  tree = AliXRDPROOFtoolkit::MakeChain(filteredList, "dEdx", 0, nChunks);
  treeEv = AliXRDPROOFtoolkit::MakeChain(filteredList, "eventInfoV0", 0, nChunks);
  treePt = AliXRDPROOFtoolkit::MakeChain(filteredList, "highPt", 0, nChunks);
  treeV0=AliXRDPROOFtoolkit::MakeChain(filteredList, "V0s", 0, nChunks);
  treeV01=AliXRDPROOFtoolkit::MakeChain(filteredList, "V0s", 0, nChunks);
  treeCosmic = AliXRDPROOFtoolkit::MakeChain(filteredList, "CosmicPairs", 0, nChunks);
  treeCosmic1 = AliXRDPROOFtoolkit::MakeChain(filteredList, "CosmicPairs", 0, nChunks);
  //
  treeEv->BuildIndex("gid");
  tree->BuildIndex("gid");
  tree->AddFriend(treeEv,"Ev");
  //
  AliAnalysisTaskFilteredTree::SetDefaultAliasesV0(treeV0);
  AliAnalysisTaskFilteredTree::SetDefaultAliasesV0(treeV01);
  SetAliases(tree);
  SetAliases(treeV0,0);   // v0 positive tracks
  SetAliases(treeV01,1);  // v0 negative tracks
  SetAliases(treePt);
  SetAliases(treeCosmic,-1,0);   // cosmic+ secondary
  SetAliases(treeCosmic1,-1,1);  // cosmic+secondries
  //
  //tree->GetBranch("esdTrack.")->SetAddress(&param);
  //tree->GetBranch("esdTrack.fTPCdEdxInfo")->SetAddress(&dEdxInfo);
}



/// TODO - for high pt tracks use combined momentum  instead of the TPC momentum
/// TODO - mean parameterization for period is not valid per run - rejected outlier runs from merging - to be done automatically
///      - each run to be QAed
///      - possible explanation:
///            - n-primary tracks versus occupancy
///            - calibration for MIP - divergence for the high dEdx particles ...
///\brief Make/fill histograms
/// \param nEvents - for test purposes number of events can be adjusted
void MakeHistograms(Int_t nEvents){
  // 2. Define histograms
  TString hisString = "";
  hisString += "logRatioTot03:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisRatioTot03(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioTot13:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisRatioTot13(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioTot23:atgl:multA:pileUpZ:sdEdxMax:#isSelected&&abs(qPt)<2.5>>hisRatioTot23(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioMax03:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisRatioMax03(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioMax13:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisRatioMax13(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioMax23:atgl:multA:pileUpZ:sdEdxMax:#isSelected&&abs(qPt)<2.5>>hisRatioMax23(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  //
  hisString += "logRatioTot01:atgl:multA:pileUpZ:sdEdxMax2:#isSelected>>hisRatioTot01(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioTot12:atgl:multA:pileUpZ:sdEdxMax0:#isSelected>>hisRatioTot12(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioTot02:atgl:multA:pileUpZ:sdEdxMax1:#isSelected&&abs(qPt)<2.5>>hisRatioTot02(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioMax01:atgl:multA:pileUpZ:sdEdxMax2:#isSelected>>hisRatioMax01(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioMax12:atgl:multA:pileUpZ:sdEdxMax0:#isSelected>>hisRatioMax12(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "logRatioMax02:atgl:multA:pileUpZ:sdEdxMax1:#isSelected&&abs(qPt)<2.5>>hisRatioMax02(100,-1.0,1.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  //
  hisString += "ratioTotMax0:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisRatioTotMax0(100,0.5,2.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "ratioTotMax1:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisRatioTotMax1(100,0.5,2.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "ratioTotMax2:atgl:multA:pileUpZ:sdEdxMax:#isSelected&&abs(qPt)<2.5>>hisRatioTotMax2(100,0.5,2.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "ratioTotMax3:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisRatioTotMax3(100,0.5,2.0,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  //
  hisString += "chi2TPC:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisChi2TPC(100,0.0001,6,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "nCross0:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisCross0(63,0.0,63,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "nCross1:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisCross1(63,0.0,63,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "nCross2:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisCross2(63,0.0,63,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "nFraction0:atgl:multA:pileUpZ:sdEdxMax1:#isSelected>>hisFraction0(120,0.5,1.1,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "nFraction1:atgl:multA:pileUpZ:sdEdxMax0:#isSelected>>hisFraction1(120,0.5,1.1,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "nFraction2:atgl:multA:pileUpZ:sdEdxMax1:#isSelected>>hisFraction2(120,0.5,1.1,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "nFraction3:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>hisFraction3(120,0.5,1.1,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "n3Fraction0:atgl:multA:pileUpZ:sdEdxMax1:#isSelected>>his3Fraction0(120,0.5,1.1,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "n3Fraction1:atgl:multA:pileUpZ:sdEdxMax0:#isSelected>>his3Fraction1(120,0.5,1.1,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "n3Fraction2:atgl:multA:pileUpZ:sdEdxMax1:#isSelected>>his3Fraction2(120,0.5,1.1,20,0,1,20,0,20000,5,-50,50,20,6,32);";
  hisString += "n3Fraction3:atgl:multA:pileUpZ:sdEdxMax:#isSelected>>his3Fraction3(120,0.5,1.1,20,0,1,20,0,20000,5,-50,50,20,6,32);";
 //
  hisString+="logdEdxMax:qP:atgl:multA:pileUpZ:#isSelected>>hisdEdxMax(300,3,8,300,-2,2,10,-1,1,5,0,20000,5,-50,50);";
  hisString+="logdEdxTot:qP:atgl:multA:pileUpZ:#isSelected>>hisdEdxTot(300,3,8,300,-2,2,10,-1,1,5,0,20000,5,-50,50);";
  hisString+="logdEdxMax:qP:atgl:multA:pileUpZ:#isSelected&&pileUpCut>>hisdEdxMaxP(300,3,8,300,-2,2,10,-1,1,5,0,20000,5,-50,50);";
  hisString+="logdEdxTot:qP:atgl:multA:pileUpZ:#isSelected&&pileUpCut>>hisdEdxTotP(300,3,8,300,-2,2,10,-1,1,5,0,20000,5,-50,50);";
  //
  //
  if (nEvents>tree->GetEntries()){
    nEvents=20*tree->GetEntries();
  }
  for (Int_t iTree=0; iTree<4; iTree++) {
    TTree *treeR=NULL;
    if (iTree==0) treeR=tree;
    if (iTree==1) treeR=treeV0;
    if (iTree==2) treeR=treeV01;
    if (iTree==3) treeR=treePt;
    TObjArray *hisArray = AliTreePlayer::MakeHistograms(treeR, hisString, "nclCut&&chi2Cut&&downscale&&abs(qPt)<1.25&&nclCutGold&&pileUpCut", 0, nEvents, -1, 15);
    for (Int_t i = 0; i < hisArray->GetEntries(); i++) gROOT->Add(hisArray->At(i));
    TTreeSRedirector *pcstream = new TTreeSRedirector(TString::Format("hisOutputs%d.root",iTree).Data(), "recreate");
    pcstream->GetFile()->cd();
    for (Int_t iHis = 0; iHis < hisArray->GetEntries(); iHis++) {
      hisArray->At(iHis)->Write();
    }
    delete pcstream;
  }
  // clean electron sample - to do  - make rough fits for electron position +-5 sigma cut min 20 %
  Int_t entries=treeV0->Draw("track0.fTPCsignal","isSelected&&(track0.P()>0.3&&track0.P()<0.9)&&abs(track0.fTPCsignal-80)<20&&abs(track1.fTPCsignal-80)<20&&ELike0>0.2","",100000);
  Double_t medianEl, rmsEl;
  AliMathBase::EvaluateUni(entries,treeV0->GetV1(),medianEl, rmsEl,0.8*entries);
  treeV0->SetAlias("medianEl",TString::Format("(%f+0)",medianEl).Data());
  treeV0->SetAlias("rmsEl",TString::Format("(%f+0)",rmsEl).Data());
  treeV0->SetAlias("cleanGamma","ELike>0.4&&abs(track1.fTPCsignal-medianEl)<2*rmsEl&&acos(v0.fPointAngle)<0.006");
  treeV01->SetAlias("cleanGamma","ELike>0.4&&&&abs(track0.fTPCsignal-medianEl)<(5.*rmsEl)&&acos(v0.fPointAngle)<0.006");
  treeV0->SetAlias("cleanEl","abs(track0.fTPCsignal-medianEl)<(5*rmsEl)");
  treeV01->SetAlias("cleanEl","abs(track0.fTPCsignal-medianEl)<(5*rmsEl)");
  entries=treeV0->Draw("track0.fTPCdEdxInfo.GetSignalTot(0):track0.fTPCdEdxInfo.GetSignalTot(1):track0.fTPCdEdxInfo.GetSignalTot(2)","isSelected&&(track0.P()>0.3&&track0.P()<0.9)&&cleanGamma&&cleanEl","goff",100000);
  treeV0->SetAlias("medianEl0",TString::Format("(%f+0)",TMath::Median(entries,treeV0->GetVal(0))).Data());
  treeV0->SetAlias("medianEl1",TString::Format("(%f+0)",TMath::Median(entries,treeV0->GetVal(1))).Data());
  treeV0->SetAlias("medianEl2",TString::Format("(%f+0)",TMath::Median(entries,treeV0->GetVal(2))).Data());
  ::Info("MakeHistograms","Electron position:\t%f",medianEl);
  ::Info("MakeHistograms","Electron rms:\t%f", rmsEl);
  //
  treeV0->SetAlias("cleanEl0","abs(track0.fTPCdEdxInfo.GetSignalTot(1)-medianEl1)<(2.*rmsEl)&&abs(track0.fTPCdEdxInfo.GetSignalTot(2)-medianEl2)<(2.0*rmsEl)");
  treeV0->SetAlias("cleanEl1","abs(track0.fTPCdEdxInfo.GetSignalTot(0)-medianEl0)<(2.*rmsEl)&&abs(track0.fTPCdEdxInfo.GetSignalTot(2)-medianEl2)<(2.0*rmsEl)");
  treeV0->SetAlias("cleanEl2","abs(track0.fTPCdEdxInfo.GetSignalTot(0)-medianEl0)<(2.*rmsEl)&&abs(track0.fTPCdEdxInfo.GetSignalTot(1)-medianEl1)<(2.0*rmsEl)");
  //
  treeV0->SetAlias("logSignal","log(track0.fTPCsignal+0)");
  treeV0->SetAlias("logMax0","log(track0.fTPCdEdxInfo.GetSignalMax(0)+0)");
  treeV0->SetAlias("logMax1","log(track0.fTPCdEdxInfo.GetSignalMax(1)+0)");
  treeV0->SetAlias("logMax2","log(track0.fTPCdEdxInfo.GetSignalMax(2)+0)");
  treeV0->SetAlias("logMax3","log(track0.fTPCdEdxInfo.GetSignalMax(3)+0)");
  treeV0->SetAlias("logTot0","log(track0.fTPCdEdxInfo.GetSignalTot(0)+0)");
  treeV0->SetAlias("logTot1","log(track0.fTPCdEdxInfo.GetSignalTot(1)+0)");
  treeV0->SetAlias("logTot2","log(track0.fTPCdEdxInfo.GetSignalTot(2)+0)");
  treeV0->SetAlias("logTot3","log(track0.fTPCdEdxInfo.GetSignalTot(2)+0)");
  //
  hisString="";
  hisString += "logSignal:tgl:multA:#cleanEl0>>hiscElectronSignal(100,3.5,5,40,-1,1,20,0,20000);";
  hisString += "logMax0:tgl:multA:#cleanEl0>>hiscElectronMax0(100,3.5,5,40,-1,1,20,0,20000);";
  hisString += "logMax1:tgl:multA:#cleanEl1>>hiscElectronMax1(100,3.5,5,40,-1,1,20,0,20000);";
  hisString += "logMax2:tgl:multA:#cleanEl2>>hiscElectronMax2(100,3.5,5,40,-1,1,20,0,20000);";
  hisString += "logTot0:tgl:multA:#cleanEl0>>hiscElectronTot0(100,3.5,5,40,-1,1,20,0,20000);";
  hisString += "logTot1:tgl:multA:#cleanEl1>>hiscElectronTot1(100,3.5,5,40,-1,1,20,0,20000);";
  hisString += "logTot2:tgl:multA:#cleanEl2>>hiscElectronTot2(100,3.5,5,40,-1,1,20,0,20000);";

  TObjArray *hisArray = AliTreePlayer::MakeHistograms(treeV0, hisString, "isSelected&&abs(qP)>0.4&&abs(qP)<2&&cleanGamma&&cleanEl", 0, nEvents, -1, 15);
  for (Int_t i = 0; i < hisArray->GetEntries(); i++) gROOT->Add(hisArray->At(i));
  TTreeSRedirector *pcstream = new TTreeSRedirector("hisOutputsEl.root", "recreate");
  pcstream->GetFile()->cd();
  for (Int_t iHis = 0; iHis < hisArray->GetEntries(); iHis++) {
      hisArray->At(iHis)->Write();
  }
  delete pcstream;
  delete hisArray;
  //
  // clean pion sample
  //
  treeV0->SetAlias("dEdxPi","AliMathBase::BetheBlochAleph(track0.fIp.P()/0.13957)");
  treeV0->SetAlias("logSignalPi","log(track0.fTPCsignal/dEdxPi+0)");
  treeV0->SetAlias("cleanPion","abs(logSignalPi-log(50.))<0.2");
  treeV0->SetAlias("logMaxPi0","log(track0.fTPCdEdxInfo.GetSignalMax(0)/dEdxPi+0)");
  treeV0->SetAlias("logMaxPi1","log(track0.fTPCdEdxInfo.GetSignalMax(1)/dEdxPi+0)");
  treeV0->SetAlias("logMaxPi2","log(track0.fTPCdEdxInfo.GetSignalMax(2)/dEdxPi+0)");
  treeV0->SetAlias("logMaxPi3","log(track0.fTPCdEdxInfo.GetSignalMax(3)/dEdxPi+0)");
  treeV0->SetAlias("logTotPi0","log(track0.fTPCdEdxInfo.GetSignalTot(0)/dEdxPi+0)");
  treeV0->SetAlias("logTotPi1","log(track0.fTPCdEdxInfo.GetSignalTot(1)/dEdxPi+0)");
  treeV0->SetAlias("logTotPi2","log(track0.fTPCdEdxInfo.GetSignalTot(2)/dEdxPi+0)");
  treeV0->SetAlias("logTotPi3","log(track0.fTPCdEdxInfo.GetSignalTot(2)/dEdxPi+0)");
  hisString="";
  hisString += "logSignalPi:tgl:multA:#cleanPion>>hiscPionSignal(100,3.3,4.5,40,-1,1,20,0,20000);";
  hisString += "logMaxPi0:tgl:multA:#cleanPion>>hiscPionMax0(100,3.3,4.5,40,-1,1,20,0,20000);";
  hisString += "logMaxPi1:tgl:multA:#cleanPion>>hiscPionMax1(100,3.3,4.5,40,-1,1,20,0,20000);";
  hisString += "logMaxPi2:tgl:multA:#cleanPion>>hiscPionMax2(100,3.3,4.5,40,-1,1,20,0,20000);";
  hisString += "logTotPi0:tgl:multA:#cleanPion>>hiscPionTot0(100,3.3,4.5,40,-1,1,20,0,20000);";
  hisString += "logTotPi1:tgl:multA:#cleanPion>>hiscPionTot1(100,3.3,4.5,40,-1,1,20,0,20000);";
  hisString += "logTotPi2:tgl:multA:#cleanPion>>hiscPionTot2(100,3.3,4.5,40,-1,1,20,0,20000);";

  hisArray = AliTreePlayer::MakeHistograms(treeV0, hisString, "abs(qP)>1&&abs(qP)<2&&isSelected", 0, nEvents, -1, 15);
  for (Int_t i = 0; i < hisArray->GetEntries(); i++) gROOT->Add(hisArray->At(i));
  pcstream = new TTreeSRedirector("hisOutputsPion.root", "recreate");
  pcstream->GetFile()->cd();
  for (Int_t iHis = 0; iHis < hisArray->GetEntries(); iHis++) {
      hisArray->At(iHis)->Write();
  }
  delete pcstream;
  gSystem->Exec("hadd -f hisOutput.root hisOutputs*.root");

}
/// TODO - add 30 bins into tgl binning



/// make histogram  maps
void MakeMaps(){
  timer.Start();
  ::Info("MakeMaps","Begin");
  TTreeSRedirector* pcstreamIn= new TTreeSRedirector("hisOutput.root", "");
  TObjArray *hisArray= new TObjArray();
  TObjArray *keyList=(TObjArray*)pcstreamIn->GetFile()->GetListOfKeys();
  for (Int_t iHis=0; iHis<keyList->GetEntries();iHis++){
    THn* his = (THn*)pcstreamIn->GetFile()->Get(keyList->At(iHis)->GetName());
    hisArray->AddLast(his);
  }
  //
  TTreeSRedirector* pcstream= new TTreeSRedirector("mapOutput.root", "recreate");
  for (Int_t i = 0; i < hisArray->GetEntries(); i++) gROOT->Add(hisArray->At(i));
  TMatrixD projectionInfo(5,5);
  projectionInfo(0,0)=0;  projectionInfo(0,1)=1;  projectionInfo(0,2)=0;
  projectionInfo(1,0)=1;  projectionInfo(1,1)=2;  projectionInfo(1,2)=0;
  projectionInfo(2,0)=2;  projectionInfo(2,1)=3;  projectionInfo(2,2)=0;  //
  projectionInfo(3,0)=3;  projectionInfo(3,1)=1;  projectionInfo(3,2)=0;
  projectionInfo(4,0)=4;  projectionInfo(4,1)=0;  projectionInfo(4,2)=0;  //
  for (Int_t iHis=0; iHis<hisArray->GetEntries(); iHis++) {
    Int_t proj[6] = {0, 1, 2, 3, 4, 5};
    THn *hisInput = (THn *) hisArray->At(iHis);
    if (hisInput==0) continue;
//    if (TString(hisInput->GetName()).Contains("hisRatio")){
//      projectionInfo(1,1)=1; projectionInfo(2,1)=2;
//    }else{
//      projectionInfo(1,1)=0; projectionInfo(2,1)=3;  ///
//    }
    if (TString(hisInput->GetName()).Contains("hisdEdx")>0) continue;
    ::Info("MakeMaps", "%s\t%d\t%d\t%f", hisInput->GetName(), hisInput->GetNdimensions(), hisInput->GetNbins(), (hisInput->GetEntries()));
    TStatToolkit::MakeDistortionMapFast(hisInput,pcstream,projectionInfo,0,0.05);
  }
  delete pcstream;
  ::Info("MakeMaps","End");
  timer.Print();
}
void LoadSummaryTrees(){
  treeMapID = LoadSummaryTree("hisc.*");
  treeMap = LoadSummaryTree("(hisR.*|hisC.*|hisF.*)");
}
///
/// \param sregexp  - input tree selection
/// \param addName  - switch to add directory name into alias name
/// \return          - tree with friend tress and aliases
TTree* LoadSummaryTree(TString sregexp, Bool_t addName) {
  //
  TTree * treeLoad =NULL;
  TPRegexp regexp(sregexp);
  TObjArray *fileList = gSystem->GetFromPipe("cat performance.list").Tokenize("\n");
  for (Int_t iFile = 0; iFile < fileList->GetEntries(); iFile++) {
    TString dName = fileList->At(iFile)->GetName();
    dName = gSystem->DirName(dName.Data());
    dName = gSystem->BaseName(dName.Data());
    TFile *fInput = TFile::Open(fileList->At(iFile)->GetName());
    TList *keys = fInput->GetListOfKeys();
    for (Int_t iKey = 0; iKey < keys->GetEntries(); iKey++) {
      printf("%s.%s\n", dName.Data(), keys->At(iKey)->GetName());
      if (regexp.Match(keys->At(iKey)->GetName()) > 0) {
        TTree *treeC = (TTree *) fInput->Get(keys->At(iKey)->GetName());
        if (treeLoad == NULL) {
          TFile *fInput2 = TFile::Open(fileList->At(iFile)->GetName());
          treeLoad = (TTree *) fInput2->Get(keys->At(iKey)->GetName());
        }
        if (addName)
          treeLoad->AddFriend(treeC, TString::Format("%s.%s", dName.Data(), keys->At(iKey)->GetName()).Data());
        if (addName == kFALSE) treeLoad->AddFriend(treeC, TString::Format("%s", keys->At(iKey)->GetName()).Data());
      }
    }
  }

  treeLoad->SetAlias("ratioMaxOK","abs(hisRatioMax03Dist.binMedian-hisRatioMax03Dist.meanG)<0.02&&abs(hisRatioMax13Dist.binMedian-hisRatioMax13Dist.meanG)<0.02&&abs(hisRatioMax23Dist.binMedian-hisRatioMax23Dist.meanG)<0.02");
  treeLoad->SetAlias("ratioTotOK","abs(hisRatioTot03Dist.binMedian-hisRatioTot03Dist.meanG)<0.02&&abs(hisRatioTot13Dist.binMedian-hisRatioTot13Dist.meanG)<0.02&&abs(hisRatioTot23Dist.binMedian-hisRatioTot23Dist.meanG)<0.02");
  treeLoad->SetAlias("ratioTotMaxOK","abs(hisRatioTotMax0Dist.meanG/hisRatioTotMax0Dist.binMedian-1)<0.03&&abs(hisRatioTotMax1Dist.meanG/hisRatioTotMax1Dist.binMedian-1)<0.03&&abs(hisRatioTotMax2Dist.meanG/hisRatioTotMax2Dist.binMedian-1)<0.03");
  treeLoad->SetAlias("isOK","entriesG>5&&ratioTotOK&&ratioMaxOK&&ratioTotMaxOK");
  treeLoad->SetAlias("isOK1","entriesG>5&&ratioTotOK&&ratioMaxOK&&ratioTotMaxOK&&pileUpZBin==3");
  treeLoad->SetAlias("isOKChi2","isOK&&abs(hisChi2TPCDist.vecLTM.fElements[1]-hisChi2TPCDist.meanG)<0.1&&hisChi2TPCDist.rmsG<0.7");
  //
  treeLoad->SetAlias("rmsMax0F","sqrt((hisRatioMax01Dist.rmsGFit**2+hisRatioMax02Dist.rmsGFit**2-hisRatioMax12Dist.rmsGFit**2)/2.)");
  treeLoad->SetAlias("rmsMax1F","sqrt((hisRatioMax01Dist.rmsGFit**2+hisRatioMax12Dist.rmsGFit**2-hisRatioMax02Dist.rmsGFit**2)/2.)");
  treeLoad->SetAlias("rmsMax2F","sqrt((hisRatioMax02Dist.rmsGFit**2+hisRatioMax12Dist.rmsGFit**2-hisRatioMax01Dist.rmsGFit**2)/2.)");
  treeLoad->SetAlias("rmsTot0F","sqrt((hisRatioTot01Dist.rmsGFit**2+hisRatioTot02Dist.rmsGFit**2-hisRatioTot12Dist.rmsGFit**2)/2.)");
  treeLoad->SetAlias("rmsTot1F","sqrt((hisRatioTot01Dist.rmsGFit**2+hisRatioTot12Dist.rmsGFit**2-hisRatioTot02Dist.rmsGFit**2)/2.)");
  treeLoad->SetAlias("rmsTot2F","sqrt((hisRatioTot02Dist.rmsGFit**2+hisRatioTot12Dist.rmsGFit**2-hisRatioTot01Dist.rmsGFit**2)/2.)");
  //
  treeLoad->SetAlias("rmsMax0","sqrt((hisRatioMax01Dist.rmsG**2+hisRatioMax02Dist.rmsG**2-hisRatioMax12Dist.rmsG**2)/2.)");
  treeLoad->SetAlias("rmsMax1","sqrt((hisRatioMax01Dist.rmsG**2+hisRatioMax12Dist.rmsG**2-hisRatioMax02Dist.rmsG**2)/2.)");
  treeLoad->SetAlias("rmsMax2","sqrt((hisRatioMax02Dist.rmsG**2+hisRatioMax12Dist.rmsG**2-hisRatioMax01Dist.rmsG**2)/2.)");
  treeLoad->SetAlias("rmsTot0","sqrt((hisRatioTot01Dist.rmsG**2+hisRatioTot02Dist.rmsG**2-hisRatioTot12Dist.rmsG**2)/2.)");
  treeLoad->SetAlias("rmsTot1","sqrt((hisRatioTot01Dist.rmsG**2+hisRatioTot12Dist.rmsG**2-hisRatioTot02Dist.rmsG**2)/2.)");
  treeLoad->SetAlias("rmsTot2","sqrt((hisRatioTot02Dist.rmsG**2+hisRatioTot12Dist.rmsG**2-hisRatioTot01Dist.rmsG**2)/2.)");
//  treeLoad->SetAlias("rmsMaxTot0","hisRatioTotMax0Dist.rmsG/hisRatioTotMax0Dist.meanG");
//  treeLoad->SetAlias("rmsMaxTot1","hisRatioTotMax1Dist.rmsG/hisRatioTotMax1Dist.meanG");
//  treeLoad->SetAlias("rmsMaxTot2","hisRatioTotMax2Dist.rmsG/hisRatioTotMax2Dist.meanG");
//  treeLoad->SetAlias("rmsMaxTot0LF","hisRatioTotMax0Dist.rmsGFit/hisRatioTotMax0Dist.meanGFit");
//  treeLoad->SetAlias("rmsMaxTot1LF","hisRatioTotMax1Dist.rmsGFit/hisRatioTotMax1Dist.meanGFit");
//  treeLoad->SetAlias("rmsMaxTot2LF","hisRatioTotMax2Dist.rmsGFit/hisRatioTotMax2Dist.meanGFit");
  //
  // ncl normalized resolution and other aliases
  for (Int_t iPad=0; iPad<3; iPad++) {
    treeLoad->SetAlias(TString::Format("rmsMaxTot%d",iPad).Data(), TString::Format("hisRatioTotMax0Dist.rmsG/hisRatioTotMax0Dist.meanG",iPad,iPad).Data());
    treeLoad->SetAlias(TString::Format("rmsMaxTot%dLF",iPad).Data(), TString::Format("hisRatioTotMax0Dist.rmsGFit/hisRatioTotMax0Dist.meanGFit",iPad,iPad).Data());
    //
    treeLoad->SetAlias(TString::Format("rmsMax%dNorm",iPad).Data(), TString::Format("rmsMax%d*sqrt(hisCross%dDist.mean0)",iPad,iPad).Data());
    treeLoad->SetAlias(TString::Format("rmsTot%dNorm",iPad).Data(), TString::Format("rmsTot%d*sqrt(hisCross%dDist.mean0)",iPad,iPad).Data());
    treeLoad->SetAlias(TString::Format("rmsMax%dNormLF",iPad).Data(), TString::Format("rmsMax%dF*sqrt(hisCross%dDist.mean0Fit)",iPad,iPad).Data());
    treeLoad->SetAlias(TString::Format("rmsTot%dNormLF",iPad).Data(), TString::Format("rmsTot%dF*sqrt(hisCross%dDist.mean0Fit)",iPad,iPad).Data());
    treeLoad->SetAlias(TString::Format("missingFrac%d",iPad).Data(), TString::Format("1-hisFraction%dDist.mean0",iPad).Data());
    treeLoad->SetAlias(TString::Format("missingFrac%dLF",iPad).Data(), TString::Format("1-hisFraction%dDist.mean0Fit",iPad).Data());
    treeLoad->SetAlias(TString::Format("missingFracRMS%d",iPad).Data(), TString::Format("hisFraction%dDist.rms0",iPad).Data());
    treeLoad->SetAlias(TString::Format("missingFracRMS%dLF",iPad).Data(), TString::Format("hisFraction%dDist.rms0Fit",iPad).Data());
    treeLoad->SetAlias(TString::Format("missingFracRMS%dExp",iPad).Data(), TString::Format("sqrt(missingFrac%dLF/hisCross%dDist.mean0)",iPad,iPad).Data());
  }
  treeLoad->SetAlias("mdEdx","(50)/sdEdxMaxCenter**2");
  treeLoad->SetAlias("multACenter1000","multACenter/1000.");
  TStatToolkit::AddMetadata(treeLoad,"mdEdx.AxisTitle","1/dEdx (1/MIP) ");
  TStatToolkit::AddMetadata(treeLoad,"atglCenter.AxisTitle","tan(#Theta)");
  TStatToolkit::AddMetadata(treeLoad,"sdEdxMaxCenter.AxisTitle","#sqrt{dEdx_{Max}}");
  TStatToolkit::AddMetadata(treeLoad,"multACenter.AxisTitle","occupancy estimator (N_{ESDtracks})");
  TStatToolkit::AddMetadata(treeLoad,"multACenter1000.AxisTitle","occupancy estimator (N_{ESDtracks}/1000)");

  treeLoad->SetMarkerStyle(21);
  treeLoad->SetMarkerSize(0.8);
  return treeLoad;
}

///
void LoadNDLocal(){
  TFile *fND= TFile::Open("ndFitAll.root");
  if (fND==NULL) {::Error("LoadNDLocal","Invalid or not existing file ndFitAll.root");return;}
  TList* listND = fND->GetListOfKeys();
  if (listND==NULL) {::Error("LoadNDLocal","Invalid or empty file ndFitAll.root");return;}
  for (Int_t iND=0; iND<listND->GetEntries();iND++){
    AliNDLocalRegression *regression=(AliNDLocalRegression *)fND->Get(listND->At(iND)->GetName());
    Int_t hashIndex=regression->GetVisualCorrectionIndex();
    AliNDLocalRegression::AddVisualCorrection(regression, hashIndex);
    if (treeMap) {
      treeMap->SetAlias(TString::Format("%sFit", regression->GetName()).Data(), TString::Format(
              "AliNDLocalRegression::GetCorrND(%d,atglMean,multAMean,sdEdxMaxMean,pileUpZCenter+0)", hashIndex).Data());
      treeMap->SetAlias(TString::Format("%sFit80", regression->GetName()).Data(), TString::Format(
              "AliNDLocalRegression::GetCorrND(%d,atglMean,multAMean,sqrt(80.),pileUpZCenter+0)", hashIndex).Data());
    }
    if (tree){
      tree->SetAlias(TString::Format("%sFit",regression->GetName()).Data(),TString::Format("AliNDLocalRegression::GetCorrND(%d,tgl,multA,dEdxMax,pileUpZ+0)",hashIndex).Data());
    }
  }
  for (Int_t iND=0; iND<12; iND++){
    regMeanG[iND]=AliNDLocalRegression::GetVisualCorrection(TString::Format("%s.meanG",regNames[iND]).Data());
    regRMS[iND]=AliNDLocalRegression::GetVisualCorrection(TString::Format("%s.rmsG",regNames[iND]).Data());
  }

}



/// \brief  Make dEdx ration local ND fits in 4D
/// \param varName
/// \param customCut
/// \param errorVar
void makeNDLocalFit(TString varName, TString customCut, TString errorVar){
  ::Info("makeNDLocalFit","makeNDLocalFit(\"%s\",\"%s\",\"%s\")",varName.Data(),customCut.Data(),errorVar.Data());
  if (treeMap==0) {
    treeMap = LoadSummaryTree();
  }
  // 1. test variable and cut string
  Int_t entriesTestVar = treeMap->Draw(varName.Data(),"1","goff",1000);
  Int_t entriesTestError = treeMap->Draw(errorVar.Data(),"1","goff",1000);
  Int_t entriesTestCut = treeMap->Draw(customCut.Data(),"1","goff",1000);
  if (entriesTestVar<=0){ ::Error("makeNDLocalFit","Invalid variable: %s",varName.Data()); return; }
  if (entriesTestCut<=0){ ::Error("makeNDLocalFit","Invalid selection: %s",customCut.Data()); return; }
  if (entriesTestError<=0){ ::Error("makeNDLocalFit","Invalid error expression: %s",errorVar.Data()); return; }
  treeMap->SetAlias("value",varName.Data());
  treeMap->SetAlias("error",errorVar.Data());
  treeMap->SetAlias("customCut",customCut.Data());
  //
  // 2.) Make NDlocal regression
  //    2.1) get histogram layout
  TStopwatch timer;
  Int_t     ndim=4;
  Int_t     nbins[4]= {10,  10,   20     , 5};  // {tgl,mult, dEdx, pileUpZCenter}
  Double_t  xmin[4] = {0,  0,     5    , -40 };
  Double_t  xmax[4] = {1,  18000, 32,  40};
  THnF* hN= new THnF("exampleFit","exampleFit", ndim, nbins, xmin,xmax);
  AliNDLocalRegression *regression= new AliNDLocalRegression;
  regression= new  AliNDLocalRegression();
  regression->SetName(varName.Data());
  Int_t hashIndex=regression->GetVisualCorrectionIndex();
  TCut cut="1";  ///TODO  in the fit currently we assumed the entrylist was applied - otherwise cut has to be specified
  regression->SetHistogram((THn*)(hN->Clone()));
  regression->MakeFit(treeMap,"value:error", "atglMean:multAMean:sdEdxMaxMean:pileUpZCenter","customCut","0.10:4000:5:50","2:2:2:2",0.01);
  AliNDLocalRegression::AddVisualCorrection(regression, hashIndex);
  treeMap->SetAlias("fit",TString::Format("AliNDLocalRegression::GetCorrND(%d,atglMean,multAMean,sdEdxMaxMean,pileUpZCenter+0)",hashIndex).Data());
  treeMap->SetAlias(TString::Format("%sFit",varName.Data()).Data(),TString::Format("AliNDLocalRegression::GetCorrND(%d,atglMean,multAMean,sdEdxMaxMean,pileUpZCenter+0)",hashIndex).Data());
  if (tree){
    tree->SetAlias(TString::Format("%sFit",varName.Data()).Data(),TString::Format("AliNDLocalRegression::GetCorrND(%d,tgl,multA,dEdxMax,pileUpZ+0)",hashIndex).Data());
  }
  TFile *fout= TFile::Open(TString::Format("ndFit%s.root",varName.Data()),"recreate");
  regression->Write(varName.Data());
  delete fout;
  timer.Print();
  treeMap->GetListOfAliases()->ls();
}
///\brief fit all ratios and RMS
void MakeNDLocalRegressions(){
  makeNDLocalFit("hisRatioTotMax0Dist.binMedian","isOK","sqrt(0.01**2+hisRatioTotMax0Dist.rmsG**2/hisRatioTotMax0Dist.entries)");
  makeNDLocalFit("hisRatioTotMax1Dist.binMedian","isOK","sqrt(0.01**2+hisRatioTotMax1Dist.rmsG**2/hisRatioTotMax1Dist.entries)");
  makeNDLocalFit("hisRatioTotMax2Dist.binMedian","isOK","sqrt(0.01**2+hisRatioTotMax2Dist.rmsG**2/hisRatioTotMax2Dist.entries)");
  makeNDLocalFit("hisRatioTotMax0Dist.meanG","isOK","sqrt(0.01**2+hisRatioTotMax0Dist.rmsG**2/hisRatioTotMax0Dist.entries)");
  makeNDLocalFit("hisRatioTotMax1Dist.meanG","isOK","sqrt(0.01**2+hisRatioTotMax1Dist.rmsG**2/hisRatioTotMax1Dist.entries)");
  makeNDLocalFit("hisRatioTotMax2Dist.meanG","isOK","sqrt(0.01**2+hisRatioTotMax2Dist.rmsG**2/hisRatioTotMax2Dist.entries)");
  makeNDLocalFit("hisRatioTotMax0Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTotMax0Dist.rmsG**2/hisRatioTotMax0Dist.entries)");
  makeNDLocalFit("hisRatioTotMax1Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTotMax1Dist.rmsG**2/hisRatioTotMax1Dist.entries)");
  makeNDLocalFit("hisRatioTotMax2Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTotMax2Dist.rmsG**2/hisRatioTotMax2Dist.entries)");
  //
  makeNDLocalFit("hisRatioTot03Dist.meanG","isOK","sqrt(0.01**2+hisRatioTot03Dist.rmsG**2/hisRatioTot03Dist.entries)");
  makeNDLocalFit("hisRatioTot13Dist.meanG","isOK","sqrt(0.01**2+hisRatioTot13Dist.rmsG**2/hisRatioTot03Dist.entries)");
  makeNDLocalFit("hisRatioTot23Dist.meanG","isOK","sqrt(0.01**2+hisRatioTot13Dist.rmsG**2/hisRatioTot03Dist.entries)");
  makeNDLocalFit("hisRatioMax03Dist.meanG","isOK","sqrt(0.01**2+hisRatioMax03Dist.rmsG**2/hisRatioMax03Dist.entries)");
  makeNDLocalFit("hisRatioMax13Dist.meanG","isOK","sqrt(0.01**2+hisRatioMax13Dist.rmsG**2/hisRatioMax03Dist.entries)");
  makeNDLocalFit("hisRatioMax23Dist.meanG","isOK","sqrt(0.01**2+hisRatioMax13Dist.rmsG**2/hisRatioMax03Dist.entries)");
  makeNDLocalFit("hisRatioTot03Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTot03Dist.rmsG**2/hisRatioTot03Dist.entries)");
  makeNDLocalFit("hisRatioTot13Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTot13Dist.rmsG**2/hisRatioTot03Dist.entries)");
  makeNDLocalFit("hisRatioTot23Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTot13Dist.rmsG**2/hisRatioTot03Dist.entries)");
  makeNDLocalFit("hisRatioMax03Dist.rmsG","isOK","sqrt(0.01**2+hisRatioMax03Dist.rmsG**2/hisRatioMax03Dist.entries)");
  makeNDLocalFit("hisRatioMax13Dist.rmsG","isOK","sqrt(0.01**2+hisRatioMax13Dist.rmsG**2/hisRatioMax03Dist.entries)");
  makeNDLocalFit("hisRatioMax23Dist.rmsG","isOK","sqrt(0.01**2+hisRatioMax13Dist.rmsG**2/hisRatioMax03Dist.entries)");
  //
  makeNDLocalFit("hisRatioTot01Dist.meanG","isOK","sqrt(0.01**2+hisRatioTot01Dist.rmsG**2/hisRatioTot01Dist.entries)");
  makeNDLocalFit("hisRatioTot12Dist.meanG","isOK","sqrt(0.01**2+hisRatioTot12Dist.rmsG**2/hisRatioTot12Dist.entries)");
  makeNDLocalFit("hisRatioTot02Dist.meanG","isOK","sqrt(0.01**2+hisRatioTot02Dist.rmsG**2/hisRatioTot02Dist.entries)");
  makeNDLocalFit("hisRatioMax01Dist.meanG","isOK","sqrt(0.01**2+hisRatioMax01Dist.rmsG**2/hisRatioMax01Dist.entries)");
  makeNDLocalFit("hisRatioMax12Dist.meanG","isOK","sqrt(0.01**2+hisRatioMax12Dist.rmsG**2/hisRatioMax12Dist.entries)");
  makeNDLocalFit("hisRatioMax02Dist.meanG","isOK","sqrt(0.01**2+hisRatioMax02Dist.rmsG**2/hisRatioMax02Dist.entries)");
  //
  makeNDLocalFit("hisRatioTot01Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTot01Dist.rmsG**2/hisRatioTot01Dist.entries)");
  makeNDLocalFit("hisRatioTot12Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTot12Dist.rmsG**2/hisRatioTot12Dist.entries)");
  makeNDLocalFit("hisRatioTot02Dist.rmsG","isOK","sqrt(0.01**2+hisRatioTot02Dist.rmsG**2/hisRatioTot02Dist.entries)");
  makeNDLocalFit("hisRatioMax01Dist.rmsG","isOK","sqrt(0.01**2+hisRatioMax01Dist.rmsG**2/hisRatioMax01Dist.entries)");
  makeNDLocalFit("hisRatioMax12Dist.rmsG","isOK","sqrt(0.01**2+hisRatioMax12Dist.rmsG**2/hisRatioMax12Dist.entries)");
  makeNDLocalFit("hisRatioMax02Dist.rmsG","isOK","sqrt(0.01**2+hisRatioMax02Dist.rmsG**2/hisRatioMax02Dist.entries)");
  //
  makeNDLocalFit("hisCross0Dist.mean0","isOK","sqrt(0.5**2+hisCross0Dist.rms0**2/hisCross0Dist.entries)");
  makeNDLocalFit("hisCross1Dist.mean0","isOK","sqrt(0.5**2+hisCross1Dist.rms0**2/hisCross1Dist.entries)");
  makeNDLocalFit("hisCross2Dist.mean0","isOK","sqrt(0.5**2+hisCross2Dist.rms0**2/hisCross2Dist.entries)");
  makeNDLocalFit("hisFraction0Dist.mean0","isOK","sqrt(0.01**2+hisFraction0Dist.rms0**2/hisFraction0Dist.entries)");
  makeNDLocalFit("hisFraction1Dist.mean0","isOK","sqrt(0.01**2+hisFraction1Dist.rms0**2/hisFraction1Dist.entries)");
  makeNDLocalFit("hisFraction2Dist.mean0","isOK","sqrt(0.01**2+hisFraction2Dist.rms0**2/hisFraction2Dist.entries)");

  makeNDLocalFit("hisChi2TPCDist.meanG","isOKChi2","sqrt(0.01**2+hisChi2TPCDist.rms0**2/hisChi2TPCDist.entries)");
  makeNDLocalFit("hisChi2TPCDist.binMedian","isOKChi2","sqrt(0.01**2+hisChi2TPCDist.rms0**2/hisChi2TPCDist.entries)");
  makeNDLocalFit("hisChi2TPCDist.rmsG","isOKChi2","sqrt(0.01**2+hisChi2TPCDist.rms0**2/hisChi2TPCDist.entries)");
  //
  treeMap->SetAlias("rmsMax0","sqrt((hisRatioMax01Dist.rmsGFit**2+hisRatioMax02Dist.rmsGFit**2-hisRatioMax12Dist.rmsGFit**2)/2.)");
  treeMap->SetAlias("rmsMax1","sqrt((hisRatioMax01Dist.rmsGFit**2+hisRatioMax12Dist.rmsGFit**2-hisRatioMax02Dist.rmsGFit**2)/2.)");
  treeMap->SetAlias("rmsMax2","sqrt((hisRatioMax02Dist.rmsGFit**2+hisRatioMax12Dist.rmsGFit**2-hisRatioMax01Dist.rmsGFit**2)/2.)");
  treeMap->SetAlias("rmsTot0","sqrt((hisRatioTot01Dist.rmsGFit**2+hisRatioTot02Dist.rmsGFit**2-hisRatioTot12Dist.rmsGFit**2)/2.)");
  treeMap->SetAlias("rmsTot1","sqrt((hisRatioTot01Dist.rmsGFit**2+hisRatioTot12Dist.rmsGFit**2-hisRatioTot02Dist.rmsGFit**2)/2.)");
  treeMap->SetAlias("rmsTot2","sqrt((hisRatioTot02Dist.rmsGFit**2+hisRatioTot12Dist.rmsGFit**2-hisRatioTot01Dist.rmsGFit**2)/2.)");
  treeMap->SetAlias("rmsMaxTot0","hisRatioTotMax0Dist.rmsG/hisRatioTotMax0Dist.meanG");
  treeMap->SetAlias("rmsMaxTot1","hisRatioTotMax1Dist.rmsG/hisRatioTotMax1Dist.meanG");
  treeMap->SetAlias("rmsMaxTot2","hisRatioTotMax2Dist.rmsG/hisRatioTotMax2Dist.meanG");
  treeMap->SetAlias("rmsMaxTot0LF","hisRatioTotMax0Dist.rmsGFit/hisRatioTotMax0Dist.meanGFit");
  treeMap->SetAlias("rmsMaxTot1LF","hisRatioTotMax1Dist.rmsGFit/hisRatioTotMax1Dist.meanGFit");
  treeMap->SetAlias("rmsMaxTot2LF","hisRatioTotMax2Dist.rmsGFit/hisRatioTotMax2Dist.meanGFit");

}


void MakeFitResolLowMult() {
  likeGausCachy->SetParameters(0.95, 1);
  treeMap->SetAlias("fitCut0", "entries>200&&isOK&&pileUpZBin==3&&multACenter<1000");
  treeMap->SetAlias("fitCut1", "entries>200&&isOK&&pileUpZBin==3");
  //
  // Create ncl normalized resolution fitter
  AliTMinuitToolkit *fitterResol0 = new AliTMinuitToolkit("AliTMinuitToolkitFitterResol0.root");
  TFormula *formulaResol0 = new TFormula("formulaResol", "[0]*(x[0]^abs([1]))*(x[1]^abs([2]))");
  fitterResol0->SetFitFunction((TF1 *) formulaResol0, kTRUE);
  fitterResol0->SetVerbose(0x1);
  fitterResol0->SetLogLikelihoodFunction(likeGausCachy);
  fitterResol0->SetInitialParam(new TMatrixD(3, 4));
  TMatrixD &initParam = (*fitterResol0->GetInitialParam());
  //
  TObjArray fitRMSNorm0(9);
  TObjArray fitRMSNorm0Error(9);
  treeMap->SetAlias("outlierCut", "(1+0)");
  for (Int_t iType = 0; iType < 3; iType++) {
    for (Int_t iPad = 0; iPad < 3; iPad++) {
      initParam(0, 0) = 1.6;initParam(0, 1) = 2;initParam(0, 2) = 0.;initParam(0, 3) = 10;
      initParam(1, 0) = 0.25;initParam(1, 1) = 2;initParam(1, 2) = 0.;initParam(1, 3) = 10;
      initParam(2, 0) = 0.5;initParam(2, 1) = 2;initParam(2, 2) = 0.;initParam(2, 3) = 10;
      TString fitValue = "";
      if (iType == 0) {
        fitValue = TString::Format("rmsMax%dNorm", iPad);
        TStatToolkit::AddMetadata(treeMap, Form("rmsMax%dNormFit0.AxisTitle", iPad), Form("RMS_{dEdx_{Max0%d}}  #sqrt{N_{CR}} - fit 0", iPad));
        TStatToolkit::AddMetadata(treeMap, Form("rmsMax%dNormLF.AxisTitle", iPad), Form("RMS_{dEdx_{Max0%d}}  #sqrt{N_{CR}} - local ND fit", iPad));
      }
      if (iType == 1) {
        fitValue = TString::Format("rmsTot%dNorm", iPad);
        TStatToolkit::AddMetadata(treeMap, Form("rmsMax%dTotFit0.AxisTitle", iPad), Form("RMS_{dEdx_{Tot%d}}  #sqrt{N_{CR}}", iPad));
        TStatToolkit::AddMetadata(treeMap, Form("rmsMax%dNormLF.AxisTitle", iPad),
                                  Form("RMS_{dEdx_{Max0%d}}  #sqrt{N_{CR}} - local ND fit", iPad));
      }
      if (iType == 2) {
        fitValue = TString::Format("rmsMaxTot%d", iPad);
        //TStatToolkit::AddMetadata(treeMap,Form("rmsMax%dTotFit0.AxisTitle",iPad),Form("RMS_{dEdx_{Tot%d}}  #sqrt{N_{CR}}",iPad));
        //TStatToolkit::AddMetadata(treeMap,Form("rmsMax%dNormLF.AxisTitle",iPad),Form("RMS_{dEdx_{Max0%d}}  #sqrt{N_{CR}} - local ND fit",iPad));
      }

      if (iPad != 2)
        fitterResol0->FillFitter(treeMap, (fitValue + ":10").Data(), "(1./sdEdxMaxMean**2):(1./sqrt(1+atglMean**2))",
                                 "fitCut0&&outlierCut", 0, 10000000);
      if (iPad == 2 && iType < 2)
        fitterResol0->FillFitter(treeMap, (fitValue + ":10").Data(), "(1./sdEdxMaxMean**2):(1./sqrt(1+atglMean**2))",
                                 "fitCut0&&atglCenter<0.7&&outlierCut", 0, 10000000);
      fitterResol0->Bootstrap(20, "report0");
      treeMap->SetAlias((fitValue + "Fit0").Data(), fitterResol0->GetFitFunctionAsAlias().Data());
      const TVectorD &fit = *(fitterResol0->GetParameters());
      const TVectorD &fitErr = *(fitterResol0->GetRMSEstimator());
      fitRMSNorm0.AddLast(fit.Clone());
      fitRMSNorm0Error.AddLast(fitErr.Clone());
      ::Info("MakeFitResol", "%s\t%2.2f\t%3.3f\t%3.3f", fitValue.Data(), fit[0], fit[1], fit[2]);

    }
  }
}
//

void MakeFitResolHighMult() {
  //
  AliTMinuitToolkit *fitterResolMaxTot = new AliTMinuitToolkit("AliTMinuitToolkitFitterResolMaxTot.root");
  TFormula *formulaResolMaxTot = new TFormula("formulaResolMaxTot", "sqrt( x[0]**2 + ([0])*(x[1]^([1]))*(x[2]^([2])) )");
  fitterResolMaxTot->SetFitFunction((TF1 *) formulaResolMaxTot, kTRUE);
  fitterResolMaxTot->SetVerbose(0x1);
  fitterResolMaxTot->SetLogLikelihoodFunction(likeGausCachy);
  fitterResolMaxTot->SetInitialParam(new TMatrixD(3, 4));
   AliTMinuitToolkit::SetPredefinedFitter("fitterResolMaxTotHighMult",fitterResolMaxTot);
  //
  for (Int_t iType = 2; iType < 3; iType++) {
    for (Int_t iPad = 0; iPad < 3; iPad++) {
      TString fitValue = "";
      TString fitData = "";
      if (iType == 0) {
        fitValue = TString::Format("rmsMax%dNorm", iPad);
        fitData = TString::Format("rmsMax%dNormFit0:(multAMean/10000.):(1/((sdEdxMaxCenter**2)*sqrt(1+atglCenter**2)))", iPad);
        TStatToolkit::AddMetadata(treeMap, Form("rmsMax%dNormFitM.AxisTitle", iPad), Form("RMS_{dEdx_{Max0%d}}  #sqrt{N_{CR}} - fit 0", iPad));
      }
      if (iType == 1) {
        fitValue = TString::Format("rmsTot%dNorm", iPad);
        fitData = TString::Format("rmsTot%dNormFit0:(multAMean/10000.):(1/((sdEdxMaxCenter**2)*sqrt(1+atglCenter**2)))", iPad);
        TStatToolkit::AddMetadata(treeMap, Form("rmsMax%dTotFitM.AxisTitle", iPad), Form("RMS_{dEdx_{Tot%d}}  #sqrt{N_{CR}}", iPad));
      }
      if (iType == 2) {
        fitValue = TString::Format("rmsMaxTot%d", iPad);
        fitData = TString::Format("rmsMaxTot%dFit0:(multAMean/10000.):(1/((sdEdxMaxCenter**2)*sqrt(1+atglCenter**2)))", iPad);
        //TStatToolkit::AddMetadata(treeMap,Form("rmsMax%dTotFit0.AxisTitle",iPad),Form("RMS_{dEdx_{Tot%d}}  #sqrt{N_{CR}}",iPad));
        //TStatToolkit::AddMetadata(treeMap,Form("rmsMax%dNormLF.AxisTitle",iPad),Form("RMS_{dEdx_{Max0%d}}  #sqrt{N_{CR}} - local ND fit",iPad));
      }
      TMatrixD &initParamMaxTot = (*fitterResolMaxTot->GetInitialParam());
      initParamMaxTot(0, 0) = 0.1;
      initParamMaxTot(0, 1) = 2;
      initParamMaxTot(0, 2) = 0.0;
      initParamMaxTot(0, 3) = 0.0;
      initParamMaxTot(1, 0) = 1;
      initParamMaxTot(1, 1) = 2;
      initParamMaxTot(1, 2) = 0.0;
      initParamMaxTot(1, 3) = 0.0;
      initParamMaxTot(2, 0) = 1;
      initParamMaxTot(2, 1) = 2;
      initParamMaxTot(2, 2) = 0.0;
      initParamMaxTot(2, 3) = 0.0;
      fitterResolMaxTot->FillFitter(treeMap, (fitValue+":1").Data(), fitData.Data(), "fitCut1&&outlierCut", 0, 10000000);
      fitterResolMaxTot->Bootstrap(20, "report1");
      treeMap->SetAlias("fit", fitterResolMaxTot->GetFitFunctionAsAlias().Data());
      treeMap->SetAlias((fitValue + "FitM").Data(), fitterResolMaxTot->GetFitFunctionAsAlias().Data());
    }
  }

  fitterResolMaxTot->FillFitter(treeMap, "rmsMaxTot0:1/(0.005+0)", "rmsMaxTot0Fit0:(multAMean/10000.):(1/((sdEdxMaxCenter**2)*sqrt(1+atglCenter**2)))", "fitCut1&&outlierCut", 0,
                                10000000);
  fitterResolMaxTot->Bootstrap(20, "report1");
  treeMap->SetAlias("fit", fitterResolMaxTot->GetFitFunctionAsAlias().Data());


}

///
void MakeFitNclFraction() {
  //
  treeMap->SetAlias("fitCutNcl","isOK&&pileUpZBin==3&&entries>100");
  AliTMinuitToolkit *fitterNclFraction = new AliTMinuitToolkit("AliTMinuitToolkitFitterNclFraction.root");
  TFormula *formulaNclFraction = new TFormula("formulaNclFraction", "[0]+[1]*((x[0]^[2])*([4]+x[1]^[3]))");
  fitterNclFraction->SetFitFunction((TF1 *) formulaNclFraction, kTRUE);
  fitterNclFraction->SetVerbose(0x1);
  fitterNclFraction->SetLogLikelihoodFunction(likeGausCachy);
  fitterNclFraction->SetInitialParam(new TMatrixD(5, 4));
  TMatrixD &initParamNcl = (*fitterNclFraction->GetInitialParam());
  TObjArray fitParam(3), fitParamError(3);
  AliTMinuitToolkit::SetPredefinedFitter("fitterNclFraction",fitterNclFraction);
  for (Int_t iPad = 0; iPad < 4; iPad++) {
    TString fitValue = "";
    initParamNcl(0, 0) = 0.01;
    initParamNcl(0, 1) = 3;
    initParamNcl(1, 0) = 0.1;
    initParamNcl(1, 1) = 3;
    initParamNcl(2, 0) = 1;
    initParamNcl(2, 1) = 5;
    initParamNcl(3, 0) = 1;
    initParamNcl(3, 1) = 5;
    initParamNcl(4, 0) = 0.01;
    initParamNcl(4, 1) = 5;
    fitValue = TString::Format("missingFrac%d", iPad);
    fitterNclFraction->FillFitter(treeMap, (fitValue + ":1/0.01").Data(), "(50./((sdEdxMaxCenter**2)*sqrt(1+atglCenter**2))):(multACenter/10000)", "fitCutNcl&&Entry$%3==1", 0, 10000000);
    fitterNclFraction->Bootstrap(20, "report1");
    treeMap->SetAlias((fitValue + "FitM").Data(), fitterNclFraction->GetFitFunctionAsAlias().Data());
    TStatToolkit::AddMetadata(treeMap, Form("missingFrac%d.AxisTitle", iPad), Form("cluster fraction : p(cli| (i-#delta)+..+(i+#delta))", iPad));
    TStatToolkit::AddMetadata(treeMap, Form("missingFrac%dLF.AxisTitle", iPad), Form("cluster fraction local fit: p(cli| (i-#delta)+..+(i+#delta))", iPad));
    TStatToolkit::AddMetadata(treeMap, Form("missingFrac%dFitM.AxisTitle", iPad), Form("F_{clGlobalFit} =  p(cli| (i-#delta)+..+(i+#delta))", iPad));
    ::Info("MakeFitNclFraction", "Pad%d", iPad);
    fitterNclFraction->GetParameters()->Print();
    fitterNclFraction->GetRMSEstimator()->Print();
    fitParam.AddAt(fitterNclFraction->GetParameters()->Clone(),iPad);
    fitParamError.AddAt(fitterNclFraction->GetRMSEstimator()->Clone(),iPad);
  }
  for (Int_t iPad=0; iPad<3; iPad++) {
    ::Info("MakeFitNclFraction","Pad%d",iPad);
    fitParam.At(iPad)->Print();
    fitParamError.At(iPad)->Print();
  }
}

///
void MakeFitTPCChi2() {
  //
  treeMap->SetAlias("fitCutChi2","isOKChi2&&pileUpZBin==3&&entries>100");
  AliTMinuitToolkit *fitterTPCChi2 = new AliTMinuitToolkit("AliTMinuitToolkitFitterTPCChi2.root");
  TFormula *formulaTPCChi2 = new TFormula("formulaTPCChi2", "[0]+[1]*((x[0]^[2])*([4]+x[1]^[3])*(x[2]^[5]))");
  fitterTPCChi2->SetFitFunction((TF1 *) formulaTPCChi2, kTRUE);
  fitterTPCChi2->SetVerbose(0x1);
  fitterTPCChi2->SetLogLikelihoodFunction(likeGausCachy);
  fitterTPCChi2->SetInitialParam(new TMatrixD(6, 4));
  TMatrixD &initParamNcl = (*fitterTPCChi2->GetInitialParam());
  TObjArray fitParam(4), fitParamError(4);
  AliTMinuitToolkit::SetPredefinedFitter("fitterTPCChi2",fitterTPCChi2);
  for (Int_t iType = 0; iType < 4; iType++) {
    TString fitValue = "";
    initParamNcl(0, 0) = 0.3;
    initParamNcl(0, 1) = 3;
    initParamNcl(1, 0) = 0.1;
    initParamNcl(1, 1) = 3;
    initParamNcl(2, 0) = 1;
    initParamNcl(2, 1) = 5;
    initParamNcl(3, 0) = 1;
    initParamNcl(3, 1) = 5;
    initParamNcl(4, 0) = 0.01;
    initParamNcl(4, 1) = 5;
    initParamNcl(5, 0) = 0.01;
    initParamNcl(5, 1) = 5;
    if (iType==0) fitValue = "hisChi2TPCDist.binMedian";
    if (iType==1) fitValue = "hisChi2TPCDist.mean";
    if (iType==2) fitValue = "hisChi2TPCDist.meanG";
    if (iType==3) fitValue = "hisChi2TPCDist.rmsG";
    fitterTPCChi2->FillFitter(treeMap, (fitValue + ":1/0.01").Data(), "(50./((sdEdxMaxCenter**2))):(multACenter/10000):sqrt(1+atglCenter**2)", "fitCutChi2&&Entry$%3==1", 0, 10000000);
    fitterTPCChi2->Bootstrap(20, "report1");
    treeMap->SetAlias((fitValue + "FitM").Data(), fitterTPCChi2->GetFitFunctionAsAlias().Data());
//    TStatToolkit::AddMetadata(treeMap, Form("missingFrac%d.AxisTitle", iPad), Form("cluster fraction : p(cli| (i-#delta)+..+(i+#delta))", iPad));
    fitterTPCChi2->GetParameters()->Print();
    fitterTPCChi2->GetRMSEstimator()->Print();
    fitParam.AddAt(fitterTPCChi2->GetParameters()->Clone(),iType);
    fitParamError.AddAt(fitterTPCChi2->GetRMSEstimator()->Clone(),iType);
  }
  for (Int_t iPad=0; iPad<4; iPad++) {
    ::Info("MakeFitTPCChi2","Pad%d",iPad);
    fitParam.At(iPad)->Print();
    fitParamError.At(iPad)->Print();
  }

}

