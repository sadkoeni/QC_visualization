/*
  //.L $NOTES/JIRA/ATO-432/code/drawHighdEdx.C+
  // LoadChain("filtered.list",10);
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);

  .L $NOTES/JIRA/ATO-432/code/tofStudies.C+
  InitTree()
  makeActiveMap();

 */
#include "TChain.h"
#include "TFile.h"
#include "AliESDtrack.h"
#include "AliXRDPROOFtoolkit.h"
#include "TVectorD.h"
#include "AliTRDgeometry.h"
#include "TH1F.h"
#include "TMath.h"
#include "TPad.h"

TChain * treeV0=NULL;
TChain * tree=NULL;
AliESDtrack *track=0;
const Int_t kNDetectors=640; // number of TRD+TOF detectors
const Int_t occuCut=8;       // occupancy cut
const Double_t kMarginR=5;
const Double_t kMarginZ=3;
TVectorD activeMap(kNDetectors);
TH1F* histoDetector=NULL;

TVectorD vecSec(7), vecZ(7), vecDet(7), vecdSec(7), vecdEdge(7), vecActive(7), vecStatus(7),vecDeadZ(7),vecDeadR(7), vecDeadDet(7);
Double_t nFindable=0;
Double_t nFound=0;
AliTRDgeometry geom;
TMatrixD zBoundary0(7,6);
TMatrixD zBoundary1(7,6);

/// \brief cache TRD geometry
void CacheTRDGeom(){
  for (Int_t iLayer=0; iLayer<6; iLayer++){
    for (Int_t iStack=0;iStack<5; iStack++){
      zBoundary0(iLayer,iStack)=geom.GetRow0(iLayer,iStack,0);
      zBoundary1(iLayer,iStack)=-geom.GetRow0(iLayer,4-iStack,0);
    }
  }
  for (Int_t iLayer=0; iLayer<6; iLayer++) zBoundary1(iLayer,5)=-zBoundary0(iLayer,0);
}

/// \brief build active detector map
void makeActiveMap(Int_t nPoints=10000, Double_t threshold=0.3){
  // Int_t nPoints=2000;
  histoDetector=new TH1F("histoDetector","histoDetector",kNDetectors,0,kNDetectors);
  for (Int_t iLayer=0; iLayer<6; iLayer++){
    TString what=TString::Format("GetDet(%d)>>+histoDetector",iLayer);
    TString where=TString::Format("isOK0&&track0.GetTRDtrkltOccupancy(%d)>10",iLayer);
    treeV0->Draw(what.Data(),where.Data(),"goff",nPoints);
  }
  treeV0->Draw("GetDet(6)>>+histoDetector","isOK0&&isTOFIn0","goff",nPoints);
  Double_t entriesMedian=TMath::Median(10,histoDetector->GetArray());
  for (Int_t iDet=0;iDet<kNDetectors; iDet++){
    activeMap(iDet) = (histoDetector->GetBinContent(iDet+1)>threshold*entriesMedian);
  }
}

/// \brief cache per track information
/// status code:
///        0 - tracklet  found
///        1 - tracklet in active zone but not found
///        2 - track in dead zone
///        3 - track in not active detector
Double_t chacheInfo(Int_t entry){
  treeV0->GetEntry(entry);
  const Double_t tanPhi=TMath::Tan(TMath::Pi()/18.);
  const Double_t edgeF0=TMath::Sqrt(1+tanPhi*tanPhi); // sqrt(1+TMath::Tan(TMath::Pi()/18.)**2)
  nFindable=0;
  nFound=0;
  for (Int_t iLayer=0; iLayer<7; iLayer++){
    Double_t radius=geom.GetTime0(iLayer%6)-3;  // to add constant
    if (iLayer==6) radius=375;
    Double_t phi=track->GetInnerParam()->GetParameterAtRadius(radius*edgeF0,5,7);
    Double_t gZ=track->GetInnerParam()->GetParameterAtRadius(radius*edgeF0,5,2);
    if (phi<0) phi+=TMath::TwoPi();
    Double_t sector=9.*phi/TMath::Pi();
    Int_t iStack=0;
    for (iStack=0; iStack<5; iStack++) {
      if (gZ > zBoundary1(iLayer, iStack)) break;
    }

    if (iStack>4) iStack=4;
    Int_t iDet=0;
    if (iLayer<6)  iDet=geom.GetDetector(iLayer,iStack,sector);    // TRD numbering
    if (iLayer==6) iDet=540+sector*5+iStack;                       // TOF numbering
    vecDet[iLayer]=iDet;
    vecSec[iLayer]=sector;
    vecZ[iLayer]=gZ;
    Double_t dSec=sector-int(sector);
    Double_t dEdge=(dSec<0.5) ? dSec*radius*tanPhi:(1-dSec)*radius*tanPhi;
    vecdSec[iLayer]=dSec;
    vecdEdge[iLayer]=dEdge;
    vecStatus[iLayer]=0;
    vecActive[iLayer]=0;
    //
    Bool_t isActive=0;
    Bool_t isFindable = kFALSE;
    Bool_t isActiveZ=(gZ<zBoundary0(iLayer,iStack) && gZ>zBoundary1(iLayer, iStack));
    Bool_t isDeadR=TMath::Abs(dEdge)<kMarginR;
    Bool_t isDeadZ=(!(gZ<zBoundary0(iLayer,iStack)-kMarginZ && gZ>zBoundary1(iLayer, iStack)+kMarginZ));   /// add savety margin 1 cm because of misalignemnt for TOF not sure of z bounaries
    Bool_t isDeadDet=activeMap(iDet)<0.5;
    Int_t status=0;
    if (iLayer<6) {
      vecActive[iLayer]=(track->GetTRDtrkltOccupancy(iLayer)>occuCut);
      if ( vecActive[iLayer] > 0) {
        nFindable++;
        nFound++;
        isFindable=kTRUE;
        status=0;
        isActive=kTRUE;
      } else {
        isFindable=! (isDeadR || isDeadZ || isDeadDet);
        if (nFindable) nFindable++;
      }
    }else{
      vecActive[iLayer]=track->IsOn(0x2000);
      isActive=track->IsOn(0x2000);
      isFindable=! (isDeadR || isDeadZ || isDeadDet);
    }
    vecDeadZ[iLayer]=isDeadZ;
    vecDeadR[iLayer]=isDeadR;
    vecDeadDet[iLayer]=isDeadDet;
    if (vecActive[iLayer]>0.5) continue;
    if (isFindable){
      vecStatus[iLayer]=1;   // not found  but active
      continue;
    }
    if (isDeadDet) {vecStatus[iLayer]=3; continue;}   // not found  but active
    vecStatus[iLayer]=2; // dead zones
  }
  return 1;
}

Double_t GetDet(Int_t iLayer){ return vecDet[iLayer%7];}
Double_t GetSec(Int_t iLayer){ return vecSec[iLayer%7];}
Double_t GetZ(Int_t iLayer){ return vecZ[iLayer%7];}
Double_t GetdSec(Int_t iLayer){ return vecdSec[iLayer%7];}
Double_t GetdEdge(Int_t iLayer){ return vecdEdge[iLayer%7];}
Double_t GetStatus(Int_t iLayer){ return vecStatus[iLayer%7];}
Double_t isActive(Int_t iLayer){ return vecActive[iLayer%7];}
Double_t isNotActive(Int_t iLayer) {return (vecDeadR[iLayer%7]+vecDeadZ[iLayer%7]+vecDeadDet[iLayer%7])>0;}
Double_t isDeadR(Int_t iLayer){ return vecDeadR[iLayer%7];}
Double_t isDeadZ(Int_t iLayer){ return vecDeadZ[iLayer%7];}
Double_t isDeadDet(Int_t iLayer){ return vecDeadDet[iLayer%7];}
Double_t TRDFound(){return nFound;}
Double_t TRDFindable(){return nFindable;}

Double_t GetStatus2(Int_t iLayer){ // joined status of layer iLayer and iLayer+1  / used to reduce dimensionality of lookup table
  if ((vecActive[iLayer%7]+vecActive[(iLayer+1)%7])>0)     return 0;                                  // at minimum one layer active
  if (((vecStatus[iLayer%7]==1)+(vecStatus[(iLayer+1)%7])==1)>0)    return 1;                         // both layer not active
  if ((vecDeadDet[iLayer%7]+vecDeadDet[(iLayer+1)%7])>0)   return 3;                                  // not responding detector
  if (vecDeadR[iLayer%7]+vecDeadZ[iLayer%7] +vecDeadR[(iLayer+1)%7]+vecDeadZ[(iLayer+1)%7]) return 2; // dead zone cross
  return 4; // undefined - should not happen
}

void InitTreeTrack(){
  tree = AliXRDPROOFtoolkit::MakeChain("filtered.list","highPt",0,1000);
  tree->SetBranchAddress("esdTrack.",&track);
  //
  tree->SetAlias("isOK0","esdTrack.fITSncls>3&&abs(esdTrack.fP[4])<2.5&&abs(esdTrack.fP[3])<1&&chacheInfo(Entry$)&&tofLength>300");
  tree->SetAlias("tofLength","tofClInfo.fElements[5]");
  tree->SetAlias("tofTime","tofClInfo.fElements[0]");
  tree->SetAlias("tofOK","abs(tofClInfo.fElements[3])<6&&abs(tofClInfo.fElements[4])<6&&tofTime>0");
}

void InitTree() {
  CacheTRDGeom();
  treeV0 = AliXRDPROOFtoolkit::MakeChain("filtered.list","V0s",0,1000);
  treeV0->SetAlias("isITSOn0","(track0.fFlags&0x1)>0");
  treeV0->SetAlias("isITSOn1","(track1.fFlags&0x1)>0");
  treeV0->SetAlias("idProton0", "type==4&&abs(tpcNsigma0.fElements[4])<2"); // stong criteria to suppress backround
  treeV0->SetAlias("idProton1", "type==2&&abs(tpcNsigma1.fElements[4])<2"); // stong criteria to suppress backround
  treeV0->SetAlias("isTOFIn0","abs(tofClInfo0.fElements[3])<20");
  treeV0->SetAlias("isTOFGold0","abs(tofClInfo0.fElements[3])<2");
  treeV0->SetAlias("isTOFProton0","abs(tofNsigma0.fElements[4])<6");
  treeV0->SetAlias("isTOFIn1","abs(tofClInfo1.fElements[3])<20");
  treeV0->SetAlias("isTOFProton1","abs(tofNsigma1.fElements[4])<6");
  //
  treeV0->SetAlias("dTRD0","(9*track0.fIp.GetParameterAtRadius(305,5,7)/pi+18)-int(9*track0.fIp.GetParameterAtRadius(305,5,7)/pi+18)");
  treeV0->SetAlias("dTRD1","(9*track0.fIp.GetParameterAtRadius(318,5,7)/pi+18)-int(9*track0.fIp.GetParameterAtRadius(318,5,7)/pi+18)");
  //treeV0->Draw("(9*track0.fIp.GetParameterAtRadius(300,5,7)/pi+18)-int(9*track0.fIp.GetParameterAtRadius(300,5,7)/pi+18)","idProton0&&track0.Pt()>0.5","",100000);
  treeV0->SetBranchAddress("track0.",&track);
  //
  treeV0->SetAlias("isOK0","abs(track0.fP[4])<2.5&&abs(track0.fP[3])<1&&chacheInfo(Entry$)");
  //treeV0->SetAlias("notActive")
}

/// code to copy paste to test variables
void DrawExample() {
  //
  treeV0->Draw("track0.GetTRDntracklets()","isITSOn0&&idProton0&&track0.Pt()>0.5&&isTOFIn0","",200000);
  treeV0->Draw("GetTRDtrkltOccupancy(0)>0:dTRD0>>his(50,0,1)","idProton0&&track0.Pt()>1","prof",100000);
  //
  treeV0->Draw("track0.GetTRDtrkltOccupancy(0)>10:GetdEdge(0)","chacheInfo(Entry$)&&track0.fTRDncls>20","prof",10000);
  //
  treeV0->Draw("!isActive(0):GetdSec(0)","isOK0&&isTOFIn0&&!(isDeadZ(0)||isDeadDet(0)||isDeadR(0))","prof",50000);
  //

  TString hisString="";
  hisString+="GetZ(0):#isActive(0)>>hisTRDZ0(800,-400,400);";
  hisString+="GetZ(1):#isActive(1)>>hisTRDZ1(800,-400,400);";
  hisString+="GetZ(2):#isActive(2)>>hisTRDZ2(800,-400,400);";
  hisString+="GetZ(3):#isActive(3)>>hisTRDZ3(800,-400,400);";


}



void makeNucleiPlots(){
  tree->Draw("tofLength/tofTime:esdTrack.P()>>hisTOF_TRD2(300,0.5,3,300,0,0.04)","tofOK&&selectionPIDMask>0&&esdTrack.GetTRDntracklets()>=2","colz",5000000);
  gPad->SaveAs("hisTOF_TRD2.pdf");
  tree->Draw("tofLength/tofTime:esdTrack.P()>>hisTOF_TRD4(300,0.5,3,300,0,0.04)","tofOK&&selectionPIDMask>0&&esdTrack.GetTRDntracklets()>=4","colz",5000000);
  gPad->SaveAs("hisTOF_TRD4.pdf");
  tree->Draw("tofLength/tofTime:esdTrack.P()>>hisTOF_TRD6(300,0.5,3,300,0,0.04)","tofOK&&selectionPIDMask>0&&esdTrack.GetTRDntracklets()>=6","colz",5000000);
  gPad->SaveAs("hisTOF_TRD6.pdf");
  //
  /*
  hisTOF_TRD2->SetMarkerStyle(21); hisTOF_TRD4->SetMarkerStyle(21); hisTOF_TRD6->SetMarkerStyle(21);
  hisTOF_TRD2->SetMarkerSize(0.5); hisTOF_TRD4->SetMarkerSize(0.5); hisTOF_TRD6->SetMarkerSize(0.5);

  x=hisTOF_TRD2->ProjectionY("TRD2_p1",100,110); x->SetLineColor(1); x.SetMarkerSize(0.4); x->Draw(""); x->SetTitle("TRD2");
  x=hisTOF_TRD4->ProjectionY("TRD4_p1",100,110); x->SetLineColor(2); x.SetMarkerSize(0.4);x->Draw("same"); x->SetTitle("TRD4");
  x=hisTOF_TRD6->ProjectionY("TRD6_p1",100,110); x->SetLineColor(4); x.SetMarkerSize(0.4);x->Draw("same"); x->SetTitle("TRD6");
  gPad->BuildLegend(0.15,0.55,0.4,0.85)
  gPad->SaveAs("hisTOF_TRDSliceExample.pdf");

  x=hisTOF_TRD2->ProjectionY("TRD2_p1",100,105); x->SetLineColor(1); x.SetMarkerSize(0.4); x->Draw(""); x->SetTitle("TRD2");
  x=hisTOF_TRD4->ProjectionY("TRD4_p1",100,105); x->SetLineColor(2); x.SetMarkerSize(0.4);x->Draw("same"); x->SetTitle("TRD4");
  x=hisTOF_TRD6->ProjectionY("TRD6_p1",100,105); x->SetLineColor(4); x.SetMarkerSize(0.4);x->Draw("same"); x->SetTitle("TRD6");
  gPad->BuildLegend(0.15,0.55,0.4,0.85)
  gPad->SaveAs("hisTOF_TRDSliceExample.pdf");
*/

}