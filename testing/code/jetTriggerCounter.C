/*
   AliDrawStyle::SetDefaults();
   AliDrawStyle::ApplyStyle("figTemplate");
   gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
   .L $NOTES/JIRA/ATO-432/code/jetTriggerCounter.C+
   LoadChains(10);

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


TTree *treeEvent=0;
TTree *treeV0=0;
TObjArray *triggerMask;
TObjArray *triggerMeaning;
void LoadChains(Int_t nFiles);
void MakeHistogramsMult(Int_t nEvents);
void MakeHistogramsV0(Int_t nEvents);
TStopwatch timer;
///
/// \param nFiles
void LoadChains(Int_t nFiles) {
  // nFiles=10
  TObjArray * triggerList = ( gSystem->GetFromPipe("cat trigger.list | grep -v \"#\"")).Tokenize("\n");
  triggerMask=new TObjArray(triggerList->GetEntries());
  triggerMeaning=new TObjArray(triggerList->GetEntries());
  for (Int_t i=0; i<triggerList->GetEntries(); i++){
    TObjArray *arr = TString(triggerList->At(i)->GetName()).Tokenize("//");
    TString name = arr->At(0)->GetName();
    name.ReplaceAll("-","_");
    name.ReplaceAll(" ","");name.ReplaceAll("\t","");
    triggerMask->AddAt(new TObjString(name),i);
    triggerMeaning->AddAt(new TObjString(arr->At(1)->GetName()),i);
    delete arr;
  }
  //
  treeEvent = AliXRDPROOFtoolkit::MakeChainRandom("filtered.list", "eventInfoV0", 0, nFiles);
  treeV0 = AliXRDPROOFtoolkit::MakeChainRandom("filtered.list", "V0s", 0, nFiles);
  AliAnalysisTaskFilteredTree::SetDefaultAliasesV0(treeV0);
  treeEvent->SetAlias("hmTrigger", "triggerClass.String().Contains(\"CVHM\")!=0");
  treeEvent->SetAlias("trdJetTrigger", "triggerClass.String().Contains(\"CINT7HJT\")!=0||triggerClass.String().Contains(\"CINT8HJT\")!=0");
  treeEvent->SetAlias("caloTrigger","triggerClass.String().Contains(\"CDMC7\")!=0||triggerClass.String().Contains(\"CEMC7\")!=0");
  treeV0->SetAlias("hmTrigger", "triggerClass.String().Contains(\"CVHM\")!=0");
  treeV0->SetAlias("trdJetTrigger", "triggerClass.String().Contains(\"CINT7HJT\")!=0||triggerClass.String().Contains(\"CINT8HJT\")!=0");
  treeV0->SetAlias("caloTrigger","triggerClass.String().Contains(\"CDMC7\")!=0||triggerClass.String().Contains(\"CEMC7\")!=0");
  treeV0->SetAlias("radius","v0.fRr");
  //
  for (Int_t i=0; i<triggerList->GetEntries(); i++) {
    treeEvent->SetAlias(TString::Format("trigger%s",triggerMask->At(i)->GetName()).Data(),TString::Format("triggerClass.String().Contains(\"%s\")!=0",triggerMask->At(i)->GetName()).Data());
    treeV0->SetAlias(TString::Format("trigger%s",triggerMask->At(i)->GetName()).Data(),TString::Format("triggerClass.String().Contains(\"%s\")!=0",triggerMask->At(i)->GetName()).Data());
  }
  delete triggerList;
}

/// \brief make multiplity nad V0 inv mass histograms
/// \param nEvents
void MakeHistogramsMult(Int_t nEvents) {
  //Int_t nEvents = 100000
  TString hisString = "";
  hisString += "mult:#!hmTrigger>>hisMultNotHM(200,0,200);";
  hisString += "mult:#hmTrigger>>hisMultHM(200,0,200);";
  hisString += "mult:#trdJetTrigger>>hisMultTRDJetTrigger(200,0,200);";
  hisString += "mult:#caloTrigger>>hisMultCaloTrigger(200,0,200);";
  for (Int_t i=0; i<triggerMask->GetEntries(); i++){
    hisString+=TString::Format("mult:#trigger%s>>hisMult%s(200,0,200);",triggerMask->At(i)->GetName(),triggerMask->At(i)->GetName());
    hisString+=TString::Format("ntracks:#trigger%s>>hisNTracks%s(200,0,200);",triggerMask->At(i)->GetName(),triggerMask->At(i)->GetName());
  }
  TObjArray *hisArray = AliTreePlayer::MakeHistograms(treeEvent, hisString, "mult>0", 0, nEvents, -1, 15);
  for (Int_t i = 0; i < hisArray->GetEntries(); i++) gROOT->Add(hisArray->At(i));
  TFile *f = TFile::Open("hisOutputMult.root","recreate");
  hisArray->Write();
  delete f;
}

/// \brief make multiplity nad V0 inv mass histograms
/// \param nEvents
void MakeHistogramsV0(Int_t nEvents) {
  //Int_t nEvents = 100000
  TString hisString = "";
  hisString+=TString::Format("K0Delta:mpt:radius:#type==8>>hisK0All(100,-0.1,0.1,50,0,1,5,0,20);");
  hisString+=TString::Format("LDelta:mpt:radius:#type==4>>hisLambdaAll(100,-0.03,0.03,50,0,1,5,0,20);");
  hisString+=TString::Format("ALDelta:mpt:radius:#type==2>>hisALambdaAll(100,-0.03,0.03,50,0,1,5,0,20);");

  for (Int_t i=0; i<triggerMask->GetEntries(); i++){
    hisString+=TString::Format("K0Delta:mpt:radius:#trigger%s&&type==8>>hisK0%s(100,-0.1,0.1,50,0,1,5,0,20);",triggerMask->At(i)->GetName(),triggerMask->At(i)->GetName());
    hisString+=TString::Format("LDelta:mpt:radius:#trigger%s&&type==4>>hisLambda%s(100,-0.03,0.03,50,0,1,5,0,20);",triggerMask->At(i)->GetName(),triggerMask->At(i)->GetName());
    hisString+=TString::Format("ALDelta:mpt:radius:#trigger%s&&type==2>>hisALambda%s(100,-0.03,0.03,50,0,1),5,0,20;",triggerMask->At(i)->GetName(),triggerMask->At(i)->GetName());
  }
  TObjArray *hisArray = AliTreePlayer::MakeHistograms(treeV0, hisString, "cutV0&&radius<20&&sqrt(kf.fChi2)<3", 0, nEvents, -1, 15);
  for (Int_t i = 0; i < hisArray->GetEntries(); i++) gROOT->Add(hisArray->At(i));
  TFile *f = TFile::Open("hisOutputV0.root","recreate");
  hisArray->Write();
  delete f;
}

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
    ::Info("MakeMaps", "%s\t%d\t%d\t%f", hisInput->GetName(), hisInput->GetNdimensions(), hisInput->GetNbins(), (hisInput->GetEntries()));
    TStatToolkit::MakeDistortionMapFast(hisInput,pcstream,projectionInfo,0,0.05);
  }
  delete pcstream;
  ::Info("MakeMaps","End");
  timer.Print();
}

///
/// \param sregexp  - input tree selection
/// \param addName  - switch to add directory name into alias name
/// \return          - tree with friend tress and aliases
TTree* LoadSummaryTree(TString sregexp, Bool_t addName) {
  //
  TTree *treeLoad = NULL;
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
  return treeLoad;
}



void MakeDraw() {
  TTree *treeMap = LoadSummaryTree("his.*", 0);
  treeMap->GetListOfFriends()->ls();
  TLegend *legend = new TLegend(0.16, 0.15, 0.8, 0.4, "K0 yield (2017)");
  legend->SetNColumns(1);

  TMultiGraph *gr = TStatToolkit::MakeMultGraph(treeMap, "",
                                                "hisK0AllDist.entries;hisK0CVHMDist.entries;hisLambdaCINT7HJTDist.entries:mptCenter",
                                                "radiusCenter<6&&mptCenter<0.2", "25;21;22;23", "1;2;4;6", 0, 1, 3, legend);
  TStatToolkit::DrawMultiGraph(gr, "ap");
  legend->Draw();
  TLegend *legendRes = new TLegend(0.16, 0.15, 0.8, 0.4, "K0 inv mass peak width (2017)");
  TMultiGraph *grRes = TStatToolkit::MakeMultGraph(treeMap, "",
                                                "hisK0AllDist.rmsG;hisK0CVHMDist.rmsG;hisLambdaCINT7HJTDist.rmsG:mptCenter",
                                                "radiusCenter<6&&mptCenter<0.2", "25;21;22;23", "1;2;4;6", 0, 1, 3, legendRes);
  TStatToolkit::DrawMultiGraph(grRes, "ap");
  legendRes->Draw();
}



/*
Int_t nEvents=5000000;
  tree->Draw("mult>>hisNotHM(200,0,200)","mult>0&&!hmTrigger","",nEvents);
  tree->Draw("mult>>hisHM(200,0,200)","mult>0&&hmTrigger","",nEvents);
  tree->Draw("mult>>hisTRD(200,0,200)","mult>0&&trdJetTrigger","",nEvents);
  tree->Draw("mult>>hisCalo(200,0,200)","mult>0&&caloTrigger","",nEvents);
*/
