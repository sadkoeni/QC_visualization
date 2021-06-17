///
/// get local count statistics
///
/*
gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
gSystem->AddIncludePath("-I\$NOTES/");
.L $NOTES/JIRA/ATO-432/code/tpcLocClusterLoop.C++g
processTree("clusterDump700HLT.root")

processTree("merged_clusterDumpFullHLT.root")
*/

//gInterpreter->AddIncludePath("\$ALICE_ROOT/include/")

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVectorF.h"
#include "TTreeStream.h"
#include "TTreePlayer.h"
#include "AliSysInfo.h"
#include <TH2F.h>
#include <TH2F.h>
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>

//ClassImp(AliSysInfo)


Int_t sectorLast=-1;
TTree * tree =0;
TTreeSRedirector *locCntStream=0;
TTreeSRedirector *QAStream=0;
TObjString fPath=0;
TFile* f;
Bool_t locBadRows[72][96];
Bool_t gloBadRows[72][96];

Double_t clPadCnt=0;
Double_t clPadSum=0;

Double_t firstPadCl=0;
Double_t lastPadCl=0;
int ilastPadCl=0;
int firstPadRowCl=0;
int lastPadRowCl=0;

std::vector<Double_t>  firstPad[96];     //vectors of cluster starts and ends in a padrow
std::vector<Double_t>  lastPad[96];

void search(int padRow, int pad);

void getBadRows(Bool_t gl=kFALSE, Bool_t loadGl=kTRUE);

void InitTree(TObjString fName){
    AliSysInfo::AddStamp("Init_start",0,0,0);
    f= TFile::Open(fName.String());
    tree = (TTree*)f->Get("clusterDump700");
//    tree = (TTree*)f->Get("clusterDumpFull");
    locCntStream = new TTreeSRedirector(TString::Format("locTrk700HLT.root"),"recreate");
    char cwd[100];
    getcwd(cwd, sizeof(cwd));
    fPath=TObjString(cwd);  
}

void processTree(TObjString fClFile){

  
  TH2F *histo = new TH2F("histo","",100,0,0,96,0,0);
  TH1F *histoPadRows = new TH1F("histo","",96,0,96);

  InitTree(fClFile);
  ULong64_t gid=-1;
  ULong64_t gidCurr=-1;
  std::vector<Int_t> padLengthSum;
  std::vector<Int_t> padRows;
  std::vector<Int_t> clusterSize;
  
  int binCnt=-1;
  int padRowCnt=-1;
  int padRowFirst=-1;
  int binxtemp=-1;
  Int_t entriesCluster  = tree->GetEntries();
  tree->SetBranchAddress("gid",&gid);
  
//  getBadRows(kFALSE,kFALSE);       //run local pad row QA

  TH1::AddDirectory(kTRUE);
  
  //loop over events and get number of local counts
  for (Int_t iCl=0; iCl<entriesCluster; iCl++){
    tree->GetEntry(iCl);
    if(gid!=gidCurr){    //if new event or last entry, then write info about previous event
      std::cout<<"prev gid: "<<gidCurr<<" next: "<<gid<<std::endl;
//      for (Int_t sec=0; sec<72; sec++){     
        cout<<"\nSECTOR: "<<sec<<" "<<"fRow:fPad>>histo1(100,,,96,0,96)"<<TString::Format("gid==%llu ",gid)<<endl;
        int cnt=tree->Draw("fRow:fPad>>histo1(100,,,96,0,96)",TString::Format("gid==%llu ",gid),"goff");
        tree->GetEntry(iCl);
        histo = (TH2F*)gDirectory->Get("histo1");
        histo->Draw();
        Int_t ent=histo->GetEntries();
        
        //find pad clusters starts/ends in each padrow
//        for(int binx=1; binx<histo->GetNbinsX()+1; binx++){       //loop over time/z
          for(int biny=1; biny<96; biny++) {  //loop over pad rows
            
            if(locBadRows[sec][biny-1] || gloBadRows[sec][biny-1]) continue;  //mask out bad rows from local chunk and from global                                                   //sum of HLT clusters contributing to track cluster
            for(int binx=1; binx<histo->GetNbinsX()+1; binx++){       //loop over pads in pad row
//              if(histo->GetBinContent(binx,biny))std::cout<<"pad: "<<histo->GetXaxis()->GetBinCenter(binx)<<" padrow: "<<histo->GetYaxis()->GetBinCenter(biny)<<" histo->GetBinContent(binx,biny): "<<histo->GetBinContent(binx,biny)<<std::endl;
              if(histo->GetBinContent(binx,biny)>0){
                for(binxtemp=binx; binxtemp<histo->GetNbinsY()+1; binxtemp++){       //loop over pads in pad row  in cluster to find its end
                  if(histo->GetBinContent(binxtemp,biny)==0){ 
                    break;
                  }
                }
                firstPad[biny].push_back(histo->GetXaxis()->GetBinCenter(binx));
                lastPad[biny].push_back(histo->GetXaxis()->GetBinCenter(binxtemp));
              }          
            }      
          }
//        }  

        for(int i=0; i<96; i++){
          for(int j=0; j < firstPad[i].size(); j++) {
//            std::cout<<"sector: "<<sec<<" first size: "<<firstPad[i].size()<<" last size: "<<lastPad[i].size()<<std::endl;
//            std::cout<<"Pad row: "<<i<<" first: "<<firstPad[i][j]<<" last: "<<lastPad[i][j]<<std::endl;
          }
        }
        
        //find (meta-) clusters along pads and padrows
        for(int i=0; i<96-1; i++){                            //loop over padrows
//        for(int i=0; i<5; i++){                            //loop over padrows
          for(int j=0; j < firstPad[i].size(); j++) {      //loop over clusters in padrow
              firstPadCl=firstPad[i][j];
//              cout<<firstPad[i].size()<<" "<<lastPad[i].size()<<endl;
              firstPadRowCl=i;
//              cout<<"search: for pad row: "<<i<<", cluster# "<<j<<": first padrow: "<<firstPadRowCl<<" first pad: "<<firstPadCl<<endl;
              search(i,j);   
              //write info for each sector
              gidCurr=gid;
              (*locCntStream)<< "locTrk" <<
                    "fPath.="<<&fPath<<                
                    "gid="<<gidCurr<<
//                    "sector="<<sec<<
                    "firstPadCl="<<firstPadCl<<
                    "lastPadCl="<<lastPadCl<<
                    "firstPadRowCl="<<firstPadRowCl<<
                    "lastPadRowCl="<<lastPadRowCl<<
                    "clPadPadRowArea="<<clPadCnt<<  
                    "\n"; 
          }
        }
    
  clPadCnt=0;
  clPadSum=0;

  firstPadCl=0;
  lastPadCl=0;
  ilastPadCl=0;
  firstPadRowCl=0;
  lastPadRowCl=0;  
  
 // reset vectors   
      for(int i=0; i<96; i++){
        firstPad[i].clear();
        lastPad[i].clear();
      }       
//  }       

        
      

  }
  if(iCl==entriesCluster-1) break;    
  }
  delete locCntStream;
  
//  TFile* fLoc= TFile::Open("locCount700HLT.root");
//  TTree* locCnt =  (TTree*)fLoc->Get("locCnt");
////  TString sel="Iteration$==4&&binCntThreshV.fElements>20";
//  TString sel="Iteration$==3&&binCntThreshV.fElements>6";   // loose cut for local testing
//  
//  //write selected gids to text file
//  ((TTreePlayer*)(locCnt->GetPlayer()))->SetScanRedirect(true);
//  ((TTreePlayer*)(locCnt->GetPlayer()))->SetScanFileName("gid.List");
//  locCnt->Scan("gid:sector:binCntThreshV.fElements",sel.Data(),"colsize=30");
//  
//  //write selected gids to root tree
//  TFile* fLocOut= TFile::Open("selectedGIDs.root","RECREATE");
//  TTree* selGid = locCnt->CopyTree(sel.Data());
//  selGid->Write();
  
//  AliSysInfo::AddStamp("Process_stop",0,0,0);
//  TFile* fProf= TFile::Open("tpcLocClusterLoopProf.root","RECREATE");
//  TTree * trProf = AliSysInfo::MakeTree("syswatch.log");
//  rename("syswatch.log","syswatchloc.log");       //TODO: add filen to name!
//  trProf->Write();
  
  TH1::AddDirectory(kFALSE);
}


void search(int i, int j){                //i ... padrow, j ... pad
  //find last pad and pad row in this cluster
//  cout<<"Start search for: pad row index: "<<i<<" pad index:"<<firstPad[i][j]<<endl;
  clPadCnt=lastPad[i][j]-1;
  lastPadRowCl=i;
  ilastPadCl=j;
  lastPadCl=firstPad[i][j];
  for(int m=i+1; m<96-1; m++){                    //loop over next padrows
    if(firstPad[m].size()==0){
//      if(lastPadRowCl-i>1)  std::cout<<"FOUND cl size:"<<clPadCnt<<" START: padrow "<<i<<" padfirst: "<<firstPad[i][j]<<" END padrow "<<lastPadRowCl<<" padfirst "<<lastPadCl<<std::endl;
      return;
    }
    for(int k=0; k < firstPad[m].size(); k++) {    //loop over clusters in next padrow
      if(abs(lastPadCl-firstPad[m][k])<5){
//        std::cout<<"checking pr:"<<m<<" pad: "<<k<<" first pad:"<<firstPad[m][k]<<"  abs(firstPad[lastPadRowCl][ilastPadCl]-firstPad[m][k]): "<< firstPad[lastPadRowCl][ilastPadCl]<<"-"<<firstPad[m][k]<<std::endl;
//        std::cout<<"compatible: padrow "<<lastPadRowCl<<" pad: "<<lastPadCl<<" and padrow "<<m<<" pad "<<firstPad[m][k]<<std::endl;
        lastPadRowCl=m;                             //last pad row of cluster 
        lastPadCl=firstPad[m][k];                   //first pad in last pad row pf cluster
        ilastPadCl=k;                               //index of first pad in last pad row pf cluster
        clPadCnt+=abs(lastPad[m][k]-firstPad[m][k]);
        if(k < firstPad[m].size()-1 && abs(firstPad[m][k+1]-firstPad[m][k+1])<3){   //check if next cluster is very close and delete it if it is
          clPadCnt+=abs(lastPad[m][k]-firstPad[m][k]);
          firstPad[m].erase(firstPad[m].begin() + k+1);   //erase pad clusters that were associated to a cluster, otherwise they will be used as seeds again
          lastPad[m].erase(lastPad[m].begin() + k+1);
        }
        firstPad[m].erase(firstPad[m].begin() + k);   //erase pad clusters that were associated to a cluster, otherwise they will be used as seeds again
        lastPad[m].erase(lastPad[m].begin() + k);  

        break;  //look no further in this pad row - go to next padrow
      }
      else if(k == firstPad[m].size()-1){
//        if(lastPadRowCl-i>1)  std::cout<<"FOUND cl size:"<<clPadCnt<<" START: padrow "<<i<<" padfirst: "<<firstPad[i][j]<<" END padrow "<<lastPadRowCl<<" padfirst "<<lastPadCl<<std::endl;
        return;  // nothing in whole pad row -> quit search
      }
    }
  }
//  if(lastPadRowCl-i>1)std::cout<<"FOUND Cross cl size:"<<clPadCnt<<" compatible: padrow "<<i<<" padfirst: "<<firstPad[i][j]<<" and padrow "<<lastPadRowCl<<" padfirst "<<lastPadCl<<std::endl;
  return;
}

//Find bad pad rows from local chunk (gl=kFALSE) or from global data set (gl=kTRUE). In the later case can either create QA info or load it (loadGl).
/*
gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
gSystem->AddIncludePath("-I\$NOTES/");
.L $NOTES/JIRA/ATO-432/code/tpcLocClusterLoop.C++g
getBadRows(kTRUE,kFALSE)   //to create QA tree and histos for global pad row QA
 */
void getBadRows(Bool_t gl, Bool_t loadGl){
  
  TH1::AddDirectory(kTRUE);
  
  TH2F *hRowOROC = new TH2F("hRowOROC","",36,36,72,96,0,96);
  TH2F *hRowIROC = new TH2F("hRowIROC","",36,0,36,63,0,63);
  
  Double_t cutIROC;
  Double_t cutOROC;
  TFile* fm, *fqa;
  TTree * t; 
  Int_t sec, row;
  Double_t cntRat, mean;
  
  if(gl){
    std::cout<<"\nRunning Pad Row QA on whole data set!\n"<<std::endl;
    fm= TFile::Open("merged_clusterDump700HLT.root","READ");
    if(fm==0x0) {
      std::cout<<"Error: merged_clusterDump700HLT.root not found!!"<<std::endl;
      return;
    }
    if(loadGl) fqa=TFile::Open("QA700HLTGlobal.root","READ");
//    t =(TTree*)fm->Get("clusterDump700");
    t =(TTree*)fm->Get("clusterDumpFull");
    if(!loadGl) QAStream = new TTreeSRedirector("QA700HLTGlobal.root","recreate");
    cutIROC=1.7;
    cutOROC=4;
  }
  else{
    std::cout<<"\nRunning Pad Row QA on local data chunk!\n"<<std::endl;
//    fm= TFile::Open("clusterDump700HLT.root","READ");
    fm= f;
    if(fm==0x0) {
      std::cout<<"Error: clusterDump700HLT.root not found!! File needs init from InitTree()"<<std::endl;
      return;
    }
    t =(TTree*)fm->Get("clusterDumpFull");
    QAStream = new TTreeSRedirector("QA700HLTLocal.root","recreate");
    cutIROC=6;
    cutOROC=10;
  }

  //Loop over OROCS
  if(gl && loadGl){
    hRowOROC = (TH2F*)fqa->Get("hRowOROC");
  }
  else {
    t->Draw("fRow:fDetector>>hRowOROC(36,36,72,96,0,96)","fDetector>35","colz goff");
    hRowOROC = (TH2F*)gDirectory->Get("hRowOROC");//->Clone();
  }
//  cout<<"matching tree entries: "<<t->Draw("fRow:fDetector>>hRowOROC(36,36,72,96,0,96)","fDetector>35","colz goff")<<endl;
  mean=hRowOROC->Integral()/(hRowOROC->GetNbinsX()*hRowOROC->GetNbinsY());
  std::cout<<"\nmean counts in OROC per pad row and sector : "<<mean<< ", cutOROC: "<<cutOROC<<std::endl;
  for(Int_t binx=1; binx<hRowOROC->GetNbinsX()+1; binx++){
    for(Int_t biny=1; biny<hRowOROC->GetNbinsY()+1; biny++){
      sec=36+binx-1;
      row=biny-1;
      cntRat=hRowOROC->GetBinContent(binx,biny)/mean;
      if(!gl || !loadGl) {
        (*QAStream)<< "QA" <<
        "sector="<<sec<<
        "row="<<row<<
        "Counts700overMean="<<cntRat<<
        "\n";
      }
      if(cntRat>cutOROC)  {
        if(!gl) locBadRows[36+binx-1][biny-1]=kTRUE;
        else    gloBadRows[36+binx-1][biny-1]=kTRUE;
        std::cout<<"mask out sector: "<<36+binx-1<<", pad row: "<<biny-1<<" content/mean = "<<cntRat<<std::endl;
      }
      else {
        if(!gl) locBadRows[36+binx-1][biny-1]=kFALSE;
        else    gloBadRows[36+binx-1][biny-1]=kFALSE;
      }
      hRowOROC->SetBinContent(binx,biny, cntRat);
    }
  }
  if(!gl || !loadGl){
    hRowOROC->Write();
  }
  
  //Loop over IROCS  
  if(gl && loadGl){
    hRowIROC = (TH2F*)fqa->Get("hRowIROC");
  }
  else {
    t->Draw("fRow:fDetector>>hRowIROC(36,0,36,63,0,63)","fDetector<35","colz");
    hRowIROC = (TH2F*)gDirectory->Get("hRowIROC");
  }
  mean=hRowIROC->Integral()/(hRowIROC->GetNbinsX()*hRowIROC->GetNbinsY());
  std::cout<<"\nmean counts in IROC per pad row and sector : "<<mean<< ", cutIROC: "<<cutIROC<<std::endl;
  for(int binx=1; binx<hRowIROC->GetNbinsX()+1; binx++){
    for(int biny=1; biny<hRowIROC->GetNbinsY()+1; biny++){
      sec=binx-1;
      row=biny-1;
      cntRat=hRowIROC->GetBinContent(binx,biny)/mean;
      if(!gl || !loadGl){
        (*QAStream)<< "QA" <<
        "sector="<<sec<<
        "row="<<row<<
        "Counts700overMean="<<cntRat<<
        "\n";    
      }
      if(cntRat>cutIROC)  {
        if(!gl) locBadRows[binx-1][biny-1]=kTRUE;
        else    gloBadRows[binx-1][biny-1]=kTRUE;
        std::cout<<"mask out sector: "<<binx-1<<", pad row: "<<biny-1<<" content/mean = "<<cntRat<<std::endl;
      }
      else{
        if(!gl) locBadRows[binx-1][biny-1]=kFALSE;
        else    gloBadRows[binx-1][biny-1]=kFALSE;
        }
      hRowIROC->SetBinContent(binx,biny, cntRat);
    }
  }

  if(!gl || !loadGl){
    hRowIROC->Write();
  }
  
  TH1::AddDirectory(kFALSE);
  return;
}
