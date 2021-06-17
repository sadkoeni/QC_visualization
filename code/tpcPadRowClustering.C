///
/// get local count statistics
///
/*
gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
gSystem->AddIncludePath("-I\$NOTES/");
.L $NOTES/JIRA/ATO-432/code/tpcPadRowClustering.C++g
processTree("merged_clusterDumpFullHLT.root")
processTree("clusterDump700HLT.root")
*/

//gInterpreter->AddIncludePath("\$ALICE_ROOT/include/")

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TVectorF.h"
#include "TTreeStream.h"
#include "TTreePlayer.h"
#include "AliSysInfo.h"
#include <TH3F.h>
#include <TH2F.h>
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>

//ClassImp(AliSysInfo)

ULong64_t gidCurr=-1;
Int_t sectorLast=-1;
TTree * tree =0;
TTreeSRedirector *locCntStream=0;
TTreeSRedirector *QAStream=0;
TObjString fPath=0;
TFile* f;
Bool_t locBadRows[72][96];
Bool_t gloBadRows[72][96];


const Int_t nTOT=2;     //define sector selection criteria: require nPads pad that have a time over threshold of at least 10*nTOT HLT Clusters with fMax>700
const Int_t nPRs=20;
const Int_t nTOT1=4;     //define sector selection criteria: require nPads pad that have a time over threshold of at least 10*nTOT HLT Clusters with fMax>700
const Int_t nPRs1=10;

//const Int_t nTOT=1;     //loose cuts for local testing
//const Int_t nPRs=12;

void getBadRows(Bool_t gl=kFALSE, Bool_t loadGl=kTRUE);

void InitTree(TObjString fName){
    AliSysInfo::AddStamp("Init_start",0,0,0);
    f= TFile::Open(fName.String());
    tree = (TTree*)f->Get("clusterDump700");    //use this for normal grid running
    Int_t entriesCluster  = tree->GetEntries();
//    tree = (TTree*)f->Get("clusterDumpFull"); //was used to retrieve the selectedGIDs.root trees after grid runs
    locCntStream = new TTreeSRedirector(TString::Format("locCount700HLT.root"),"recreate");
    char cwd[100];
    getcwd(cwd, sizeof(cwd));
    fPath=TObjString(cwd);  
}

void processTree(TObjString fClFile){

  
  TH3F *histo = new TH3F("histo","",4,0,1024,40,-40,40,96,0,96);
//  TH1F *histoPR = new TH1F("histo","",1024,0,1024,96,0,96);
  TH2F *histo2D = new TH2F("histo2","",40,-40,40,96,0,96);

  InitTree(fClFile);
  ULong64_t gid=-1;
  Double_t locCnt700;
  Int_t locCnt700Nbins;
  TVectorF binCntV(11);
  TVectorF binCntSumV(11);
  int binCnt=-1;
  int padRowCnt=-1;
  int padRowFirst=-1;
  int first=0;
  int last=0;

  TVectorF maxConsPR(11);
  TVectorF found(11);
  TVectorF firstPR(11);
  TVectorF lastPR(11);
  TVectorF firstPad(11);
  TVectorF lastPad(11);
  TVectorF firstPRmax(11);
  TVectorF lastPRmax(11);
  TVectorF firstPadmax(11);
  TVectorF maxConsPad(11);


  
  Bool_t isSel=kFALSE;

  Int_t binxLow;
  Int_t binxHi;
  
  Int_t entriesCluster  = tree->GetEntries();
  tree->SetBranchAddress("gid",&gid);
  
  getBadRows(kFALSE,kFALSE);       //run local pad row QA
//  getBadRows(kTRUE,kTRUE);  //load global pad row QA
  TH1::AddDirectory(kTRUE);
  //loop over events and get number of local counts
  for (Int_t iCl=0; iCl<entriesCluster; iCl++){
	  
    
//    if(iCl%(int(entriesCluster/10000))==0)   std::cout<<"Done "<<iCl<<" out of "<<entriesCluster<<std::endl;
    
    tree->GetEntry(iCl);
    if(gid!=gidCurr){    //if new event or last entry, then write info about previous event
      std::cout<<"prev gid: "<<gidCurr<<" next: "<<gid<<std::endl;
      for (Int_t sec=0; sec<72; sec++){     
        int cnt=tree->Draw("fRow:fY:fTimeBin>>histo1(4,0,1024,40,-40,40,96,0,96)",TString::Format("gid==%llu && fDetector==%d",gid,sec),"goff");
        tree->GetEntry(iCl);		//needed otherwise loop over tree entries is not working
        histo = (TH3F*)gDirectory->Get("histo1");
        histo->Draw();
        		      
        Int_t ent=histo->GetEntries();
        locCnt700=0;
        locCnt700Nbins=0;
        for(int binx=1; binx<histo->GetNbinsX()+1; binx++){       //loop over time/z
              binxLow= histo->GetXaxis()->GetBinLowEdge(binx);
              binxHi= histo->GetXaxis()->GetBinUpEdge(binx);
              for(int exp=0; exp<=10; exp++){
                maxConsPR[exp]=-1;
                found[exp]=-1;
                firstPR[exp]=-1;
                lastPR[exp]=-1;
                firstPad[exp]=-1;
                lastPad[exp]=-1;
                firstPRmax[exp]=-1;
                lastPRmax[exp]=-1;
                firstPadmax[exp]=-1;
                maxConsPad[exp]=-1;
              }
              for(int binz=2; binz<histo->GetNbinsZ(); binz++) {  //loop over pad rows -> find largest cluster in pad row direction
                  for(int exp=0; exp<=10; exp++){       
                    found[exp]=-1;	
                }
              if(locBadRows[sec][binz-1] || gloBadRows[sec][binz-1]) continue;             //mask out bad rows from local chunk and from global

              for(int biny=0; biny<histo->GetNbinsY()+1; biny++){     //loop over pads in pad row
                  binCnt= TMath::Max(histo->GetBinContent(binx,biny,binz-1),TMath::Max(histo->GetBinContent(binx,biny,binz+1),histo->GetBinContent(binx,biny,binz)));
                  for(int exp=0; exp<=10; exp++){                   //loop over TOT cuts      
                    if(binCnt>=TMath::Nint(TMath::Power(1.5,exp))){
                      found[exp]=1;
                    }      
                  }
              }
              //handle pad row clusters
              for(int exp=0; exp<=10; exp++){
                if(found[exp]==1 && lastPR[exp]>=firstPR[exp]  ){ 	 //start of cluster
					firstPR[exp]=binz;            
				}       
                else if( (found[exp]==-1 && firstPR[exp]>lastPR[exp]) || (found[exp]==1 && (firstPR[exp]>lastPR[exp] && binz==histo->GetNbinsZ()))){ //end of cluster 
                  lastPR[exp]=binz;
                  if(lastPR[exp]-firstPR[exp] >= maxConsPR[exp]){
                    maxConsPR[exp]=lastPR[exp]-firstPR[exp];
                    firstPRmax[exp]=firstPR[exp];
                    lastPRmax[exp]=lastPR[exp];			
                    }
                  }
                   
              }        
          }
              
          for(int biny=2; biny<histo->GetNbinsY(); biny++) {  //loop over pad positions -> find largest cluster in pad direction
                for(int exp=0; exp<=10; exp++){       
                    found[exp]=-1;	
                }    
                  for(int binz=0; binz<histo->GetNbinsZ()+1; binz++){     //loop over pads row at pad position
                    if(locBadRows[sec][binz] || gloBadRows[sec][binz]) continue;             //mask out bad rows from local chunk and from global
                    binCnt= TMath::Max(histo->GetBinContent(binx,biny+1,binz),TMath::Max(histo->GetBinContent(binx,biny-1,binz),histo->GetBinContent(binx,biny,binz)));
                    for(int exp=0; exp<=10; exp++){                   //loop over TOT cuts      
                      if(binCnt>=TMath::Nint(TMath::Power(1.5,exp))){
                        found[exp]=1;
                      }      
                    }
                }
//              //handle pad clusters
              for(int exp=0; exp<=10; exp++){
                if(found[exp]==1 && lastPad[exp]>=firstPad[exp]  ){ 	 //start of cluster
					firstPad[exp]=biny;            
				}       
                else if( (found[exp]==-1 && firstPad[exp]>lastPad[exp]) || (found[exp]==1 && (firstPad[exp]>lastPad[exp] && biny==histo->GetNbinsY()))){ //end of cluster 
                  lastPad[exp]=biny;
                  if(lastPad[exp]-firstPad[exp] >= maxConsPad[exp]){
                    maxConsPad[exp]=lastPad[exp]-firstPad[exp];
                    firstPadmax[exp]=histo->GetYaxis()->GetBinCenter(firstPad[exp]);
                    }
                  }
                   
              }                 
          
          }
              
        //write info for each sector
        gidCurr=gid;
        if(maxConsPR[nTOT]>nPRs && maxConsPR[nTOT1]>nPRs1 && firstPR[nTOT]<=firstPR[nTOT1] && lastPR[nTOT]>=lastPR[nTOT1] ) isSel=kTRUE;
        else isSel=kFALSE;
        
//        if(sel){      //check time spread of leading HLT clusters to reject high eta tracks
//        tree->Draw("fRow:fTimeBin>>histo2(20,,,96,20,96)",TString::Format("gid==%llu && fDetector==%d && fTimeBin>=%d && fTimeBin<%d",gid,sec,timeBinLow,timeBinHigh),"goff");
//        
//          histoPR = (TH3F*)gDirectory->Get("histo2");
//          histoPR->Draw();
//          for(int binPR=firstPRmax[nTOT]; binPR<firstPRmax[nTOT]+maxPRcons[nTOT]; binPR++) {  //loop over pad rows -> find largest cluster in pad row direction
//            if(locBadRows[sec][binPR-1] || gloBadRows[sec][binPR-1]) continue;             //mask out bad rows from local chunk and from global
//            for(int binT=0; binT<256; binT++){
//              if(histoPR->GetBinContent(binT,binPR)>0 && wait>5) {  //start time cluster
//                
//              }
//          
//              }
//            }
//          }
        if(iCl%1000==0 || isSel){
        (*locCntStream)<< "locCnt" <<
                    "fPath.="<<&fPath<<                
                    "gid="<<gidCurr<<
                    "sector="<<sec<<
                    //~ "timeBin="<<binxCent<<
                    "timeBinLow="<<binxLow<<
                    "timeBinHigh="<<binxHi<<
                    "NBins700="<<cnt<<
                    "maxPRcons.="<<&maxConsPR<<
                    "firstPR.="<<&firstPRmax<<
//                    "lastPR.="<<&lastPRmax<<
                    "firstPad.="<<&firstPadmax<<
                    "maxPadcons.="<<&maxConsPad<<
//                    "clusterPadCnt.="<<&cntPadmax<<
                    "isSelected="<<isSel<<
                    "\n"; 
        }
//        if(isSel){
//          //get histo of first time bin for visual validity check of found clusters:
//          tree->Draw("fRow:fY>>histo2(40,-40,40,96,0,96)",TString::Format("gid==%llu && fDetector==%d && fTimeBin<256",gid,sec),"goff");
//          tree->GetEntry(iCl);	//needed otherwise loop over tree entries is not working
//          histo2D = (TH2F*)gDirectory->Get("histo2");
//          TCanvas *c1 = new TCanvas("myCanvasName1","The Canvas Title1",1600,1200);
//          histo2D->Draw("colz"); 	
//          c1->Draw();
//          c1->Print("./histo2d_1stTimeBin_"+TString::Format("GID_%llu_sec_%d.png",gid,sec)); 
//          }
        }        
      }
      if(iCl==entriesCluster-1) break;
    }
  }
  delete locCntStream;
  
  TFile* fLoc= TFile::Open("locCount700HLT.root");
  TTree* locCnt =  (TTree*)fLoc->Get("locCnt");
  TString sel=TString::Format("Iteration$ == %d && maxPRcons.fElements > %d",nTOT,nPRs); 

  
  //write selected gids to root tree
  TFile* fLocOut= TFile::Open("selectedGIDs.root","RECREATE");
  TTree* selGid = locCnt->CopyTree("");
  selGid->Write();
  
  AliSysInfo::AddStamp("Process_stop",0,0,0);
  TFile* fProf= TFile::Open("tpcPadRowClusteringProf.root","RECREATE");
  TTree * trProf = AliSysInfo::MakeTree("syswatch.log");
  rename("syswatch.log","syswatchloc.log");       //TODO: add filen to name!
  trProf->Write();
  
  TH1::AddDirectory(kFALSE);
}




//Find bad pad rows from local chunk (gl=kFALSE) or from global data set (gl=kTRUE). In the later case can either create QA info or load it (loadGl).
/*
gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
gSystem->AddIncludePath("-I\$NOTES/");
.L $NOTES/JIRA/ATO-432/code/tpcPadRowClustering.C++g
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
    t =(TTree*)fm->Get("clusterDump700");
//    t =(TTree*)fm->Get("clusterDumpFull");
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
    t =tree;
//    t =(TTree*)fm->Get("clusterDump700");
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
