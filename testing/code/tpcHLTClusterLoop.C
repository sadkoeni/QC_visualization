///
/// analyze HLT clusters tree
///
/*
gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
gSystem->AddIncludePath("-I\$NOTES/");
.L $NOTES/JIRA/ATO-432/code/tpcHLTClusterLoop.C++g

processTree(0,kTRUE)    
processTree(0,kFALSE)    
 */

#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "AliTPCclusterMI.h"
#include "TVectorF.h"
#include <TH3D.h>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>
#include "TTreeStream.h"
#include "AliSysInfo.h"
#include <TH2D.h>
#include <iostream>
#include "TMath.h"

ULong64_t gidCurr=-1;
Int_t sectorLast=-1;
TTree * tree =0;
TTree * gidTree =0;
TObjString fName=0;
TObjString fPath=0;
TTreeSRedirector *clFullstream=0;
TTreeSRedirector *pcstream=0;
//TTreeSRedirector *clstream=0;
TTreeSRedirector *cl700stream=0;


TFile* f, *fGids;
Bool_t dump=kFALSE;   //is current event dumped
Bool_t hasDumped=kFALSE; //has an event been dumped so far

void InitTree(Int_t filen, Bool_t readGids){

//    TFile f("/lustre/nyx/alice/users/slehner/LHC15oHLTLoop/245954/53/TPCsignal.root");
    AliSysInfo::AddStamp("Init_start",0,0,0);
    if(readGids) {
      clFullstream = new TTreeSRedirector("clusterDumpFullHLT.root","recreate");
    }
    fName ="TPCsignal.root";
    f= TFile::Open(fName.String(),"read");
    if(!readGids){
    pcstream = new TTreeSRedirector(TString::Format("eventDumpHLT.root"),"recreate");
    cl700stream = new TTreeSRedirector("clusterDump700HLT.root","recreate");
    }
    tree = (TTree*)f->Get("Clusters");
//    tree->SetDirectory(0);
    char cwd[100];
    getcwd(cwd, sizeof(cwd));
    fPath=TObjString(cwd);  
    AliSysInfo::AddStamp("Init_stop",0,0,0);
}

/// TODO  -  add sector Q counters = sum cluster->GetQ()
/// TODO  -  add TOT counter equivalent - sum of the GetQ() under condition -  cluster->GetMax()>thr
/// TODO  -  Higher level information - X,Y,Z   COG + SUM of Q  and RMS (in 3D)  per sector with Q>thr
/// TODO  -  Do the same calculation but in nsigma window - less sensitive to back-ground
///
void processTree(Int_t filen=0, Bool_t readGids=kFALSE){
  
  
  AliSysInfo::AddStamp("Process_start",0,0,0);
  InitTree(filen,readGids);
  if(readGids) {
    fGids = TFile::Open("selectedGIDs.root","READ");
    gidTree = (TTree*)fGids->Get("locCnt");
  }
  
  AliTPCclusterMI *cluster=0;
  ULong64_t gid=0;
  //
  Double_t sectorQSum;
  Double_t sectorQSum700;
  Double_t sectorQSum500;  
  Double_t sectorTOTSum700;
  Double_t sectorTOTSum500;
  Int_t sec;
  Double_t Q;
  Double_t Qmax;
  Float_t fx;
  Float_t fy;
  Float_t fz; 
  Float_t gx;
  Float_t gy;
  Float_t gz;
  
  TH3D* hQdis500= new TH3D("pos500","" , 100, -260, 260, 100,-260, 260, 100, -220, 220);      
  TH3D* hQdis700= new TH3D("pos700","" , 100, -260, 260, 100,-260, 260, 100, -220, 220);
  //Centre Of Gravity
  Float_t xCOG500;
  Float_t yCOG500;
  Float_t zCOG500; 
  Float_t xCOG700;
  Float_t yCOG700;
  Float_t zCOG700; 
  //StandardDeviation
  Float_t xSig500;
  Float_t ySig500;
  Float_t zSig500; 
  Float_t xSig700;
  Float_t ySig700;
  Float_t zSig700;
  //Correlation
  Float_t xyCorr500;
  Float_t xzCorr500;
  Float_t yzCorr500; 
  Float_t xyCorr700;
  Float_t xzCorr700;
  Float_t yzCorr700;
     
  int first=0;
  int last=0;
  
  Int_t entriesCluster  = tree->GetEntries();
  tree->SetBranchAddress("Cl.",&cluster);
  tree->SetBranchAddress("gid",&gid);
  tree->SetBranchAddress("gx",&gx);
  tree->SetBranchAddress("gy",&gy);
  tree->SetBranchAddress("gz",&gz);      
  
  for (Int_t iCl=0; iCl<entriesCluster; iCl++){
    tree->GetEntry(iCl);
//    cout<<tree->GetEntries()<<endl;
//    cout<<gid<<" "<<gx<<endl;
    sec=cluster->GetDetector();
    if(iCl==0){
      gidCurr=gid;
      sectorLast=sec;
    }
   
    if ((sec!=sectorLast || iCl==entriesCluster-1) && !readGids){         //if new sector or last entry, then get info on previous sector
      
      xCOG500=hQdis500->GetMean(1);
      yCOG500=hQdis500->GetMean(2);  
      zCOG500=hQdis500->GetMean(3);
      xCOG700=hQdis700->GetMean(1);
      yCOG700=hQdis700->GetMean(2);  
      zCOG700=hQdis700->GetMean(3);
      
      xSig500=hQdis500->GetStdDev(1);
      ySig500=hQdis500->GetStdDev(2);  
      zSig500=hQdis500->GetStdDev(3);
      xSig700=hQdis700->GetStdDev(1);
      ySig700=hQdis700->GetStdDev(2);  
      zSig700=hQdis700->GetStdDev(3);
      
      xyCorr500=hQdis500->GetCorrelationFactor(1,2);
      xzCorr500=hQdis500->GetCorrelationFactor(1,3);  
      yzCorr500=hQdis500->GetCorrelationFactor(2,3);
      xyCorr700=hQdis700->GetCorrelationFactor(1,2);
      xzCorr700=hQdis700->GetCorrelationFactor(1,3);  
      yzCorr700=hQdis700->GetCorrelationFactor(2,3); 
      
      hQdis700->Reset("ICESM");
      hQdis500->Reset("ICESM");
      sectorLast=sec;
      
if(!readGids){
      (*pcstream)<< "eventDump" <<
                  "gid="<< gidCurr<<           
                  "sector="<< sec<<           
                  "fPath.="<<&fPath<<
                  "sectorQ="<< sectorQSum<<
                  "sectorTOTSum500="<< sectorTOTSum500<<
                  "sectorTOTSum700="<< sectorTOTSum700<<        
                  "sectorQSum500="<< sectorQSum500<<
                  "sectorQSum700="<< sectorQSum700<<    
                  "xCOG500="<< xCOG500<<
                  "yCOG500="<< yCOG500<<
                  "zCOG500="<< zCOG500<<
                  "xCOG700="<< xCOG700<<
                  "yCOG700="<< yCOG700<<
                  "zCOG700="<< zCOG700<<    
                  "xSig500="<< xSig500<<
                  "ySig500="<< ySig500<<
                  "zSig500="<< zSig500<<
                  "xSig700="<< xSig700<<
                  "ySig700="<< ySig700<<
                  "zSig700="<< zSig700<<
                  "xyCorr500="<< xyCorr500<<
                  "xzCorr500="<< xzCorr500<<
                  "yzCorr500="<< yzCorr500<<
                  "xyCorr700="<< xyCorr700<<
                  "xzCorr700="<< xzCorr700<<
                  "yzCorr700="<< yzCorr700<<
                  "\n"; 
      
      }      
      
    }  
    
    if (gid!=gidCurr || iCl==entriesCluster-1){    //if new event or last entry, write info on previous event
         
      last=iCl-1;
      cout<<"current gid: "<<gidCurr<<" new: "<<gid<<" first: "<<first<< " last: "<<last<<endl;         
    
      if(readGids){
        if(gidTree->Draw("gid",TString::Format("gid==%llu && isSelected",gid),"goff")){     //select clusters of events that are in the selectedGIDs file and which were selected
          cout<<gidTree->Draw("gid",TString::Format("gid==%llu",gid),"goff")<<" for: "<<gid<<endl;
          dump=kTRUE;
        }
        else {
          cout<<"not selected"<<endl;
          dump=kFALSE;
        }
      }
    
      gidCurr=gid;
 
      first=iCl;  
    }
    if(readGids && !dump) continue;

    if(readGids && dump){ 
      (*clFullstream)<< "clusterDumpFull" <<  
          "gid="<< gid<<
          "gx="<< gx<<   
          "gy="<< gy<<
          "gz="<< gz<<   
          "Cluster.="<< cluster<<
          "\n";   
      hasDumped=kTRUE;
      continue;
    }
    
    fx=cluster->GetX();
    fy=cluster->GetY();
    fz=cluster->GetZ();   
    
    Q=cluster->GetQ();
    Qmax=cluster->GetMax();
    sectorQSum+=Q;
    if(Qmax>500){
      sectorQSum500+=Q;
      sectorTOTSum500+=1; 
      hQdis500->Fill(fx,fy,fz,Q);
      if(Qmax>700){
        sectorQSum700+=Q;
        sectorTOTSum700+=1;
        hQdis700->Fill(fx,fy,fz,Q);
        
        //Fill clusterDump700 tree
        (*cl700stream)<< "clusterDump700" <<  
                "gid="<< gidCurr<<
                "fPath.="<<&fPath<<    
                "gx="<< gx<<   
                "gy="<< gy<<
                "gz="<< gz<<   
                "Cluster.="<< cluster<<
                "\n";
      }
    }
    
  }
  delete pcstream;
//  delete clstream;
  delete cl700stream;
  delete clFullstream;
  
  if(readGids){
    remove("TPCsignal.root");
//    AliSysInfo::AddStamp("Process_stop",0,0,0);
//    TFile* fProf= TFile::Open("tpcHLTClusterLoopFullProf.root","RECREATE");
//    TTree * trProf = AliSysInfo::MakeTree("syswatch.log");
//    rename("syswatch.log","syswatchtpcfull.log");       //TODO: add filen to name!
//    trProf->Write();
  }
//  if(!hasDumped) remove("clusterDumpFullHLT.root")
//  else{
//  AliSysInfo::AddStamp("Process_stop",0,0,0);
//  TFile* fProf= TFile::Open(TString::Format("tpcHLTClusterLoopProf_%d.root",filen),"RECREATE");
//  TTree * trProf = AliSysInfo::MakeTree("syswatch.log");
//  rename("syswatch.log","syswatchtpc.log");       //TODO: add filen to name!
//  trProf->Write();
//  }
}


  /*
   *determines the cut limits for IROC and OROC sectors based on the percentile
gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
gSystem->AddIncludePath("-I\$NOTES/");
.L $NOTES/JIRA/ATO-432/code/tpcHLTClusterLoop.C
getLimits(10)
 */
void getLimits(Int_t perc){ 
  AliSysInfo::AddStamp("getLimits_start",0,0,0);
//  fName ="eventDumpHLT.root";
  fName ="merged_eventDump.root";    
  f= TFile::Open(fName.String());
  
  tree = (TTree*)f->Get("eventDump");
  TH1 *hHisto=0;
  Double_t norm;
  Double_t limits[2]={0,0};   //IROC and OROC cut limit for the specified percentile
  
  tree->Draw("Sum$((Iteration$<36)*(sectorTOTSum700.fElements/sectorQ.Sum())*(sectorTOTSum700.fElements>100))>>htemp(100,0,0)","Iteration$<36&&sectorTOTSum700.fElements>0","");
  hHisto = (TH1*)gPad->GetPrimitive("htemp");
  norm=hHisto->Integral();
  hHisto = hHisto->GetCumulative();
  hHisto->Scale(1./norm);
  hHisto->Draw();
  for(int bin=0; bin<hHisto->GetNbinsX(); bin++){
    if(hHisto->GetBinContent(bin)>perc*1./100.){
        limits[0]=hHisto->GetBinCenter(bin); 
        break;
    }
  }  
  
  tree->Draw("Sum$((Iteration$>=36)*(sectorTOTSum700.fElements/sectorQ.Sum())*(sectorTOTSum700.fElements>100))>>htemp(1000,0,0)","Iteration$>=36&&sectorTOTSum700.fElements>0","");
  hHisto = (TH1*)gPad->GetPrimitive("htemp");
  norm=hHisto->Integral();
  hHisto = hHisto->GetCumulative();
  hHisto->Scale(1./norm);
  hHisto->Draw();
  for(int bin=0; bin<hHisto->GetNbinsX(); bin++){
    if(hHisto->GetBinContent(bin)>perc*1./100.){
      limits[1]=hHisto->GetBinCenter(bin); 
      break;
    }
  } 
  cout<<perc<<" % percentile limit IROC: "<<limits[0]<<" OROC: "<<limits[1]<<endl;
  
  AliSysInfo::AddStamp("getLimits_stop",0,0,0);
  return;
}