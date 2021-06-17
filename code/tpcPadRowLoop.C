///
/// get pad row statistics
///
/*
gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
gSystem->AddIncludePath("-I\$NOTES/");
.L $NOTES/JIRA/ATO-432/code/tpcPadRowLoop.C
processTree("clusterDumpHLT.root")
 */

#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include <TH3D.h>
#include <stdio.h>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>

ULong64_t gidCurr=-1;
Int_t sectorLast=-1;
TTree * tree =0;
TTreeSRedirector *prstream=0;
TObjString fPath=0;
TFile* f;

void InitTree(TObjString fName){
    AliSysInfo::AddStamp("Init_start",0,0,0);
    f= TFile::Open(fName.String());
    tree = (TTree*)f->Get("clusterDump");
    prstream = new TTreeSRedirector("padRowDumpHLT.root","recreate");
    char cwd[100];
    getcwd(cwd, sizeof(cwd));
    fPath=TObjString(cwd);  
}

void processTree(TObjString fClFile){
  InitTree(fClFile);

  AliTPCclusterMI *cluster=0;
  ULong64_t gid=0;

  Int_t sec;
  Double_t Q;
  Double_t Qmax;
  
  TVectorF nPadRows700(72);
  TVectorF nPadRows500(72);
  
  TVectorF firstPadRow700(72);
  TVectorF firstPadRow500(72);
  
  TVectorF maxConsecPadRow700(72);    //max number of consecutive pad rows with high signal
  TVectorF maxConsecPadRow500(72);
    
  int first=0;
  int last=0;
  
  Int_t entriesCluster  = tree->GetEntries();
  tree->SetBranchAddress("Cluster.",&cluster);
  tree->SetBranchAddress("gid",&gid);
    
  for (Int_t iCl=0; iCl<entriesCluster; iCl++){
    tree->GetEntry(iCl);
    sec=cluster->GetDetector();
    if(iCl==0) gidCurr=gid;

    if(sec!=sectorLast || iCl==entriesCluster-1){         //if new sector or last entry
     
      if(gid!=gidCurr || iCl==entriesCluster-1){    //if new event or last entry, then write info about previous event
        cout<<"current gid: "<<gidCurr<<" next: "<<gid<<endl;      

        (*prstream)<< "padRowDump" <<
                    "fPath.="<<&fPath<<                
                    "gid="<< gidCurr<<           
//                    "nPadRows500.="<< &nPadRows500<<
//                    "nPadRows700.="<< &nPadRows700<< 
                    "firstPadRow700.="<< &firstPadRow700<<
//                    "firstPadRow500.="<< &firstPadRow500<< 
                    "maxConsecPadRow700.="<< &maxConsecPadRow700<<
//                    "maxConsecPadRow500.="<< &maxConsecPadRow500<<
                    "maxConsecPadRow700.="<< &maxConsecPadRow700<<
//                    "maxConsecPadRow500.="<< &maxConsecPadRow500<<
                    "\n"; 

        gidCurr=gid;
        for(int i =0; i<72; i++) {      //reset sector statistics
          nPadRows500[i]=0;
          nPadRows700[i]=0;
          firstPadRow700[i] = -1;
          firstPadRow500[i] = -1;
          maxConsecPadRow700[i] = 0;
          maxConsecPadRow500[i] = 0;
        }  
      }
      if(iCl==entriesCluster-1) break;
      //Get sector info
      Bool_t found500=kFALSE;
      Bool_t found700=kFALSE;
      
      Bool_t prev500=kFALSE;
      Bool_t prev700=kFALSE;

      Int_t tempmaxConsec500=0;
      Int_t tempmaxConsec700=0;

//Get pad row statistics      
      for(int pr=0; pr<96; pr++){  
  
        cout<<pr<<" "<<sec<<" "<<gid<<endl;
        //get pad row info for signal>700
        if(tree->Draw("Cluster.fQ",TString::Format("gid==%ld && Cluster.fRow==%d && Cluster.fMax>700 && Cluster.fDetector==%d ",gid,pr,sec),"goff")){
          tree->GetEntry(iCl);  //read tree again after TTree::Draw() was called          
          nPadRows700[sec]+=1;
          if(!found700) {
            firstPadRow700[sec]=pr;
            found700=kTRUE;
          }
          maxConsecPadRow700[sec]+=1;
        
        //get pad row info for signal>700
//        if(tree->Draw("Cluster.fQ",TString::Format("gid==%ld && Cluster.fRow==%d && Cluster.fMax>700 && Cluster.fDetector==%d ",gid,pr,sec),"goff")){
//          tree->GetEntry(iCl);  //read tree again after TTree::Draw() was called 
//          nPadRows700[sec]+=1;
//          if(!found700) {
//            firstPadRow700[sec]=pr;
//            found700=kTRUE;
//          }
//          maxConsecPadRow700[sec]+=1;
//        }  
//        else{         //if no high signal in this pad row
//          tree->GetEntry(iCl);  //read tree again after TTree::Draw() was called
//          if(maxConsecPadRow700[sec]>tempmaxConsec700)  tempmaxConsec700=maxConsecPadRow700[sec];
//          maxConsecPadRow700[sec]=0;
//        }
        }  
        else{         //if no high signal in this pad row
          tree->GetEntry(iCl);  //read tree again after TTree::Draw() was called
          if(maxConsecPadRow700[sec]>tempmaxConsec700)  tempmaxConsec700=maxConsecPadRow700[sec];
          maxConsecPadRow700[sec]=0;
        }
    }
    maxConsecPadRow500[sec]=tempmaxConsec500;
    maxConsecPadRow700[sec]=tempmaxConsec700;
      
      sectorLast=cluster->GetDetector();  
    }
    
   
  }
  delete prstream;
  AliSysInfo::AddStamp("Process_stop",0,0,0);
  TFile* fProf= TFile::Open("tpcLocCntLoopProf.root","RECREATE");
  TTree * trProf = AliSysInfo::MakeTree("syswatch.log");
  trProf->Write();
}
