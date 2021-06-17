/*
  gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
  ocdbpath="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/";
  run=244340;
  .L $AliPhysics_SRC/PWGPP/CalibMacros/CPass0/ConfigCalibTrain.C
  ConfigCalibTrain(run,ocdbpath)

  filename="/lustre/nyx/alice/DETdata/triggeredRaw/alice/cern.ch/user/p/pwg_pp/triggeredRaw/alice/data/2015/LHC15n/000244340/0010/rawSelected0.root"

  .L $NOTES/JIRA/ATO-432/code/TPCHighdEdxRawLoop.C+
  makeFileLoop(filename)
 */

#include "AliTPCAltroMapping.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliRawReaderRoot.h"
#include "AliTPCRawStreamV3.h"
#include "AliLog.h"
#include "AliTPCROC.h"
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <map>

std::map<unsigned long long, float> mapQ;
std::map<unsigned long long, float> mapTOT500;

float getQMax(unsigned long long gid);
float getTOT500Max(unsigned long long gid);

void makeFileLoop(const char *filename) {
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("rawDump.root","new");
  AliTPCAltroMapping **mapping = AliTPCcalibDB::Instance()->GetMapping();
  AliTPCParam *fParam = AliTPCcalibDB::Instance()->GetParameters();
  // TODO - use only pads from  active channel map   AliTPCcalibDB::Instance()->GetActiveChannelMap()
  //
  const Int_t kNIS = fParam->GetNInnerSector();	//number of inner sectors
  const Int_t kNOS = fParam->GetNOuterSector();	//number of outer sectors
  const Int_t kNS = kNIS + kNOS;	//all sectors
  AliRawReader *rawReader = new AliRawReaderRoot(filename);
  Int_t evtNr = 0;
  Float_t pos[2];
  TVectorF sectorQSum(72);
  TVectorF sectorTOTSum100(72);
  TVectorF sectorTOTSum300(72);  
  TVectorF sectorTOTSum500(72);
  TVectorF sectorTOTSum700(72);
  TVectorF sectorQSum100(72);
  TVectorF sectorQSum300(72);
  TVectorF sectorQSum500(72);
  TVectorF sectorQSum700(72);
  //
  TObjString fName(filename);
  while (rawReader->NextEvent()) {
    evtNr++;
    printf("Event\t%d\n", evtNr);
    AliTPCRawStreamV3 *rawStreamV3 = new AliTPCRawStreamV3(rawReader, (AliAltroMapping **) mapping);
    /// reset counters
    sectorQSum.Zero(); sectorTOTSum700.Zero(); sectorTOTSum500.Zero();
    //
    ULong64_t gid=rawReader->GetEventIdAsLong();
    UInt_t    timeStamp=rawReader->GetTimestamp();
    for (Int_t iSec = 0; iSec < kNS; iSec++) {
      // printf("Sec%d\n", iSec);
      Int_t nRows = 0; //number of rows in sector
      Int_t nDDLs = 0; //number of DDLs
      Int_t indexDDL = 0; //DDL index
      if (iSec < kNIS) {
        nRows = fParam->GetNRowLow();
        nDDLs = 2;
        indexDDL = iSec * 2;
      } else {
        nRows = fParam->GetNRowUp();
        nDDLs = 4;
        indexDDL = (iSec - kNIS) * 4 + kNIS * 2;
      }

      //
      // Load the raw data for corresponding DDLs
      //
      rawReader->Reset();
      rawReader->Select("TPC", indexDDL, indexDDL + nDDLs - 1);
      //
      while (rawStreamV3->NextDDL()) {
        if (AliLog::GetGlobalDebugLevel() > 2) rawStreamV3->PrintRCUTrailer();
        while (rawStreamV3->NextChannel()) {
          Int_t iSector = rawStreamV3->GetSector();                       //  current sector
          Int_t iRow = rawStreamV3->GetRow();                          //  current row
          Int_t iPad = rawStreamV3->GetPad();                          //  current pad
          AliTPCROC::Instance()->GetPositionGlobal(iSector,iRow,iPad,pos);
          while (rawStreamV3->NextBunch()) {
            UInt_t startTbin = rawStreamV3->GetStartTimeBin();
            Int_t bunchLength = rawStreamV3->GetBunchLength();
            const UShort_t *sig = rawStreamV3->GetSignals();
            //ProcessBunch(isector, iRow, iPad, bunchLength, startTbin, sig);
            for (Int_t iTimeBin = 0; iTimeBin < bunchLength; iTimeBin++) {
              Float_t signal = (Float_t) sig[iTimeBin];
              //printf("%d\t%d\t%d\t\t%d\n",isector,iRow,iPad,startTbin--,signal);
              Int_t timeBin=  startTbin+  iTimeBin;
              sectorQSum[iSector]+=signal;
              if (signal>100) {sectorTOTSum100[iSector]+=1; 
                               sectorQSum100[iSector]+=signal;}                            
              if (signal>300) {sectorTOTSum300[iSector]+=1; 
                               sectorQSum300[iSector]+=signal;}
              if (signal>500) {sectorTOTSum500[iSector]+=1; 
                               sectorQSum500[iSector]+=signal;}                            
              if (signal>700) {sectorTOTSum700[iSector]+=1; 
                               sectorQSum700[iSector]+=signal;}

              (*pcstream) << "rawDump" <<
                          "event=" << evtNr <<                    /// event within file
                          "gid="<<gid<<                           /// global unique ID
                          "sector=" << iSector <<                /// sector id
                          "row=" << iRow <<                      /// pad row
                          "pad=" << iPad <<                      /// pad number
                          "x=" <<pos[0] <<
                          "y=" <<pos[1] <<
                          "time=" << timeBin <<
                          "signal=" << signal <<
                          "\n";
            }
          }
        }
      }
    }
    (*pcstream) << "eventDump" <<
                "event=" << evtNr <<                            ///
                "gid="<< gid<<                                  ///
                "timeStamp="<<timeStamp<<
                "fName.="<<&fName<<
                "sectorQ.="<< &sectorQSum<<
                "sectorTOTSum100.="<< &sectorTOTSum100<<
                "sectorTOTSum300.="<< &sectorTOTSum300<<            
                "sectorTOTSum500.="<< &sectorTOTSum500<<
                "sectorTOTSum700.="<< &sectorTOTSum700<<
                "sectorQSum100.="<< &sectorQSum100<<
                "sectorQSum300.="<< &sectorQSum300<<            
                "sectorQSum500.="<< &sectorQSum500<<
                "sectorQSum700.="<< &sectorQSum700<<            
                "\n";
  }
  delete pcstream;
}

Double_t GetQMax(Double_t gid){
  return  mapQ[gid];
}

void makeDrawExample(Char_t* fName){
/*
	gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
	.L TPCHighdEdxRawLoop.C+
	makeDrawExample("244359.root")
*/ 	
	TFile* f = TFile::Open(fName);
	TCanvas *c1 = new TCanvas("c1","c1",1400,1000);

//	TTreeSRedirector *pcstream = new TTreeSRedirector(TString::Format("%s_corrTree.root",fName),"new");
        
	std::map<unsigned long long, float>::iterator it;
	
	unsigned long long gid;
        TVectorF* sectorQSum=0;
        TVectorF* sectorTOT500Sum=0;
	
	TTree* rawTree= (TTree*)f->Get("rawDump");
	TTree* evTree= (TTree*)f->Get("eventDump");
	
	auto nevent = evTree->GetEntries();
	cout<<nevent<<endl;
	
	evTree->SetBranchAddress("gid",&gid);
	evTree->SetBranchAddress("sectorTOTSum500.",&sectorTOT500Sum);
	evTree->SetBranchAddress("sectorQ.",&sectorQSum);
	
	for(int i =0;i<nevent;i++){
		evTree->GetEntry(i);
		if(sectorQSum->Max()>0)	mapQ.insert(std::pair<unsigned long long,float>(gid,sectorQSum->Max()));	
                if(sectorTOT500Sum->Max()>0) mapTOT500.insert(std::pair<unsigned long long,float>(gid,sectorTOT500Sum->Max()));
	}

//        rawTree->Draw("x:y:signal","getQMax(gid)>0","colz");
//        c1->Print("./QMax.png");			

        rawTree->Draw("getTOT500Max(gid):getQMax(gid)","getTOT500Max(gid)>0","prof");
        c1->Print("./corrTotQ_prof.png");


}

float getQMax(unsigned long long gid){
    return (mapQ.find(gid))->second;
}

float getTOT500Max(unsigned long long gid){
    return (mapTOT500.find(gid))->second;
}
