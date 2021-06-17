/*
Local Test:
gSystem->AddIncludePath("-I\$NOTES/");
.L recHLT.C
recHLT("rawSelected0.root")

Grid Test:
TGrid::Connect("alien://",0,0,"t");
.L recHLT.C
recHLT()
 */

//gInterpreter->AddIncludePath("\$ALICE_ROOT/include/")

#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
//#include "AliTPCclusterMI.h"
//#include "AliCDBManager.h"
//#include "AliReconstruction.h"
//#include "AliTPCReconstructor.h"
//#include "AliTPCtracker.h"
#include "TROOT.h"
#include "TVectorF.h"
#include <TH3D.h>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>
//#include "TTreeStream.h"
//#include "AliSysInfo.h"
#include <TH2D.h>
#include <iostream>
//#include "TMath.h"
//#include "TString.h"


//void recHLT(const char *filename="raw.root", int nevents=-1, const char* ocdb="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB/", TString additionalRecOptions=TString("")) {
void recHLT(const char *filename="/alien/alice/data/2015/LHC15o/000245064/raw/15000245064037.914.root", int filen=0, int nevents=-1, const char* ocdb="raw://", TString additionalRecOptions=TString("")) {
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for 2010 RAW data
  //
  /////////////////////////////////////////////////////////////////////////////////////////
  // Set the CDB storage location
//  gSystem->RedirectOutput(TString::Format("recHLT_%d.log",filen));
  AliSysInfo::AddStamp("Init_start",0,0,0);
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb);
  if(!(man->InitFromSnapshot("../../../OCDB.root", kFALSE))) printf("Failed to Init From Snapshot!!");
  TObjString* os;
  AliReconstruction rec;
  rec.SetCDBSnapshotMode("../../../OCDB.root");
  // Set reconstruction flags (skip detectors here if neded with -<detector name>
  rec.SetRunLocalReconstruction("ALL"); /// Check - we need only rec points
  rec.SetRunTracking("");               /// Check - we need only rec points
  rec.SetFillESD("");
  rec.SetRunVertexFinder(kFALSE);
  // QA options
  //  rec.SetRunQA("Global:ESDs") ;
  //  rec.SetRunQA(":") ;
  //  rec.SetRunQA("ALL:ALL") ;
  rec.SetRunQA("") ; // disable all

  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kFALSE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);

  // Set additional reconstruction options.
  // They are in the form "Detector:value;Detector2:value" in a single string.
  // For instance: additionalRecOptions="TPC:useHLTorRAW"
  {
    TIter nexttok( additionalRecOptions.Tokenize(";") );
    while (( os = (TObjString *)nexttok() )) {
      TString detOpt = os->String();
      Ssiz_t idx = detOpt.Index(":");
      if (idx < 0) continue;
      TString detOptKey = detOpt(0,idx);
      TString detOptVal = detOpt(idx+1,999);
      if (detOptKey.IsNull() || detOptVal.IsNull()) continue;
      printf("Setting additional reconstruction option: %s --> %s\n", detOptKey.Data(),
                                                                      detOptVal.Data());
      rec.SetOption(detOptKey.Data(), detOptVal.Data());
    }
  }

  // Set protection against too many events in a chunk (should not happen)
  if (nevents>0) rec.SetEventRange(0,nevents);

  // Upload CDB entries from the snapshot (local root file) if snapshot exist
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    rec.SetCDBSnapshotMode("OCDB.root");
  }

  // Specific AD storage, see https://alice.its.cern.ch/jira/browse/ALIROOT-6056
  //  rec.SetSpecificStorage("AD/Calib/TimeSlewing", "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal");

  //load any specific local storages
  if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {
    gROOT->LoadMacro("localOCDBaccessConfig.C");
    localOCDBaccessConfig();
  }

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);
  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);
  // Delete recpoints
  rec.SetDeleteRecPoints("TPC TRD");
  // Set 100% of friends
  // rec.SetFractionFriends(2.0);
  AliLog::Flush();
  //

  AliTPCReconstructor::SetStreamLevel(AliTPCtracker::kStreamClDumpLocal);  // option to DUMP HLT clusters


  rec.Run();
  
  AliSysInfo::AddStamp("Process_stop",0,0,0);
  TFile* fProf= TFile::Open(TString::Format("recHLTProf_%d.root",filen),"RECREATE");
  TTree * trProf = AliSysInfo::MakeTree("syswatch.log");
  rename("syswatch.log","syswatchrec.log");   //TODO: add filen to name!
  trProf->Write();
}
