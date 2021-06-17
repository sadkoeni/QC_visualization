


void highdEdxFilter(const char *filename, int filen){
  
  gInterpreter->AddIncludePath("\$ALICE_ROOT/include/");
  
  gSystem->RedirectOutput(TString::Format("recHLT_%d.log",filen)); 
  std::cout<<"\nRunning recHLT on "<<filename<<std::endl;
  gROOT->LoadMacro("recHLT.C");
  recHLT(filename,filen);

  gSystem->RedirectOutput(TString::Format("tpcHLTClusterLoop_%d.log",filen));  
  std::cout<<"\nRunning tpcHLTClusterLoop on output of recHLT("<<filename<<")"<<std::endl;
  gROOT->ProcessLine(".L tpcHLTClusterLoop.C++");
  gROOT->ProcessLine(TString::Format("processTree(%d)",filen));

  gSystem->RedirectOutput(TString::Format("tpcPadRowClustering_%d.log",filen)); 
  std::cout<<"\nRunning tpcPadRowClustering.C"<<std::endl;
  gROOT->ProcessLine(".L tpcPadRowClustering.C");
  gROOT->ProcessLine("processTree(\"clusterDump700HLT.root\")");
  
  gSystem->RedirectOutput(TString::Format("tpcHLTFullClusterLoop_%d.log",filen)); 
  std::cout<<"\nRunning tpcHLTClusterLoop for full events"<<std::endl;
  gROOT->ProcessLine(".L tpcHLTClusterLoop.C++");
  gROOT->ProcessLine("processTree(0,kTRUE)");  

  std::cout<<"\nDone"<<std::endl;

  return;
  
}