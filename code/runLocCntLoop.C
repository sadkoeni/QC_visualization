void runLocCntLoop(){
  
  gInterpreter->AddIncludePath("\$ALICE_ROOT/include/");
  gSystem->RedirectOutput(TString::Format("locCntLoop_%d.log")); 
  std::cout<<"\nRunning locCntLoop.C"<<std::endl;
  gROOT->ProcessLine(".L tpcLocCntLoop.C++");
  processTree("merged_clusterDump700HLT.root");
  
  gROOT->ProcessLine(".L tpcHLTClusterLoop.C++");
  processTree(0,kTRUE);   

  std::cout<<"\nDone"<<std::endl;

  return;
  
}