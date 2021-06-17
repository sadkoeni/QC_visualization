/*
  gSystem->AddIncludePath("-I\$ALICE_ROOT/include/");
  .L $NOTES/JIRA/ATO-432/code/monoAnalysis.C+g
  loadSeeds();
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
 */

void drawTrackNclPlotQA(){
  //
  treeSeeds1->Draw("fNClustersInner:fNClustersOuter");   // Not filled yet - test agreement between counter, improving dEdx correction better match expected
  // observation - Many time Outer loop lost, Inner lost in fewer cases
  //
  treeSeeds1->Draw("fNClustersInner<fNClustersOuter*0.5");  // for inner in 0.79 % of cases
  treeSeeds1->Draw("fNClustersOuter<fNClustersInner*0.5");
  // for outer happen in 25 % of cases = bug fixed 18.03 ->  0.0060469557
  //
  treeSeeds1->Draw("fNClustersOuter==0:clM.fX","fNClustersInner>40","prof"); // for outer happen in 25 % of cases
 treeSeeds1->Draw("fNClustersOuter<fNClustersInner*0.5:clM.fX","fNClustersInner>40","prof"); // for outer happen in 25 % of cases
  // chi2 not filled
  // parameters along trajectory not filled
}

void makeDefaultPlots(){
  treeSeeds1->Draw("seed.fRmsYMedian0");   // mean ~ 0,18 cm
  treeSeeds1->Draw("seed.fRmsZMedian0");  // mean ~ 0.27 cm
  // energy loss correction vissible <100 MeV
  treeSeeds1->Draw("log(seed.fParamInner.P()/seed.fParamOuter.P()):sqrt(seed.fParamInner.P()*seed.fParamOuter.P())>>his(50,0,0.5)","sqrt(seed.fParamInner.P()*seed.fParamOuter.P())<1","prof")
  gPad->SaveAs("momentumLoss.png");
  //
  treeSeeds1->Draw("pullP2:seed.fParamArrayIn.fData[].fX>>hisPullP2(20,80,250,100,-6,6)","nclMin>100&&trackletOK","colz");
  ((TH2F*)treeSeeds1->GetHistogram())->FitSlicesY();
  hisPullP2_2->SetMarkerStyle(21)
  hisPullP2_2->Draw("same");
  gPad->SaveAs("pullP2vsX.png");
  //  ~ rms ` 1 as it should be
  treeSeeds1->Draw("pullP3:seed.fParamArrayIn.fData[].fX>>hisPullP3(20,80,250,100,-6,6)","nclMin>100&&trackletOK","colz");
  ((TH2F*)treeSeeds1->GetHistogram())->FitSlicesY();
  hisPullP3_2->SetMarkerStyle(21)
  hisPullP3_2->Draw("same");
  gPad->SaveAs("pullP3vsX.png");
  //  ~ rms ` 1 as it should be


}