/*
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
 .L $NOTES/JIRA/ATO-432/code/drawHighdEdx.C+
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  LoadSummaryTrees();
  LoadNDLocal();
  MakeFitResolLowMult();
  MakeFitResolHighMult();
  MakeFitNclFraction()
 */

/// TODO plots/fits
///  number of missing clusters - demonstrate poisson characteristic
/*

  treeMap->Draw("hisFraction0Dist.rms0:sqrt((1-hisFraction0Dist.mean0Fit)/hisCross0Dist.binMedian+0.012**2)","fitCut1","colz")
  htemp->Fit("pol1");
 */




void Init(){
  treeMap->SetAlias("rmsScaling","((dEdxMaxCenter**0.25)*sqrt(1+tglCenter**2)**0.25)");
  //
  treeMap->SetMarkerStyle(21);
  treeMap->SetMarkerSize(0.6);
  gPad->SetRightMargin(0.15);
  //
  treeMap->Draw("rmsMax0F:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==2&&multACenter<14000","colz");
  gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.15);
  AliDrawStyle::TPadApplyStyle("figTemplateTRD",gPad);

}


void DrawintrinsicResolution(){
  treeMap->Draw("rmsMax0F:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==5&&multACenter<14000","colz");
  gPad->SaveAs("rmsMax0F_multACenter_dEdxMaxCenter.png");
  treeMap->Draw("rmsMax0F*rmsScaling:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==5&&multACenter<14000","colz");
  gPad->SaveAs("rmsMax0FrmsScaling_multACenter_dEdxMaxCenter.png");
  //
  treeMap->Draw("rmsTot0F:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==5&&multACenter<14000","colz");
  gPad->SaveAs("rmsTot0F_multACenter_dEdxMaxCenter.png");
  treeMap->Draw("rmsTot0F*rmsScaling:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==5&&multACenter<14000","colz");
  gPad->SaveAs("rmsTot0FrmsScaling_multACenter_dEdxMaxCenter.png");

  treeMap->Draw("rmsMax1F:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==5&&multACenter<14000","colz");
  gPad->SaveAs("rmsMax1F_multACenter_dEdxMaxCenter.png");
  treeMap->Draw("rmsMax1F*rmsScaling:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==5&&multACenter<14000","colz");
  gPad->SaveAs("rmsMax1FrmsScaling_multACenter_dEdxMaxCenter.png");

  treeMap->Draw("rmsTot1F:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==5&&multACenter<14000","colz");
  gPad->SaveAs("rmsTot1F_multACenter_dEdxMaxCenter.png");
  treeMap->Draw("rmsTot1F*rmsScaling:multACenter:dEdxMaxCenter","isOK&&hisRatioMax01Dist.entries>100&&pileUpZBin==3&tglBin==5&&multACenter<14000","colz");
  gPad->SaveAs("rmsTot1FrmsScaling_multACenter_dEdxMaxCenter.png");

}

void drawTransfer(){
  //
  tree->Draw("ratioTotMax0:hisRatioTotMax0Dist.meanGFit:(pullTotMax)<6","nclCut&&chi2Cut","colz",10000);
  gPad->SaveAs("ratioTotMax0_RatioTotMax0Dist.meanGFit_Out6.png");
  tree->Draw("ratioTotMax1:hisRatioTotMax1Dist.meanGFit:(pullTotMax)<6","nclCut&&chi2Cut","colz",10000);
  gPad->SaveAs("ratioTotMax1_RatioTotMax1Dist.meanGFit_Out6.png");
}


void DrawResol(){
  //gPad->SetRightMargin(0.15);
  {
    treeMap->Draw("rmsMax0NormFit0:sdEdxMaxCenter:atglCenter", "fitCut0&&multABin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax0NormFit0_sdEdxMaxCenter_atglCenter.gif");
    treeMap->Draw("rmsMax0NormLF:rmsMax0NormFit0:atglCenter", "fitCut0&&multABin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax0NormLF_rmsMax0NormFit0_atglCenter.gif");
  }
  {
    treeMap->Draw("rmsMax1NormFit0:sdEdxMaxCenter:atglCenter", "fitCut0&&multABin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax1NormFit0_sdEdxMaxCenter_atglCenter.gif");
    treeMap->Draw("rmsMax1NormLF:rmsMax1NormFit0:atglCenter", "fitCut0&&multABin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax1NormLF_rmsMax1NormFit0_atglCenter.gif");
  }
    {
    treeMap->Draw("rmsMax2NormFit0:sdEdxMaxCenter:atglCenter", "fitCut0&&multABin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax2NormFit0_sdEdxMaxCenter_atglCenter.gif");
    treeMap->Draw("rmsMax2NormLF:rmsMax2NormFit0:atglCenter", "fitCut0&&multABin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax2NormLF_rmsMax2NormFit0_atglCenter.gif");
      treeMap->Draw("rmsMax2NormLF/rmsMax2NormFit0:atglCenter:sdEdxMaxCenter", "fitCut0&&multABin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
      gPad->SaveAs("rmsMax2NormLF_rmsMax2NormFit0_atglCenter_sdEdxMaxCenter.gif");
  }

  {
    treeMap->Draw("rmsMax0NormLF/rmsMax0NormFit0:sdEdxMaxCenter:multACenter", "fitCut1&&atglBin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax0NormLF_rmsMax0NormFit0_sdEdxMaxCenter_multACenter.gif");

    treeMap->Draw("rmsMax1NormLF/rmsMax1NormFit0:sdEdxMaxCenter:multACenter", "fitCut1&&atglBin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax1NormLF_rmsMax1NormFit0_sdEdxMaxCenter_multACenter.gif");
  }

  {
    treeMap->Draw("rmsMax1NormLF/rmsMax1NormFit1:sdEdxMaxCenter:multACenter","fitCut1&&atglBin==1","colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("rmsMax1NormLF_rmsMax1NormFit1_sdEdxMaxCenter_multACenter.gif");
  }

}

void DrawMissing() {
  {
    treeMap->Draw("missingFrac0FitM:mdEdx:multACenter1000", "fitCutNcl&&atglBin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("missingFrac0FitM_sdEdxMaxCenter_multACenter.gif");
  }
  {
    treeMap->Draw("missingFrac0FitM:multACenter1000:sdEdxMaxCenter", "fitCutNcl&&atglBin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("missingFrac0FitM_multACenter_sdEdxMaxCenter.gif");
  }
  {
    treeMap->Draw("missingFrac0FitM:missingFrac0LF:multACenter1000", "fitCutNcl&&atglBin==1&&abs(sdEdxMaxBin-3)<=3", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("missingFrac0FitM_missingFrac0LF_multACenter1000.gif");
  }
  //
  {
    treeMap->Draw("missingFrac1FitM:mdEdx:multACenter1000", "fitCutNcl&&atglBin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("missingFrac1FitM_sdEdxMaxCenter_multACenter.gif");
  }
  {
    treeMap->Draw("missingFrac1FitM:multACenter1000:sdEdxMaxCenter", "fitCutNcl&&atglBin==1", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("missingFrac1FitM_multACenter_sdEdxMaxCenter.gif");
  }
  {
    treeMap->Draw("missingFrac1FitM:missingFrac1LF:multACenter1000", "fitCutNcl&&atglBin==1&&abs(sdEdxMaxBin-3)<=3", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("missingFrac1FitM_missingFrac1LF_multACenter1000.gif");
  }

  {
    treeMap->Draw("missingFracRMS0:missingFracRMS0Exp", "fitCutNcl", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("missingFracRMS0_missingFracRMS0Exp.gif");
  }
  {
    treeMap->Draw("missingFracRMS1:missingFracRMS1Exp", "fitCutNcl", "colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("missingFracRMS1_missingFracRMS1Exp.gif");
  }
}

void DrawChi2(){

   {
     treeMap->Draw("hisChi2TPCDist.meanGFitM:multACenter1000:mdEdx","fitCutChi2&&atglBin==1","colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("hisChi2TPCDist.meanGFitM_multACenter1000_mdEdx.gif");
  }
     {
     treeMap->Draw("hisChi2TPCDist.meanGFitM:mdEdx:multACenter1000","fitCutChi2&&atglBin==1","colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("hisChi2TPCDist.meanGFitM_mdEdx_multACenter1000.gif");
  }
   {
     treeMap->Draw("hisChi2TPCDist.rmsGFitM/hisChi2TPCDist.meanGFitM:multACenter1000:sdEdxMaxCenter","fitCutChi2&&atglBin==1","colz");
    TStatToolkit::AdaptHistoMetadata(treeMap, 0, "colz");
    gPad->SaveAs("hisChi2TPCDist.rmsGFitM_hisChi2TPCDist.meanGFitM_multACenter1000_mdEdx.gif");
  }


}


void CheckInvariant(){

  treeMap->Draw("(exp(hisRatioMax01Dist.meanGFit)/exp(hisRatioTot01Dist.meanGFit))/(hisRatioTotMax1Dist.meanGFit/hisRatioTotMax0Dist.meanGFit)","isOK&&entries>200");
  treeMap->Draw("(exp(hisRatioMax02Dist.meanGFit)/exp(hisRatioTot02Dist.meanGFit))/(hisRatioTotMax2Dist.meanGFit/hisRatioTotMax0Dist.meanGFit)","isOK");
  treeMap->Draw("(exp(hisRatioMax12Dist.meanGFit)/exp(hisRatioTot12Dist.meanGFit))/(hisRatioTotMax2Dist.meanGFit/hisRatioTotMax1Dist.meanGFit)","isOK");
}


TString GetLatexFormula(){
  fitter= AliTMinuitToolkit::GetPredefinedFitter("fitterNclFraction");
  TVectorD *param  = fitter->GetParameters();
  TVectorD *errors = fitter->GetRMSEstimator();

  TString format="([2.2p0]#pm[2.2e0])";
}


void DrawHisto(){
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  hisdEdxMax->Projection(0,1)->Draw("colz");

}