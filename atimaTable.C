///
TTree  tree =0;
void initTree(){
  Double_t mip=1.52; ///MeV/g cm−2
  Double_t density=1.661; /// #mg/cm^2
  Double_t mass0=TDatabasePDG::Instance()->GetParticle("proton")->Mass()*1000;   /// mass in MeV

  tree.ReadFile("elossAr.txt");
  tree.SetAlias("mass",Form("A*%f",mass0));
  tree.SetAlias("p","sqrt((E+mass)**2-mass**2)");
  tree.SetAlias("bg","p/mass");
  tree.SetAlias("dEdxAleph",Form("Z*Z*AliExternalTrackParam::BetheBlochAleph(bg)*%f*0.001",mip));
}

voiad MakeGraphs(){
   tree->Draw("dEdx/dEdxAleph:bg:Z","stopping>1&&bg>0.05","colz");
}