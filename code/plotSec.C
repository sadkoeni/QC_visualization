#include <stdio.h>
#include <iomanip>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>


void plotSec(){

TCanvas *c1 = new TCanvas("c1","c1",1400,1000);
TCanvas *c2 = new TCanvas("c2","c2",1400,1000);

gStyle->SetOptStat(0);

TChain* chainSel  = new TChain("locCnt");
TString selName("selectedGIDs.root");
if(chainSel->Add(selName,0)==0) return;

TChain* chainCl   = new TChain("clusterDumpFull");
TString clName("merged_clusterDumpFullHLT.root");
//~ TString clName("merged_clusterDumpFullHLT.root");
chainCl->Add(clName,0);

TFile* out;


ULong_t gid;
ULong_t gidTemp=-1;
Int_t sector;
Int_t sectorTemp=-1;
Int_t NBins200;
Bool_t isSel;
Int_t timlo;
Int_t timhi;
TVectorF* maxConsPR=0;

TBranch *brgid  = chainSel->GetBranch("gid");
TBranch *brse  = chainSel->GetBranch("sector");
TBranch *brsel  = chainSel->GetBranch("isSelected");
TBranch *brnbins  = chainSel->GetBranch("NBins700");
TBranch *brtimlo  = chainSel->GetBranch("timeBinLow");
TBranch *brtimhi  = chainSel->GetBranch("timeBinHigh");

TBranch *brmaxConsPR = 0;
chainSel->SetBranchAddress("maxPRcons.",&maxConsPR,&brmaxConsPR);

brgid->SetAddress(&gid);
brse->SetAddress(&sector);
brnbins->SetAddress(&NBins200);
brsel->SetAddress(&isSel);
brtimlo->SetAddress(&timlo);
brtimhi->SetAddress(&timhi);


int entries=chainSel->GetEntries();
for(int it=0; it<entries;it++){	
    chainSel->GetEntry(it);
	if((gidTemp==gid && sector==sectorTemp) || !isSel) continue;
	cout<<"new gid and sector found: "<<gid<<" "<<sector<<" (*maxConsPR)[4]: "<<(*maxConsPR)[4]<<", (*maxConsPR)[3]"<<(*maxConsPR)[3]<<endl;
	gidTemp=gid;
	sectorTemp=sector;
	out= TFile::Open(TString::Format("%llu_Histos_fMax200_wSelSec.root",gid),"RECREATE");
	cout<<it<<" out of "<<entries<<endl;

	c1->cd();
	gPad->SetLogz(0);	

	chainCl->Draw("gy:gx>>htemp(200,,,200,,)",TString::Format("abs(fQ)*(fMax>200) *(gid==%llu && fDetector==%d && fTimeBin>=%d && fTimeBin<%d) ",gid, sector,timlo,timhi),"colz");
	TH2 *myh = (TH2*)gPad->GetPrimitive("htemp"); 
	myh->SetName(TString::Format("gid_%llu_Sec_%d_YX",gid,sector));
	myh->SetTitle(TString::Format("YX fQ*(fMax>200) *(gid==%llu) && sector == %d, time: %d-%d",gid,sector,timlo,timhi));
	out->cd();
	myh->Write();
	c1->Print(TString::Format("Gid_%llu_YX_200_Sec_%d_time_%d-%d",gid,sector,timlo,timhi)+"_fMax.png");
			
	c1->cd();
	gPad->SetLogz();	
	chainCl->Draw("gy:gz>>htemp(200,,,200,,)",TString::Format("abs(fQ)*(fMax>200) *(gid==%llu && fDetector==%d && fTimeBin>=%d && fTimeBin<%d )",gid, sector,timlo,timhi),"colz");
	myh = (TH2*)gPad->GetPrimitive("htemp"); 
	myh->SetName(TString::Format("gid_%llu_Sec_%d_YZ",gid,sector));
	myh->SetTitle(TString::Format("YZ fQ*(fMax>200) *(gid==%llu)&& sector == %d, time: %d-%d",gid,sector,timlo,timhi));
	out->cd();
	myh->Write();
	c1->Print(TString::Format("Gid_%llu_YZ_200_Sec_%d_time_%d-%d",gid,sector,timlo,timhi)+"_fMax.png");

	c1->cd();
	gPad->SetLogz();
	chainCl->Draw("gy:gx>>htemp(4000,,,4000,,)",TString::Format("abs(fQ)*(fMax>200) *(gid==%llu)",gid),"colz");
	TH2 *myh = (TH2*)gPad->GetPrimitive("htemp"); 
	myh->SetName(TString::Format("gid_%llu_Sec_ALL_Time_ALL_YX",gid));
	myh->SetTitle(TString::Format("YX fQ*(fMax>200) *(gid==%llu), time: %d-%d",gid,timlo,timhi));
	out->cd();
	myh->Write();
	c1->Print(TString::Format("Gid_%llu_YX_200_Sec_ALL",gid)+"_fMax.png");
			
	c1->cd();
	gPad->SetLogz();	
	chainCl->Draw("gy:gz>>htemp(4000,,,4000,,)",TString::Format("abs(fQ)*(fMax>200) *(gid==%llu)",gid),"colz");
	myh = (TH2*)gPad->GetPrimitive("htemp"); 
	myh->SetName(TString::Format("gid_%llu_Sec_ALL_Time_ALL_YZ",gid));
	myh->SetTitle(TString::Format("YZ fQ*(fMax>200) *(gid==%llu), time: %d-%d",gid,timlo,timhi));
	out->cd();
	myh->Write();
	c1->Print(TString::Format("Gid_%llu_YZ_200_Sec_ALL",gid)+"_fMax.png");
	
	//for verification of pad row cluster finding also plot the sectors in fRow and fY coordinates for fMax>200		
	c1->cd();
	gPad->SetLogz(0);	
	chainCl->Draw("fRow:fY>>htemp(40,-40,40,96,0,96)",TString::Format("(fMax>200) *(gid==%llu&& fDetector==%d&& fTimeBin>=%d && fTimeBin<%d) ",gid, sector,timlo,timhi),"colz");
	myh = (TH2*)gPad->GetPrimitive("htemp"); 
	myh->SetName(TString::Format("gid_%llu_Sec_%d_YZ",gid));
	myh->SetTitle(TString::Format("fRow:fY, fMax>200, gid=%llu, sector=%d, time: %d-%d",gid,sector,timlo,timhi));
	out->cd();
	myh->Write();
	c1->Print(TString::Format("Gid_%llu_fRow_fY_200_Sec_%d_time_%d-%d",gid,sector,timlo,timhi)+"_fMax.png");
	c1->cd();
	gPad->SetLogz(0);	
	
	//~ chainCl->Draw("fRow:Time>>htemp(40,-40,40,96,0,96)",TString::Format("(fMax>200) *(gid==%llu&& fDetector==%d&& fTimeBin>=%d && fTimeBin<%d) ",gid, sector,timlo,timhi),"colz");
	//~ myh = (TH2*)gPad->GetPrimitive("htemp"); 
	//~ myh->SetName(TString::Format("gid_%llu_Sec_%d_YZ",gid));
	//~ myh->SetTitle(TString::Format("fRow:fY, fMax>200, gid=%llu, sector=%d, time: %d-%d",gid,sector,timlo,timhi));
	//~ out->cd();
	//~ myh->Write();
	//~ c1->Print(TString::Format("Gid_%llu_fRow_fY_200_Sec_%d_time_%d-%d",gid,sector,timlo,timhi)+"_fMax.png");

	//~ chainCl->Draw("fRow:fPad>>htemp(140,,,96,,)",TString::Format("fQ*(fMax>200) *(gid==%llu && fDetector==%d)",gid, sector),"colz");
	//~ myh = (TH2*)gPad->GetPrimitive("htemp"); 
	//~ myh->SetName(TString::Format("Gid_%llu_Sec_%d_RowPad_200",gid,sector));
	//~ myh->SetTitle(TString::Format("Row:Pad fQ*(fMax>200) *(gid==%llu) && sector == %d",gid,sector));
	//~ out->cd();
	//~ myh->Write();
	//~ c1->Print(TString::Format("Gid_%llu_Sec_%d_RowPad_200",gid,sector)+"_fMax.png");
			
	//~ c1->cd();
	//~ gPad->SetLogz();	
	//~ chainCl->Draw("fRow:fTimeBin>>htemp(1000,,,96,,)",TString::Format("fQ*(fMax>200) *(gid==%llu && fDetector==%d)",gid, sector),"colz");
	//~ myh = (TH2*)gPad->GetPrimitive("htemp"); 
	//~ myh->SetName(TString::Format("Gid_%llu_Sec_%d_RowTime_200",gid,sector));
	//~ myh->SetTitle(TString::Format("Row:Time fQ*(fMax>200) *(gid==%llu && sector==%d)",gid,sector));
	//~ out->cd();
	//~ myh->Write();
	//~ c1->Print(TString::Format("Gid_%llu_Sec_%d_RowTime_200",gid,sector)+"_fMax.png");
	
	//~ chainCl->Draw("gz:gy:gx>>htemp3(200,,,200,,,200,,)",TString::Format("fQ*(fMax>200) *(gid==%llu)",gid),"colz");
	//~ TH3 *myh3 = (TH3*)gPad->GetPrimitive("htemp3"); 
	//~ myh3->SetName(TString::Format("gid_%llu_YXZ",gid));
	//~ myh3->SetTitle(TString::Format("ZYX fQ*(fMax>200) *(gid==%llu)",gid));
	//~ out->cd();
	//~ myh3->Write();
	}

}	
