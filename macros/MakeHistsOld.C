#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TRandom3.h>

Bool_t Accept( Double_t colDist, Double_t colRad, TH2* beamFunc );
TRandom3 *randoms = 0;

void MakeHists(){
  
  //TString amoDir = "/scratch/simong/testBad3.root";
  TString amoDir = "/w/work5/home/simong/CBremGen/Adapt/Copper3/Repeat_0/Events_0_*";
  //TString polDir = "/w/work5/home/simong/CBremGen/Adapt/DiamondBadTest/Repeat_0/Events_0_*";
  TString polDir = "/w/work5/home/simong/CBremGen/Adapt/Diamond5/Repeat_103/Events_103_*.root";
  //TString polDir = "/w/work5/home/simong/CBremGen/Adapt/Diamond6/Repeat_0/Events_0_*.root";
  //TString polDir = "/scratch/simong/Events_102_*";

  TChain* amoChain = new TChain("phono");

  amoChain->Add(amoDir);

  TFile* outFile = new TFile("outPlotsTemp2.root","RECREATE");

  //Define Histograms
  Int_t    ebins = 400;
  Double_t emin  = 0;
  Double_t emax  = 1.6;

  TH1F* amoHist    = new TH1F("amohist",   "amohist",   ebins,emin,emax);
  TH1F* amoHistPol = new TH1F("amohistpol","amohistpol",ebins,emin,emax);
  TH1F* amoHistPro0 = new TH1F("amohistpro0","amohistpro0",ebins,emin,emax);
  TH1F* amoHistPro1 = new TH1F("amohistpro1","amohistpro1",ebins,emin,emax);
  TH1F* amoHistPro0Pol = new TH1F("amohistpro0pol","amohistpro0pol",ebins,emin,emax);
  TH1F* amoHistPro1Pol = new TH1F("amohistpro1pol","amohistpro1pol",ebins,emin,emax);

  TH2F* amoCircHist2D  = new TH2F("amoCircHist2D", "amoCircHist2D", ebins,emin,emax,400,0,1);
  TH1F* amoCircHistPol = new TH1F("amoCirchistpol","amoCirchistpol",ebins,emin,emax);

  TH1F* polHist    = new TH1F("polhist",   "polhist",   ebins,emin,emax);
  TH1F* polHistPol = new TH1F("polhistpol","polhistpol",ebins,emin,emax);
  TH1F* polHistPro0 = new TH1F("polhistpro0","polhistpro0",ebins,emin,emax);
  TH1F* polHistPro1 = new TH1F("polhistpro1","polhistpro1",ebins,emin,emax);
  TH1F* polHistPro0Pol = new TH1F("polhistpro0pol","polhistpro0pol",ebins,emin,emax);
  TH1F* polHistPro1Pol = new TH1F("polhistpro1pol","polhistpro1pol",ebins,emin,emax);

  TH2F* circHist2D  = new TH2F("circhist2D", "circhist2D", ebins,emin,emax,400,0,1);
  TH1F* circHistPol = new TH1F("circhistpol","circhistpol",ebins,emin,emax);

//   TH1F* polHistOv    = new TH1F("polhistOv",   "polhistOv",   ebins,emin,emax);
//   TH1F* polHistOvPol = new TH1F("polhistOvpol","polhistOvpol",ebins,emin,emax);

//   TH1F* polHistOvA    = new TH1F("polhistOvA",   "polhistOvA",   ebins,emin,emax);
//   TH1F* polHistOvAPol = new TH1F("polhistOvApol","polhistOvApol",ebins,emin,emax);

//   TH1F* polHistOff    = new TH1F("polhistOff",   "polhistOff",   ebins,emin,emax);
//   TH1F* polHistOffPol = new TH1F("polhistOffpol","polhistOffpol",ebins,emin,emax);

//   TH1F* polHistOffA    = new TH1F("polhistOffA",   "polhistOffA",   ebins,emin,emax);
//   TH1F* polHistOffAPol = new TH1F("polhistOffApol","polhistOffApol",ebins,emin,emax);

  TTreeReader     fReader;  //!the tree reader

  auto beamFunc = new TF2("bigaus","bigaus");
  
  Double_t xBeamRad = 0.0005; //m
  Double_t yBeamRad = 0.0005; //m

  beamfunc->SetParameters(1,0,xBeamRad,0,yBeamRad);
  //Could initiate 

  Double_t thetaLim = 0.0001;
  Int_t analN = 1000000;

  //  Double_t thetaLim = 0.001;
  //Int_t analN = 100000000;

  //Double_t thetaLim = 0.00006;
  //Int_t analN = 10000000000;
  TTreeReaderValue<Double_t>       diffXS    = {fReader, "diffXS"};
  TTreeReaderValue<Double_t>       wgt       = {fReader, "wgt"};
  TTreeReaderValue<Double_t>       polar     = {fReader, "polar"};
  TTreeReaderValue<Double_t>       circ1     = {fReader, "circular"};
  TTreeReaderArray<Double_t>       k         = {fReader, "k"};

  TTreeReaderArray<Double_t>       qphoton   = {fReader, "qphoton"};
  TTreeReaderValue<Int_t>          process   = {fReader, "process"};
    
  fReader.SetTree(amoChain);
  cout << amoChain->GetEntries() << " " << analN << endl;
  fReader.SetEntriesRange(0,analN);
  Int_t j=0; 
  while(fReader.Next()){
//     cout << j++ << endl;
//     cout << "HI" << endl;

    TVector3 vec(k[0],k[1],k[2]);
    Double_t qpho = qphoton[2];
    //if(qpho>-0.000001){
    if(vec.Theta()<thetaLim){
      Double_t e = k[2];
      Double_t weight = (*wgt);
      Double_t Xsec   = (*diffXS);
      Double_t pol    = (*polar);
      Double_t circ   = (*circ1);
      Int_t    pro    = (*process);
      
      if(pro){
	amoHistPro1->Fill(e,weight*Xsec);
	amoHistPro1Pol->Fill(e,weight*Xsec*((pol*2)-1));
      }
      else {
	amoHistPro0->Fill(e,weight*Xsec);
	amoHistPro0Pol->Fill(e,weight*Xsec*((pol*2)-1));
      }

      amoHist->Fill(e,weight*Xsec);

      amoHistPol->Fill(e,weight*Xsec*((pol*2)-1));
      amoCircHist2D->Fill(e,circ,weight*Xsec);
      amoCircHistPol->Fill(e,weight*Xsec*circ);
    }
  }

  TChain* polChain = new TChain("phono");//,"phono2");
  //polChain->Reset();
  polChain->Add(polDir);
  
  TTreeReader     fReader2;  //!the tree reader

  TTreeReaderValue<Double_t>       diffXS2    = {fReader2, "diffXS"};
  TTreeReaderValue<Double_t>       wgt2       = {fReader2, "wgt"};
  TTreeReaderValue<Double_t>       polar2     = {fReader2, "polar"};
  TTreeReaderValue<Double_t>       circ2      = {fReader2, "circular"};
  TTreeReaderArray<Double_t>       k2         = {fReader2, "k"};

  TTreeReaderArray<Double_t>       qphoton2   = {fReader2, "qphoton"};
  TTreeReaderValue<Int_t>          process2   = {fReader2, "process"};

  //fReader2.Restart();
  //Int_t analN = 100000000;
  cout << polChain->GetEntries() << " " << analN << endl;
  fReader2.SetTree(polChain);
  fReader2.SetEntriesRange(0,analN);
  Int_t i = 0;
  while(fReader2.Next()){
//     cout << qphoton2[2] << endl;
//     cout << i++ << endl;
//     cout << "HI" << endl;
    Double_t qpho = qphoton2[2];

    TVector3 vec(k2[0],k2[1],k2[2]);
    //cout << vec.Theta() << endl;

    Double_t e = k2[2];
    Double_t weight = (*wgt2);
    Double_t Xsec   = (*diffXS2);
    Double_t pol    = (*polar2);
    Double_t circ   = (*circ2);
    Int_t    pro    = (*process2);
    //if(pro) continue;

    //if(qpho>-0.000001){
    //if(qpho>-0.0000001){
    if(vec.Theta()<thetaLim){



      //if(vec.Theta()<0.0005){
      polHist->Fill(e,weight*Xsec);
      polHistPol->Fill(e,weight*Xsec*((pol*2)-1));
      circHist2D->Fill(e,circ,weight*Xsec);
      circHistPol->Fill(e,weight*Xsec*circ);

      if(pro){
	polHistPro1->Fill(e,weight*Xsec);
	polHistPro1Pol->Fill(e,weight*Xsec*((pol*2)-1));
      }
      else{
	polHistPro0->Fill(e,weight*Xsec);
	polHistPro0Pol->Fill(e,weight*Xsec*((pol*2)-1));
      }
    }

    // vec[0]/=2;
    
//     if(vec.Theta()<0.0005){
//       polHistOv->Fill(e,weight*Xsec);
//       polHistOvPol->Fill(e,weight*Xsec*((pol*2)-1));
//     }

//     vec[0]*=2;
//     vec[1]/=2;
    
//     if(vec.Theta()<0.0005){
//       polHistOvA->Fill(e,weight*Xsec);
//       polHistOvAPol->Fill(e,weight*Xsec*((pol*2)-1));
//     }

//     vec[1]*=2;

//     vec.RotateX(0.001);
    
//     if(vec.Theta()<0.001){
//       polHistOff->Fill(e,weight*Xsec);
//       polHistOffPol->Fill(e,weight*Xsec*((pol*2)-1));
//     }

//     vec.RotateX(-0.001);
//     vec.RotateY(0.001);
    
//     if(vec.Theta()<0.001){
//       polHistOffA->Fill(e,weight*Xsec);
//       polHistOffAPol->Fill(e,weight*Xsec*((pol*2)-1));
//     }

  }
  TString eMaxS = "1.6";
  TF1* circFunc = new TF1("circ","(4*x/"+eMaxS+"-x*x/("+eMaxS+"*"+eMaxS+"))/(4-4*x/"+eMaxS+"+3*x*x/("+eMaxS+"*"+eMaxS+"))",0,1.6);

  amoHistPol->Divide(amoHist);
  amoCircHistPol->Divide(amoHist);
  polHistPol->Divide(polHist);
  polHistPro0Pol->Divide(polHistPro0);
  polHistPro1Pol->Divide(polHistPro1);
  circHistPol->Divide(polHist);

//   polHistOvPol->Divide(polHistOv);
//   polHistOvAPol->Divide(polHistOvA);

//   polHistOffPol->Divide(polHistOff);
//   polHistOffAPol->Divide(polHistOffA);
  amoHist->Write();
  amoHistPol->Write();
  amoHistPro0->Write();
  amoHistPro1->Write();
  amoHistPro0Pol->Write();
  amoHistPro1Pol->Write();
  amoCircHist2D->Write();
  amoCircHistPol->Write();
  polHist->Write();
  polHistPol->Write();
  polHistPro0->Write();
  polHistPro1->Write();
  polHistPro0Pol->Write();
  polHistPro1Pol->Write();
  circHist2D->Write();
  circHistPol->Write();
  circFunc->Write();

  circFunc->SetNpx(1000);
  circFunc->SetLineWidth(1);
  circFunc->SetLineColor(2);

  TCanvas* can = new TCanvas("can","can",1800,1800);

  can->SetLogz();
  can->Divide(2,2);

  can->cd(1);
  gPad->SetLogz();
  amoCircHist2D->Draw("colz");
  amoCircHistPol->Draw("same hist");
  circFunc->Draw("same hist");

  can->cd(2);
  gPad->SetLogz();
  amoCircHist2D->Draw("colz");
  circHistPol->Draw("same hist");
  circFunc->Draw("same hist");

  can->cd(3);

  TH1* tempHist = (TH1*)circHistPol->Clone("divided");

  tempHist->Divide(circFunc);
  tempHist->SetMaximum(1.1);
  tempHist->SetMinimum(0.9);
  tempHist->Draw("hist");
  
  can->cd(4);

  TH1* tempHist2 = (TH1*)amoCircHistPol->Clone("divided2");

  tempHist2->Divide(circFunc);
  tempHist2->SetMaximum(1.1);
  tempHist2->SetMinimum(0.9);
  tempHist2->SetLineColor(2);
  tempHist2->Draw("hist same");

  

  can->SaveAs("CircularPol.pdf");
//   polHistOff->Write();
//   polHistOffPol->Write();
//   polHistOffA->Write();
//   polHistOffAPol->Write();

//   polHistOv->Write();
//   polHistOvPol->Write();
//   polHistOvA->Write();
//   polHistOvAPol->Write();


  outFile->Close();  

}

//--------------------------------------------------
// Acceptance calculation
//--------------------------------------------------

Bool_t Accept( Double_t colDist, Double_t colRad, TH2* beamFunc ){ 

  Bool_t accept = 0;

  return accept;  

}
