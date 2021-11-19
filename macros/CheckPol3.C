#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TF1.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TSystem.h>

using namespace std;

Double_t GausOnBase(Double_t *, Double_t *);
void parFromHuman(Double_t beamMeV = 1508.0, Double_t edgeMeV = 650.0, Double_t spreadMeV = 10.0, 
		  Double_t colliDist_m = 2.5, Double_t colliRad_mm = 3, Int_t nVec = 2, Double_t *par=NULL);
void enhFromParams(Double_t *par=NULL);
Double_t efit(const Double_t *);


TCanvas *genCanvas=NULL;
TF1 *gausFit;

//Some enumerators and names
enum {
  THETA,  // [0] theta      main angle responsible for coherent edge cutoffs
  SIGMA,  // [1] sigma      smearing of theta
  THETAR, // [2] thetar     relative angle resonsible for colli cutoffs
  SIGMAR, // [3] sigmar     smearing of colli cutoff angle
  //  POWER,
  E0MEV,  // [4] beam energy
  NVEC,   // [5] nvec       no of vectors contributing
  IVEC};  // [6] ivec[]     array of intensities of vectors up to nvec.

//Approx Form factor is F(g^2) = (q^2 + b^(-2)) ^ -2
//Where b= 111 x Z^(-1/3) (x 925 to get into units of crystal lattice)
const Double_t B = 0.247892436;  //where did I get that ? Timm ?
const Double_t A=0.03;           //made up for now, need to get the actual no for this later
const Double_t k=26.5601;        //put in formula for k later (my own stonehenge paper)

//Double_t beamMeV = 1508.0;
Double_t beamMeV = 1508.0;
Double_t edgeGuess = 645.0;
//Double_t colliDist_m = 0.0001;
//Double_t colliRad_mm = 10;
Double_t colliDist_m = 2.5;
Double_t colliRad_mm = 3;
Int_t nVec=2;

const Int_t VECTORS[]={2,4,6,8,10};    //list of the vectors to be included (022,044);

Int_t THETASTEPS = 201;          //no of steps in convoluting with gaussian
Double_t LOWFIT = 308.0;         //how far below the peak in MeV to start fitting

Double_t fitMinEnergy;
Int_t fitMinBin;
Double_t fitMaxEnergy;
Int_t fitMaxBin;
Int_t verbose=0;
Double_t bestPar[11];
Double_t bestChisq;

Int_t counter = 1;
 
//these are really just convenient arrays - they don't ever get plotted.
TH1F *weightHist  = NULL;
TH2F *thetaWeight = NULL;
TH2F *thetaPol    = NULL;
TH2F *thetaTtot   = NULL;
TH2F *thetaItot   = NULL;

//All histograms and related are globals
TH1F *histP   = NULL;  //pol
TH1F *histE   = NULL;  //enh from calculation
TH1D *tempRes =NULL;

TH1F *histD = NULL;  //enh from data to be fitted
TH1D *edgeProj =NULL;
TH1D *edgeProjGaus =NULL;
TH1D *hMaxPol =NULL;
TH1D *hChi2 =NULL;

const Double_t *Errors;
Double_t lowError;
Double_t highError;

//fit of a gaussian on a baseline, for gausFit,
Double_t GausOnBase(Double_t *x, Double_t *par) {
  Double_t arg = 0;
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];
   Double_t fitval = par[3] + par[0]*TMath::Exp(-0.5*arg*arg);
   return fitval;
}

//void CheckPol3( TString inputFile = "/scratch/simong/FileCheckPol2.root", TString outFile = "/scratch/simong/FileCheckPol3.root" ){
//void CheckPol3( TString inputFile = "/scratch/simong/FileCheckPol2Alt106.root", TString outFile = "/scratch/simong/FileCheckPol3Alt106-6-NEW.root", Double_t Baseline = 1.06, Double_t PeakRatio = 1.75){
void CheckPol3( TString inputFile = "/w/work14/simong/PolOut/All-Enh-CHEAT.root", TString outFile = "/w/work14/simong/PolOut/All-PolRANDOM.root", Double_t Baseline = 1.06, Double_t PeakRatio = 1.75){

  //Setup Canvas for checking stuff
  if(genCanvas) delete genCanvas;
  genCanvas = new TCanvas("genCanvas","genCanvas",5,5,1005,655);
  genCanvas->Divide(4,2);   
  genCanvas->GetPad(1)->SetGridx(1);
  genCanvas->GetPad(1)->SetGridy(1);
  genCanvas->GetPad(2)->SetGridx(1);
  genCanvas->GetPad(2)->SetGridy(1);
  genCanvas->GetPad(3)->SetGridx(1);
  genCanvas->GetPad(3)->SetGridy(1);
  genCanvas->GetPad(3)->SetGridx(1);
  genCanvas->GetPad(3)->SetGridy(1);

  //Open input files
  TFile*  iFile   = new TFile(inputFile,"read"    );
  TFile*  oFile   = new TFile(outFile,  "recreate");
  
  TH2F*    hFileEG     = (TH2F*)iFile->Get("FileEG");  
  TH2F*    hFileEGP    = (TH2F*)iFile->Get("FileEGP");  
  TH2F*    hFileEGR    = (TH2F*)iFile->Get("FileEGR");
  
  //  TH1F*    hFileEGamor = (TH1F*)iFile->Get("FileEGamor");
  TH2F*    hFileEGamor  = (TH2F*)iFile->Get("FileEGamor");
  TH2F*    hFileEGamorP = (TH2F*)iFile->Get("FileEGamorP");
  TH2F*    hFileEGamorR = (TH2F*)iFile->Get("FileEGamorR");
  
  TH2F*    hFileEGenh   = (TH2F*)iFile->Get("FileEGenh");
  TH2F*    hFileEGenhP  = (TH2F*)iFile->Get("FileEGenhP");
  TH2F*    hFileEGenhR  = (TH2F*)iFile->Get("FileEGenhR");
  
  TH2F*    hFileEGenhs = (TH2F*)hFileEGenh->Clone("FileEGenhs");

  TH2F*    hFileEGenhcal  = (TH2F*)hFileEGenh->Clone("FileEGenhcal");

  TH2F*    hFileEGpol     = (TH2F*)hFileEGenh->Clone("FileEGpol");

  TH2F*    hFileEGres     = (TH2F*)hFileEGenh->Clone("FileEGres");

  TH1D*    tempProj       = ((TH2F*)hFileEGenh)->ProjectionX("projX");

  TH2F*    hPolEdge       = new TH2F("PolEdge","PolEdge",15,630,660,20,0.4,0.5);

  tempRes      = ((TH2F*)hFileEGenh)->ProjectionX("tempRes");

  edgeProj     = ((TH2F*)hFileEGenh)->ProjectionY("edgeProj");
  edgeProjGaus = ((TH2F*)hFileEGenh)->ProjectionY("edgeProjGaus");

  hMaxPol      = ((TH2F*)hFileEGenh)->ProjectionY("MaxPol");
  hChi2        = ((TH2F*)hFileEGenh)->ProjectionY("Chi2");

  hChi2->Reset();
  hChi2->SetMinimum(0.9);

  tempRes->Reset();
  hMaxPol->Reset();
  edgeProj->Reset();
  edgeProjGaus->Reset();
  
  histD    = (TH1F*)tempProj->Clone("projX2");
  histD->Reset();
  hFileEGres->Reset(); 
  hFileEGpol->Reset(); 
  hFileEGenhcal->Reset();

  UInt_t nBins = hFileEGenhs->GetNbinsX();
  Double_t diff1, diff2, lowmean, highmean1, highmean2,fitedge,scalefac,scale;
  Double_t par[11];
  ROOT::Math::Minimizer* min;
  Char_t name[30];

  fitedge = 0;

  gausFit=new TF1("gausFit",GausOnBase,0,100,4);

  //Get polarisation
  //  for(Int_t i=195; i<210;i++){
  //  for(Int_t i=1; i<hFileEGenh->GetNbinsY(); i+=10){
  for(Int_t i=425; i<474; i++){

    histD->Reset();

    //Get rid of zeros
    for(UInt_t n=1;n<=nBins-1;n++){
      if(hFileEGenhs->GetBinContent(n,i)<0.1){
	hFileEGenhs->SetBinContent(n,i,hFileEGenhs->GetBinContent(n+1,i));
	hFileEGenhs->SetBinError(n,i,hFileEGenhs->GetBinError(n+1,i));
      }
    }
    //Get rid of zeros 2nd pass
    for(UInt_t n=1;n<=nBins-1;n++){
      if(hFileEGenhs->GetBinContent(n,i)<0.1){
	hFileEGenhs->SetBinContent(n,i,hFileEGenhs->GetBinContent(n+1,i));
	hFileEGenhs->SetBinError(n,i,hFileEGenhs->GetBinError(n+1,i));
      }
    }
    //  Get rid of spikes up and down
    for(UInt_t n=2;n<=nBins-1;n++){
      diff1=(hFileEGenhs->GetBinContent(n,i)-hFileEGenhs->GetBinContent(n-1,i))/hFileEGenhs->GetBinContent(n-1,i);
      diff2=(hFileEGenhs->GetBinContent(n,i)-hFileEGenhs->GetBinContent(n+1,i))/hFileEGenhs->GetBinContent(n+1,i);

      if (((fabs(diff1)>0.03)&&(fabs(diff2)>0.03))&&(fabs(diff1-diff2)<0.1)){
	hFileEGenhs->SetBinContent(n,i,0.5*(hFileEGenhs->GetBinContent(n-1,i)+hFileEGenhs->GetBinContent(n+1,i)));
	hFileEGenhs->SetBinError(n,i,0.5*(hFileEGenhs->GetBinError(n-1,i)+hFileEGenhs->GetBinError(n+1,i)));
      }
    }
  

    hFileEGenhs->GetYaxis()->SetRange(i,i);

    //    hFileEGenhs->Scale(1/(Double_t)hFileEGenh->GetNbinsY());

    if(hFileEGenhs->GetEntries()==0)  continue;

    cout << hFileEGenhs->GetEntries() << endl;

//     lowmean=10000000.0;
//     for(UInt_t n=125;n<=160;n++){
//       hFileEGenhs->GetXaxis()->SetRange(n-1,n+6);
//       cout << hFileEGenhs->Integral() << endl;
//       if((hFileEGenhs->Integral()<lowmean)){
// 	lowmean=hFileEGenhs->Integral();
//       }
//     }
//     lowmean/=8;
//     if(lowmean==0)  continue;

//     cout << lowmean << endl;

//     if(lowmean !=10000000.0){
//       hFileEGenhs->GetXaxis()->SetRange();
//       for(UInt_t n=0;n<nBins;n++){
// 	hFileEGenhs->SetBinContent(n,i,hFileEGenhs->GetBinContent(n,i)*(Baseline*6.0/(lowmean)));
// 	hFileEGenhs->SetBinError(n,i,hFileEGenhs->GetBinError(n,i)*(Baseline*6.0/(lowmean)));
// 	//	hFileEGenhs->SetBinContent(n,i,hFileEGenhs->GetBinContent(n,i)*(6.0/(lowmean)));
//       }
//     }

//     highmean1=0;
//     for(UInt_t n=115;n<=135;n++){
//       hFileEGenhs->GetXaxis()->SetRange(n-6,n+1);
//       //      cout << hFileEGenhs->Integral() << endl;
//       if((hFileEGenhs->Integral()>highmean1)){
// 	highmean1=hFileEGenhs->Integral();
// 	highmean1/=8;
//       }
//     }

//     highmean2=0;
//     for(UInt_t n=180;n<=210;n++){
//       hFileEGenhs->GetXaxis()->SetRange(n-6,n+1);
//       //      cout << hFileEGenhs->Integral() << endl;
//       if((hFileEGenhs->Integral()>highmean2)){
// 	highmean2=hFileEGenhs->Integral();
// 	highmean2/=8;
//       }
//     }

//     if(highmean1!=0 && highmean2!=0){
//       hFileEGenhs->GetXaxis()->SetRange();
//       for(UInt_t n=0;n<nBins;n++){
// 	hFileEGenhs->SetBinContent(n,i,hFileEGenhs->GetBinContent(n,i)*((highmean1)/((highmean2)*1.6)));
// 	hFileEGenhs->SetBinError(n,i,hFileEGenhs->GetBinError(n,i)*((highmean1)/((highmean2)*1.6)));
// 	//	hFileEGenhs->SetBinContent(n,i,hFileEGenhs->GetBinContent(n,i)*(6.0/(lowmean)));
//       }
//     }

//     scale = Baseline/lowmean;
//     //    scale = (((highmean1)/((highmean2)*PeakRatio))+(Baseline/(lowmean)))/2;

//     if(highmean1!=0 && highmean2!=0 && lowmean !=1000000000.0){
//       hFileEGenhs->GetXaxis()->SetRange();
//       for(UInt_t n=0;n<nBins;n++){
// 	hFileEGenhs->SetBinContent(n,i,hFileEGenhs->GetBinContent(n,i)*scale);
// 	hFileEGenhs->SetBinError(n,i,hFileEGenhs->GetBinError(n,i)*scale);
// 	//	hFileEGenhs->SetBinContent(n,i,hFileEGenhs->GetBinContent(n,i)*(6.0/(lowmean)));
//       }
//     }

//     cout << "SCALE = " << scale << endl;
//     cout << (Baseline/(lowmean)) << endl;
//     cout << ((highmean1)/((highmean2)*PeakRatio)) << endl;
    

    histD->Add(((TH2F*)hFileEGenhs)->ProjectionX(""));

    histD->GetXaxis()->SetRange(125,160);
    histD->SetMaximum(1.2*histD->GetMaximum());
    histD->SetMinimum(0.0);
    histD->GetXaxis()->SetRange();

    histE = (TH1F*)histD->Clone("histE");
    histE->Sumw2(0);
    histE->Reset();
    histE->GetXaxis()->SetRange(125,160);
    histE->SetMaximum(1.2*histD->GetMaximum());
    histE->SetMinimum(0.0);
    histE->GetXaxis()->SetRange();
    histP = (TH1F*)histD->Clone("histP");
    histP->Reset();
    histP->SetMaximum(1);

    edgeProj->SetMaximum(700);
    edgeProj->SetMinimum(600);

    genCanvas->cd(1);
    histD->Draw("P");
    //Fit
    //Now try to make some guesses at the initial parameters
    histD->GetXaxis()->SetRange(185,210);

    Int_t maxBin2 = histD->GetMaximumBin();

    histD->GetXaxis()->SetRange(125,160);
    gausFit->SetRange(histD->GetBinCenter(histD->GetMaximumBin()),histD->GetBinCenter(histD->GetMaximumBin())+50.0);
    gausFit->SetParameter(1,histD->GetBinCenter(histD->GetMaximumBin()));
    gausFit->SetParameter(2,10.0);
    gausFit->SetParameter(3,1.0);
    histD->Fit(gausFit,"rN");

    lowmean=0.0;

    Int_t maxBin = histD->GetMaximumBin();

    //Get the edge from the derivative
    for(Float_t d = histD->GetBinCenter(maxBin);d < histD->GetBinCenter(maxBin+90.0);d+=0.1){ 
      //	cout << gausFit->Derivative(d) << endl;
      if(gausFit->Derivative(d)<lowmean){
	lowmean=gausFit->Derivative(d);
	fitedge=d;
      }
    }

    cout << gausFit->GetParameter(1) << endl;
    cout << gausFit->GetParError(1) << endl;
    cout << gausFit->GetParameter(2) << endl;
    cout << gausFit->GetParError(2) << endl;

    edgeProjGaus->SetBinContent(i,fitedge);
    edgeProjGaus->SetBinError(i,gausFit->GetParError(1)+(gausFit->GetParError(2)*TMath::Abs(fitedge-gausFit->GetParameter(1)))/(gausFit->GetParameter(2)));


//     cout << beamMeV << endl;
//     cout << fitedge << endl;
//     cout << edgeGuess << endl;
//     cout << gausFit->GetParameter(2) << endl;

    histD->GetXaxis()->SetRange(125,160);

    if(abs(fitedge-edgeGuess)>50) continue;

    //Now we have enough information to set the basic parameters
    parFromHuman(beamMeV,fitedge,gausFit->GetParameter(2),colliDist_m,colliRad_mm,nVec,par);

    //set the intensities
    for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
      par[IVEC+v] = histD->GetBinContent(maxBin)*2.0/((Double_t)VECTORS[v]*(Double_t)VECTORS[v]);      //tailing off as 1/VECTORS[v]^2
    }
    
    cout << histD->GetBinContent(maxBin) << " " << histD->GetBinContent(maxBin2) << endl;
    cout << (histD->GetBinContent(maxBin)-1)/(histD->GetBinContent(maxBin2)-1) << endl;
    if(histD->GetEntries()==0) continue;
    cout <<  par[IVEC+0]/par[IVEC+1] << endl;

    enhFromParams(par);

    //Redo the intensities according to a the calc / data ration
    scalefac=histD->GetMaximum()/histE->GetMaximum();
    for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
      par[IVEC+v]*=scalefac;
    }

    histE->GetXaxis()->SetRange(50,200);
    histD->GetXaxis()->SetRange(50,200);
//     histE->GetXaxis()->SetRange();
//     histD->GetXaxis()->SetRange();
    enhFromParams(par);
    histE->SetLineColor(2);
    histE->Draw("same");    
    genCanvas->cd(2);
    histP->SetLineColor(2);
    histP->Draw();
    genCanvas->cd(3);
    edgeProj->SetMarkerSize(0.3);
    edgeProj->SetMarkerStyle(3);
    edgeProj->Draw("P");
    edgeProjGaus->SetMarkerSize(0.3);
    edgeProjGaus->SetMarkerStyle(3);
    edgeProjGaus->SetMarkerColor(2);
    edgeProjGaus->Draw("same");
    genCanvas->Update();
    gSystem->ProcessEvents();
    
    //Set the range of the fit to be some sensible amount below peak and just past the 2nd peak.
    fitMinBin=histE->FindBin(histD->GetBinCenter(histD->GetMaximumBin())-LOWFIT);
    fitMaxBin    = histE->FindBin(par[E0MEV]/(((2.0/4.25)*((par[E0MEV]/histD->GetBinCenter(histD->GetMaximumBin()))-1.0))+1));

    histE->GetXaxis()->SetRange(fitMinBin,fitMaxBin);
    histD->GetXaxis()->SetRange(fitMinBin,fitMaxBin);
    histP->GetXaxis()->SetRange(fitMinBin,fitMaxBin);

    cout << "fitMaxBin " << fitMaxBin << endl;
    min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simple");
    
    // set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
    min->SetMaxIterations(10000);  // for GSL 
    min->SetTolerance(0.01);
    min->SetPrintLevel(1);
    

    ROOT::Math::Functor ft(&efit,IVEC+nVec);     
    min->SetFunction(ft);
    
  
    //Now set the variables
//     min->SetLimitedVariable(THETA,   "Theta",   par[THETA],      par[THETA]/1000.0,  0.9*par[THETA], 1.1*par[THETA]);
//     min->SetLimitedVariable(SIGMA,   "Sigma",   1.5*par[SIGMA],      par[SIGMA]/100.0,  0.5*par[SIGMA],  80.0*par[SIGMA]);
//     min->SetLimitedVariable(THETAR,  "Thetar",  15*par[THETAR],     par[THETAR]/100.0, 0.5*par[THETAR], 80.0*par[THETAR]);
//     min->SetLimitedVariable(SIGMAR,  "Sigmar",  15*par[SIGMAR],     par[SIGMAR]/100.0, 0.5*par[SIGMAR], 80.0*par[SIGMAR]);
    min->SetLimitedVariable(THETA,   "Theta",   par[THETA],      par[THETA]/100.0,  0.9*par[THETA], 1.1*par[THETA]);
    min->SetLimitedVariable(SIGMA,   "Sigma",   2.5*par[SIGMA],  par[SIGMA]/100.0,  par[SIGMA],  5.0*par[SIGMA]);
    min->SetLimitedVariable(THETAR,  "Thetar",  2.0*par[THETAR],     par[THETAR]/100.0, 0.2*par[THETAR], 5.0*par[THETAR]);
    //min->SetFixedVariable(THETAR,  "Thetar",  40.0*par[THETAR]);
    min->SetLimitedVariable(SIGMAR,  "Sigmar",  0.5*par[SIGMAR], par[SIGMAR]/100.0, 0.1*par[SIGMAR], 20*par[SIGMAR]);
    //    min->SetFixedVariable(SIGMAR,  "Sigmar",  0.01*par[SIGMAR]);
    //    min->SetLimitedVariable(POWER,  "Power",  par[POWER], par[POWER]/100.0, 0.1*par[POWER], 20*par[POWER]);
    min->SetFixedVariable(E0MEV,     "E0MeV",   par[E0MEV]);  //no of vectors
    min->SetFixedVariable(NVEC,      "Nvec",    par[NVEC]);  //no of vectors
    cout << "C" << endl;
    for(int n=0;n<nVec;n++){
      sprintf(name,"Vec0%d%d", VECTORS[n],VECTORS[n]);
      //      min->SetFixedVariable(n+IVEC, name, 2*par[n+IVEC]); 
      min->SetLimitedVariable(n+IVEC, name, 0.9*par[n+IVEC], par[n+IVEC]/100.0, 0.5*par[n+IVEC], 15*par[n+IVEC]); 
    }
    
    verbose=1;             //make it show the fit as it's improving
    bestChisq=100000.00;   //set this high for starters
    
    min->Minimize();

    //    if(bestChisq>1.02) continue;

//     Errors = min->Par();
    cout << "PARAMETERS" << endl;
    cout << par[0] << endl;
    cout << par[1] << endl;
    cout << par[2] << endl;
    cout << par[3] << endl;
    cout << par[4] << endl;
    cout << par[5] << endl;
    cout << par[6] << endl;
    cout << par[7] << endl;
    cout << par[8] << endl;
  
    enhFromParams(bestPar);

    histE->GetXaxis()->SetRange(125,160);

    hPolEdge->Fill(fitedge,histP->GetBinContent(histP->GetMaximumBin()));

    hMaxPol->SetBinContent(i,histP->GetBinContent(histP->GetMaximumBin()));

    gausFit->SetRange(histE->GetBinCenter(histE->GetMaximumBin()),histE->GetBinCenter(histE->GetMaximumBin())+50.0);
    gausFit->SetParameter(1,histE->GetBinCenter(histE->GetMaximumBin()));
    gausFit->SetParameter(2,10.0);
    gausFit->SetParameter(3,1.0);
    histE->Fit(gausFit,"rN");

    maxBin = histE->GetMaximumBin();

    lowmean=0.0;

    //Get the edge from the derivative
    for(Float_t d = histE->GetBinCenter(maxBin);d < histE->GetBinCenter(maxBin+90.0);d+=0.1){ 
      //	cout << gausFit->Derivative(d) << endl;
      if(gausFit->Derivative(d)<lowmean){
	lowmean=gausFit->Derivative(d);
	fitedge=d;
      }
    }

    edgeProj->SetBinContent(i,fitedge);
    //    edgeProj->SetBinError(i,gausFit->GetParError(1)+(gausFit->GetParError(2)*TMath::Abs(fitedge-gausFit->GetParameter(1)))/(gausFit->GetParameter(2)));

    for(Int_t n=0; n<tempRes->GetNbinsX(); n++){
      hFileEGres->SetBinContent(n,i,tempRes->GetBinContent(n));
    }

    cout << fitedge << endl;
    cout << gausFit->GetParameter(1) << endl;
    cout << gausFit->GetParError(1) << endl;
    cout << gausFit->GetParameter(2) << endl;
    cout << gausFit->GetParError(2) << endl;

    histE->GetXaxis()->SetRange();

    hChi2->SetBinContent(i,bestChisq);

    histE->GetXaxis()->SetRange();
    histD->GetXaxis()->SetRange();
    histP->GetXaxis()->SetRange();

    genCanvas->cd(1);
    histD->Draw("P");
    histE->SetLineColor(2);
    histE->Draw("same");    
    genCanvas->cd(2);
    histP->SetLineColor(2);
    histP->Draw();
    genCanvas->cd(3);
    edgeProj->SetMarkerSize(0.3);
    edgeProj->SetMarkerStyle(3);
    edgeProj->Draw("P");
    edgeProjGaus->SetMarkerSize(0.3);
    edgeProjGaus->SetMarkerStyle(3);
    edgeProjGaus->SetMarkerColor(2);
    edgeProjGaus->Draw("same");
    genCanvas->cd(4);
    hMaxPol->SetMarkerSize(0.3);
    hMaxPol->SetMarkerStyle(3);
    hMaxPol->Draw("P");
    genCanvas->cd(5);
    hPolEdge->Draw("colz");
    genCanvas->cd(6);
    hChi2->SetMarkerSize(0.3);
    hChi2->SetMarkerStyle(3);
    hChi2->Draw("P");

    cout << "LOOK HERE " <<  histP->GetBinContent(histP->GetMaximumBin()) << endl;

    gSystem->Sleep(500);

    histD->GetXaxis()->SetRange();

    for(UInt_t n=1;n<=nBins-1;n++){

      hFileEGenhcal->SetBinContent(n,i,histE->GetBinContent(n));
      hFileEGpol   ->SetBinContent(n,i,histP->GetBinContent(n));

      cout << TMath::Sqrt((938.27+histP->GetBinLowEdge(n))*(938.27+histP->GetBinLowEdge(n))-(histP->GetBinLowEdge(n)*histP->GetBinLowEdge(n))) << " " << histP->GetBinContent(n) << endl;

    }



    //    break;

  }
  
  hFileEGenhs->GetYaxis()->SetRange();


  genCanvas->Write();
  edgeProj->Write();
  edgeProjGaus->Write();
  hFileEG->Write();
  hFileEGP->Write();
  hFileEGR->Write();
  hFileEGamor->Write();
  hFileEGamorP->Write();
  hFileEGamorR->Write();
  hFileEGenh->Write();
  hFileEGenhP->Write();
  hFileEGenhR->Write();
  hFileEGenhs->Write();
  hFileEGenhcal->Write();
  hFileEGpol->Write();
  hFileEGres->Write();
  hPolEdge->Write();

  oFile->Close();
  iFile->Close();

  cout << endl;
  cout << inputFile << endl;
  cout << outFile << endl << endl;

}
  
void parFromHuman(Double_t beamMeV, Double_t edgeMeV, Double_t spreadMeV, Double_t colliDist_m, Double_t colliRad_mm, Int_t nVec, Double_t *par){

  //takes some physical quantities and makes them into parameters, then calls the 
  //enhFromParams function.
  
  
  //  Double_t par[10];                                                           //array of parameters
  Int_t g = 2;                                                                //variables used in CLAS note
  Double_t E0 = beamMeV;
  Double_t Eg = edgeMeV;
  
  
  par[THETA]  = k/(g*E0*E0*((1/Eg)-(1/E0)));                                  //theta from edge and beam energy
  par[SIGMA]  = (par[THETA]-(k/(g*E0*E0*((1/(Eg-spreadMeV))-(1/E0)))))/3.0;   //spread in theta from spread in edge 
  par[THETAR] = E0*0.001*5.0*colliRad_mm/colliDist_m;                         //cut from collimator
  //par[THETAR] = E0*0.0005*5.0*colliRad_mm/colliDist_m;                         //cut from collimator
  par[SIGMAR] = par[THETAR]*par[SIGMA]/par[THETA];                            //smear in above same fractional sigma as above
  //par[SIGMAR] = 0.5*par[THETAR]*par[SIGMA]/par[THETA];                            //smear in above same fractional sigma as above
  //  par[POWER]  = -0.5;
  par[E0MEV]  = E0;                                                           //beam energy
  par[NVEC]   = (Double_t)nVec;                                                         //no of harmonics

  for(int v=0;v<par[NVEC];v++){                                               //give the vectors intensities
    par[IVEC+v] = 2.0/(Double_t)VECTORS[v];                                   //tailing off as 1/VECTORS[v]
    //cout << IVEC+v << "  v   " << par[IVEC+v] << endl; 
  }

  enhFromParams(par);

}  


//The main customized fitting function which gets called by MINUIT
Double_t efit(const Double_t *parms){
  
  Double_t chisq = 1.0;
  Double_t delta;
  Double_t b1,b2;
  Double_t err;
  Double_t *par = (Double_t*)parms;

  histE->Reset("ICE"); //reset the histogram
  //  cout << "par[0]= " << par[0] << endl;

  //call the function to make the enhancement and polarization
  enhFromParams(par);

  chisq = 1.0;
  //loop over all the required bins in the histogram to work out a chisq
  for(int n=fitMinBin;n<=fitMaxBin;n++){
    b1=histE->GetBinContent(n);
    b2=histD->GetBinContent(n);
    err=1;//histE->GetBinError(n)/histD->GetBinError(n);
    delta=(b1-b2)/err;
    chisq+=(delta*delta);

    tempRes->SetBinContent(n,delta/((b1+b2)/2));
    //note - not a proper chisq because its an enhancement

  }
   
  fprintf(stderr,"Chisq: \t%6.2f\t\r",chisq);


  if(chisq<bestChisq){
    bestChisq=chisq;
    for(int n=0;n<10;n++){
      bestPar[n]=par[n];
    }
    if(verbose){
      if(10%(counter++)){
	//if verbose, draw this on the canvas for every iteration to see how it's going
	histE->GetXaxis()->SetRange(50,200);
	histD->GetXaxis()->SetRange(50,200);
// 	histE->GetXaxis()->SetRange();
// 	histD->GetXaxis()->SetRange();
	genCanvas->cd(1);
	histD->SetLineColor(4);
	histD->Draw("P");
	histD->SetMinimum(0.0);
	histD->SetMaximum(3.0);
	genCanvas->cd(1);
	histE->Draw("same");
	
	genCanvas->cd(2);
	histP->Draw();
	histP->SetMinimum(0);
	histP->SetMaximum(1);

	genCanvas->cd(7);
	tempRes->SetMarkerSize(0.3);
	tempRes->SetMarkerStyle(3);
	tempRes->SetMarkerColor(4);
	tempRes->Draw("PL");
	
	genCanvas->Draw();   
	
	genCanvas->Update();
	gSystem->ProcessEvents();
	counter=1;
      }
    }
  }
  return chisq;
  
}

void enhFromParams(Double_t *par){
  //make an enhancement and corresponding polarization from some the parameters as defined in the CLAS note.
  //this function is can be called stand alone, but will also ba called many times from the fitting function

  Double_t xd[10];
  Double_t xc[10];
  Double_t Q[10];
  Double_t cohContrib;
  Double_t cohTotal;
  Double_t phiTotal;
  Double_t etotal;
  Double_t ptotal;
  Double_t x=0.0;
  Int_t    g=0;
  Double_t weight=0.0;
  Double_t weightError=0.0;
  Double_t weightSum=0.0;
  Double_t polSum=0.0;
  Double_t weightSumError=0.0;
  Double_t polSumError=0.0;
  Double_t phi,chi,cd;
  Double_t amo;
  Int_t jbin=0;
 
  //loop over sigma
  // for(int p=0;p<10;p++){
  //  cout << p << ": " << par[p] << ", ";
  //}
  //cout << endl;

  // if needed, make some hists
//   if(!histE){
//     histE        = new TH1F("Enhancement", "Enhancement",  1000, 0, par[E0MEV]);
//     histP        = new TH1F("Polarization", "Polarization",1000, 0, par[E0MEV]);
//     histE->SetMinimum(0);
//     histP->SetMinimum(0);
//     histP->SetMaximum(1);
//   }    
  if(!thetaPol){
    weightHist   = new TH1F("weightHist",  "weightHist", THETASTEPS+1, 0, THETASTEPS+1 );
    thetaWeight  = new TH2F("thetaWeight", "thetaWeight",histE->GetNbinsX(), histE->GetXaxis()->GetXbins()->GetArray(), THETASTEPS+1,0, THETASTEPS+1);
    thetaPol     = new TH2F("thetaPol",    "thetaPol",   histE->GetNbinsX(), histE->GetXaxis()->GetXbins()->GetArray(), THETASTEPS+1,0, THETASTEPS+1);
    thetaItot    = new TH2F("thetaItot",   "thetaItot",  histE->GetNbinsX(), histE->GetXaxis()->GetXbins()->GetArray(), THETASTEPS+1,0, THETASTEPS+1);
    weightHist ->Sumw2();
    thetaWeight->Sumw2();
    thetaPol   ->Sumw2();
    thetaItot  ->Sumw2();
  }

  //reset them all for fresh filling
  histE->Reset("ICE");
  histP->Reset("ICE");
  thetaPol->Reset("ICE");
  thetaItot->Reset("ICE");
  weightHist->Reset("ICE");
  thetaWeight->Reset("ICE");

  for(Double_t j=par[THETA]-3.0*par[SIGMA];j<=par[THETA]+3.001*par[SIGMA];j+=(6.0*par[SIGMA])/THETASTEPS){
    
    weight=TMath::Gaus(j,par[THETA],par[SIGMA]);   //get the weight from the gaussian
    weightSum+=weight;                             //add to sum      
    
    //find the discontinuity for each vector
    for(int v=0;v<par[NVEC];v++){
      g=VECTORS[v];
      xd[v]=1.0/((k/(g*par[E0MEV]*j))+1.0);
      Q[v]=(1.0-xd[v])/xd[v];
      xc[v]=xd[v]/(1+((par[THETAR]*par[THETAR])*(1-xd[v])));
    }

    //loop over all bins in the histogram
    for(int bin=1;bin<=histE->GetNbinsX();bin++){
      x=histE->GetBinCenter(bin)/par[E0MEV];            //find the value of the bin
      //amo=1/x;                                    //assume amo = inc = 1/x over region of interest
      //      amo=TMath::Power(x,par[POWER]);                                    //assume amo = inc = 1/x over region of interest
      //      amo=TMath::Power(x,-0.5);                                    //assume amo = inc = 1/x over region of interest
      amo=TMath::Power(x,-1);                                    //assume amo = inc = 1/x over region of interest
      
      cohTotal=0.0;
      phiTotal=0.0;
      
      //loop over all the vectors
      for(int v=0;v<par[NVEC];v++){
	if(x>xd[v]) continue;           //only do up to x_dg
	 
	//work out chi and phi
	phi=(2*Q[v]*Q[v]*x*x)/((1-x)*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1))));
	chi=((Q[v]*Q[v])/(1-x))*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1)));
	//	chi=((Q[v]*Q[v]*x)/(1-x))*(1+((1-x)*(1-x))-((4*Q[v]*Q[v]*x*x/(1-x))*(((1-x)/(Q[v]*x))-1))); //RETURN
	//	cout << j  << "  " << chi << endl;
	cd=0.5*(1+TMath::Erf((x-xc[v])/(TMath::Sqrt(2)*par[SIGMAR])));

	//get coherent contrib for the vector
	cohContrib=cd*par[IVEC+v]*chi;

	//add to the total and update the phi total
	cohTotal+=cohContrib;
	phiTotal+=cohContrib*phi;
	
      }
      if(cohTotal>0.0) {
	phiTotal/=cohTotal;   //divide by the cohTotal to get the weighted dmean phi
	//cout << x << " " << phiTotal << " " << cohTotal << " " << weight << endl;	 
      }

      //      cout << j << " " << cohTotal << " " << phiTotal << endl << endl;

      //enhancement = coherent total + inc (or amo).
      etotal=1+(cohTotal)/(amo);
      //      etotal=(amo+cohTotal)/amo;
      //and pol like this
      //      ptotal=phiTotal*cohTotal/(cohTotal + amo);
      ptotal=phiTotal*cohTotal;
      //ptotal=-phiTotal*(1-(1/etotal));

      //add the weighted contribution to the enhancement
      histE->Fill(x*par[E0MEV],weight*etotal);
//       histE->Fill(x*par[E0MEV],cohTotal);
//       genCanvas->cd(4);
//       histE->Draw();
      //histE->Fill(x*par[E0MEV],etotal);

      //keep the pol for this x,theta coord
      thetaPol->Fill(x*par[E0MEV],jbin,ptotal);

      //keep the total intensity for this x,theta coord
      thetaItot->Fill(x*par[E0MEV],jbin,cohTotal+amo);
    }
    
    //save the weight for this theta point
    weightHist->Fill(jbin,weight);
    jbin++;

  }
  //normalize the sum of the weighted enhancements
  histE->Scale(1.0/weightSum);
  
  
  //loop over each x bin, adding the weighted contribs from each theta pos
  for(int bin=1; bin<=histP->GetNbinsX(); bin++){
    weightSum=0.0;
    polSum=0.0;
    
    for(int jb=1;jb<=weightHist->GetNbinsX();jb++){
      weight=weightHist->GetBinContent(jb);
      weightError=weightHist->GetBinError(jb);

      //      polSum+=thetaPol->GetBinContent(bin,jb)*thetaItot->GetBinContent(bin,jb)*weight;
      polSum+=thetaPol->GetBinContent(bin,jb)*weight;
      weightSum+=thetaItot->GetBinContent(bin,jb)*weight;

      polSumError+=thetaPol->GetBinError(bin,jb)*weightError;
      weightSumError+=thetaItot->GetBinError(bin,jb)*weightError;
      //polSum+=thetaPol->GetBinContent(bin,jb)*weight;
      //weightSum+=weight;
    }
    polSum/=weightSum;
    polSumError/=weightSumError;
    //    histP->Fill(histP->GetBinCenter(bin),polSum);
    histP->SetBinContent(bin,polSum);
    histP->SetBinError(bin,polSumError);
  } 
}

//------------------------------------------


void CheckAll(){

  //  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt100-2.root",1.0,1.75);
  //  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt101-2.root",1.01,1.75);
  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt102-2.root",1.02,1.75);
  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt103-2.root",1.03,1.75);
  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt104-2.root",1.04,1.75);
  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt105-2.root",1.05,1.75);
  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt106-2.root",1.06,1.75);
  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt107-2.root",1.07,1.75);
  CheckPol3("/scratch/simong/FileCheckPol2Alt100.root","/scratch/simong/FileCheckPol3Alt110-2.root",1.10,1.75);



}
