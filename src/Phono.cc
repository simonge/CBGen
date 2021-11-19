//
// phonoCorr.C - computes the coherent bremsstrahlung energy and
//                polarization spectrum with phonon corrections.
//
// author: richard.t.jones at uconn.edu
// version: march 25, 2017

//#define TEST_SAMPLE_DQ2 1

#include <iostream>
#include <math.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TRandom3.h>

#include "Complex.h"
#include "TPhoton.h"
#include "TLepton.h"
#include "TCrossSection.h"
#include "TThreeVectorReal.h"
#include "TFourVectorReal.h"
#include "TLorentzBoost.h"
#include "constants.h"
#include "Phono.hh"

using namespace std;

//----------------------------------------------------------------------------------------------------------------------

Phono::Phono(Int_t nEvents = 1000, Int_t seed=0, TString fileName="testB.root", TString BeamFile="setup/MainzSettings.txt", TString TargetFile="setup/Diamond.txt", TString AdapterFile = "Adapter.astate", TString AdapterOut = "Adapter.astate", Int_t outputROOT = 0 )
  : statefileCorr(AdapterFile), fSeed(seed), statefileOut(AdapterOut), fROOTout(outputROOT), fFileName(fileName)
{

  if(fROOTout){
    fFile       = new TFile(fFileName, "recreate");
    MakeTree();
    MakeEventTree();
    MakeG4Tree();
  }

  SetBeam(BeamFile);
  SetRadiator(TargetFile);

  nsamples = nEvents;

  Qmax        = 2 * M_PI * hbarc * pow(3 / (4 * M_PI) * Ncell / Vcell, 1/3.);
  kT          = kBoltzmann * temperature;
  fPBeam      = sqrt(sqr(fBeamE) - sqr(mElectron));
  beamVec     = new ROOT::Math::PxPyPzMVector(0,0,fPBeam,mElectron); // CHECK and CORRECT!!!
  electronVec = new ROOT::Math::PxPyPzMVector(0,0,0,0);
  gammaVec    = new ROOT::Math::PxPyPzMVector(0,0,0,0);
  polVec      = new TVector3();
  qlong_220   = fEdge / (fBeamE - fEdge) * sqr(mElectron) / (2 * fBeamE);
  qLong0      = kmin / (fBeamE - kmin) * sqr(mElectron) / (2 * fBeamE);

  fRandom = new TRandom3(0);
  hL1dual = 0;
  L1dual0 = L1dual(0,1000);
  hRDdual = 0;

}

//----------------------------------------------------------------------------------------------------------------------

Phono::Phono(Int_t seed=0, TString fileName="testB.root", Double_t beamE=1.6, Double_t edge=0.4, Double_t eMin = 0.05, Int_t nEvents = 1000 )
  :  fFileName(fileName)
{

  if(fROOTout){
    fFile       = new TFile(fFileName, "recreate");
    MakeTree();
    MakeEventTree();
    MakeG4Tree();
  }

  nsamples = nEvents;
 
  fBeamE        = beamE;
  kmin          = eMin;
  fEdge         = edge;

  Qmax        = 2 * M_PI * hbarc * pow(3 / (4 * M_PI) * Ncell / Vcell, 1/3.);
  kT          = kBoltzmann * temperature;
  fPBeam      = sqrt(sqr(fBeamE) - sqr(mElectron));
  beamVec     = new ROOT::Math::PxPyPzMVector(0,0,fPBeam,mElectron); // CHECK and CORRECT!!!
  electronVec = new ROOT::Math::PxPyPzMVector(0,0,0,0);
  gammaVec    = new ROOT::Math::PxPyPzMVector(0,0,0,0);
  polVec      = new TVector3();
  qlong_220   = fEdge / (fBeamE - fEdge) * sqr(mElectron) / (2 * fBeamE);
  qLong0      = kmin  / (fBeamE - kmin)  * sqr(mElectron) / (2 * fBeamE);

  fRandom = new TRandom3(seed);
  hL1dual = 0;
  L1dual0 = L1dual(0,1000);
  hRDdual = 0;


}

//----------------------------------------------------------------------------------------------------------------------

Phono::Phono(TString setupFile)
{

  TString Key, Value;

  ifstream paramFile(setupFile);

  TString beamFile     = "";
  TString radiatorFile = "";
  fFileName = "";

  while( paramFile >> Key >> Value ){
    if(Key.EqualTo("Repeats"))           nRepeats      = Value.Atoi();
    if(Key.EqualTo("Seed"))              fRandom       = new TRandom3(Value.Atoi());
    if(Key.EqualTo("EventN"))            nsamples      = Value.Atoi();
    if(Key.EqualTo("OutDirectory"))      fFileName     = Value+fFileName;
    if(Key.EqualTo("OutFile"))           fFileName     = fFileName+Value;
    if(Key.EqualTo("BeamSettings"))      beamFile      = Value;
    if(Key.EqualTo("RadiatorSettings"))  radiatorFile  = Value;
    if(Key.EqualTo("AdapterFile"))       statefileCorr = Value;
  }
  
  if(fROOTout){
    fFile       = new TFile(fFileName, "recreate");
    MakeTree();
    MakeEventTree();
    MakeG4Tree();
  }

  if(!beamFile.EqualTo(""))     SetBeam(beamFile);
  if(!radiatorFile.EqualTo("")) SetRadiator(radiatorFile);

  Qmax        = 2 * M_PI * hbarc * pow(3 / (4 * M_PI) * Ncell / Vcell, 1/3.);
  kT          = kBoltzmann * temperature;
  fPBeam      = sqrt(sqr(fBeamE) - sqr(mElectron));
  beamVec     = new ROOT::Math::PxPyPzMVector(0,0,fPBeam,mElectron); // CHECK and CORRECT!!!
  electronVec = new ROOT::Math::PxPyPzMVector(0,0,0,0);
  gammaVec    = new ROOT::Math::PxPyPzMVector(0,0,0,0);
  polVec      = new TVector3();
  qlong_220   = fEdge / (fBeamE - fEdge) * sqr(mElectron) / (2 * fBeamE);
  qLong0      = kmin / (fBeamE - kmin) * sqr(mElectron) / (2 * fBeamE);

  fRandom = new TRandom3(0);
  hL1dual = 0;
  L1dual0 = L1dual(0,1000);
  hRDdual = 0;

}

//----------------------------------------------------------------------------------------------------------------------

Phono::Phono()
{

  fFileName   = "testGlueX.root";
  fFile       = new TFile(fFileName, "recreate");
  MakeTree();
  MakeEventTree();
  MakeG4Tree();

  Qmax        = 2 * M_PI * hbarc * pow(3 / (4 * M_PI) * Ncell / Vcell, 1/3.);
  kT          = kBoltzmann * temperature;
  fPBeam      = sqrt(sqr(fBeamE) - sqr(mElectron));
  beamVec     = new ROOT::Math::PxPyPzMVector(0,0,fPBeam,mElectron); // CHECK and CORRECT!!!
  electronVec = new ROOT::Math::PxPyPzMVector(0,0,0,0);
  gammaVec    = new ROOT::Math::PxPyPzMVector(0,0,0,0);
  polVec      = new TVector3();
  qlong_220   = fEdge / (fBeamE - fEdge) * sqr(mElectron) / (2 * fBeamE);
  qLong0      = kmin / (fBeamE - kmin) * sqr(mElectron) / (2 * fBeamE);

  fRandom = new TRandom3(0);
  hL1dual = 0;
  L1dual0 = L1dual(0,1000);
  hRDdual = 0;

}

//----------------------------------------------------------------------------------------------------------------------

Double_t Phono::nBoseEinstein(Double_t x)
{
   Double_t xfactor = x * kBoltzmann * TDebye / kT;
   if (xfactor > 1e-8)
      return 1 / (exp(xfactor) - 1);
   else
      return 1 / (xfactor + xfactor * xfactor / 2);
}

//----------------------------------------------------------------------------------------------------------------------

Double_t Phono::L1dual(Double_t R, Int_t ndiv=1000)
{
   Double_t RQmax = (R + 1e-18) * Qmax / hbarc;
   Double_t dx = 1. / ndiv;
   Double_t sum = 0;
   for (Int_t i=0; i < ndiv; ++i) {
      Double_t x = (i + 0.5) * dx;
      Double_t term = (2 * nBoseEinstein(x) + 1) * sin(RQmax * x);
      sum += term;
   }
   sum /= RQmax * ndiv;
   return sum;
}

//----------------------------------------------------------------------------------------------------------------------

TH1D* Phono::L1dual_hist(Int_t nbins, Double_t Rmax, Double_t Rmin=1e-12, Int_t ndiv=1000)
{
   Double_t log10Rmin = log(Rmin) / log_10;
   Double_t log10Rmax = log(Rmax) / log_10;
   TH1D *h = new TH1D("hL1dual", "L1dual(R) vs log10(R)",
                      nbins, log10Rmin, log10Rmax);
   for (Int_t i=0; i < nbins; ++i) {
      Double_t R = pow(10, h->GetXaxis()->GetBinCenter(i+1));
      h->SetBinContent(i+1, L1dual(R, ndiv * 5));
   }
   return h;
}

//----------------------------------------------------------------------------------------------------------------------

void Phono::make_hL1dual(Double_t R, Int_t nbins=1000)
{
   Double_t Rmin;
   Double_t Rmax;
   if (hL1dual) {
      nbins = hL1dual->GetNbinsX();
      Rmax = pow(10, hL1dual->GetXaxis()->GetBinUpEdge(nbins));
      if (R <= Rmax)
         return;
      else
         delete hL1dual;
   }
   Rmax = (R > 1e-9)? R : 1e-9;
   Rmin = (R < 1e-12)? R : 1e-12;
   hL1dual = L1dual_hist(nbins, Rmax, Rmin, nbins);
}

//----------------------------------------------------------------------------------------------------------------------

Double_t Phono::RDdual(Double_t q2, Double_t R)
{
   make_hL1dual(R);
   Double_t log10R = log(R) / log_10;
   Double_t L1d = hL1dual->Interpolate(log10R);
   Double_t q2factor = 3 * q2 / (2 * Matom * kBoltzmann * TDebye);
   return R * exp(-q2factor * L1dual0) * (exp(q2factor * L1d) - 1);
}

//----------------------------------------------------------------------------------------------------------------------

TH1D* Phono::RDdual_hist(Double_t q2, Double_t Rmax, Double_t Rmin=1e-12, Int_t ndiv=1000)
{
   Double_t log10Rmin = log(Rmin) / log_10;
   Double_t log10Rmax = log(Rmax) / log_10;
   TH1D *h = new TH1D("hRDdual", "R Ddual(q2,R) vs log10(R)",
                      ndiv, log10Rmin, log10Rmax);
   for (Int_t i=0; i < ndiv; ++i) {
      Double_t R = pow(10, h->GetXaxis()->GetBinCenter(i+1));
      h->SetBinContent(i+1, RDdual(q2, R));
   }
   return h;
}

//----------------------------------------------------------------------------------------------------------------------

void Phono::make_hRDdual(Double_t q2, Double_t R, Int_t ndiv=1000)
{
   Double_t Rmin;
   Double_t Rmax;
   if (hRDdual) {
      Int_t nbins = hRDdual->GetNbinsX();
      Rmax = hRDdual->GetXaxis()->GetBinUpEdge(nbins);
      static Double_t lastq2 = 0;
      if (lastq2 == q2 && R <= Rmax)
         return;
      else
         delete hRDdual;
      lastq2 = q2;
   }
   Rmax = (R > 1e-9)? R : 1e-9;
   Rmin = (R < 1e-12)? R : 1e-12;
   hRDdual = RDdual_hist(q2, Rmax, Rmin, ndiv);
}

//----------------------------------------------------------------------------------------------------------------------

TH1D* Phono::htest_hist(TH1D *clone)
{
   TH1D *h = (TH1D*)gROOT->FindObject("htest");
   if (h == 0) {
      h = (TH1D*)clone->Clone("htest");
      h->Reset();
   }
   return h;
}

//----------------------------------------------------------------------------------------------------------------------

Double_t Phono::DQ2(Double_t q2, Double_t Q, Int_t ndiv=1000)
{
   Double_t Rconverged = 1e-6;
   make_hL1dual(Rconverged);
   make_hRDdual(q2, Rconverged, ndiv);
   TH1D *htest = htest_hist(hRDdual);
   Int_t nbins = hRDdual->GetNbinsX();
   Double_t Q_hbarc = (Q + 1e-50) / hbarc;
   Double_t Rmin = pow(10, hRDdual->GetXaxis()->GetBinLowEdge(1));
   Double_t RD_ = RDdual(q2, Rmin);
   Double_t sum = RD_ * sin(Rmin * Q_hbarc) / (Rmin * Q_hbarc);
   for (Int_t i=0; i < nbins; ++i) {
      Double_t R = pow(10, hRDdual->GetXaxis()->GetBinCenter(i+1));
      Double_t RD = hRDdual->GetBinContent(i+1);
      Double_t term = (RD - RD_) * cos(R * Q_hbarc);
      htest->SetBinContent(i+1, term);
      sum += term;
      RD_ = RD;
   }
   Double_t q2factor = 3 * q2 / (2 * Matom * kBoltzmann * TDebye);
   return sum / (2 * sqr(M_PI) * hbarc);
}

//----------------------------------------------------------------------------------------------------------------------

TH1D* Phono::DQ2_hist(Double_t q2, Double_t Qmax2, Int_t nbins=1000)
{
   TH1D *h = new TH1D("hDQ2", "DQ2(q2,Q)", nbins, 0, Qmax2);
   for (Int_t i=0; i < nbins; ++i) {
      Double_t Q = Qmax2 * (i + 0.5) / nbins;
      Double_t thisDQ2 = DQ2(q2, Q, nbins * 100);
      h->SetBinContent(i+1, thisDQ2);
   }
   return h;
}

//----------------------------------------------------------------------------------------------------------------------

Double_t Phono::sample_hkl(Int_t hkl[3], const Double_t u2[3])
{
   const Double_t hexp = 80.;
   const Double_t kexp = 80.;
   const Double_t lexp = 4.;
   const Double_t log_hexp_factor = log(1 - exp(-1/hexp));
   const Double_t log_kexp_factor = log(1 - exp(-1/kexp));
   const Double_t log_lexp_factor = log(1 - exp(-1/lexp));
   hkl[0] = floor(-hexp * log(fabs(2 * u2[0] - 1)));
   hkl[1] = floor(-kexp * log(fabs(2 * u2[1] - 1)));
   hkl[2] = floor(-lexp * log(fabs(2 * u2[2] - 1)));
   Double_t logwgt = log(8);
   logwgt += hkl[0] / hexp - log_hexp_factor;
   logwgt += hkl[1] / kexp - log_kexp_factor;
   logwgt += hkl[2] / lexp - log_lexp_factor;
   if (u2[0] < 0.5) {
      hkl[0] = -hkl[0] - 1;
   }
   if (u2[1] < 0.5) {
      hkl[1] = -hkl[1] - 1;
   }
   if (u2[2] < 0.5) {
      hkl[2] = -hkl[2] - 1;
   }

   //Diamond FCC
   if(crystalStructure.EqualTo("Diamond")){
     if (hkl[0] / 2 * 2 == hkl[0]) {
       Int_t hklsum = hkl[0] + hkl[1] + hkl[2];
       if (hkl[1] / 2 * 2 != hkl[1] ||
	   hkl[2] / 2 * 2 != hkl[2] ||
	   hklsum / 4 * 4 != hklsum)
	 {
	   return -999;
	 }
     }
     else if (hkl[1] / 2 * 2 == hkl[1] || 
	      hkl[2] / 2 * 2 == hkl[2])
       {
	 return -999;
       }
   }
   else{
     //Face-centered cubic (Copper)
     if(!((hkl[0]/2*2 != hkl[0] &&
	   hkl[1]/2*2 != hkl[1] &&
	   hkl[2]/2*2 != hkl[2]) ||
	  (hkl[0]/2*2 == hkl[0] &&
	   hkl[1]/2*2 == hkl[1] &&
	   hkl[2]/2*2 == hkl[2]) 
	  )
	)
       {
	 return -999;
       }
   }
   return logwgt;
}

//----------------------------------------------------------------------------------------------------------------------

void Phono::prepare_sampling()
{
   TFile roothistos("/w/work1/home/simong/CBremGen/phonoCorr.root");
   hDQ2q_high_q2 = (TH2D*)roothistos.Get("hDQ2q_high_q2");
   hDQ2_low_q2 = (TH2D*)roothistos.Get("hDQ2_low_q2");
   hDQ2q_high_q2->SetDirectory(0);
   hDQ2_low_q2->SetDirectory(0);
}

//----------------------------------------------------------------------------------------------------------------------
// Define the vectors that orient the crystal as follows:
//  1. (2,2,0) lies in the x-z plane with a small negative z component
//  2. (-2,2,0) points close to the y axis with a larger -z component
//  3. (0,0,1) points close the the z axis (beam direction)
//----------------------------------------------------------------------------------------------------------------------

void Phono::SetAxis(Bool_t perp, Double_t rand1, Double_t rand2){

   Double_t tiltV = TMath::ASin(qlong_220 / qtotal_220);
   Double_t tiltH = fHAngle*TMath::DegToRad(); // radians

   if(fRockingAngle!=0){
     tiltV += (1-2*rand1)*fRockingAngle*TMath::DegToRad();
     tiltH += (1-2*rand2)*fRockingAngle*TMath::DegToRad();
   }

   Double_t qUnit = 2*M_PI * hbarc / aLattice;

   qAxis[0].SetPolar(qUnit, M_PI/2, 0);
   qAxis[1].SetPolar(qUnit, M_PI/2, M_PI/2);
   qAxis[2].SetPolar(qUnit, 0, 0);
   for (Int_t i=0; i < 3; ++i) {
     qAxis[i].Rotate(posZhat, +M_PI/4); 
      if(perp){
	qAxis[i].Rotate(posYhat, -tiltH);
	qAxis[i].Rotate(posXhat, +tiltV);
      }
      else{
	qAxis[i].Rotate(posXhat, +tiltH);
	qAxis[i].Rotate(posYhat, -tiltV);
      }
   }
}


//----------------------------------------------------------------------------------------------------------------------
// Setup output tree
//----------------------------------------------------------------------------------------------------------------------
void Phono::MakeTree(){

   fTree = new TTree("phono", "coherent bremsstrahlung one-phonon process");
   fTree->Branch("wgt", &row.wgt, "wgt/D");
   fTree->Branch("logwgts", &row.logwgts[0], "logwgts[8]/D");
   fTree->Branch("logwgts", &row.logwgts[0], "logwgts[5]/D");
   fTree->Branch("qphoton", &row.qphoton[0], "qphoton[3]/D");
   fTree->Branch("Qphonons", &row.Qphonons[0], "Qphonons[3]/D");
   fTree->Branch("diffXS", &row.diffXS, "diffXS/D");
   fTree->Branch("diffXSnocrystal", &row.diffXSnocrystal, "diffXSnocrystal/D");
   fTree->Branch("polar", &row.polar, "polar/D");
   fTree->Branch("circular", &row.circular, "circular/D");
   fTree->Branch("k", &row.k[0], "k[3]/D");
   fTree->Branch("costhetastar", &row.costhetastar, "costhetastar/D");
   fTree->Branch("hkl", &row.hkl[0], "hkl[3]/I");
   fTree->Branch("process", &row.process, "process/I");
   fTree->Branch("u", &u[0], "u[9]/D");
   //fTree->Branch("u", &u[0], "u[11]/D");
   
}


//----------------------------------------------------------------------------------------------------------------------
// Setup output tree
//----------------------------------------------------------------------------------------------------------------------
void Phono::MakeEventTree(){

   fEventTree = new TTree("Events", "Geant4 polarised Bremsstrahlung events");

   fEventTree->Branch("BeamE",        &fBeamE,     "BeamE/D"        );
   fEventTree->Branch("Polarisation", &row.polar,  "Polarisation/D" );
   fEventTree->Branch("PolVec",       &polVec                       );
   fEventTree->Branch("DiffXS",       &row.diffXS, "DiffXS/D"       );
   fEventTree->Branch("Weight",       &Weight,     "Weight/D"       );
   fEventTree->Branch("Electron",     &electronVec                  );
   fEventTree->Branch("Gamma",        &gammaVec                     );
   
}


//----------------------------------------------------------------------------------------------------------------------
// Setup output tree
//----------------------------------------------------------------------------------------------------------------------
void Phono::MakeG4Tree(){

   fG4Tree = new TTree("h1", "Geant4 tagger Bremsstrahlung events");

   fG4Tree->Branch("X_vtx",        &fVx,    "X_vtx/F"           );
   fG4Tree->Branch("Y_vtx",        &fVy,    "Y_vtx/F"           );
   fG4Tree->Branch("Z_vtx",        &fVz,    "Z_vtx/F"           );

   fG4Tree->Branch("Px_e",         &fPxe,   "Px_e/F"            );
   fG4Tree->Branch("Py_e",         &fPye,   "Py_e/F"            );
   fG4Tree->Branch("Pz_e",         &fPze,   "Pz_e/F"            );
   fG4Tree->Branch("En_e",         &fEe,    "En_e/F"            );

   fG4Tree->Branch("Px_p",         &fPxp,   "Px_p/F"            );
   fG4Tree->Branch("Py_p",         &fPyp,   "Py_p/F"            );
   fG4Tree->Branch("Pz_p",         &fPzp,   "Pz_p/F"            );
   fG4Tree->Branch("En_p",         &fEp,    "En_p/F"            );

   fG4Tree->Branch("Polx",         &polX,   "Polx/F"            );
   fG4Tree->Branch("Poly",         &polY,   "Poly/F"            );
   fG4Tree->Branch("Polz",         &polZ,   "Polz/F"            );
   fG4Tree->Branch("Weight",       &Weight, "Weight/D"          );
   
}

//----------------------------------------------------------------------------------------------------------------------
// Generate coherent bremsstrahlung photons from a sum over
// reciprocal lattice vectors, where the gamma radiation is
// accompanied by a phonons radiated or absorbed in the
// radiator crystal assumed to be in equilibrium at a fixed
// temperature T. The results are saved in a root tree.
//----------------------------------------------------------------------------------------------------------------------

Int_t Phono::gentree()
{
   statefileCorr.ReplaceAll(".root", "Corr.astate");

   if (hDQ2_low_q2 == 0 || hDQ2q_high_q2 == 0) {
      prepare_sampling();
   }
   if(fROOTout) fFile->cd();

   //Define the crystal orientation
   SetAxis(1);
   TThreeVectorReal q220(2 * qAxis[0] + 2 * qAxis[1]);
   TThreeVectorReal q_220(-2 * qAxis[0] + 2 * qAxis[1]);
   std::cout << "direction (2,2,0) is "; (q220 / q220.Length()).Print();
   std::cout << "direction (-2,2,0) is "; (q_220 / q_220.Length()).Print();

   // Create the patricles
   TLepton eIn(mElectron), eOut(mElectron);
   TPhoton gOut;

   AdaptiveSampler samplerCorr(fNDims,fRandom);  

   Int_t trial  = 0;
   Int_t trial0 = 0;
   Int_t trial1 = 0;
   Int_t trial2 = 0;
   Int_t trial3 = 0;
   Int_t trial4 = 0;
   Int_t trial5 = 0;
   Int_t trial6 = 0;

   // Repeats to test sampler
   for (Int_t repeat=0; repeat < nRepeats; repeat++) {
     trial  = 0;
     trial0 = 0;
     trial1 = 0;
     trial2 = 0;
     trial3 = 0;
     trial4 = 0;
     trial5 = 0;
     trial6 = 0;

     if (FILE *file = fopen(statefileCorr, "r")) {
       samplerCorr.restoreState(statefileCorr.Data());
//        samplerCorr.reset_stats();
//        samplerCorr.check_subsets();
     }
     else {
       samplerCorr.saveState(statefileCorr.Data());
     } 
     samplerCorr.reset_stats();
  
     if(nRepeats>1){
       TString fFileNameRep(fFileName);
       fFileNameRep.ReplaceAll(".root", TString::Itoa(repeat,10)+".root");
	 
       fFile       = new TFile(fFileNameRep, "recreate");
       fFile->cd();
       MakeTree();
       MakeEventTree();
       MakeG4Tree();
     }

     //samplerCorr.display_tree();

     // Main event gen loop
     for (Int_t sample=0; sample < nsamples; ++sample) {
       
       
       trial++;
            
       //if(sample%1000==0) cout << "Sample: " << sample << "/" << nsamples << " Trials: " << trial << "\r" << flush;

   //     std::cout << trial << " " << trial0 << " " << trial1 << " " << trial2 << " " << trial3 << " " << trial4 << std::endl;

       Double_t logwgts[10] = {};
       logwgts[0] = log(samplerCorr.sample(u));
       if(fRockingAngle>0) SetAxis(1,u[9],u[10]);
       logwgts[1] = sample_hkl(row.hkl, &u[0]);
       Double_t logwgt = logwgts[0] + logwgts[1];
       if (logwgt < -99) {
	 samplerCorr.feedback(u,0);
	 sample--;
	 trial0++;
	 continue;
       }
       TThreeVectorReal qLattice;
       qLattice = row.hkl[0] * qAxis[0] + 
	 row.hkl[1] * qAxis[1] +
	 row.hkl[2] * qAxis[2];

       //       cout << qLattice[0] << " " << qLattice[1] << " " << qLattice[2] << endl;
       
       TFourVectorReal qphoton;


       if(u[8]<0.5){
	 row.process = 1;

	 // The virtual photon momentum is generated by the following steps.
	 // 1. Throw the longitudinal component q[3] according to the density
	 //    f(q[3]) = 1 / [q[3] (1 + q[3] / qLong0)] where qLong0 is
	 //    a constant that sets the minimum q[3]. This corresponds to
	 //    q[3] = qLong0 / (2^u - 1) where u ~ Unif(0,1), and
	 //    w = |d q[3] / du| = log(2) q[3] (1 + q[3] / qLong0)
	
	 Double_t qLong = qLong0 / (pow(2, u[3]) - 1 + 1e-99);
	 logwgts[2] = log(log(2)) + log(qLong) + log(1 + qLong / qLong0);
	
	 // 2. Get total phonon longitudinal momentum Q[3] from q[3], as
	 //    Q[3] = qLattice[3] - q[3]
	
	 TThreeVectorReal Qphonons;
	 Qphonons[3] = qLattice[3] + qLong;
	 Double_t Qlong2 = sqr(Qphonons[3]);
	
	 // 3. Generate total phonon transverse momentum Q[1],Q[2] as
	 //    Q[1] = QT * cos(QTphi)
	 //    Q[2] = QT * sin(QTphi)
	 //    where QTphi = 2 pi u0, u0 ~ Unif(0,1), and
	 //    QT^2 = Q[3]^2 [(1 + QTmax^2 / Q[3]^2)^u1 - 1], u1 ~ Unif(0,1),
	 //    Jacobian{Q[1],Q[2]; u0,u1} = pi Q^2 log(1 + QTmax^2 /Q[3]^2).
	 //    However the Q^2 gets factored Int_to the D(Q,q2) function below
	 //    to keep it from blowing up as Q->0, so the correct weight is
	 //    w = pi log(1 + QTmax^2 / Q[3]^2)
	 //    The value of QTmax is not critical because the function D*Q^2
	 //    converges very rapidly at large Q^2. See the code below for
	 //    how a value for QTmax is obtained.
	
	 Double_t QTmax2;
	 Double_t q2est = qLattice.LengthSqr();
	 Double_t q2_split = hDQ2q_high_q2->GetYaxis()->GetBinLowEdge(1);
	 if (q2est > q2_split) {
	   Int_t imax = hDQ2q_high_q2->GetNbinsX();
	   Double_t xmax = hDQ2q_high_q2->GetXaxis()->GetBinUpEdge(imax);
	   QTmax2 = q2est * sqr(xmax);
	 }
	 else {
	   Int_t imax = hDQ2_low_q2->GetNbinsX();
	   Double_t xmax = hDQ2_low_q2->GetXaxis()->GetBinUpEdge(imax);
	   QTmax2 = sqr(xmax);
	 }
	 Double_t QT = sqrt(Qlong2 * (pow(1 + QTmax2 / Qlong2, u[4]) - 1));
	 Double_t QTphi = 2 * M_PI * u[5];
	 Qphonons[1] = QT * cos(QTphi);
	 Qphonons[2] = QT * sin(QTphi);
	 logwgts[3] = log(M_PI) + log(log(1 + QTmax2 / Qlong2));
	
	 // 4. Generate the transverse components of the photon momentum as
	 //    q[1] = qLattice[1] - Q[1]
	 //    q[2] = qLattice[2] - Q[2]
	
	 qphoton = TFourVectorReal(0, qLattice - Qphonons);
	
	 // 5. Multiply the produce of w factors listed above by the D*Q^2
	 //    function to generate a final weight for this event. The total
	 //    cross section is estimated by the average value of this weight
	 //    multiplied by the differential cross section.
	
	 Double_t q2 = qphoton.LengthSqr();
	 if (q2 > q2_split) {
	   Int_t imax = hDQ2q_high_q2->GetNbinsX();
	   Int_t jmax = hDQ2q_high_q2->GetNbinsY();
	   Double_t xmax = hDQ2q_high_q2->GetXaxis()->GetBinLowEdge(imax);
	   Double_t xmin = hDQ2q_high_q2->GetXaxis()->GetBinLowEdge(1);
	   Double_t ymax = hDQ2q_high_q2->GetYaxis()->GetBinLowEdge(jmax);
	   Double_t ymin = hDQ2q_high_q2->GetYaxis()->GetBinLowEdge(1);
	   Double_t x = Qphonons.Length() / sqrt(q2);
	   x = (x < xmin)? xmin : (x > xmax)? xmax : x;
	   Double_t y = (q2 < ymin)? ymin : (q2 > ymax)? ymax : q2;
	   Double_t DQ2q = hDQ2q_high_q2->Interpolate(x,y);
	   if (DQ2q > 0)
	     logwgts[4] = log(DQ2q) - log(q2) / 2;
	   else
	     logwgts[4] = -999;
	 }
	 else {
	   Int_t imax = hDQ2_low_q2->GetNbinsX();
	   Int_t jmax = hDQ2_low_q2->GetNbinsY();
	   Double_t xmax = hDQ2_low_q2->GetXaxis()->GetBinLowEdge(imax);
	   Double_t xmin = hDQ2_low_q2->GetXaxis()->GetBinLowEdge(1);
	   Double_t ymax = hDQ2_low_q2->GetYaxis()->GetBinLowEdge(jmax);
	   Double_t ymin = hDQ2_low_q2->GetYaxis()->GetBinLowEdge(1);
	   Double_t x = Qphonons.Length();
	   x = (x < xmin)? xmin : (x > xmax)? xmax : x;
	   Double_t y = (q2 < ymin)? ymin : (q2 > ymax)? ymax : q2;
	   Double_t DQ2a = hDQ2_low_q2->Interpolate(x,y);
	   if (DQ2a > 0)
	     logwgts[4] = log(DQ2a) - log(y / q2);
	   else
	     logwgts[4] = -999;
	 }
	 row.Qphonons[0] = Qphonons[1];
	 row.Qphonons[1] = Qphonons[2];
	 row.Qphonons[2] = Qphonons[3];

       }
       // Zero phonon exchange
       else{

	 row.process = 0;
	 qphoton = TFourVectorReal(0, qLattice);
	 row.Qphonons[0] = 0;
	 row.Qphonons[1] = 0;
	 row.Qphonons[2] = 0;

	 logwgts[2] = 0;
	 logwgts[3] = 0;
	 Double_t q2factor = -qphoton.InvariantSqr();
	 q2factor   *= 3/(2*Matom*kBoltzmann*TDebye);
	 logwgts[4]  = -q2factor*L1dual0;

       }

       //cout << qLattice[0] << " " << qLattice[1] << " " << qLattice[2] << endl;
       
       logwgt += logwgts[2] + logwgts[3] + logwgts[4];
       if (logwgt < -99) {
	 samplerCorr.feedback(u,0);
	 sample--;
	 trial1++;
	 continue;
       }
      
       row.qphoton[0] = qphoton[1];
       row.qphoton[1] = qphoton[2];
       row.qphoton[2] = qphoton[3];


       // So far we have generated the total 4-momentum of the final
       // photon-electron system. Now we need to generate the angles
       // for its "decay", which is best done in the rest frame.
 
       TFourVectorReal pin(fBeamE, 0, 0, fPBeam);
       TFourVectorReal ptot(pin + qphoton);
       Double_t M2 = ptot.InvariantSqr();
       if (M2 < sqr(mElectron)) {
	 samplerCorr.feedback(u,0);
	 sample--;
	 trial2++;
	 continue;
       }
       Double_t kstar = (M2 - sqr(mElectron)) / (2 * sqrt(M2));


       // In the comments in TCrossSection::Bremsstralhung I read the following:
       //   "To get a simple expression for the density of final states,
       //    I redefined the solid angle for the outgoing photon around
       //    the momentum axis of the final electron+photon, rather than
       //    the incoming electron direction."
       // To be consistent with this, I need to make sure the azimuthal angle
       // I generate is with respect to the ptot axis, not the beam axis. The
       // jacobian for this is dk / dcosthetastar = gammastar betastar kstar.
 
       Double_t phistar = u[7] * 2*M_PI;
       Double_t thetastar = acos(u[6] * 2 - 1);
      
       TThreeVectorReal khat = ptot / ptot.Length();
       TThreeVectorReal ihat(-khat[3], 0, khat[1]);
       TThreeVectorReal jhat(khat);
       jhat.Cross(ihat);
       TThreeVectorReal uhat(sin(thetastar) * cos(phistar) * ihat +
			     sin(thetastar) * sin(phistar) * jhat +
			     cos(thetastar) * khat);
       TFourVectorReal kout(kstar, kstar * uhat);
       TThreeVectorReal beta(-(TThreeVector)ptot / ptot[0]);
       kout.Boost(beta);

       //cout << row.process << " " << kout[0] <<  " " << kmin << endl;
       if (kout[0] < kmin) {
	 samplerCorr.feedback(u,0);
	 sample--;
	 trial3++;
	 continue;
       }
       Double_t jacob = ptot.Length() * kstar / sqrt(M2);
       logwgts[5] = log(jacob * 4*M_PI);
       logwgt += logwgts[5];
       if (logwgt < -99) {
	 samplerCorr.feedback(u,0);
	 sample--;
	 trial4++;
	 continue;
       }
       
       row.k[0] = kout[1];
       row.k[1] = kout[2];
       row.k[2] = kout[3];
       row.costhetastar = cos(thetastar);

       // compute the differential cross section, polarization
 
       TFourVectorReal pout(ptot - kout);
       eIn.SetMom(pin);
       eOut.SetMom(pout);
       gOut.SetMom(kout);
       eIn.SetPol(zeroVector);
       eOut.AllPol();
       gOut.AllPol();
       row.diffXSnocrystal = TCrossSection::Bremsstrahlung(eIn, eOut, gOut);

       //ADDED
       eIn.SetPol(posZhat);
       gOut.SetPol(posZhat);
       Double_t sameXS = TCrossSection::Bremsstrahlung(eIn, eOut, gOut);
       gOut.SetPol(negZhat);
       Double_t oppXS =  TCrossSection::Bremsstrahlung(eIn, eOut, gOut);
       row.circular = (sameXS-oppXS)/row.diffXSnocrystal;

       gOut.SetPol(uhat.SetPolar(1, M_PI/2., -kout.Phi()));
       row.polar = TCrossSection::Bremsstrahlung(eIn, eOut, gOut);
       row.polar /= row.diffXSnocrystal;

       
       
       // include the atomic and crystal form factors
       Double_t Sff = (row.hkl[0] / 2 * 2 == row.hkl[0])? 8 : 8 / sqrt(2);
       Double_t betacut = 111 * pow(Zatom, -1/3.) / mElectron;
       Double_t Fff = 1 / (1 + sqr(betacut) * qphoton.LengthSqr());
       Double_t XffSqr_d3q = pow(2 * PI_ * hbarc, 3) / Vcell;
       row.diffXS = row.diffXSnocrystal*sqr(Sff) * sqr(Zatom) * sqr(1-Fff) * XffSqr_d3q;

       //cout << row.diffXS << " " << row.process << endl;


       Double_t wgt = exp(logwgt);
       for (Int_t i=0; i < 10; i++)
	 row.logwgts[i] = logwgts[i];
       row.wgt = wgt;
       if (!std::isfinite(wgt)) {
	 std::cerr << "phonoCorr warning: found infinite or undefined wgt,"
		   << " discarding the following event:" << std::endl
		   << "  wgt = " << row.wgt << std::endl
		   << "  logwgts = " 
		   << row.logwgts[0] << ", "
		   << row.logwgts[1] << ", "
		   << row.logwgts[2] << ", "
		   << row.logwgts[3] << ", "
		   << row.logwgts[4] << ", "
		   << row.logwgts[5] << ", "
		   << row.logwgts[6] << ", "
		   << row.logwgts[7] << ", "
		   << row.logwgts[8] << ", "
		   << row.logwgts[9] << std::endl
		   << "  hkl = "
		   << row.hkl[0] << ", " << row.hkl[1] << ", " << row.hkl[2]
		   << "  qphoton = " 
		   << row.qphoton[0] << ", " 
		   << row.qphoton[1] << ", "
		   << row.qphoton[2] << std::endl
		   << "  Qphonons = " 
		   << row.Qphonons[0] << ", " 
		   << row.Qphonons[1] << ", " 
		   << row.Qphonons[2] << std::endl
		   << "  diffXS = "
		   << row.diffXS << ", polar = "
		   << row.polar << std::endl
		   << "  k = "
		   << row.k[0] << ", "
		   << row.k[1] << ", "
		   << row.k[2] << std::endl
		   << std::endl
		   << "  costhetastar = " << row.costhetastar
		   << std::endl;
	 wgt = 0;
       }


       // feedback to the sampler
       samplerCorr.feedback(u, row.diffXS * wgt);

       // apply a minimum cut on wgt * diffXS for saving this event
 
       Double_t minPfactor = 1e7;
       //if((!onlyType && trial%2) || (onlyType && onlyInCoh) ) minPfactor = 1e9;
       Double_t thisPfactor = wgt * row.diffXS;
       //       cout << wgt << " " << row.diffXS << " " << thisPfactor << " " << thisPfactor / minPfactor <<endl;
       if (thisPfactor == 0) {
	 sample--;
	 trial5++;
	 continue;
       }
       else if (thisPfactor < minPfactor) {
	 Double_t uu = fRandom->Uniform();
	 if (uu > thisPfactor / minPfactor){
	   sample--;
	   trial6++;
	   continue;
	 }
	 row.wgt = minPfactor / row.diffXS;
       }

       if((sample % 5000 == 0) && nRepeats!=1) samplerCorr.adapt();

       //std::cout << trial0 << " " << trial1 << " " << trial2 << " " << trial3 << " " << trial4 << " " << trial5 << " " << trial6 << endl << std::endl;
       trial0 = 0;
       trial1 = 0;
       trial2 = 0;
       trial3 = 0;
       trial4 = 0;
       trial5 = 0;
       trial6 = 0;

       TVector3 temp = TVector3(row.k[0],row.k[1],row.k[2]).Cross(TVector3(row.k[0],row.k[1],row.k[2]).Cross(TVector3(row.qphoton[0],row.qphoton[1],row.qphoton[2])));
       polVec->SetXYZ(temp.X(),temp.Y(),temp.Z());

       polVec->SetMag(row.polar);
       Weight = row.wgt*row.diffXS;
       gammaVec->SetCoordinates(row.k[0],row.k[1],row.k[2],0);
       electronVec->SetCoordinates(-row.k[0]-row.qphoton[0]-row.Qphonons[0],-row.k[1]-row.qphoton[1]-row.Qphonons[1],fBeamE-mElectron-row.k[2]-row.qphoton[2]-row.Qphonons[2],mElectron);

       // G4 output variables
       fVx  = 0;
       fVy  = 0;
       fVz  = 0;

       fPxe = electronVec->Px();
       fPye = electronVec->Py();
       fPze = electronVec->Pz();
       fEe  = electronVec->E();

       fPxp = gammaVec->Px();
       fPyp = gammaVec->Py();
       fPzp = gammaVec->Pz();
       fEp  = gammaVec->E();
      
       polX = polVec->X();
       polY = polVec->Y();
       polZ = polVec->Z();
       
       if(fROOTout){
	 fTree->Fill();
	 fEventTree->Fill();
	 fG4Tree->Fill();
       }
     }
     if(fROOTout){
       fTree->Write();
       fEventTree->Write();
       fG4Tree->Write();
       fFile->Close();
     }
//      cout << trial << endl;
//      cout << nsamples << endl;
     std::cout << "sampler-corr reports efficiency " << samplerCorr.getEfficiency() << std::endl;
//      samplerCorr.adapt();
     
//      TString outState(statefileCorr); 
//      outState.ReplaceAll(".astate", TString::Itoa(fSeed,10)+".astate");
//      statefileCorr = outState;
//      samplerCorr.saveState(outState.Data(),1);
     samplerCorr.saveState(statefileOut.Data(),1);

     //samplerCorr.reset_stats();
   }

   return nsamples;

}

//----------------------------------------------------------------------------------------------------------------------

void Phono::SetBeam(TString bFile)
{

  TString Key, Value;

  ifstream paramFile(bFile);


  while( paramFile >> Key ){
    if(Key.EqualTo("BeamEnergy"))        paramFile >> fBeamE;
    if(Key.EqualTo("MinEnergy"))         paramFile >> kmin;
    if(Key.EqualTo("PeakEnergy"))        paramFile >> fEdge;
    if(Key.EqualTo("HAngle"))            paramFile >> fHAngle;
  }
  
}


//----------------------------------------------------------------------------------------------------------------------

void Phono::SetRadiator(TString rFile)
{

  TString Key, Value;

  ifstream paramFile(rFile);

  while( paramFile >> Key ){
    if(Key.EqualTo("Atomic"))          paramFile >> Zatom;
    if(Key.EqualTo("AtomicMass"))      paramFile >> Matom;
    if(Key.EqualTo("CellN"))           paramFile >> Ncell;
    if(Key.EqualTo("CellVol"))         paramFile >> Vcell;
    if(Key.EqualTo("LatticeA"))        paramFile >> aLattice;
    if(Key.EqualTo("DebyeTemp"))       paramFile >> TDebye;
    if(Key.EqualTo("Temperature"))     paramFile >> temperature;
    if(Key.EqualTo("RockingAngle"))    paramFile >> fRockingAngle;
    if(Key.EqualTo("Structure"))       paramFile >> crystalStructure;
  }

}


//----------------------------------------------------------------------------------------------------------------------

