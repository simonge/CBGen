//
// phonoCorr.C - computes the coherent bremsstrahlung energy and
//                polarization spectrum with phonon corrections.
//
// author: richard.t.jones at uconn.edu
// version: march 25, 2017

//#define TEST_SAMPLE_DQ2 1

#ifndef Phono_h
#define Phono_h 1

#include <iostream>
#include <math.h>
#include "TMath.h"
#include "Math/Vector4D.h"

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
#include "TVector3.h"

#include "AdaptiveSampler.hh"


struct treerow {
  Double_t wgt;
  Double_t logwgts[10];
  Double_t qphoton[3];
  Double_t Qphonons[3];
  Double_t diffXSnocrystal;
  Double_t diffXS;
  Double_t polar;
  Double_t circular;
  Double_t k[3];
  Double_t costhetastar;
  Int_t hkl[3];
  Int_t process;
};

const TThreeVectorReal zeroVector(0,0,0);
const TThreeVectorReal posXhat(1,0,0);
const TThreeVectorReal negXhat(-1,0,0);
const TThreeVectorReal posYhat(0,1,0);
const TThreeVectorReal negYhat(0,-1,0);
const TThreeVectorReal posZhat(0,0,1);
const TThreeVectorReal negZhat(0,0,-1);

//---------------------------------------------------------------------------------------

class Phono {
  
public:
  Phono(Int_t nEvents, Int_t seed, TString fileName, TString BeamFile, TString TargetFile, TString AdapterFile, TString AdapterOut, Int_t outputROOT );
  Phono(Int_t seed, TString fileName, Double_t beamE, Double_t edge, Double_t eMin, Int_t nEvents );
  Phono(TString settingsFile);
  Phono();
  
  Int_t    gentree();


private:

  TString  fFileName;
  TFile*   fFile;
  TTree*   fTree;
  TTree*   fEventTree;
  TTree*   fG4Tree;
  treerow  row;
  TString  statefileCorr = "Adapter.astate";
  TString  statefileOut  = "Adapter.astate";
  Int_t    fSeed;

  Int_t    fROOTout = 1;
  Int_t    nRepeats = 1;
  Int_t    nsamples = 10000;

  //Constants
  const Double_t hbarc       = 1.97327e-16;    // GeV.m  (hbar*c constant)
  const Double_t kBoltzmann  = 8.617333262e-14; // GeV/K (Boltzmann constant)
  const Double_t log_10 = log(10.);

  // parameters for the beamline
  Double_t fBeamE = 12.0;      // GeV
  Double_t fPBeam;      // Beam momentum


  // physical parameters for the target crystal
  Double_t fEdge = 9.0;      // GeV
  Double_t fHAngle = 60.0000;

  Double_t qtotal_220 = 9.8e-6; // GeV
  Double_t qlong_220;
  Double_t aLattice    = 3.567e-10;  // m
  Double_t TDebye      = 2200;       // K
  Double_t temperature = 300;        // K
  Double_t kT;                       // GeV
  
  Double_t Ncell  = 8;
  Double_t Vcell  = 45.5e-30;       // m^3
  Double_t Zatom  = 6;              // GeV/c^2
  Double_t Matom  = 12 * 0.932;     // GeV/c^2
  Double_t Qmax; // GeV/c
  Double_t fRockingAngle = 0;
  TString  crystalStructure = "Diamond";

  Double_t kmin = 1; // low-energy cutoff in bremsstrahlung spectrum (GeV)
  Double_t qLong0;

  TThreeVectorReal qAxis[3];

  const static Int_t    fNDims = 11;
  //const static Int_t    fNDims = 9;
  Double_t u[fNDims];
  
  Double_t L1dual0;
  
  Bool_t onlyType  = 0;
  Bool_t onlyInCoh = 0;


  TH2D* hDQ2q_high_q2;
  TH2D* hDQ2_low_q2;
  TH1D* hL1dual;
  TH1D* hRDdual;

  TRandom3* fRandom;

  // Event Tree Variables 
  Double_t  Weight;
  TVector3* polVec;
  ROOT::Math::PxPyPzMVector* beamVec;
  ROOT::Math::PxPyPzMVector* electronVec;
  ROOT::Math::PxPyPzMVector* gammaVec;

  // G4 Tree Variables
  Float_t fVx;
  Float_t fVy;
  Float_t fVz;

  Float_t fPxe;
  Float_t fPye;
  Float_t fPze;
  Float_t fEe;

  Float_t fPxp;
  Float_t fPyp;
  Float_t fPzp;
  Float_t fEp;

  Float_t polX;
  Float_t polY;
  Float_t polZ;

  

//   void     unif01(Int_t n, Double_t *u){fRandom->RndmArray(n,u);}
  void     SetAxis(Bool_t perp, Double_t rand1=0, Double_t rand2=0);
  void     MakeTree();
  void     MakeEventTree();
  void     MakeG4Tree();

  Double_t nBoseEinstein(Double_t x);
  Double_t L1dual(Double_t R, Int_t ndiv);

  TH1D*    L1dual_hist(Int_t nbins, Double_t Rmax, Double_t Rmin, Int_t ndiv);
  void     make_hL1dual(Double_t R, Int_t nbins);
  void     set_temperature(Double_t T) { kT = kBoltzmann * T;}
  Double_t RDdual(Double_t q2, Double_t R);
  TH1D*    RDdual_hist(Double_t q2, Double_t Rmax, Double_t Rmin, Int_t ndiv);
  void     make_hRDdual(Double_t q2, Double_t R, Int_t ndiv);
  TH1D*    htest_hist(TH1D *clone);
  Double_t DQ2(Double_t q2, Double_t Q, Int_t ndiv);
  TH1D*    DQ2_hist(Double_t q2, Double_t Qmax, Int_t nbins);
  Double_t sample_hkl(Int_t hkl[3], const Double_t u[3]);
  void     prepare_sampling();

  void     SetBeam(TString);
  void     SetRadiator(TString);
  
  inline Int_t sqr(Int_t x) { return x*x; }
  inline Float_t sqr(Float_t x) { return x*x; }
  inline Double_t sqr(Double_t x) { return x*x; }
  inline Complex_t sqr(Complex_t x) { return x*x; }


};
//---------------------------------------------------------------------------------------

#endif
