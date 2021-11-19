
#include "TMath.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TF1.h"  
#include "TH1.h"   
#include "TH2.h"   
#include "TH3.h"   
#include "TRandom3.h"   
#include "TLorentzVector.h"   
#include "TVector3.h"   
#include <iostream>

#include "Phono.hh"
//#include "Amor.hh"

using namespace std;

int main( int argc, char** argv ){

  Phono* generator;
  //Amor*  generatorA;
  Int_t nEvents = 20000;

  if(argc==1){
    generator = new Phono();
  }
  else{
    TString  outDir       = "/scratch/simong/";
    TString  outFile      = "test.root";
    Int_t    seed         = 1234567890;
    Double_t beamE        = 1.6;
    Double_t peakE        = 0.5;
    TString  settingsFile = "";
    TString  adapterFile  = "";
    TString  adapterOut   = "";
    TString  beamFile     = "";
    TString  targetFile   = "";
    Int_t    writeROOT    = 0;

    for ( Int_t i=1; i<argc; i=i+2 ) {    
      if      ( TString(argv[i]) == "-r" ) settingsFile = argv[i+1];
      else if ( TString(argv[i]) == "-d" ) outDir       = argv[i+1];    
      else if ( TString(argv[i]) == "-f" ) outFile      = argv[i+1];    
      else if ( TString(argv[i]) == "-e" ) beamE        = atof(argv[i+1]);
      else if ( TString(argv[i]) == "-p" ) peakE        = atof(argv[i+1]);
      else if ( TString(argv[i]) == "-s" ) seed         = atoi(argv[i+1]);      
      else if ( TString(argv[i]) == "-n" ) nEvents      = atoi(argv[i+1]);     
      else if ( TString(argv[i]) == "-b" ) beamFile     = argv[i+1];   
      else if ( TString(argv[i]) == "-t" ) targetFile   = argv[i+1];   
      else if ( TString(argv[i]) == "-a" ) adapterFile  = argv[i+1]; 
      else if ( TString(argv[i]) == "-o" ) adapterOut   = argv[i+1];
      else if ( TString(argv[i]) == "-root" ) writeROOT = atoi(argv[i+1]);
    }
    
    if(!settingsFile.EqualTo("")){
      generator = new Phono(settingsFile);
    }
    else{
      generator = new Phono(nEvents,seed,outDir+outFile,beamFile,targetFile,outDir+adapterFile,outDir+adapterOut,writeROOT);
      //       generator = new Phono(seed,outDir+outFile,beamE,peakE,0.05,nEvents);
    }
    

    
  }

  generator->gentree();

}
