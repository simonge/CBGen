
#include <TROOT.h>   
#include <iostream>
#include <TRandom3.h>
#include <TSystem.h>
#include "AdaptiveSampler.hh"

using namespace std;

int main( int argc, char** argv ){

  TString inDir   = "";
  TString outFile = "";

  for ( Int_t i=1; i<argc; i=i+2 ) {    
    if      ( TString(argv[i]) == "-d" ) inDir   = argv[i+1];
    else if ( TString(argv[i]) == "-o" ) outFile = argv[i+1];    
  }
    
  if(inDir.EqualTo("") || outFile.EqualTo("")){
    cout << "Input directory and output file must be defined" << endl;
    exit(0);
  }

  TRandom3 *fRandom = new TRandom3(0);
  AdaptiveSampler sampler(9,fRandom);

  void* dirp = gSystem->OpenDirectory(inDir);
  TString fileName;
  Int_t n = 0;


  while(fileName = gSystem->GetDirEntry(dirp)){
    if(fileName.EqualTo("")) break;
    if(!fileName.EndsWith(".astate")) continue;
    if(n==0){
      sampler.restoreState((inDir+fileName).Data());
    }
    else{
      sampler.mergeState((inDir+fileName).Data());
    }
    n++;
  }

  sampler.adapt();
  sampler.saveState(outFile.Data(),1);

}
