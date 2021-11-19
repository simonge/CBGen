#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChainElement.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;

void SinglePol( Int_t fileNumber,
		  TString fileType ){

  //Open input files
  TString DataDir = "/scratch/simong/";
  TString DataAmo = "/scratch/simong/";
  TString outFile = "/scratch/simong/" + TString::Itoa(fileNumber,10) + "Pol.root";
  TString suffix  = "Out.root";
  TString inputFile = DataDir + TString::Itoa(fileNumber,10) + suffix;

  TChain oTree("PID","datachain");

  oTree.Add(inputFile);
 
  TFile*  oFile   = new TFile(outFile,"recreate");
 
  Int_t fMaxTagged = 100;

  //Input Tree
  Double_t       fTGamma[fMaxTagged];
  Int_t          fTaggedN;
  Int_t          fPrompt[fMaxTagged];
  TClonesArray   *fTagged;
  
  UInt_t bins      = 352;
  Double_t *fEBins = new Double_t[bins+1];
  //Beam Energy Bins
  ifstream f;
  Double_t tempE;
  f.open("/home/simong/A2Edam/data/EgBins.txt.ALL");
  for(UInt_t i = 0; i <= bins; i++){
    
    f >> tempE;
    fEBins[i] = tempE;
    //    fEBins[i] = ToCME(tempE);

  }
  f.close();

  fTagged  = new TClonesArray("TLorentzVector",fMaxTagged);

  oTree.SetCacheSize(10000000);

  oTree.SetBranchAddress("TagTime",fTGamma);
  oTree.SetBranchAddress("FNTagged",&fTaggedN);
  oTree.SetBranchAddress("Prompt",fPrompt);
  oTree.SetBranchAddress("Tagged",&fTagged);

  oTree.SetBranchStatus("*",0);
  oTree.SetBranchStatus("FNTagged",1);
  oTree.SetBranchStatus("TagTime",1);
  oTree.SetBranchStatus("Prompt",1);
  oTree.SetBranchStatus("Tagged",1);

  Int_t LowFiles  = 18700;
  Int_t HighFiles = 19500;
  UInt_t NFiles    = HighFiles - LowFiles;

  Double_t *FileBins = new Double_t[NFiles+1];

  for(UInt_t i = 1; i <= NFiles; i++){

    FileBins[i] = LowFiles+i;
    
  }
  FileBins[0] = FileBins[1]-1;

  TH2F*    hFileEGPERPP        = new TH2F("FileEGperpP",
				      "FileEGperpP",
				      bins, fEBins,
				      NFiles, FileBins);

  TH2F*    hFileEGPERPR        = new TH2F("FileEGperpR",
				      "FileEGperpR",
				      bins, fEBins,
				      NFiles, FileBins);

  TH2F*    hFileEGPARAP        = new TH2F("FileEGparaP",
				      "FileEGparaP",
				      bins, fEBins,
				      NFiles, FileBins);

  TH2F*    hFileEGPARAR        = new TH2F("FileEGparaR",
				      "FileEGparaR",
				      bins, fEBins,
				      NFiles, FileBins);

  TH2F*    hFileEGamorP    = new TH2F("FileEGamorP",
				      "FileEGamorP",
				      bins, fEBins,
				      NFiles, FileBins);
  
  TH2F*    hFileEGamorR    = new TH2F("FileEGamorR",
				      "FileEGamorR",
				      bins, fEBins,
				      NFiles, FileBins);

  tempE = 0;

  for(Int_t i = 0; i<oTree.GetEntries(); i+=1){
    oTree.GetEntry(i);

    for(Int_t j=0;j<fTaggedN;j++){

      tempE = ((TLorentzVector*)fTagged->At(j))->E();

      if(fileType == "PARA"){
	if(fTGamma[j]<0 || fTGamma[j]>200) continue;
	if(fPrompt[j]==0) hFileEGPARAP->Fill(tempE,fileNumber);
	if(fPrompt[j]==1) hFileEGPARAR->Fill(tempE,fileNumber); 
      }
      else if(fileType == "PERP"){
	if(fTGamma[j]<0 || fTGamma[j]>200) continue;
	if(fPrompt[j]==0) hFileEGPERPP->Fill(tempE,fileNumber);
	if(fPrompt[j]==1) hFileEGPERPR->Fill(tempE,fileNumber); 
      }      
      else{
	if(fTGamma[j]<0 || fTGamma[j]>200) continue;
	if(fPrompt[j]==0) hFileEGamorP->Fill(tempE,fileNumber);
	if(fPrompt[j]==1) hFileEGamorR->Fill(tempE,fileNumber); 
      }      

    }  
  }

  TH2F*    hFileEGPARA = (TH2F*)(hFileEGPARAP->Clone("FileEGpara"));
  hFileEGPARA->Add(hFileEGPARAR,-0.125);

  TH2F*    hFileEGPERP = (TH2F*)(hFileEGPERPP->Clone("FileEGperp"));
  hFileEGPERP->Add(hFileEGPERPR,-0.125);
  
  TH2F*    hFileEGamor = (TH2F*)(hFileEGamorP->Clone("FileEGamor"));
  hFileEGamor->Add(hFileEGamorR,-0.125);

  hFileEGPARA->Write();
  hFileEGPARAP->Write();
  hFileEGPARAR->Write();
  hFileEGPERP->Write();
  hFileEGPERPP->Write();
  hFileEGPERPR->Write();
  hFileEGamor->Write();
  hFileEGamorP->Write();
  hFileEGamorR->Write();

  oFile->Close();

}
