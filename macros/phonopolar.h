//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr  1 13:02:12 2017 by ROOT version 6.06/08
// from TTree phono/coherent bremsstrahlung one-phonon process
// found on file: phonopolar.root
//////////////////////////////////////////////////////////

#ifndef phonopolar_h
#define phonopolar_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>

// Headers needed by this particular selector


class phonopolar : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<Double_t> qphoton = {fReader, "qphoton"};
   TTreeReaderArray<Double_t> Qphonons = {fReader, "Qphonons"};
   TTreeReaderValue<Double_t> diffXS = {fReader, "diffXS"};
   TTreeReaderValue<Double_t> polar = {fReader, "polar"};
   TTreeReaderArray<Double_t> k = {fReader, "k"};
   TTreeReaderValue<Double_t> wgt = {fReader, "wgt"};
   TTreeReaderArray<Double_t> logwgts = {fReader, "logwgts"};
   TTreeReaderArray<Int_t>    hkl = {fReader, "hkl"};
   TTreeReaderValue<Double_t> costhetastar = {fReader, "costhetastar"};
   TTreeReaderValue<Int_t>    process = {fReader, "process"};


   phonopolar(TTree * /*tree*/ =0) { }
   virtual ~phonopolar() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(phonopolar,0);

 protected:
   TH1D *hcolspectrum;
   TH1D *hcolspectrum2;
   TProfile *hcolpolarization;
   TH1D *hcolqmag;
   TH1D *lhist;
   TH2D *hkhist;
   TH2D *qThist;
   TH1D *q3hist;
   TH2D *hintens[9];
   TH2D *h1polar[9];
   TH2D *h2polar[9];
   TH2D *h3polar[9];
   TH2D *h4polar[9];
   TH2D *hintensMod[9];
   TH2D *h1polarMod[9];
   TH2D *h2polarMod[9];
   TH2D *h3polarMod[9];
   TH2D *h4polarMod[9];


   TH2D *hDQ2_low_q2[2];
   TH2D *hDQ2q_high_q2[2];
};

#endif

#ifdef phonopolar_cxx
void phonopolar::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t phonopolar::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef phonopolar_cxx
