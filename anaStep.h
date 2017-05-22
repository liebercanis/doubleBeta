//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 22 13:17:26 2017 by ROOT version 6.08/00
// from TTree ntStep/ step variables 
// found on file: legendTree-2017-4-22-13-16-17.root
//////////////////////////////////////////////////////////

#ifndef anaStep_h
#define anaStep_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
enum TrackStatus { active=1, hitPMT=2, absorbed=4, boundaryAbsorbed=8,
                      absorbedLAr=16, inactive=32, hitWLS = 64, totalInternal=128, backScatter=256, notBoundary=512,
                      scint=2*notBoundary, 
                      eIoni=2*scint, 
                      hIoni=2*eIoni, 
                      ionIoni=2*hIoni,
                      compton= 2*ionIoni,
                      hitGe=2*compton, 
                      isBad=2*hitGe};


class anaStep {
public :

   TH1F* hScintTime;
   TH1F* hWlsTime;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         ev;
   Float_t         parent;
   Float_t         pdg;
   Float_t         status;
   Float_t         microsec;
   Float_t         length;
   Float_t         energy;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_parent;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_status;   //!
   TBranch        *b_microsec;   //!
   TBranch        *b_length;   //!
   TBranch        *b_energy;   //!

   anaStep(TTree *tree=0);
   virtual ~anaStep();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef anaStep_cxx

Int_t anaStep::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t anaStep::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void anaStep::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev, &b_ev);
   fChain->SetBranchAddress("parent", &parent, &b_parent);
   fChain->SetBranchAddress("pdg", &pdg, &b_pdg);
   fChain->SetBranchAddress("status", &status, &b_status);
   fChain->SetBranchAddress("microsec", &microsec, &b_microsec);
   fChain->SetBranchAddress("length", &length, &b_length);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   Notify();
}

Bool_t anaStep::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void anaStep::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t anaStep::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef anaStep_cxx
