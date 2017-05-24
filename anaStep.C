#define anaStep_cxx
#include "anaStep.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
anaStep::anaStep()  
{
  char tag[80] = "2017-4-22-15-32-50";
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  TString fileName = TString("legendTree-") + TString(tag) + TString(".root");
  printf(" looking for file %s\n",fileName.Data());
  TFile *f = new TFile(fileName,"readonly");
  if(!f) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }
  fChain = (TChain*) f->Get("ntStep");
  fChain->ls();

   Init();
  
   // open ouput file and make some histograms
  TString outputFileName = TString("anaStep.root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  outfile->cd();

   printf(" opening output file %s \n",outputFileName.Data());

   hScintTime = new TH1F("ScintTime","  scint time microsec ",1200,0,12);
   hWlsTime = new TH1F("WlsTime"," WLS time microsec ",1200,0,12);
   gDirectory->pwd();

   Loop();

   TF1 *f1 = new TF1("f1", "expo", .01, 12);
   hScintTime->Fit("f1", "R");
   TF1 *f2 = new TF1("f2", "expo", .01, 12);
   hWlsTime->Fit("f2", "R");
   double c1=f1->GetParameter(0);
   double c2=f2->GetParameter(0);
   double r1=f1->GetParameter(1);
   double r2=f2->GetParameter(1);
   

   printf(" scint tau = %f WLS tau = %f \n",-1/r1,-1/r2);
  
   
   outfile->Write();
}

anaStep::~anaStep()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void anaStep::Loop()
{
//   In a ROOT session, you can do:
//      root> .L anaStep.C
//      root> anaStep t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int istatus = int(status);
      bool isScint =   istatus&TrackStatus::scint ;
      bool isWls   =   istatus&TrackStatus::hitWLS;
      if ( isScint&&!isWls )  hScintTime->Fill(tglobal);
      if ( isWls )  hWlsTime->Fill(tglobal);
   }
}
