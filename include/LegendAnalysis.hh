#ifndef LegendAnalysis_h
#define LegendAnalysis_h 1

#include "globals.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"

//#include "g4root.hh"
#include "TTree.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "LTEvent.hxx"
// singleton class for root file handling
// M.G. 
// .. is it multi-thread safe?
class LegendAnalysis
{
  private:
    LegendAnalysis() { Initialize(); };
    void Initialize();
    TFile   *fHistFile;
    TFile   *fTreeFile;
    TTree  *fTree;
    LTEvent *fEvent;
    static LegendAnalysis* fLegendAnalysis; 
    // Disabled (not implemented) copy constructor and asignment.
    LegendAnalysis(const LegendAnalysis&);
    LegendAnalysis& operator=(const LegendAnalysis&);
    TH1F *hOptical;
    TH1F *hWls;
    TH1F *hPmtHits;
    TH1F *hEElectron;
    TH1F *hEGamma;
    TH1F *hScintYield;
    TH1F *hWlsYield;
    TH1F *hPmtHitCount;
    
 
  public:
    ~LegendAnalysis() {
      printSummary();
      // normalize histogram
      bool ok = fHistFile->cd("OpBoundary");
      if(ok) {
        TH1F* hOptBoundary = (TH1F*) gDirectory->Get("OptBoundary"); 
        TH1F* hOptBoundaryWeight = (TH1F*) gDirectory->Get("OptBoundaryWeight"); 
        TH1F* hOptBoundaryAve = (TH1F*) gDirectory->Get("OptBoundaryAve"); 
        for(unsigned ib=0; ib<hOptBoundary->GetNbinsX(); ++ib) {
          G4int nbin = hOptBoundary->GetBinContent(ib);
          if(nbin<1) continue;
          G4double wbin = hOptBoundaryWeight->GetBinContent(ib);
          hOptBoundaryAve->SetBinContent( ib,wbin/G4double(nbin) );
        }
      }
      
      // Be careful when writing the final Tree header to the file!
      fTreeFile->cd();
      G4cout<<" LegendAnalysis:  working root directory  is  "; fTreeFile->pwd();G4cout <<  G4endl; 
      fTreeFile->ls();
      //fTree->Print();
      G4cout<<" LegendAnalysis: fTreeFile->Write(); " << G4endl; 
      fTreeFile->Write();
      G4cout<<" LegendAnalysis: fTreeFile->Close(); " << G4endl; 
      fTreeFile->Close();
      G4cout<<" LegendAnalysis: fHistFile->Write(); " << G4endl; 
      fHistFile->Write();
      G4cout<<" LegendAnalysis: fHistFile->Close(); " << G4endl; 
      fHistFile->Close();
      
    } 
    static LegendAnalysis* Instance();
    
    TDirectory *topTreeDir() { return fTreeFile->GetDirectory(NULL); }
    TDirectory *topHistDir() { return fHistFile->GetDirectory(NULL); }
    TTree *getTree() { return fTree; }
    TFile *getTreeFile() { return fTreeFile; }
    TFile *getHistFile() { return fHistFile; }
    void anaEvent(const G4Event* anEvent);
    void anaTrajectories(G4TrajectoryContainer* trajectoryContainer);
    LTEvent* getEvent() { return fEvent;}
    void printSummary();
}
      
    
;
#endif
