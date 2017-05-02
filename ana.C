#include <vector>
#include <TLorentzVector.h>
#include <TVector3.h>


void ana(char *chtag = "1493395802")//Feb9_50V_5000V")
{

  TString tag(chtag);
  TString inputFileName = TString("legend-")+tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
  TTree *ltree=NULL;
  unsigned aSize = 0;

  // tree has to be in file
  ltree = (TTree*) infile->Get("LTree");
  if(ltree) aSize=ltree->GetEntriesFast();
  printf(" ltree with %i entries \n",int(aSize));

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("ana-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());
  TH1F* hTrajZ = new TH1F("TrajZ"," trajectory z ",320, -210, 110 );
  hTrajZ->SetXTitle(" absolute position Z (mm)");
  TH1F* hTrajHitZ = new TH1F("TrajHitZ"," pmt hit trajectory z ",320, -210, 110 );
  hTrajHitZ->SetXTitle(" absolute position Z (mm)");
  TH1F* hOpticalYield = new TH1F("OpticalYield"," # Ar Scint / e- energy  ",1000, 0,50000);
  hOpticalYield->SetXTitle(" photons/energy(MeV) ");

  TH1F* hArYield = new TH1F("ArYield"," # Ar Scint / e- energy  ",1000, 0,50000);
  hArYield->SetXTitle(" photons/energy(MeV) ");


  TH2F* hPhotonsVsEnergy = new TH2F("PhotonsVsEnergy"," # detected photons versus e- energy  ",100,0,1.0,100,0,1000);
  hPhotonsVsEnergy->SetXTitle(" primary energy(MeV) ");
  hPhotonsVsEnergy->SetYTitle(" detected photons ");

  
  LTEvent *ev = new LTEvent();
  ltree->SetBranchAddress("event",&ev);

  for(unsigned entry =0; entry < aSize; ++entry ) {
    if(entry%100==0) printf("\t entry %i \n",entry);
    ltree->GetEntry(entry);

    //
    hOpticalYield->Fill( float(ev->nOptPhotons)/ev->ePrimary);
    hArYield->Fill( float(ev->nArScint)/ev->ePrimary);
    //float evertex  =  ev->pvertex[0].particle[0].KEnergy;
    hPhotonsVsEnergy->Fill( ev->ePrimary , float(ev->nPmtHits) );


    // loop over trajectories
    for(unsigned itr=0; itr<ev->traject.size() ; ++itr) {
      LTTraject trj = ev->traject[itr];
      hTrajZ->Fill(trj.Position[0].Z());
      if(trj.Type==4) hTrajHitZ->Fill(trj.Position[0].Z());
    }
  }

  // end of ana 
  outfile->Write();
}
