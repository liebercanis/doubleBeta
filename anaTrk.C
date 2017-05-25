#include <vector>
#include <TLorentzVector.h>
#include <TVector3.h>

enum TrackStatus { active=1, hitPMT=2, absorbed=4, boundaryAbsorbed=8,
                      absorbedLAr=16, inactive=32, hitWLS = 64, totalInternal=128, backScatter=256, notBoundary=512,
                      scint=2*notBoundary, 
                      eIoni=2*scint, 
                      hIoni=2*eIoni, 
                      ionIoni=2*hIoni,
                      compton= 2*ionIoni,
                      hitGe=2*compton, 
                      isBad=2*hitGe};

void anaTrk(TString tag = "2017-4-25-11-29-5")
{
  TString inputFileName = TString("legendTree-")+tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
  TTree *trkTree=NULL;
  unsigned aSize = 0;

  // tree has to be in file
  trkTree = (TTree*) infile->Get("trkTree");
  if(!trkTree) {
    printf(" trkTree not found! \n");
    return;
  }
  if(trkTree) aSize=trkTree->GetEntriesFast();
  printf(" trkTree with %i entries \n",int(aSize));

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("anaTrk-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());
  // add histograms here
  double rmin=55;
  double rmax=355;
  TH1::SetDefaultSumw2(true); // turn on error bars 	
  TH1F *hRtrk = new TH1F("Rtrk"," radial position ",rmax*5,0,rmax);
  TH1F *hRtrkWls = new TH1F("RtrkWls"," radial position ",rmax*5,0,rmax);

  TH1F *hRScaled = new TH1F("RScaled"," (r/rmax)^2 ",1000,0,1);
  TH1F *hRScaledWls = new TH1F("RScaledWls"," (r/rmax)^2 ",1000,0,1);
  TH1F *hZtrk = new TH1F("Ztrk"," z position ",1700,-210,110);
  TH1F *hZtrkWls = new TH1F("ZtrkWls"," z position ",1700,-210,110);

  TH1F *hElectronStepLength =  new TH1F("ElectronStepLength"," e- step length (mm) ",1000,0,0.2);
  TH1F *hElectronStepEloss =  new TH1F("ElectronStepEloss"," e- step deposited energy (ev) ",200,0,20);
  hElectronStepEloss->SetXTitle(" electron energy loss per step (eV) ");
  TH1F *hEStepElectron =  new TH1F("EStepElectron"," e- step weighted by ke ",1000,0,20); 
  hEStepElectron->SetXTitle(" length (mm) "); 
  hEStepElectron->SetYTitle(" dE/dlength eV/mm "); 
  TH1F *hEStepOptical  =  new TH1F("EStepOptical","  optical step weighted by ke ",1000,0,1000);
  hEStepOptical->SetXTitle(" length (mm) "); 
  hEStepOptical->SetYTitle(" dE/dlength eV/mm "); 

  vector<double> errorStepElectron(hEStepElectron->GetNbinsX()+1);// zeroth bin is underflow.
  vector<double> errorStepOptical(hEStepOptical->GetNbinsX()+1);
  for(unsigned i1=0; i1<errorStepElectron.size(); ++i1) errorStepElectron[i1]=0;
  for(unsigned i2=0; i2<errorStepOptical.size(); ++i2) errorStepOptical[i2]=0;

  // create pointer to branch and set branch address
  LTTrack *trk = new LTTrack();
  trkTree->SetBranchAddress("track",&trk);
  double sumwElectron=0;
  double sumwOptical=0;
  double fanoLAr = 0.1;

  //TString theName("");
  for(unsigned entry =0; entry < aSize ; ++entry ) {
    trkTree->GetEntry(entry);
    if(entry%100000==0) printf("\t entry %i \n",entry);
    TVector3 position = trk->position;
    TString name = trk->traject.Name;
    /* for printing some names
    if(theName!=name)  {
      printf(" trajectory name is %s step %i \n",name.Data(),trk->nstep);
      theName=name;
    }
    */
    // default error is  sqrt(sum of squares of weight)  we want  sqrt(sum of weight) assuming error on ke is sqrt(ke)
    double ystep = trk->length; // appears to be accumulate length of track to this step.
    if( (name == TString("opticalphoton")) && trk->edep >0 ) {
      int obin = hEStepOptical->FindBin(ystep);
      hEStepOptical->SetBinContent(obin,hEStepOptical->GetBinContent(obin)+trk->edep); // weight by energy
      errorStepOptical[obin+1] += trk->edep;
      hElectronStepEloss->Fill(trk->edep);
    }
    if ((name == TString("e-" )) && trk->edep > 0) {
      int ebin = hEStepElectron->FindBin(ystep);
      hEStepElectron->SetBinContent(ebin,hEStepElectron->GetBinContent(ebin)+trk->edep); // weight by energy
      hElectronStepLength->Fill(trk->stepLength);
      if(ebin+1<errorStepElectron.size()) errorStepElectron[ebin+1] += fanoLAr*trk->edep;
      else printf(" ebin %u is bigger than vector size %u \n",ebin+1, errorStepElectron.size());
      bool isEIoni = trk->status & eIoni;
      /*printf("step %i length %f steplength %f length-steplength %f ke %E edep %E ke+edep status %X eIoni %i \n",
         trk->nstep,trk->length, trk->stepLength,trk->length - trk->stepLength,  trk->ke,  trk->edep,  trk->ke + trk->edep,trk->status,isEIoni); 
         */
      
    }
    int status = trk->status;
    double radius = sqrt(position.X()*position.X()+position.Y()*position.Y());
    double rscale = pow( (radius-rmin)/(rmax-rmin),2);
    if(status&hitWLS) {
      hRtrkWls->Fill(radius);
      hRScaledWls->Fill(rscale) ;
      hZtrkWls->Fill( position.Z());
    } else {
      hRtrk->Fill(radius);
      hRScaled->Fill(rscale);
      hZtrk->Fill( position.Z());
    }
  }

    // take square rootsi
  for(unsigned i1=0; i1<errorStepElectron.size(); ++i1) errorStepElectron[i1]=sqrt(errorStepElectron[i1]);
  for(unsigned i2=0; i2<errorStepOptical.size(); ++i2) errorStepOptical[i2]=sqrt(errorStepOptical[i2]);


  hEStepElectron->SetError(&errorStepElectron[0]);
  hEStepOptical->SetError(&errorStepOptical[0]);

  TCanvas* cele = new TCanvas("stepElectron","stepElectron");
  hEStepElectron->Fit("expo");
  hEStepElectron->Draw();


  TCanvas* copt = new TCanvas("stepOptical","stepOptical");
  hEStepOptical->Fit("expo");
  hEStepOptical->Draw();


  // end of ana 
  outfile->Write();
}
