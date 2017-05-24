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

void anaTrk(TString tag = "2017-4-23-12-42-14")
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
  TH1::SetDefaultSumw2(true); // turn on error bars 	
  double rmax=355;
  double rmin=55;
  TH1F *hRtrk = new TH1F("Rtrk"," radial position ",rmax*5,0,rmax);
  TH1F *hRtrkWls = new TH1F("RtrkWls"," radial position ",rmax*5,0,rmax);

  TH1F *hRScaled = new TH1F("RScaled"," (r/rmax)^2 ",1000,0,1);
  TH1F *hRScaledWls = new TH1F("RScaledWls"," (r/rmax)^2 ",1000,0,1);
  TH1F *hZtrk = new TH1F("Ztrk"," z position ",1700,-210,110);
  TH1F *hZtrkWls = new TH1F("ZtrkWls"," z position ",1700,-210,110);

  int maxsteps = 250;
  TH1F *hEStepElectron =  new TH1F("EStepElectron"," e- step weightd by ke ",maxsteps,0,maxsteps);
  TH1F *hEStepOptical  =  new TH1F("EStepOptical","  optical step weightd by ke ",maxsteps,0,maxsteps);

  vector<double> errorStepElectron(maxsteps+1);// zeroth bin is underflow.
  vector<double> errorStepOptical(maxsteps+1);

 
  // create pointer to branch and set branch address
  LTTrack *trk = new LTTrack();
  trkTree->SetBranchAddress("track",&trk);
  double sumwElectron=0;
  double sumwOptical=0;

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
    if(name == TString("opticalphoton")) {
      hEStepOptical->Fill(float(trk->nstep),trk->ke );
      errorStepOptical[trk->nstep+1] += trk->ke;
    }
    if(name == TString("e-" )) {
      hEStepElectron->Fill(float(trk->nstep),trk->ke ); // weight by energy
      errorStepElectron[trk->nstep+1] += trk->ke;
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

  // take square roots
  for(unsigned i1=0; i1<errorStepElectron.size(); ++i1) errorStepElectron[i1]=sqrt(errorStepElectron[i1]);
  for(unsigned i2=0; i2<errorStepOptical.size(); ++i2) errorStepOptical[i2]=sqrt(errorStepOptical[i2]);



  hEStepElectron->SetError(&errorStepElectron[0]);
  hEStepOptical->SetError(&errorStepOptical[0]);
  // end of ana 
  outfile->Write();
}
