//#include <TLorentzVector.h>
#include <TVector3.h>

enum { SCINT,ABS,OTHER,OUT,CONVERT,BOUNDABS,ABSGEWLS,ABSGESCINT,WLS,PMT,LAST };


enum G4OpBoundaryProcessStatus {  Undefined,
                                  Transmission, FresnelRefraction,
                                  FresnelReflection, TotalInternalReflection,
                                  LambertianReflection, LobeReflection,
                                  SpikeReflection, BackScattering,
                                  Absorption, Detection, NotAtBoundary,
                                  SameMaterial, StepTooSmall, NoRINDEX,
                                  PolishedLumirrorAirReflection,
                                  PolishedLumirrorGlueReflection,
                                  PolishedAirReflection,
                                  PolishedTeflonAirReflection,
                                  PolishedTiOAirReflection,
                                  PolishedTyvekAirReflection,
                                  PolishedVM2000AirReflection,
                                  PolishedVM2000GlueReflection,
                                  EtchedLumirrorAirReflection,
                                  EtchedLumirrorGlueReflection,
                                  EtchedAirReflection,
                                  EtchedTeflonAirReflection,
                                  EtchedTiOAirReflection,
                                  EtchedTyvekAirReflection,
                                  EtchedVM2000AirReflection,
                                  EtchedVM2000GlueReflection,
                                  GroundLumirrorAirReflection,
                                  GroundLumirrorGlueReflection,
                                  GroundAirReflection,
                                  GroundTeflonAirReflection,
                                  GroundTiOAirReflection,
                                  GroundTyvekAirReflection,
                                  GroundVM2000AirReflection,
                                  GroundVM2000GlueReflection,
                                  Dichroic };

// max is 16 bits!
enum TrackStatus { active=1, scint=2, absorbed=4, boundaryAbsorbed=8,
                      absorbedLAr=16, inactive=32, hitWLS = 64, totalInternal=128, backScatter=256, fresnelRefract=512,
                      fresnelReflect = 2*fresnelRefract,
                      spikeReflect=2*fresnelReflect, 
                      hitPMT=2*spikeReflect, 
                      hitGe=2*hitPMT,
                      absGe=2*hitGe,
                      isBad=2*absGe};


double hc=1.23984193E+3; // eV-nm  KE unit is eV

TString bnames[Dichroic+1];

int setBorderNames() 
{
  int in=0;
  bnames[in++]="Undefined";   
  bnames[in++]="Transmission";  
  bnames[in++]="FresnelRefraction";
  bnames[in++]="FresnelReflection";
  bnames[in++]="TotalInternalReflection";
  bnames[in++]= "LambertianReflection";
  bnames[in++]="LobeReflection";
  bnames[in++]="SpikeReflection";
  bnames[in++]="BackScattering";
  bnames[in++]="Absorption";
  bnames[in++]="Detection";
  bnames[in++]="NotAtBoundary";
  bnames[in++]="SameMaterial";
  bnames[in++]="StepTooSmall";
  bnames[in++]="NoRINDEX";
  bnames[in++]="PolishedLumirrorAirReflection";
  bnames[in++]="PolishedLumirrorGlueReflection";
  bnames[in++]="PolishedAirReflection";
  bnames[in++]="PolishedTeflonAirReflection";
  bnames[in++]="PolishedTiOAirReflection";
  bnames[in++]="PolishedTyvekAirReflection";
  bnames[in++]="PolishedVM2000AirReflection";
  bnames[in++]="PolishedVM2000GlueReflection";
  bnames[in++]="EtchedLumirrorAirReflection";
  bnames[in++]="EtchedLumirrorGlueReflection";
  bnames[in++]="EtchedAirReflection";
  bnames[in++]="tchedTeflonAirReflection";
  bnames[in++]="EtchedTiOAirReflection";
  bnames[in++]="EtchedTyvekAirReflection";
  bnames[in++]="EtchedVM2000AirReflection";
  bnames[in++]="EtchedVM2000GlueReflection";
  bnames[in++]="GroundLumirrorAirReflection";
  bnames[in++]="GroundLumirrorGlueReflection";
  bnames[in++]="GroundAirReflection";
  bnames[in++]="GroundTeflonAirReflection";
  bnames[in++]="GroundTiOAirReflection";
  bnames[in++]="GroundTyvekAirReflection";
  bnames[in++]="GroundVM2000AirReflection";
  bnames[in++]="GroundVM2000GlueReflection";
  bnames[in++]="Dichroic";
  return in;
}

std::vector<std::string> trackStatusNames;
std::vector<int> trackStatusValue;


unsigned setTrackStatusNames()
{
  trackStatusValue.clear();
  trackStatusValue.push_back(active);
  trackStatusValue.push_back(scint);
  trackStatusValue.push_back(absorbed);
  trackStatusValue.push_back(boundaryAbsorbed);
  trackStatusValue.push_back(absorbedLAr);
  trackStatusValue.push_back(inactive);
  trackStatusValue.push_back(hitWLS);
  trackStatusValue.push_back(totalInternal);
  trackStatusValue.push_back(backScatter);
  trackStatusValue.push_back(fresnelRefract);
  trackStatusValue.push_back(fresnelReflect);
  trackStatusValue.push_back(spikeReflect);
  trackStatusValue.push_back(hitPMT);
  trackStatusValue.push_back(hitGe);
  trackStatusValue.push_back(absGe);
  trackStatusValue.push_back(isBad);

  trackStatusNames.clear();
  trackStatusNames.push_back("active");
  trackStatusNames.push_back("scint");
  trackStatusNames.push_back("absorbed");
  trackStatusNames.push_back("boundaryAbsorbed");
  trackStatusNames.push_back("absorbedLAr");
  trackStatusNames.push_back("inactive");
  trackStatusNames.push_back("hitWLS");
  trackStatusNames.push_back("totalInternal");
  trackStatusNames.push_back("backScatter");
  trackStatusNames.push_back("fresnelRefract");
  trackStatusNames.push_back("fresnelReflect");
  trackStatusNames.push_back("spikeReflect");
  trackStatusNames.push_back("hitPMT");
  trackStatusNames.push_back("hitGe");
  trackStatusNames.push_back("absGe");
  trackStatusNames.push_back("isBad");
  return  trackStatusNames.size();
}

void calcRatio(double a, double b, double &r, double &e) 
{
  r=0; e=0;
  if(b<1) return;
  r=a/b;
  if(a<1) return;
  // assume errors are root(n)A
  e = r * sqrt( 1./a + 1./b );
}

void printTrackStatus(int in,int status)
{
  printf("track status:  ");
  for (int i=0; i< in; ++i) if( status& (1<<i) )  printf("\t bit %i %s ",i,trackStatusNames[i].c_str());
  printf("\n");
}


void printTrackStatusNames()
{
  for (unsigned i=0; i< trackStatusNames.size() ; ++i) {
    printf(" bit %i hex val %X %s \n",i,trackStatusValue[i],trackStatusNames[i].c_str());
  }
}



void printBorderNames(int in)
{
  for (int i=0; i< in ; ++i) printf(" bit %i %s \n",i,bnames[i].Data());
}


void anaTrk(TString tag = "2017-5-27-11-23-26")
{
  int inames =setBorderNames();
  //printBorderNames(inames);
  unsigned itnames = setTrackStatusNames();
  printTrackStatusNames();
  TString inputFileName = TString("legendTree-")+tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName,"READONLY");
  TTree *trkTree=NULL;
  unsigned aSize = 0;

  // tree has to be in file
  infile->GetObject("trkTree",trkTree);
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
  double rmax= 160;
  double zmax= 160;

  double xoffset = 200;  // shift xrange -350,50 to -150,+150 
  double yoffset = 0;
  double zoffset = 50;

  double xvoffset = 200;  // shift xrange -350,50 to -150,+150 
  double yvoffset = 0;
  double zvoffset = 50;


  TH1::SetDefaultSumw2(true); // turn on error bars

  TH1F *hLambdaG = new TH1F("LambdaG"," absGe photon wavelength ",500,100,600);
  TH1F *hCerenkovBound = new TH1F("CerenkovBound"," last boundary of other ",Dichroic,0,Dichroic);
  
  TH2F *hPMTWLS = new TH2F("PMTWLS"," PMTWLS z versus r ",rmax*5,0,rmax,zmax*2,-zmax,zmax);
  TH2F *hgroup1Physical = new TH2F("group1Physical"," z versus r ",rmax*5,0,rmax,zmax*2,-zmax,zmax);
  TH2F *hgroup1PhysicalXYStart = new TH2F("group1PhysicalXYStart"," y versus x ",zmax*2,-zmax,zmax,zmax*2,-zmax,zmax);
  TH2F *hgroup1PhysicalXYEnd = new TH2F("group1PhysicalXYEnd"," y versus x ",zmax*2,-zmax,zmax,zmax*2,-zmax,zmax);

  TNtuple *ntele = new TNtuple("ntele","ntele","rvert:xvert:yvert:zvert:rele:zele:parent:ke");
  TNtuple *ntopt = new TNtuple("ntopt","ntopt","rvert:zvert:parent:rlog:rtrk:ztrk:lambda:length:time:1fstat:lar:status:boundary");
  TNtuple *ntpos = new TNtuple("ntpos","start-end","rstart:xstart:ystart:zstart:rend:xend:yend:zend:lambda:length:fstat:lar:boundary");

  TH1F *hLambdaLar = new TH1F("LambdaLar"," Ar scint photon absorbed in larPhysical wavelength ",500,100,600);
  TH1F *hLambdaWls = new TH1F("LambdaWls"," WLS photon wavelength ",500,100,600);
  TH1F *hLambdaCerenkov = new TH1F("LambdaCerenkov"," other photon wavelength ",500,100,600);
  TH1F *hLambdaAbsGeScint = new TH1F("LambdaAbsGeScint"," absGe photon wavelength ",500,100,600);
  TH1F *hLambdaAbsGeWls = new TH1F("LambdaAbsGeWls"," absGe photon wavelength ",500,100,600);

  TH1F *hLambdaPmt = new TH1F("LambdaPmt"," pmt detected photon wavelength ",500,100,600);
 

  

  TH2F *hXYvert = new TH2F("XYvert"," primary vertex Y versus X ",zmax/5,-zmax,zmax,zmax/5,-zmax,zmax);
  TH2F *hRZvert = new TH2F("RZvert"," primary vertex z versus r ",rmax/10,0,rmax,zmax/5,-zmax,zmax);
  TH1F *hRvert = new TH1F("Rvert"," primary vertex radial position ",rmax/10,0,rmax);
  TH1F *hZvert = new TH1F("Zvert"," primary vertex z position ",zmax/5,-zmax,zmax);
  

  TH1F *hLambdaScint = new TH1F("LambdaScint"," Ar scint photon absorbed wavelength ",500,100,600);
  TH1F *hRtrkScint = new TH1F("RtrkScint"," Ar Scint radial position ",rmax/10,0,rmax);
  TH1F *hZtrkScint = new TH1F("ZtrkScint"," Ar Scint z position ",zmax/5,-zmax,zmax);

  TH1F *hLambdaNotAt = new TH1F("LambdaNotAt"," not at boundary  ",500,100,600);
   
  
    

  TH1F *hRele = new TH1F("Rele"," electron radial position ",rmax/5,0,rmax);
  TH1F *hZele = new TH1F("Zele"," electron z position ",zmax/5,-zmax,zmax);
  TH1F *hRele0 = new TH1F("Rele0"," electron radial position ",rmax/10,0,rmax);
  TH1F *hZele0 = new TH1F("Zele0"," electron z position ",zmax/5,-zmax,zmax);
  TH1F *hRele1 = new TH1F("Rele1"," electron radial position ",rmax/10,0,rmax);
  TH1F *hZele1 = new TH1F("Zele1"," electron z position ",zmax/5,-zmax,zmax);

  TH1F *hLambdaBoundAbs = new TH1F("LambdaBoundAbs"," photon boundary absorbed ",500,100,600);
  TH1F *hRtrkBoundAbs = new TH1F("RtrkBoundAbs"," photon R boundary absorbed",rmax*5,0,rmax);
  TH1F *hZtrkBoundAbs = new TH1F("ZtrkBoundAbs","  photon Z boundary absorbed",zmax*2,-zmax,zmax);

  TH1F *hLambdaAbs = new TH1F("LambdaAbs"," absorbed wavelength ",500,100,600);
  TH1F *hRtrkAbs = new TH1F("RtrkAbs"," absorbed radial position ",rmax/10,0,rmax);
  TH1F *hZtrkAbs = new TH1F("ZtrkAbs"," absorbed  z position ",zmax/5,-zmax,zmax);
  
  TH1F *hRtrkWlsStart = new TH1F("RtrkWlsStart"," WLS starting radial position ",rmax/10,0,rmax);
  TH1F *hRtrkStart = new TH1F("RtrkStart"," radial start position  ",rmax/10,0,rmax);
  TH1F *hRtrkNotAt = new TH1F("RtrkNotAt"," not at  radial position ",rmax/10,0,rmax);
  TH1F *hZtrkWlsStart = new TH1F("ZtrkWlsStart"," WLS starting z position ",zmax/5,-zmax,zmax);
  TH1F *hZtrkStart = new TH1F("ZtrkStart"," z start position  ",zmax/5,-zmax,zmax);
  TH1F *hZtrkNotAt = new TH1F("ZtrkNotAt"," not at z position ",zmax/5,-zmax,zmax);
  
  hRtrkStart->SetMarkerColor(kBlue); hZtrkStart->SetMarkerColor(kBlue);
  hRtrkStart->SetMarkerStyle(22); hZtrkStart->SetMarkerStyle(22);

  TH1F *hRtrkEnd = new TH1F("RtrkEnd"," radial end position  ",rmax*5,0,rmax);
  TH1F *hZtrkEnd = new TH1F("ZtrkEnd"," z end position  ",zmax*2,-zmax,zmax);
  hRtrkEnd->SetMarkerColor(kRed); hZtrkEnd->SetMarkerColor(kRed);
  hRtrkEnd->SetMarkerStyle(23); hZtrkEnd->SetMarkerStyle(23);
  
  
  TH1F *hRtrkLar = new TH1F("RtrkLar"," radial position absorbed in larPhysical ",400,0,800);
  TH1F *hZtrkLar = new TH1F("ZtrkLar"," z position absorbed larPhysical ",600,-600,600);
  
  TH1F *hRtrk = new TH1F("Rtrk"," radial position ",rmax*5,0,rmax);
  TH1F *hRtrkWls = new TH1F("RtrkWls"," radial position ",rmax*5,0,rmax);

  TH1F *hRtrkPmt = new TH1F("RtrkPmt"," radial position ",rmax*5,0,rmax);

  TH1F *hRtrkCerenkov = new TH1F("RtrkCerenkov"," other radial position ",rmax*5,0,rmax);
  TH1F *hRScaledCerenkov = new TH1F("RScaledCerenkov"," (r/rmax)^2 ",1000,0,1);
  TH1F *hZtrkCerenkov = new TH1F("ZtrkCerenkov"," other z position ",zmax*2,-zmax,zmax);

  TH1F *hRScaledAbs = new TH1F("RScaledAbs"," (r/rmax)^2 ",1000,0,1);

  TH1F *hRtrkAbsGe = new TH1F("RtrkAbsGe"," other radial position ",rmax*5,0,rmax);
  TH1F *hRScaledAbsGe = new TH1F("RScaledAbsGe"," (r/rmax)^2 ",1000,0,1);
  TH1F *hZtrkAbsGe = new TH1F("ZtrkAbsGe"," other z position ",zmax*2,-zmax,zmax);

  

  TH1F *hRScaled = new TH1F("RScaled"," (r/rmax)^2 ",1000,0,1);
  TH1F *hRScaledWls = new TH1F("RScaledWls"," (r/rmax)^2 ",1000,0,1);
  TH1F *hZtrk = new TH1F("Ztrk"," z position ",zmax*2,-zmax,zmax);
  TH1F *hZtrkWls = new TH1F("ZtrkWls"," WLS z position ",zmax*2,-zmax,zmax);
  TH1F *hZtrkPmt = new TH1F("ZtrkPmt"," PMT z position ",zmax*2,-zmax,zmax);

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
  int barray[8]; // 64 bits
  int event =-1;
  int nEvents =0;
  for(unsigned entry =0; entry < aSize  ; ++entry ) {
    trkTree->GetEntry(entry);
    unsigned  nhistory = trk->positionHistory.size();
    if(trk->evId != event) {
      ++nEvents;
      event=trk->evId;
    }
    // what are the tracks with nhistory=0?
    if(nhistory<1) continue;
    std::vector<int> boundaryStatus = trk->boundaryStatus;
    int boundaryLast=0;
    if(boundaryStatus.size()>0) boundaryLast = boundaryStatus[boundaryStatus.size()-1];
    TString name = trk->particleName;
    TString volName = trk->physVolName.Data();
   
    bool printOut=false;
    if(entry%1000000==0) printOut=true;
    if(entry%1000000==0) printf("\t entry %i %s %s %s status %X boundary size  %lu \n",
        entry,trk->physVolName.Data(), trk->process.Data(), trk->particleName.Data(),trk->status,boundaryStatus.size());
    TVector3 position = trk->position;
    TVector3 vposition = trk->vertPosition;

    
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
      if(obin+1<errorStepOptical.size()) errorStepOptical[obin+1] += trk->edep;
      hElectronStepEloss->Fill(trk->edep);
    }
    if ((name == TString("e-" )) && trk->edep > 0) {
      int ebin = hEStepElectron->FindBin(ystep);
      hEStepElectron->SetBinContent(ebin,hEStepElectron->GetBinContent(ebin)+trk->edep); // weight by energy
      hElectronStepLength->Fill(trk->stepLength[trk->stepLength.size()-1]);
      if(ebin+1<errorStepElectron.size()) errorStepElectron[ebin+1] += fanoLAr*trk->edep;
      else printf(" ebin %u is bigger than vector size %u \n",ebin+1, (unsigned) errorStepElectron.size());
      /*
      printf(" trkid %i parent %i step %i length %f steplength %f length-steplength %f ke %E edep %E ke+edep status %X eIoni %i \n",
        trk->trkId,trk->parentId, trk->nstep,trk->length, trk->stepLength,trk->length - trk->stepLength,  
        trk->ke,  trk->edep,  trk->ke + trk->edep,trk->status,isEIoni); 
        */
      
    }
    int status = trk->status;
    double xshift = position.X()+xoffset;
    double yshift = position.Y()+yoffset;
    double zshift = position.Z()+zoffset;
    double radius = sqrt(xshift*xshift+yshift*yshift);
    double rscale = pow( radius/rmax,2);

    double xvert = vposition.X()+xvoffset;
    double yvert = vposition.Y()+yvoffset;
    double zvert = vposition.Z()+zvoffset;
    double vradius = sqrt(xvert*xvert+yvert*yvert);

  
    if(name == TString("e-")) {
      ntele->Fill(vradius,xvert,yvert,zvert,radius,zshift,double(trk->parentId),trk->ke);
      hRele->Fill(radius);
      hZele->Fill(zshift);
      if(trk->parentId==0) {
        hRele0->Fill(radius);
        hZele0->Fill(zshift);
      } else {
        hRele1->Fill(radius);
        hZele1->Fill(zshift);
      }
    }

    // look at optical photons only
    if(name != TString("opticalphoton"))  continue;

    bool trkback = status&backScatter;


    //if( trkback!=boundback ) printf(" %i %i boundary %X status %X\n",boundback,trkback, trk->boundaryStatus,status);

   // if(trk->boundaryStatus==11) printf(" boundary %X status %X\n",trk->boundaryStatus,status);
   // if(status&backScatter)      printf(" is backscatter boundary %X status %X\n",trk->boundaryStatus,status);
    // remove backscatter
    //if(status&backScatter) continue;

    double lambda = hc/trk->ke;
    double flar = 0;
    if(volName==TString("larPhysical")) flar=1;
    double fstat = double(status);
    /*
    double fstat= 0; 
    if(status&totalInternal) fstat = double(totalInternal);
    if(status&fresnelReflect) fstat = double(fresnelReflect);
    if(status&fresnelRefract) fstat = double(fresnelRefract);
    if(status&spikeReflect) fstat = double(spikeReflect);
    if(status&absorbedLAr) fstat = double(absorbedLAr);
    if(status&boundaryAbsorbed) fstat = double(boundaryAbsorbed);
    if(status&absorbed) fstat = double(absorbed);
    if(status&hitWLS) fstat = double(hitWLS);
    if(status&hitPMT) fstat  = double(hitPMT);
    if(status&absGe) fstat  = double(absGe);
    */
    
    if(status&absGe) hLambdaG->Fill(lambda);

    ntopt->Fill(vradius,zvert,double(trk->parentId),position.Perp(),radius,zshift,
        lambda,trk->length,trk->time,fstat,flar,double(status),double(boundaryLast));

  // start of track
    double xstart=0;
    double ystart=0;
    double zstart=0;
    double rstart=0;
    double xend=0;
    double yend=0;
    double zend=0;
    double rend=0;
    
    if(nhistory>0) {
      TVector3 pstart = trk->positionHistory[0];
      xstart = pstart.X()+xoffset;
      ystart = pstart.Y()+yoffset;
      zstart = pstart.Z()+zoffset;
      rstart = sqrt(xstart*xstart+ystart*ystart);
      hRtrkStart->Fill(rstart);
      hZtrkStart->Fill(zstart);

      TVector3 pend = trk->positionHistory[nhistory-1];
      xend = pend.X()+xoffset;
      yend = pend.Y()+yoffset;
      zend = pend.Z()+zoffset;
      rend = sqrt(xend*xend+yend*yend);
      hRtrkEnd->Fill(rend);
      hZtrkEnd->Fill(zend);
      ntpos->Fill(rstart,xstart,ystart,zstart,rend,xend,yend,zend,lambda,trk->length,fstat,flar,double(boundaryLast));
    } else
      printf(" \t >>>>>> nhistory = %u \n",nhistory);
       
    hRvert->Fill(vradius);
    hZvert->Fill(zvert);
    hRZvert->Fill(vradius,zvert);
    hXYvert->Fill(xvert,yvert);
    

    if(status&scint) {
     hRtrkScint->Fill(radius);
     hZtrkScint->Fill(zshift);
     hLambdaScint->Fill(lambda);
    }
    

    if(volName==TString("larPhysical")) {
      hRtrkLar->Fill(radius);
      hZtrkLar->Fill(zshift);
      hLambdaLar->Fill(lambda);
    }
    if(status&boundaryAbsorbed) {
      hRtrkBoundAbs->Fill(radius);
      hZtrkBoundAbs->Fill(zshift);
      hLambdaBoundAbs->Fill(lambda);  
    }

    if(trk->parentId>1&&printOut) {
      //printf(" *** Ge status %i pre/post %s / %s \n",status,trk->preName.Data(),trk->postName.Data());
      printTrackStatus(itnames,status);
      unsigned ilast =  boundaryStatus.size()-1;
      printf("\t %i %i %s  pre %s   \n",ilast, boundaryStatus[ilast],bnames[boundaryStatus[ilast]].Data(),trk->boundaryName[ilast].c_str());
      //printf("     boundary size %lu \n", boundaryStatus.size() );
      //for(unsigned ib=0; ib< boundaryStatus.size(); ++ib ) 
      printf("\n");
    }

    if(trk->boundaryName[trk->boundaryName.size()-1].find("phys_PMTWLS")!=std::string::npos) hPMTWLS->Fill(radius,zshift);
    if(trk->boundaryName[trk->boundaryName.size()-1]==TString("group1Physical")) hgroup1Physical->Fill(radius,zshift);
    if(trk->boundaryName[trk->boundaryName.size()-1]==TString("group1Physical")) hgroup1PhysicalXYStart->Fill(xstart,ystart);
    if(trk->boundaryName[trk->boundaryName.size()-1]==TString("group1Physical")) hgroup1PhysicalXYEnd->Fill(xend,yend);

    if(volName!=TString("larPhysical")) {
      //if(!(status&hitWLS) ) hLambdaScint->Fill(lambda);
      //if(!(status&scint&&status&absorbed)) {
      if(status&hitPMT) {
        hRtrkPmt->Fill(radius);
        hZtrkPmt->Fill(zshift);
        hLambdaPmt->Fill(lambda);
     } else if(status&absGe) {
      if(0) {
          printf(" *** absorbed GE status %i pre/post %s / %s \n",status,trk->preName.Data(),trk->postName.Data());
          printTrackStatus(itnames,status);
          printf("     GEGEGEGEGE boundary size %lu \n", boundaryStatus.size() );
          for(unsigned ib=0; ib< boundaryStatus.size(); ++ib ) 
            printf("\t %i %i %s  ",ib, boundaryStatus[ib],bnames[boundaryStatus[ib]].Data());
          printf("\n");
        }
        
       
       if(status&hitWLS) hLambdaAbsGeWls->Fill(lambda);
       else hLambdaAbsGeScint->Fill(lambda);
       hRtrkAbsGe->Fill(radius);
       hRScaledAbsGe->Fill(rscale);
       hZtrkAbsGe->Fill( zshift);
     } else if(status&hitWLS) {
       hRtrkWls->Fill(radius);
       hRScaledWls->Fill(rscale) ;
       hZtrkWls->Fill(zshift);
       hLambdaWls->Fill(lambda);
       hRtrkWlsStart->Fill(rstart);
       hZtrkWlsStart->Fill(zstart);
       
     } else if(status&absorbed) {  // absorbed inside cylinder
       if(0) {
          printf(" *** absorbed track status %i pre/post %s / %s \n",status,trk->preName.Data(),trk->postName.Data());
          printTrackStatus(itnames,status);
          printf("     boundary size %lu \n", boundaryStatus.size() );
          for(unsigned ib=0; ib< boundaryStatus.size(); ++ib ) 
            printf("\t %i %i %s   \n",ib, boundaryStatus[ib],bnames[boundaryStatus[ib]].Data());
          printf("\n");
        }
        
        hLambdaAbs->Fill(lambda);
        hRtrkAbs->Fill(radius);
        hRScaledAbs->Fill(rscale);
        hZtrkAbs->Fill( zshift);
      } else if(boundaryLast==NotAtBoundary) {
        hLambdaNotAt->Fill(lambda);
        hRtrkNotAt->Fill(radius);
        hZtrkNotAt->Fill( zshift);
        if(0) {
          printf(" *** reflect track process %s status %i pre/post %s / %s ",trk->process.Data(),status,trk->preName.Data(),trk->postName.Data());
          if(status&spikeReflect) printf(" spike ");
          else if(status&totalInternal) printf(" total ");
          else printf(" what? ");
          printf("\n");
          printTrackStatus(itnames,status);
         
          printf("  reflectreflectreflect boundary size %lu ", boundaryStatus.size() );
          for(unsigned ib=0; ib< boundaryStatus.size(); ++ib ) 
            printf("\t %i %i %s  %s   ",ib, boundaryStatus[ib],bnames[boundaryStatus[ib]].Data(),trk->boundaryName[ib].c_str());
          printf(" \n");
        }
      } else if(boundaryLast==Absorption&&trk->process.Contains("Cerenkov") )  {
        hLambdaCerenkov->Fill(lambda);
        hRtrkCerenkov->Fill(radius);
        hRScaledCerenkov->Fill(rscale);
        hZtrkCerenkov->Fill( zshift);
        hCerenkovBound->Fill(double( boundaryStatus[ boundaryStatus.size()-1 ] ));
     } else {
       if(1) {
          printf(" *** other track process %s status %i pre/post %s / %s ",trk->process.Data(),status,trk->preName.Data(),trk->postName.Data());
          if(status&spikeReflect) printf(" spike ");
          else if(status&totalInternal) printf(" total ");
          else printf(" what? ");
          printf("\n");
          printTrackStatus(itnames,status);
         
          printf("  otherotherother boundary size %lu ", boundaryStatus.size() );
          for(unsigned ib=0; ib< boundaryStatus.size(); ++ib ) 
            printf("\t %i %i %s  %s   ",ib, boundaryStatus[ib],bnames[boundaryStatus[ib]].Data(),trk->boundaryName[ib].c_str());
          printf(" \n");
        }
     }
    }
  }

    // take square roots
  for(unsigned i1=0; i1<errorStepElectron.size(); ++i1) errorStepElectron[i1]=sqrt(errorStepElectron[i1]);
  for(unsigned i2=0; i2<errorStepOptical.size(); ++i2) errorStepOptical[i2]=sqrt(errorStepOptical[i2]);
     
  hLambdaCerenkov->SetLineColor(kGreen);    
  hRtrkCerenkov->SetLineColor(kGreen); 
  hRScaledCerenkov->SetLineColor(kGreen); 
  hZtrkCerenkov->SetLineColor(kGreen); 

  hEStepElectron->SetError(&errorStepElectron[0]);
  hEStepOptical->SetError(&errorStepOptical[0]);

  TCanvas* cele = new TCanvas("stepElectron","stepElectron");
  //hEStepElectron->Fit("expo");
  hEStepElectron->Draw();


  TCanvas* copt = new TCanvas("stepOptical","stepOptical");
  //hEStepOptical->Fit("expo");
  hEStepOptical->Draw();

  TCanvas* crtrkLar= new TCanvas("RtrkLar","RtrkLar");
  hRtrkLar->Draw();

  TCanvas* cztrkLar= new TCanvas("ZtrkLar","ZtrkLar");
  hZtrkLar->Draw(); 
  
  hRtrkScint->SetMarkerColor(kGreen);
  hZtrkScint->SetMarkerColor(kGreen);
  hLambdaScint->SetMarkerColor(kGreen);
  hRtrkScint->SetLineColor(kGreen);
  hZtrkScint->SetLineColor(kGreen);
  hLambdaScint->SetLineColor(kGreen);
  
  hRtrkWls->SetLineColor(kRed); hRtrkWls->SetMarkerColor(kRed);
  hZtrkWls->SetLineColor(kRed); hZtrkWls->SetMarkerColor(kRed);
  hLambdaWls->SetLineColor(kRed);hLambdaWls->SetMarkerColor(kRed);

 
  hRtrkAbs->SetLineColor(kBlue); hRtrkAbs->SetMarkerColor(kBlue);
  hZtrkAbs->SetLineColor(kBlue); hZtrkAbs->SetMarkerColor(kBlue);
  hLambdaAbs->SetLineColor(kBlue);hLambdaAbs->SetMarkerColor(kBlue);

  hRtrkAbsGe->SetLineColor(kOrange); hRtrkAbs->SetMarkerColor(kOrange);
  hZtrkAbsGe->SetLineColor(kOrange); hZtrkAbs->SetMarkerColor(kOrange);
  hLambdaAbsGeScint->SetLineColor(kOrange);hLambdaAbs->SetMarkerColor(kOrange);
  
  
  gStyle->SetOptStat(0);
  gStyle->SetOptLogy();

  TCanvas* crtrk= new TCanvas("Rtrk","Rtrk");
  hRtrkScint->Draw();  
  hRtrkWls->Draw("sames"); hRtrkAbs->Draw("sames");
  hRtrkAbsGe->Draw("sames");
  crtrk->Print(".pdf");
  

  TCanvas* cztrk = new TCanvas("Ztrk","Ztrk");
  hZtrkWls->Draw(); hZtrkScint->Draw("sames"); hZtrkAbs->Draw("sames"); //hZtrkCerenkov->Draw("sames");
  hZtrkAbsGe->Draw("sames");
  cztrk->Print(".pdf");


  TCanvas* clambda = new TCanvas("lambda","lambda");
  hLambdaScint->Draw();  hLambdaWls->Draw("sames"); hLambdaAbs->Draw("sames"); //hLambdaCerenkov->Draw("sames");
  hLambdaAbsGeScint->Draw("sames");
  hLambdaAbsGeWls->Draw("sames");
  clambda->Print(".pdf");
  
 
  // report status words
  printf("\t fstat values: absorbed %i hitWLS %i hitPMT %i \n",absorbed,hitWLS,hitPMT);
  

  double count[LAST];
  vector<TString> cname;
  cname.push_back("        all scint");
  cname.push_back("         absorbed");
  cname.push_back("     abs Cerenkov");
  cname.push_back("          outside");
  cname.push_back("     WLS  convert");
  cname.push_back("     WLS absorbed");
  cname.push_back("        absGe WLS");
  cname.push_back("      absGe scint");
  cname.push_back("              WLS");
  cname.push_back("              PMT");

  count[SCINT]  = hLambdaScint->Integral();
  count[BOUNDABS] = hLambdaBoundAbs->Integral();
  count[ABS] = hLambdaAbs->Integral();
  count[WLS]  = hLambdaWls->Integral();
  count[PMT] = hLambdaPmt->Integral();
  count[ABSGEWLS] = hLambdaAbsGeWls->Integral();
  count[ABSGESCINT] = hLambdaAbsGeScint->Integral();
  count[CONVERT] = hLambdaNotAt->Integral();
  count[OUT] = hLambdaLar->Integral();
  count[OTHER] = hLambdaCerenkov->Integral();
  

  printf("\n\tphoton totals: %s  # events = %i \n ",tag.Data(),nEvents);

  printf(" %s %.0f \n",cname[SCINT].Data(),count[SCINT]);
  double notScint = 0;
  for(int ic=SCINT+1; ic < LAST; ++ic) {
    notScint += count[ic];
    double ratio=0,eratio=0;
    calcRatio(count[ic],count[SCINT],ratio,eratio);
    printf(" %s %.0f  percent scint = %0.3f +/- %0.3f  \n",cname[ic].Data(),count[ic],ratio*100,eratio*100);
  }
  double cratio=0, ecratio=0;
  double wlsTotal = count[BOUNDABS]+count[WLS]+count[ABSGEWLS];
  calcRatio(wlsTotal,count[CONVERT],cratio,ecratio);
  printf(" fraction WLS of converted = %0.3f +/- %0.3f  \n",cratio,ecratio);
  
  printf(" sum of all not scint is %i \n", int(notScint));

  TH1F* hZtrkDivide = (TH1F*) hZtrkWlsStart->Clone("ZtrkDivide");
  hZtrkDivide->Divide(hZtrkNotAt);

  // end of ana 
  outfile->Write();

}
