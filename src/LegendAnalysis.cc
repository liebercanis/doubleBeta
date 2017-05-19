/// M.G. started with exgps example 
// $Id: HistoManager.cc 83882 2014-09-22 11:09:30Z maire $
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "LegendAnalysis.hh"
#include "LegendTrajectory.hh"
#include "G4RunManager.hh"
#include "UserTrackInformation.hh"
#include "G4VProcess.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "UserEventInformation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LegendAnalysis* LegendAnalysis::fLegendAnalysis=NULL;

LegendAnalysis* LegendAnalysis::Instance() {
  if (! fLegendAnalysis) fLegendAnalysis = new LegendAnalysis();
  return fLegendAnalysis;
}

void LegendAnalysis::Initialize()
{
 // open new ouput file with time stamp.
  time_t tnow;
  time(&tnow);

  struct tm* timeinfo;
  timeinfo = localtime(&tnow);

  char chtime[80];
  sprintf(chtime,"%u-%u-%u-%u-%u-%u",timeinfo->tm_year-100+2000,timeinfo->tm_mon,timeinfo->tm_mday,
      timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
  G4String fTreeFileName = G4String("legendTree-") + G4String(chtime) + G4String(".root");
  fTreeFile=new TFile(fTreeFileName.data(),"RECREATE");
  G4String fHistFileName = G4String("legendHist-") + G4String(chtime) + G4String(".root");
  fHistFile=new TFile(fHistFileName.data(),"RECREATE");


  G4String gmess= G4String(" LegendAnalysis: ************  output tree file  ") + fTreeFileName + G4String(" hist file ") + fHistFileName + G4String(" ************ ");
  G4cout << gmess << G4endl;
  fHistFile->cd();
  hOptical = new TH1F("Optical"," optical photons (nm) ",900,100,1000);
  hOptical->GetYaxis()->SetTitle(" photons/nm ");
  hOptical->GetXaxis()->SetTitle(" wavelength (nm) ");

  hWls = new TH1F("Wls"," Wls photons (nm) ",900,100,1000);
  hWls->GetYaxis()->SetTitle(" photons/nm ");
  hWls->GetXaxis()->SetTitle(" wavelength (nm) ");
  

  hPmtHits = new TH1F("PmtHits"," pmt photons (nm) ",900,100,1000);
  hPmtHits->GetYaxis()->SetTitle(" photons/nm ");
  hPmtHits->GetXaxis()->SetTitle(" wavelength (nm) ");
  

  hEElectron = new TH1F("EElectron"," electrons ",1000,0,1000);
  hEElectron->GetYaxis()->SetTitle(" electrons/KeV ");
  hEElectron->GetXaxis()->SetTitle(" energy (KeV) ");
  
  hEGamma = new TH1F("EGamma"," gammas ",1000,0,1000);
  hEGamma->GetYaxis()->SetTitle(" gamma/KeV ");
  hEGamma->GetXaxis()->SetTitle(" energy (KeV) ");

  topHistDir()->ls();

  fTreeFile->cd();

  // make tree in output file
  fTree = new TTree("LTree","LTree");
  //fTree->SetMaxTreeSize(1000000);
  fTree->SetMaxTreeSize(25);
  fEvent = new LTEvent();
  fTree->Branch("event",&fEvent);
  topTreeDir()->ls();

}
    
void  LegendAnalysis::anaEvent( const G4Event *anEvent)
{
  //printf(" ++++++++++++++++++++++++++++++++++++++++++++++++ \n");
  //printf("      LegendAnalysis:anaEvent called \n");
  //printf(" ++++++++++++++++++++++++++++++++++++++++++++++++ \n");
  if(!anEvent) {
    G4cout << " LegendAnalysis called with NULL G4Event pointer!!! " << G4endl;
    return;
  }
  UserEventInformation* eventInformation = (UserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();
  fEvent->evId = anEvent->GetEventID();
  fEvent->nPVert = anEvent->GetNumberOfPrimaryVertex();
  for(int iv=0; iv < fEvent->nPVert; ++iv) {
    G4PrimaryVertex* pvi = anEvent->GetPrimaryVertex(iv);
    LTPVertex lpv;
    G4VPhysicalVolume *physVol = eventInformation->GetPrimaryPhysicalVolume();
    
    lpv.physVolumeName = physVol->GetName();
    // fill primary vertex info
    lpv.VertexId=iv;
    lpv.Position.SetXYZT(pvi->GetX0(),pvi->GetY0(),pvi->GetZ0(),pvi->GetT0());
    G4ThreeVector pvrel = pvi->GetPosition() - physVol->GetTranslation();
    lpv.RelPosition.SetXYZT(pvrel.x(),pvrel.y(),pvrel.z(),pvi->GetT0());
    //G4VUserPrimaryVertexInformation* pvinfo= pvi->GetUserInformation() 
    lpv.nParticles =pvi->GetNumberOfParticle();
    for(int ip=0; ip< lpv.nParticles; ++ip ) {
      G4PrimaryParticle* pvPart = pvi->GetPrimary(ip); 
      LTParticle part;
      part.TrackId= pvPart->GetTrackID();
      part.VertexId=iv;
      part.PDG=     pvPart->GetPDGcode();
      part.Mass=    pvPart->GetMass();
      part.Charge=  pvPart->GetCharge();
      part.KEnergy= pvPart->GetKineticEnergy();
      part.Momentum.SetPxPyPzE( pvPart->GetPx(),pvPart->GetPy(),pvPart->GetPz(),pvPart->GetTotalEnergy());
      //part.print(ip);
      lpv.particle.push_back(part);
      if(iv==0&&part.TrackId==1) {
        fEvent->PDG=part.PDG;
        fEvent->ePrimary = pvPart->GetKineticEnergy();
      }
    }
    //lpv.print(fEvent->evId);
    fEvent->pvertex.push_back(lpv);
  }
  anaTrajectories( anEvent->GetTrajectoryContainer());

  // and end of analysis save this event
  //fEvent->print();
  fTree->Fill();
  //printf(" +++++++++++++++++++ Leaving Legend Analysis +++++++++++++++++++++++++++++ \n");
}
   
void  LegendAnalysis::anaTrajectories(G4TrajectoryContainer* trajectoryContainer)
{
  if(trajectoryContainer==NULL) {
    G4cout << " LegendAnalysis:: anaTrajectories called with null container " << G4endl;
    return;
  }
  fEvent->nTraj = trajectoryContainer->entries();
  for(int ij=0; ij < fEvent->nTraj; ++ij) {
    LegendTrajectory *gtrj = dynamic_cast<LegendTrajectory*>((*(trajectoryContainer))[ij]);
    const G4Track* aTrack= gtrj->GetTrack();

     if(!aTrack) {
      G4cout << "  LegendAnalysis:: no track for this trajectory  " << ij << endl;
      continue;
    }

    UserTrackInformation* trackInformation=(UserTrackInformation*) aTrack->GetUserInformation();
    if(!trackInformation) {
      continue;
      G4cout << " LegendAnalysis:: track definition is  " << aTrack->GetDefinition()->GetParticleName()<< " no user info " << endl;
    }
    
    // find creator process 
    const G4VProcess* creator=aTrack->GetCreatorProcess();
    //if(!creator) 
      //G4cout << " LegendAnalysis:: track definition is  " << aTrack->GetDefinition()->GetParticleName()<< " no creator " << endl;

    LTTraject ltraj;
   // fill from trajectory 
    ltraj.TrajId = gtrj->GetTrackID() ;  
    ltraj.ParentId = gtrj->GetParentID();        
    //PrimaryId = gtrj-> ;        
    ltraj.PDG = gtrj->GetPDGEncoding();       
    //Mass = gtrj-> ;   
    ltraj.Charge = gtrj->GetCharge(); 

    // each element in std::vector corresponds to a point on the particle path
    // typdef CLHEP::Hep3Vector G4ThreeVector;
    G4ThreeVector gpos0 = gtrj->GetPoint(0)->GetPosition();
    TLorentzVector r4(gpos0.x(),gpos0.y(),gpos0.z(),0);// were is the time?
    
    // G4int npoints = gtrj->GetPointEntries();
    /*for(G4int ip=0; ip<npoints ; ++ip ) {
      G4VTrajectoryPoint* gtrjp = gtrj->GetPoint(ip);
      G4ThreeVector GetPosition() 
    }*/

    ltraj.Position.push_back(r4);
    G4ThreeVector momentum3 = gtrj->GetInitialMomentum();
    TVector3 p3(momentum3.x(),momentum3.y(),momentum3.z());
    ltraj.KE = p3.Mag();// don't have mass
    ltraj.Momentum.push_back(p3);

    // set primary 
    if(gtrj->IsPrimary()) {
      ltraj.Type = LTTrajectType::PRI;
      /*
      if(creator) 
        //G4cout << " LegendAnalysis primary creator process " << creator->GetProcessName() 
        //<< " track status " << trackInformation->GetTrackStatus() <<  " PRIMARY TRACK track definition is  " << aTrack->GetDefinition()->GetParticleName()
        //<< "  KE=" << ltraj.KE << G4endl;
      if(!creator)  
       G4cout << " LegendAnalysis primary no creator process " 
        << " track status " << trackInformation->GetTrackStatus() <<  " PRIMARY TRACK track definition is  " << aTrack->GetDefinition()->GetParticleName()
        << "  KE=" << ltraj.KE << G4endl;
        */
    }

    //G4cout <<  " ANALYSIS id " <<   ltraj.TrajId  << " parent "  << ltraj.ParentId << " PDG "  <<  ltraj.PDG << " charge " <<
    //  ltraj.Charge << G4endl;
    if(ltraj.Charge==0) {
      if(gtrj->GetParticleName()=="opticalphoton") fEvent->eOptical += ltraj.KE;
      else fEvent->eNeutral += ltraj.KE;
    } else {
      //gtrj->ShowTrajectory();
      fEvent->eCharged += ltraj.KE;
    }

    //G4cout << " ANALYSIS E charged " << fEvent->eCharged << " neutral  " <<  fEvent->eNeutral << " optical  " << fEvent->eOptical << G4endl;
    
    ltraj.Name = TString(gtrj->GetParticleName().data());
    
    //enum LTTrajectType {UNK,PRI,SCI,WLS,HIT};
    ltraj.Type = LTTrajectType::UNK;
    if(gtrj->GetParticleName()=="opticalphoton") {
      ++fEvent->nTrajOptPhotons;
      G4double waveLength =  h_Planck*c_light/ltraj.KE/nm;//700 nm
      if(gtrj->IsPmtHit()) {
        ++fEvent->nTrajPmtHits;
        fEvent->ePmt += ltraj.KE;
        hPmtHits->Fill(waveLength);
        ltraj.Type = LTTrajectType::PMTHIT;
        //G4cout << " LegendAnalysis hitPMT creator process " << creator->GetProcessName() 
        //  << " track status " << trackInformation->GetTrackStatus() << " is hit  " 
        //  << gtrj->IsHit() <<  " is wls " << gtrj->IsWLS() << "  KE=" << ltraj.KE  << G4endl;
      } else if(gtrj->IsWLS()) {
        ++fEvent->nTrajWlsScint;
        hWls->Fill(waveLength);
        ltraj.Type = LTTrajectType::WLS;
      } else {
        ++fEvent->nTrajArScint;
        hOptical->Fill(waveLength);
        ltraj.Type = LTTrajectType::SCI;
      }
      // not optical 
     } 
    
    if(gtrj->IsGeHit()) {
      ltraj.Type = LTTrajectType::GEHIT;
      ++fEvent->nTrajGeHits;
    }

    if(ltraj.PDG==11) { // electron
      hEElectron->Fill(ltraj.KE/keV);
    } else if(ltraj.PDG==22) { // gamma 
      hEGamma ->Fill(ltraj.KE/keV);
    }
    
    // save in tree
    fEvent->traject.push_back(ltraj);
  }
}
