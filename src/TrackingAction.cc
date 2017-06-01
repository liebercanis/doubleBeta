//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: TrackingAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Legend/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
#include "TrackingAction.hh"
#include "UserTrackInformation.hh"
#include "G4EventManager.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction()
{
  // create directory 
  fDir = LegendAnalysis::Instance()->topHistDir()->mkdir("track");
  fDir->cd();
  G4cout<<" TrackingAction working root directory  is  " << G4endl;  
  gDirectory->pwd();
  G4double LowE = 2.4796*eV;//500 nm
  G4double HighE = 12.3984*eV;//100 nm
  /* 
  G4double LowE = 1.7712*eV;//700 nm
  G4double HighE = 12.3984*eV;//100 nm
  hWLSPhotonE = new TH1F("StepWLSPhotonE"," photon energy from WLS",1000,LowE,HighE);
   */
  G4double LowWLS =  h_Planck*c_light/(700.0*nm);//700 nm
  G4double HighWLS = h_Planck*c_light/(200.0*nm);//200 nm
  hTrackScintE = new TH1F("TrackScintE"," scint photon energy in LAr",1000,LowE,HighE);
  hTrackPhotonE = new TH1F("TrackPhotonE"," all photon energy in LAr",1000,LowE,HighE);
  hTrackScintYield = new TH1F("TrackScintYield"," scint photon yield",1000,0,5);
  hTrackScintYield->SetXTitle(" mean scint photons / primary energy (MeV) ");
  hAbsorbedPhotonE = new TH1F("AbsorbedPhotonE"," absorbed scint photon energy in LAr",1000,LowE,HighE);
  hWLSPhotonE = new TH1F("WLSPhotonE"," WLS photon energy ",1000,LowWLS,HighWLS);
  hPMTPhotonE = new TH1F("PMTPhotonE"," WLS photon energy ",1000,LowWLS,HighWLS);
  hCherenkovPhotonE  = new TH1F("CherenkovPhotonE"," WLS photon energy ",1000,LowWLS,HighWLS);
  hTrackStatus = new TH1F("TrackStatus"," track status ", TrackBit::MaxHistogramBit+1,0, TrackBit::MaxHistogramBit+1);
  
  // must be in top directory for ChangeFile to work
  LegendAnalysis::Instance()->topTreeDir()->cd();
  // make tree in output file
  fTrackTree = new TTree("trkTree","trkTree");
  //fTree->SetMaxTreeSize(1000000);
  fTrackTree->SetMaxTreeSize(25);
  fLTTrack = new LTTrack();
  fTrackTree->Branch("track",&fLTTrack);

  G4cout << " ...  = " << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //Let this be up to the user via vis.mac
  //  fpTrackingManager->SetStoreTrajectory(true);
  fpTrackingManager->SetTrajectory(new LegendTrajectory(aTrack) );
  fpTrackingManager->SetUserTrackInformation(new UserTrackInformation);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack){
  if(!aTrack) { 
    G4cout << " WARNING TrackingAction called with NULL track!  " << G4endl;
    return;
  }

   LTEvent* ltEvent = LegendAnalysis::Instance()->getEvent();
  
  //The UI command /tracking/storeTrajectory _bool_ does the same.
  fpTrackingManager->SetStoreTrajectory(true);
  LegendTrajectory* trajectory=(LegendTrajectory*) fpTrackingManager->GimmeTrajectory();
  if(!trajectory) {
    G4cout << " WARNING  TrackingAction::PostUserTrackingAction no trajectory found so returning " << G4endl;
    return;
  }

  //trajectory->SetForceDrawTrajectory(true);
  trajectory->SetDrawTrajectory(false);
  UserTrackInformation* trackInformation=(UserTrackInformation*) aTrack->GetUserInformation();
  
  
  //LegendAnalysis::Instance()->FillTrajectory(trajectory);
  G4double kineticE = aTrack->GetKineticEnergy();//Returns energy in MeV
  
  trajectory->SetParentId(trackInformation->GetParentId());
  const G4VProcess* creator=aTrack->GetCreatorProcess();
  if(!creator&&trackInformation->GetParentId()!=0) {
    G4cout << " WARNING TrackingAction called with NULL track G4VProcess!  parent id " << trackInformation->GetParentId()   << G4endl;
    hTrackStatus->Fill(isBad); 
    return;
  } 

  if(trackInformation->IsPrimary()) {  
    trajectory->SetDrawTrajectory(true);
    trajectory->SetPrimary();
    ltEvent->ePrimary=kineticE;
    //G4cout << " TrackingAction PRIMARY TRACK track definition is  " << aTrack->GetDefinition()->GetParticleName() << " is prim? " << trajectory->IsPrimary() << G4endl;
  }

   float evertex  =  ltEvent->ePrimary;
   //G4cout << " \t TRACKINGACTION event vertex energy " << evertex << G4endl;
   G4double meanScintE = h_Planck*c_light/(128.0*nm);

  if(aTrack->GetDefinition() ==G4OpticalPhoton::OpticalPhotonDefinition()){
    ++ltEvent->nOptPhotons;
    hTrackStatus->Fill(trackInformation->GetTrackBit()); 
    
    hTrackPhotonE->Fill(kineticE);
    if(trackInformation->GetTrackStatus()&absorbed) hAbsorbedPhotonE->Fill(kineticE);
    if(creator->GetProcessName() ==  "Scintillation") {
      ++ltEvent->nArScint;
      hTrackScintE->Fill(kineticE);
      if(evertex>0) hTrackScintYield->Fill( kineticE /evertex /meanScintE);
      //if(evertex>0) G4cout << " \t TRACKINGACTION event vertex energy " << evertex << "  photon energy ratio " << kineticE /evertex << G4endl;
    }
    if(creator->GetProcessName() == "Cerenkov") hCherenkovPhotonE->Fill(kineticE);
    
    // use track status set in SteppingAction
    trajectory->SetTrackStatus(trackInformation->GetTrackStatus());
    if(trackInformation->GetTrackStatus()&hitPMT) {
      hPMTPhotonE->Fill(kineticE);
      trajectory->SetDrawTrajectory(false);
      ++ltEvent->nPmtHits;
    } else if(trackInformation->GetTrackStatus()&hitWLS) {
      trajectory->SetDrawTrajectory(false);
      hWLSPhotonE->Fill(kineticE);
      ++ltEvent->nWlsScint; 
      
    }  
  } //optical 

  if(trackInformation->GetTrackStatus()&hitGe) {
      trajectory->SetDrawTrajectory(true);
      ++ltEvent->nGeHits; 
      //G4cout << " TrackingAction hitGe ishit? " << trajectory->IsGeHit() << " nGeHits" << ltEvent->nGeHits << G4endl;
  }

 if(trackInformation->GetTrackStatus()&eIoni) {
      trajectory->SetDrawTrajectory(true);
  }

 if(trackInformation->GetTrackStatus()&hIoni) {
      trajectory->SetDrawTrajectory(true);
  }

  if(trackInformation->GetTrackStatus()&ionIoni) {
      trajectory->SetDrawTrajectory(true);
  }  
  
  // fill tree
  // Hereafter we call current volume the volume where the step has just gone through. Geometrical informations are available from preStepPoint.
  G4StepPoint* currentPoint = aTrack->GetStep()->GetPreStepPoint();
  G4StepPoint* nextPoint    = aTrack->GetStep()->GetPostStepPoint();
  //G4VTouchable and its derivates keep these geometrical informations. We retrieve a touchable by creating a handle for it:
  G4TouchableHandle ctouch = currentPoint->GetTouchableHandle();
  fLTTrack->physVolName = ctouch->GetVolume()->GetName();
  fLTTrack->copy = ctouch->GetCopyNumber();

  // To check that the particle is leaving the current volume (i.e. it is at the last step in the volume; the postStepPoint is at the boundary):
  fLTTrack->isLeaving = nextPoint->GetStepStatus() == fGeomBoundary;
  
  fLTTrack->evId = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  fLTTrack->trkId = aTrack->GetTrackID();
  fLTTrack->parentId = aTrack->GetParentID();
  fLTTrack->status=trackInformation->GetTrackStatus();
  fLTTrack->boundaryStatus = trackInformaton->GetBoundaryProcessStatus();
  fLTTrack->length=aTrack->GetTrackLength();
  fLTTrack->nstep=aTrack->GetCurrentStepNumber();
  fLTTrack->stepLength=aTrack->GetStepLength();
  G4ThreeVector trkPos = aTrack->GetPosition();
  fLTTrack->position.SetXYZ(trkPos.x(),trkPos.y(),trkPos.z());
  G4ThreeVector vertPos = aTrack->GetVertexPosition(); 
  fLTTrack->vertPosition.SetXYZ(vertPos.x(),vertPos.y(),vertPos.z());
  // Global time (time since the current event began)
  fLTTrack->time=aTrack->GetGlobalTime()/microsecond;
  // Local time (time since the current track began)
  fLTTrack->trkTime=aTrack->GetLocalTime()/microsecond;
  fLTTrack->ke=aTrack->GetKineticEnergy()/electronvolt;
  fLTTrack->edep=aTrack->GetStep()->GetTotalEnergyDeposit()/electronvolt;
  fLTTrack->particleName = aTrack->GetDefinition()->GetParticleName();
  //fLTTrack->print();
  fTrackTree->Fill();
}


