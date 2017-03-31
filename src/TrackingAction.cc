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
#include "LegendTrajectory.hh"
#include "UserTrackInformation.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction()
{
  // create directory 
  fDir = LegendAnalysis::Instance()->topDir()->mkdir("track");
  fDir->cd();
  G4cout<<" TrackingAction working root directory  is  " << G4endl;  
  gDirectory->pwd();
  G4double LowE = 2.4796*eV;//500 nm
  G4double HighE = 12.3984*eV;//100 nm
  G4double LowWLS =  h_Planck*c_light/(700.0*nm);//700 nm
  G4double HighWLS = h_Planck*c_light/(200.0*nm);//200 nm
  hTrackPhotonE = new TH1F("TrackPhotonE"," scint photon energy in LAr",1000,LowE,HighE);
  hAbsorbedPhotonE = new TH1F("AbsorbedPhotonE"," absorbed scint photon energy in LAr",1000,LowE,HighE);
  hWLSPhotonE = new TH1F("WLSPhotonE"," WLS photon energy ",1000,LowWLS,HighWLS);
  hPMTPhotonE = new TH1F("PMTPhotonE"," WLS photon energy ",1000,LowWLS,HighWLS);
  hCherenkovPhotonE  = new TH1F("CherenkovPhotonE"," WLS photon energy ",1000,LowWLS,HighWLS);
  hTrackStatus = new TH1F("TrackStatus"," track status ",7,0,7);
  
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
  //The UI command /tracking/storeTrajectory _bool_ does the same.
  fpTrackingManager->SetStoreTrajectory(true);
  LegendTrajectory* trajectory=(LegendTrajectory*) fpTrackingManager->GimmeTrajectory();
  if(!trajectory) {
    G4cout << " WARNING  TrackingAction::PostUserTrackingAction no trajectory found so returning " << G4endl;
    return;
  }

  UserTrackInformation* trackInformation=(UserTrackInformation*) aTrack->GetUserInformation();
  
  //LegendAnalysis::Instance()->FillTrajectory(trajectory);
  G4double totE = aTrack->GetTotalEnergy();//Returns energy in MeV
  //G4double KE = aTrack->GetKineticEnergy();//Returns energy in MeV
  
  const G4VProcess* creator=aTrack->GetCreatorProcess();
  if(!creator) { 
    G4cout << " WARNING TrackingAction called with NULL track G4VProcess!  energy = " << totE << G4endl;
    return;
  } 
  //else 
  //  G4cout<<" PostUser tracking action process is " << creator->GetProcessName() << " energy " << totE << G4endl;

  if(aTrack->GetDefinition() ==G4OpticalPhoton::OpticalPhotonDefinition()){
    //enum TrackStatus { active=1, hitPMT=2, absorbed=4, boundaryAbsorbed=8,hitSphere=16, inactive=14}; 
    if(trackInformation->GetTrackStatus()&active) hTrackStatus->Fill(0); 
    else if(trackInformation->GetTrackStatus()&hitPMT) hTrackStatus->Fill(1); 
    else if(trackInformation->GetTrackStatus()&absorbed) hTrackStatus->Fill(2); 
    else if(trackInformation->GetTrackStatus()&boundaryAbsorbed) hTrackStatus->Fill(3); 
    else if(trackInformation->GetTrackStatus()&absorbedLAr) hTrackStatus->Fill(4); 
    else if(trackInformation->GetTrackStatus()&inactive) hTrackStatus->Fill(5); 
    else  hTrackStatus->Fill(6); 
    
    if(trackInformation->GetTrackStatus()&hitPMT) {
      hPMTPhotonE->Fill(totE);
    }
    
    if(creator->GetProcessName() ==  "Scintillation"){ 
      //LAr_Spectrum->Fill(KE);
      hTrackPhotonE->Fill(totE);
      if(trackInformation->GetTrackStatus()&absorbed){
        hAbsorbedPhotonE->Fill(totE);
        trajectory->SetDrawTrajectory(true);
      }
    } else if(creator->GetProcessName() ==  "OpWLS"){
      hWLSPhotonE->Fill(totE);
      trajectory->SetDrawTrajectory(true);
    } else if(creator->GetProcessName() == "Cerenkov") {
      hCherenkovPhotonE->Fill(totE);
    } else G4cout << " UNKNOWN optical photon creator  " << creator->GetProcessName()<< " definition  " << aTrack->GetDefinition() << " ???? " << G4endl;
  }
}
