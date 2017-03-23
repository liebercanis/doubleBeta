//
// ********************************************************************
// * License and Disclaimer                                           *
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
:detector(det), eventaction(evt)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

  // get volume of the current step
  G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4String volumename = 			volume->GetName();
	G4Track* aTrack = 					step->GetTrack();
	G4ThreeVector pos = 				aTrack->GetPosition();
	G4ThreeVector dir = 				aTrack->GetMomentumDirection();
	G4double length = 					step->GetStepLength();


  // check if we are in scoring volume
	if (volumename(0,4) == "source"){
			G4ThreeVector dir2 = eventaction->direction;
			if (dir2.x()*dir2.x()+dir2.y()*dir2.y()+dir2.z()*dir2.z() == 0){
				eventaction->direction = dir;
				eventaction->position = pos;
				//G4cout << " d " <<dir.cosTheta() << " " << dir.phi() << G4endl;
				//G4cout << " p " <<pos.cosTheta() << " " << pos.phi() << G4endl;
			}
	}


	if (volumename.contains("P") || volumename.contains("B")){
      G4int k = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
			//G4cout <<  k << " " <<length/cm << G4endl;
			eventaction->FillDetector(k,length);
	}

}


