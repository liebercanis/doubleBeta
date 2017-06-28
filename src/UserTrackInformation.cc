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
// $Id: UserTrackInformation.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Legend/src/UserTrackInformation.cc
/// \brief Implementation of the UserTrackInformation class
//
//
#include "UserTrackInformation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UserTrackInformation::UserTrackInformation() 
  : fStatus(active),fReflections(0),fForcedraw(false),fPrimary(false),fInToGe(0),fOutOfGe(0),fSpikeReflection(0) 
{
  fProcessName="";
  fPreName="";
  fPostName="";
  fParentId=-1;
  fPostStepStatus.clear();
  fBoundaryStatus.clear();
  fBoundaryName.clear();
  fPositionHistory.clear();
  fPositionEnergy.clear();
  fStepLength.clear();
  fStepKE.clear();

  fTrackStatusNames.clear();
  fTrackStatusNames.push_back("active");
  fTrackStatusNames.push_back("scint");
  fTrackStatusNames.push_back("absorbed");
  fTrackStatusNames.push_back("boundaryAbsorbed");
  fTrackStatusNames.push_back("absorbedLAr");
  fTrackStatusNames.push_back("inactive");
  fTrackStatusNames.push_back("hitWLS");
  fTrackStatusNames.push_back("totalInternal");
  fTrackStatusNames.push_back("backScatter");
  fTrackStatusNames.push_back("fresnelRefract");
  fTrackStatusNames.push_back("fresnelReflect");
  fTrackStatusNames.push_back("spikeReflect");
  fTrackStatusNames.push_back("hitPMT");
  fTrackStatusNames.push_back("hitGe");
  fTrackStatusNames.push_back("absGe");
  fTrackStatusNames.push_back("isBad");
  // G4cout << " UserTrackInfomation Track Status bits: " << G4endl;
  //for (unsigned i=0; i< fTrackStatusNames.size() ; ++i) G4cout << " bit " << i << " value " << pow(2,i) << "  " << fTrackStatusNames[i] << G4endl;
  
  
}
UserTrackInformation::UserTrackInformation(const  UserTrackInformation& right) 
{
    //G4int 
      fStatus=right.fStatus;
    //G4int
      fParentId=right.fParentId;
    //G4bool
      fPrimary=right.fPrimary;
    //G4int  
      fReflections=right.fReflections;
    //G4int 
      fSpikeReflection=right.fSpikeReflection;
    //G4bool 
      fForcedraw=right.fForcedraw;
    //G4String
      fProcessName=right.fProcessName;
    //G4String 
      fPreName=right.fPreName;
    //G4String 
      fPostName=right.fPostName;
    //G4int 
      fInToGe=right.fInToGe;
    //G4int 
      fOutOfGe=right.fOutOfGe;
    //G4ThreeVector 
      position=right.position; // end of track 
    //std::vector<int>
      fPostStepStatus=right.fPostStepStatus; //  fGeomBoundary=1 
    //std::vector<int> 
      fBoundaryStatus=right.fBoundaryStatus;
    //std::vector<std::string> 
      fBoundaryName=right.fBoundaryName;
    //std::vector<G4ThreeVector> 
      fPositionHistory=right.fPositionHistory;
    //std::vector<G4double> 
      fPositionEnergy=right.fPositionEnergy;
    //std::vector<G4double> 
      fStepLength=right.fStepLength;
    //std::vector<G4double> 
      fStepKE=right.fStepKE;  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UserTrackInformation::~UserTrackInformation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UserTrackInformation::AddTrackStatusFlag(int s)
{
  if(s&active) //track is now active
    fStatus&=~inactive; //remove any flags indicating it is inactive
  else if(s&inactive) //track is now inactive
    fStatus&=~active; //remove any flags indicating it is active
  fStatus|=s; //add new flags
}
