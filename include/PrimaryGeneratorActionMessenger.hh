//////////////////////////////////////////////////////////////////////
// This code implementation is the intellectual property of the
// ton-scale 0vbb in Germanium collaboration. It is based on Geant4, an
// intellectual property of the RD44 GEANT4 collaboration.
//
// *********************
//
// Neither the authors of this software system, nor their employing
// institutes, nor the agencies providing financial support for this
// work make any representation or warranty, express or implied,
// regarding this software system or assume any liability for its use.
// By copying, distributing or modifying the Program (or any work based
// on the Program) you indicate your acceptance of this statement,
// and all its terms.
///////////////////////////////////////////////////////////////////////////
#ifndef PrimaryGeneratorActionMessenger_h
#define PrimaryGeneratorActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "PrimaryGeneratorAction.hh"

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorActionMessenger: public G4UImessenger
{
  public:
    PrimaryGeneratorActionMessenger(PrimaryGeneratorAction*);
   ~PrimaryGeneratorActionMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    PrimaryGeneratorAction* genAction;

    G4UIdirectory*             Dir;
    G4UIdirectory*             detDir;
	  G4UIcmdWithAString*        EnergyDistTypeCmd;
		G4UIcmdWithAString*        AngDistTypeCmd;
		G4UIcmdWithAString*        PhysicalVolumeNameCmd;
		G4UIcmdWithAString*        ParticleDefinitionCmd;
		G4UIcmdWithAString*        SourcePositionTypeCmd;
		G4UIcmdWithoutParameter*   ShowGeneratorCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
