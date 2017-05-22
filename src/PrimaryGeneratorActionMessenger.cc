///////////////////////////////////////////////////////////////////////////
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
#include "PrimaryGeneratorActionMessenger.hh"

#include "LegendAnalysis.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction* action)
: genAction(action)
{
  detDir = new G4UIdirectory("/generator/");
  detDir->SetGuidance(" setup options for particle generator ");

  ShowGeneratorCmd = new G4UIcmdWithoutParameter("/generator/show",this);
  ShowGeneratorCmd->SetGuidance(" show generator parameters ");
  ShowGeneratorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  PhysicalVolumeNameCmd = new G4UIcmdWithAString("/generator/volume",this);
  PhysicalVolumeNameCmd->SetGuidance("set generator particle ");
  PhysicalVolumeNameCmd->SetParameterName("volume",false);
  PhysicalVolumeNameCmd->SetCandidates("e- gamma alpha");
  PhysicalVolumeNameCmd->AvailableForStates(G4State_Idle);

  ParticleDefinitionCmd = new G4UIcmdWithAString("/generator/particle",this);
  ParticleDefinitionCmd->SetGuidance("set generator particle ");
	ParticleDefinitionCmd->SetParameterName("particle",false);
	ParticleDefinitionCmd->SetCandidates("e- gamma alpha");
	ParticleDefinitionCmd->AvailableForStates(G4State_Idle);


  SourceTypeCmd = new G4UIcmdWithAString("/generator/SourceType",this);
  SourceTypeCmd->SetGuidance("set particle source type");
  SourceTypeCmd->SetParameterName("source",false);
  SourceTypeCmd->SetCandidates("Ar39 Point ");
  SourceTypeCmd->AvailableForStates(G4State_Idle);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  G4cout << " PrimaryGeneratorActionMessenger::SetNewValue " << command << " value " << newValue << G4endl;
  if( command == PhysicalVolumeNameCmd) genAction->SetPhysicalVolumeByName(newValue);
  if( command == ParticleDefinitionCmd ) genAction->SetParticle(newValue);
	if( command == SourceTypeCmd ) genAction->SetSource(newValue);
  if (command == ShowGeneratorCmd) genAction->Show();	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
