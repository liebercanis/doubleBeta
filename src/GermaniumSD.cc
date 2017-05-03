//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                            MaGe Simulation                                //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
//                                                          
//
//---------------------------------------------------------------------------//
/**
 * SPECIAL NOTES:
 * 
 * Original desing by Alexander Klimenko
 */
// 
//---------------------------------------------------------------------------//
/**
 *
 * AUTHOR:  Markus Knapp
 * CONTACT: @CONTACT@
 * FIRST SUBMISSION: 
 * 
 * REVISION:
 *
 */
//---------------------------------------------------------------------------//
//



#include "GermaniumSD.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

//#include "io/MGLogger.hh"//Gerda log maker...will flush later


GermaniumSD::GermaniumSD(G4String name, G4int nCopy) 
:G4VSensitiveDetector(name),fCopy(nCopy)
{
  if(fCopy==0) {
    fDir = LegendAnalysis::Instance()->topDir()->mkdir("GeThits");
    fDir->cd();
    hTime = new TH1F("GeHitsTime"," Ge hit time  ",2000,0,4000);
    hTime->GetYaxis()->SetTitle(" hits )");
    hTime->GetXaxis()->SetTitle(" time (ns) ");

    hEnergy = new TH1F("GeEnergy"," absorbed energy ",2000,0,2);
    hEnergy->GetYaxis()->SetTitle(" hits ");
    hEnergy->GetXaxis()->SetTitle(" energy (keV) ");
  }
  
}

GermaniumSD::~GermaniumSD()
{
}

void GermaniumSD::Initialize(G4HCofThisEvent* /*HCE*/)
{
  //G4cout <<  "\tGermaniumSD Initialized   " <<  G4endl;
}

G4bool GermaniumSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  return ProcessHits_constStep(aStep,NULL);
}

G4bool GermaniumSD::ProcessHits_constStep(const G4Step* aStep, G4TouchableHistory* )
{
  if(aStep==NULL) {
    //G4cout << " ProcessHits_constStep called with null step " << G4endl;
    return false;
  }
  G4ParticleDefinition* particleType = aStep->GetTrack()->GetDefinition();
  G4String particleName = particleType->GetParticleName();
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double gtime = aStep->GetPostStepPoint()->GetGlobalTime();   //measured in nanoseconds;
  G4double time = aStep->GetTrack()->GetGlobalTime();
  if(edep>0) G4cout << "GermanimSD::ProcessHits_constStep  " << time << " global time   " << gtime << " edep " << edep << " particle name " << particleName << G4endl; 
  //const G4VPhysicalVolume* physVol = aStep->GetPostStepPoint()->GetPhysicalVolume();
  hTime->Fill( aStep->GetTrack()->GetGlobalTime()/ns ); //convert to ns
  hEnergy->Fill(edep/keV);
  return true;
}

void GermaniumSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  //if(HCID<0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  //HCE->AddHitsCollection( HCID, PMTCollection );
}

void GermaniumSD::clear()
{

} 

void GermaniumSD::DrawAll()
{
  
} 

void GermaniumSD::PrintAll()
{

} 




