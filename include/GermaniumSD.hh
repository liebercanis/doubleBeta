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
/**
 * CLASS DECLARATION:  GermaniumSD.hh
 *
 * DESCRIPTION:
 *
 * AUTHOR:
 *
 * REVISION: MM-DD-YYYY
 *
 */

#ifndef GermaniumSD_h
#define GermaniumSD_h 1

// ---------------------------------------------------------------------------
#include "LegendAnalysis.hh"
//#include "GermaniumSDHit.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDKineticEnergyFilter.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"

class GermaniumSD : public G4VSensitiveDetector
{

public:
  GermaniumSD(G4String name, G4int nCopy=0);
  ~GermaniumSD();

  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory* );
  // -- NB: Overloaded version to be used with the LAr Instrumentation.
  // It is manually called from the process manager.
  G4bool ProcessHits_constStep(const G4Step* aStep, G4TouchableHistory* );
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();

private:
  TDirectory *fDir;
  TH1F* hTime;
  TH1F* hEnergy;
  TNtuple *ntGe;
  G4int fCopy;
  G4SDKineticEnergyFilter *fKineticFilter;

};

#endif

