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
// $Id$
//
/// \file SteppingAction.hh
/// \brief Definition of the SteppingAction class

#ifndef SteppingAction_h
#define SteppingAction_h 1
#include "LegendAnalysis.hh"

#include "G4UserSteppingAction.hh"
#include "G4OpBoundaryProcess.hh"

#include "globals.hh"

class DetectorConstruction;
class EventAction;

/// Stepping action class
///
/// It holds data member fEnergy for accumulating the energy deposit
/// in a selected volume step by step.
/// The selected volume is set from  the detector construction via the
/// SetVolume() function. The accumulated energy deposit is reset for each
/// new event via the Reset() function from the event action.

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(DetectorConstruction*, EventAction*);
    virtual ~SteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);


  private:
    bool  GeDebug; 
    TDirectory* fDir;
    TNtuple* ntStep;
    TNtuple* ntGeStep;
    TH1F* hWLSPhotonE;
    TH1F* hBoundaryStatus;
    TH1F* hParticleType;

    G4bool fOneStepPrimaries;
    G4OpBoundaryProcessStatus fExpectedNextStatus;
    
    DetectorConstruction* detector;
    EventAction*          eventaction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
