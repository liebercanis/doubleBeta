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
// $Id: TrackingAction.hh 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Legend/include/TrackingAction.hh
/// \brief Definition of the TrackingAction class
//
//
#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "LegendAnalysis.hh"

#include "G4UserTrackingAction.hh"
#include "G4VTrajectory.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
class TrackingAction : public G4UserTrackingAction {

  public:

    TrackingAction();
    virtual ~TrackingAction() {};

    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
    TDirectory *fDir;
    TH1F *hTrackPhotonE;
    TH1F *hTrackScintE;
    TH1F *hAbsorbedPhotonE;
    TH1F *hWLSPhotonE;
    TH1F *hPMTPhotonE;
    TH1F *hCherenkovPhotonE;
    TH1F *hTrackStatus;
  private:
  
 
};

#endif
