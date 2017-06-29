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
//
// $Id: G4OpWLS.hh,v 1.4 2006-06-29 21:08:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// Optical Photon WaveLength Shifting (WLS) Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpWLS.hh
// Description: Discrete Process -- Wavelength Shifting of Optical Photons 
// Version:     1.0
// Created:     2003-05-13
// Author:      John Paul Archambault
//              (Adaptation of G4Scintillation and G4OpAbsorption)
// Updated:     2005-07-28 add G4ProcessType to constructor
//              2006-05-07 - add G4VWLSTimeGeneratorProfile
// mail:        gum@triumf.ca
//              jparcham@phys.ualberta.ca
//
////////////////////////////////////////////////////////////////////////
// modified by MGold NMcFadden Thu Mar  2 10:25:05 MST 2017

#ifndef OpWLS_h
#define OpWLS_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4VWLSTimeGeneratorProfile.hh"

#include "LegendAnalysis.hh"
#include "TH1F.h"


// Class Description:
// Discrete Process -- Bulk absorption of Optical Photons.
// Class inherits publicly from G4VDiscreteProcess
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class OpWLS : public G4VDiscreteProcess 
{

public: // Without description

  ////////////////////////////////
  // Constructors and Destructor
  ////////////////////////////////

  OpWLS(const G4String& processName = "OpWLS",
                   G4ProcessType type = fOptical);

  ~OpWLS();

  ////////////
  // Methods
  ////////////

public: // With description

  G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
  // Returns true -> 'is applicable' only for an optical photon.

  G4double GetMeanFreePath(const G4Track& aTrack,
			   G4double ,
			   G4ForceCondition* );
  // Returns the absorption length for bulk absorption of optical
  // photons in media with a specified attenuation length.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				  const G4Step&  aStep);
  // This is the method implementing bulk absorption of optical
  // photons.

  G4PhysicsTable* GetIntegralTable() const;
  // Returns the address of the WLS integral table.

  void DumpPhysicsTable() const;
  // Prints the WLS integral table.

  void UseTimeProfile(const G4String name);
  // Selects the time profile generator

private:

  void BuildThePhysicsTable();
  // Is the WLS integral table;

protected:

  G4VWLSTimeGeneratorProfile* WLSTimeGeneratorProfile;
  G4PhysicsTable* theIntegralTable;
  TDirectory *fDir;
  TH1F * hOpWLSIntegral;
  TH1F * hOpWLSphotons;

};

////////////////////
// Inline methods
////////////////////

inline
G4bool OpWLS::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
G4PhysicsTable* OpWLS::GetIntegralTable() const
{
  return theIntegralTable;
}

inline
void OpWLS::DumpPhysicsTable() const
{
  G4int PhysicsTableSize = theIntegralTable->entries();
  G4PhysicsOrderedFreeVector *v;
 
  for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
    {
      v = (G4PhysicsOrderedFreeVector*)(*theIntegralTable)[i];
      v->DumpValues();
    }
}

#endif /* OpWLS_h */

