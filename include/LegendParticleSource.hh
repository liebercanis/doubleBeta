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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// ParticleSource header
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////

#ifndef LegendParticleSource_h
#define LegendParticleSource_h 1

#include "LegendAnalysis.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4Navigator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"

//#include "LegendParticleSourceMessenger.hh"


class LegendParticleSource : public G4VPrimaryGenerator {

   public:
     LegendParticleSource(); 
     ~LegendParticleSource ();
     void GeneratePrimaryVertex(G4Event *evt);

   public:

     // position distribution 
     void SetPhysicalVolume( G4VPhysicalVolume* physVol) { thePhysicalVolume=physVol ; }
     void SetPhysicalVolumeByName( G4String physical_name);
     
     G4VPhysicalVolume* GetPhysicalVolume() { return thePhysicalVolume; }
     void Show(){
       G4cout <<
         " **************** LegendParticleSource ********** " << G4endl << 
         " \t physical volume " << physVolumeName <<
         " \t source type is   " << SourceType << 
         " \t source position type is   " << SourcePosType;
       if(particle_definition) G4cout << " \t particle is "  << particle_definition->GetParticleName(); 
       G4cout << G4endl;
     }
     void GeneratePointSource();
     void GeneratePointsInVolume();
     G4bool IsSourceConfined();
     G4bool IsInArgon(G4ThreeVector rp);
     void ConfineSourceToVolume();
 
     void SetPosDisType(G4String PosType) { SourcePosType = PosType;}

     void SetPosDisShape(G4String shapeType){ Shape = shapeType;}

     void SetCenterVector(G4ThreeVector center){ centerVector = center;}

     void SetHalfZ(G4double zhalf) { halfz = zhalf;}

     void SetRadius(G4double radius){ Radius = radius;}

     void SetAngDistType(G4String atype) { AngDistType = atype;}

     void SetParticleMomentumDirection (G4ParticleMomentum aDirection) { particle_momentum_direction =  aDirection.unit(); }
 
     void SetSourceType(G4String DisType) {SourceType = DisType;}

     void SetVerbosity(int vL) {
       verbosityLevel = vL;
       G4cout << " LegendParticleSource **** Verbosity Set to: " << verbosityLevel << G4endl;
     }
     
     void GenerateIsotropicFlux();

     // energy distribution 
     void SetSource(G4String type){ SourceType = type; } 
     G4String GetSourceType(){ return SourceType ;}
     void SetMonoEnergy(G4double menergy) { particle_energy = menergy;}
     
     inline G4double GetParticleEnergy() { return particle_energy;}

     // Ar39 energy
     void   GenAr39Energy() { particle_energy = hAr39Theory->GetRandom(); } // unit is MeV
     
     // particle properties
     void SetParticleDefinition(G4ParticleDefinition * aParticleDefinition);
     inline void SetParticleCharge(G4double aCharge)
        { particle_charge = aCharge; }
  
   private:

     // position distribution
     G4String physVolumeName;
     G4String SourcePosType;
     G4String Shape;
     G4double halfz;
     G4double Radius;
     G4ThreeVector centerVector;
     G4bool Confine;
     G4String AngDistType;
     G4double MinTheta, MaxTheta, MinPhi, MaxPhi;
     G4double Phi;
     G4String SourceType;

     // particle properties 
     G4int                  NumberOfParticlesToBeGenerated;
     G4ParticleDefinition*  particle_definition;
     G4ParticleMomentum     particle_momentum_direction;
     G4double               particle_energy;
     G4double               particle_charge;
     G4ThreeVector          particle_position;
     G4double               particle_time;
     G4ThreeVector          particle_polarization;

     // Verbose
     G4int verbosityLevel;
     G4VPhysicalVolume* thePhysicalVolume;

     //LegendParticleSourceMessenger *theMessenger;
     G4Navigator *gNavigator;
     TH1D *hAr39Theory;
     TH1D *hAr39Gen;
     TDirectory *fDir;

  
};


#endif

