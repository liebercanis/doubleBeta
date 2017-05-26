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
     void SetAngDistTypeByName(G4String type){AngDistType = type;}
     void SetSourcePositionByName(G4String type){SourcePosType = type;}

     void Show(){
       G4cout <<
         " **************** LegendParticleSource ********** " << G4endl << 
         " \t physical volume       " << physVolumeName <<
         " \t angular dis type      " << AngDistType << 
         " \t energy  dis type      " << EnergyDisType << 
         " \t source position type  " << SourcePosType;

       if(particle_definition) G4cout << " \t particle is "  << particle_definition->GetParticleName();
       else G4cout << " \t no particle definition ";

       G4cout << G4endl;
     }

     G4VPhysicalVolume* GetPhysicalVolume() { return thePhysicalVolume; }
     void SetPosDisType(G4String);
     void SetPosDisShape(G4String);
     void SetCenterVector(G4ThreeVector);
     void SetHalfZ(G4double);
     void SetRadius(G4double);
     void GeneratePointSource();
     void GeneratePointsInVolume();
     void GeneratePointsOnGeSurface();
     G4bool IsSourceConfined();
     G4bool IsInArgon(G4ThreeVector rp);
     void ConfineSourceToVolume();
  
     // angular distribution
     void SetParticleMomentumDirection(G4ParticleMomentum);
     void GenerateIsotropicFlux();
     void SetMomentumIntoGeSurface();

     // energy distribution 
     void SetEnergyDistTypeByName(G4String type){EnergyDisType=type;}
     void SetMonoEnergy(G4double);
     void GenerateMonoEnergetic();
     G4String GetEnergyDisType(){ return EnergyDisType ;}
     
     inline G4double GetParticleEnergy() { return particle_energy;}
     void   GenAr39Energy() { particle_energy = hAr39Theory->GetRandom(); } // unit is MeV
     void   GenAlphaEnergy() { particle_energy = 5.5*MeV; } //Random Alpha energy
     

     // verbosity
     void SetVerbosity(G4int);
  
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
     G4String EnergyDisType;
     G4double MonoEnergy;

     // particle properties 
     G4int                  NumberOfParticlesToBeGenerated;
     G4ParticleDefinition*  particle_definition;
     G4ParticleMomentum     particle_momentum_direction;
     G4double               particle_energy;
     G4double               particle_charge;
     G4ThreeVector          particle_position;
     G4double               particle_time;
     G4ThreeVector          particle_polarization;
     G4ThreeVector          rSurfNormal;

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

