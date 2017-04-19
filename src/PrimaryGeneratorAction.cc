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
// $Id: PrimaryGeneratorAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Legend/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(){
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
 
  G4String particleName;
  fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="e-"));
  //Default energy,position,momentum
  fParticleGun->SetParticleEnergy(pGun_nrg);//511.0*keV);
  //position group 1  (-205.165517,2.829276,-55.000261) rmax 147.244883  zmax 142.671828
	//position group 2  (206.799000,0.000000,-58.770646) rmax 143.629220  zmax 146.940251
  fParticleGun->SetParticlePosition(G4ThreeVector(0,0,1));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  // Ar39 spectrum
  G4String pathFile = "External_data/Ar39Theory.root";
  TFile *ar39File = new TFile(pathFile.data());
  if (!ar39File ) 
    G4cout<<" PrimaryGeneratorActon ERROR:: file " << pathFile << " not found " << G4endl;
  else
    G4cout<<" PrimaryGeneratorAction INFO:: file " << pathFile << " opened " << G4endl;
  hAr39Theory=NULL;
  ar39File->GetObject("theory",hAr39Theory);
  if (!hAr39Theory ) 
    G4cout<<" PrimaryGeneratorAction ERROR:: no theory TH1F in file " << pathFile <<G4endl;
  else 
    G4cout<<" PrimaryGeneratorAction info hAr39Theory found " <<G4endl;
  
  // create directory 
  fDir = LegendAnalysis::Instance()->topDir()->mkdir("generate");
  fDir->cd();
  G4cout<<" PrimaryGeneratorAction working root directory  is  " << fDir->GetName() << G4endl;  
  gDirectory->pwd();
  gDirectory->Append(hAr39Theory);
  hAr39Test = (TH1D*) hAr39Theory->Clone("Ar39Test");
  hAr39Test->Reset();
  //for(int itry=0; itry<10000; ++itry) hAr39Test->Fill( getAr39Energy() );
  
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  // Ar39 event
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="e-"));
  G4double energy = getAr39Energy();
  hAr39Test->Fill( energy );
  fParticleGun->SetParticleEnergy( energy );
  fParticleGun->SetParticlePosition(G4ThreeVector(7,7,0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));  
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
