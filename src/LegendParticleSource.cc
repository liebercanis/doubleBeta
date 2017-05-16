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
//
// ParticleSource program
// --------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// This particle source is a shortened version of G4GeneralParticleSource by
// C Ferguson, F Lei & P Truscott (University of Southampton / DERA), with
// some minor modifications.
//////////////////////////////////////////////////////////////////////////////
//  This version modifed by M. Gold, April 2017 
//  allows generating events inside physical volume
//

#include <cmath>

#include "LegendParticleSource.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"


LegendParticleSource::LegendParticleSource(G4String physical_name ) 
  :physVolumeName(physical_name)
{

  thePhysicalVolume = G4PhysicalVolumeStore::GetInstance()->GetVolume(physVolumeName);
  G4cout << " LegendParticleSource physical volume is " << thePhysicalVolume->GetName() << G4endl;
  
  
  G4VPhysicalVolume *larPhys = G4PhysicalVolumeStore::GetInstance()->GetVolume("larPhysical");
  centerVector = larPhys->GetTranslation();



  // Ar39 spectrum
  G4String pathFile = "External_data/Ar39Theory.root";
  TFile *ar39File = new TFile(pathFile.data());
  if (!ar39File ) 
    G4cout<<" LegendParticleSource ERROR:: file " << pathFile << " not found " << G4endl;
  else
    G4cout<<" LegendParticleSource INFO:: file " << pathFile << " opened " << G4endl;
  hAr39Theory=NULL;
  ar39File->GetObject("theory",hAr39Theory);
  if (!hAr39Theory ) 
    G4cout<<" LegendParticleSource ERROR:: no theory TH1F in file " << pathFile <<G4endl;
  else 
    G4cout<<" LegendParticleSource info hAr39Theory found " <<G4endl;
  
  // create directory 
  fDir = LegendAnalysis::Instance()->topHistDir()->mkdir("particleSource");
  fDir->cd();
  G4cout<<" LegendParticleSource working root directory  is  " << fDir->GetName() << G4endl;  
  gDirectory->pwd();
  gDirectory->Append(hAr39Theory);
  hAr39Gen = (TH1D*) hAr39Theory->Clone("Ar39Gen");
  hAr39Gen->Reset();
  
  
  NumberOfParticlesToBeGenerated = 1;
  particle_definition = NULL;
  G4ThreeVector zero(0., 0., 0.);
  particle_momentum_direction = G4ParticleMomentum(1., 0., 0.);
  particle_energy = 1.0*MeV;
  particle_position = zero;
  particle_time = 0.0;
  particle_polarization = zero;
  particle_charge = 0.0;

  SourcePosType = "Volume";
  Confine = true;
 
  AngDistType = "iso"; 
  MinTheta = 0.;
  MaxTheta = pi;
  MinPhi = 0.;
  MaxPhi = twopi;

  EnergyDisType = "Ar39";
  //EnergyDisType = "Mono";
  MonoEnergy = 1*MeV;
  G4cout<<" LegendParticleSource source type is   " << EnergyDisType;
  if(EnergyDisType = "Mono") G4cout << " energy " << MonoEnergy;
  G4cout << G4endl;
  

  verbosityLevel = 0;

  //theMessenger = new LegendParticleSourceMessenger(this);
  gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
}

LegendParticleSource::~LegendParticleSource()
{
  //delete theMessenger;
}

void LegendParticleSource::SetPosDisType(G4String PosType) 
{
  SourcePosType = PosType;
}

void LegendParticleSource::SetPosDisShape(G4String shapeType)
{
  Shape = shapeType;
}

void LegendParticleSource::SetCenterVector(G4ThreeVector coordsOfCentre)
{
  centerVector = coordsOfCentre;
}

void LegendParticleSource::SetHalfZ(G4double zhalf)
{
  halfz = zhalf;
}

void LegendParticleSource::SetRadius(G4double radius)
{
  Radius = radius;
}

void LegendParticleSource::ConfineSourceToVolume()
{
  if(verbosityLevel == 2) G4cout << physVolumeName << G4endl;
  Confine=false;
  if(thePhysicalVolume) Confine=true;

}


void LegendParticleSource::SetAngDistType(G4String atype)
{
  AngDistType = atype;
}


void LegendParticleSource::GeneratePointSource()
{
  // Generates Points given the point source.
  if(SourcePosType == "Point")
    particle_position = centerVector;
  else
    if(verbosityLevel >= 1)
      G4cout << " LegendParticleSource **** Error SourcePosType is not set to Point" << G4endl;
}


void LegendParticleSource::GeneratePointsInVolume()
{
  G4VSolid *solid = thePhysicalVolume->GetLogicalVolume()->GetSolid();
  G4ThreeVector rmin;
  G4ThreeVector rmax;
  solid->Extent(rmin,rmax);
  G4ThreeVector rpoint;
  centerVector = thePhysicalVolume->GetTranslation(); // absolute position of the center of this phyisical volume

  bool isInside=false;
  while(!isInside){
    G4double rx = ( rmax.x() - rmin.x() )*(2.*G4UniformRand()-1.);
    G4double ry = ( rmax.y() - rmin.y() )*(2.*G4UniformRand()-1.);
    G4double rz = ( rmax.z() - rmin.z() )*(2.*G4UniformRand()-1.);
    rpoint.set(rx,ry,rz);                         // relative position 
    particle_position = rpoint +  centerVector ; // absolute position
    isInside = solid->Inside(rpoint)== EInside::kInside && IsInArgon(particle_position) ;
  }
}

G4bool LegendParticleSource::IsInArgon(G4ThreeVector rp) 
{
  bool isit = true;
  // list of all volumes
  G4PhysicalVolumeStore* theStore = G4PhysicalVolumeStore::GetInstance();
  for(G4int istore = 0; istore< theStore->size() ; ++istore ){
    G4VPhysicalVolume *pvol = theStore->at(istore);
    G4String sname = pvol->GetName();

    // if this is a volume to reject, see if it is inside
    if (  (sname.find("B8") != string::npos)
       ||(sname.find("P4") != string::npos )
       ||(sname.find("WLS") != string::npos )
       ||(sname.find("Housing") != string::npos )
       ||(sname.find("Glass") != string::npos )
       ||(sname.find("VM2000") != string::npos )
       ||(sname.find("Cryo") != string::npos )) { 
      G4VSolid *solid = pvol->GetLogicalVolume()->GetSolid();
      G4ThreeVector rtran  = pvol->GetTranslation();
      G4ThreeVector rel = rp - rtran; // point relative to center of solid
      if(solid->Inside(rel)== EInside::kInside) isit = false;
      //G4cout << " vol " << sname << " isit? " << isit << endl;
    }
  }
  return isit; 
}

G4bool LegendParticleSource::IsSourceConfined()
{

  // Method to check point is within the volume specified
  if(Confine == false) {
    G4cout << "Error: Confine is false" << G4endl;
    return false;
  }
  G4ThreeVector null(0.,0.,0.);
  G4ThreeVector *ptr = &null;

  // Check particle_position is within VolName
  G4VPhysicalVolume *theVolume =gNavigator->LocateGlobalPointAndSetup(particle_position,ptr,true);
  //G4cout << " LegendParticleSource **** Particle is in volume " << physVolumeName << " " << theVolume->GetName() << G4endl;
  if( theVolume->GetName() == physVolumeName) {
    if(verbosityLevel >= 1)
    //if( 1)
      G4cout << " LegendParticleSource **** Particle is in volume " << physVolumeName << G4endl;
    return(true);
  }
  else
    return(false);
}


void LegendParticleSource::SetParticleMomentumDirection
   (G4ParticleMomentum aDirection) {

  particle_momentum_direction =  aDirection.unit();
}


void LegendParticleSource::GenerateIsotropicFlux()
{

  G4double rndm, rndm2;
  G4double px, py, pz;

  G4double sintheta, sinphi, costheta, cosphi;
  rndm = G4UniformRand();
  costheta = std::cos(MinTheta) - rndm * (std::cos(MinTheta) - std::cos(MaxTheta));
  sintheta = std::sqrt(1. - costheta*costheta);
  
  rndm2 = G4UniformRand();
  Phi = MinPhi + (MaxPhi - MinPhi) * rndm2; 
  sinphi = std::sin(Phi);
  cosphi = std::cos(Phi);

  px = -sintheta * cosphi;
  py = -sintheta * sinphi;
  pz = -costheta;

  G4double ResMag = std::sqrt((px*px) + (py*py) + (pz*pz));
  px = px/ResMag;
  py = py/ResMag;
  pz = pz/ResMag;

  particle_momentum_direction.setX(px);
  particle_momentum_direction.setY(py);
  particle_momentum_direction.setZ(pz);

  // particle_momentum_direction now holds unit momentum vector.
  if(verbosityLevel >= 2)
    G4cout << "Generating isotropic vector: " << particle_momentum_direction << G4endl;
}


void LegendParticleSource::SetEnergyDisType(G4String DisType)
{
  EnergyDisType = DisType;
}

void LegendParticleSource::SetMonoEnergy(G4double menergy)
{
  MonoEnergy = menergy;
}

void LegendParticleSource::GenerateMonoEnergetic()
{
  particle_energy = MonoEnergy;
}


void LegendParticleSource::SetVerbosity(int vL)
{
  verbosityLevel = vL;
  G4cout << " LegendParticleSource **** Verbosity Set to: " << verbosityLevel << G4endl;
}

void LegendParticleSource::SetParticleDefinition
  (G4ParticleDefinition* aParticleDefinition)
{
  particle_definition = aParticleDefinition;
  particle_charge = particle_definition->GetPDGCharge();
}


void LegendParticleSource::GeneratePrimaryVertex(G4Event *evt)
{

  if(particle_definition==NULL) {
    G4cout << " LegendParticleSource **** No particle has been defined!" << G4endl;
    return;
  }
  
  // Position
  G4bool srcconf = false;
  G4int LoopCount = 0;
  
  while(srcconf == false)  {
    if(SourcePosType == "Point")
      GeneratePointSource();
    else if(SourcePosType == "Volume")
      GeneratePointsInVolume();
    else {
      G4cout << " LegendParticleSource **** Error: SourcePosType undefined" << G4endl;
      G4cout << "Generating point source" << G4endl;
      GeneratePointSource();
    }
    if(Confine == true) {
      srcconf = IsSourceConfined();
      // if source in confined srcconf = true terminating the loop
      // if source isnt confined srcconf = false and loop continues
    }
    else if(Confine == false)
      srcconf = true; // terminate loop
    
    LoopCount++;
    if(LoopCount == 100000) {
      G4cout << " LegendParticleSource *************************************" << G4endl;
        G4cout << "LoopCount = 100000" << G4endl;
        G4cout << "Either the source distribution >> confinement" << G4endl;
        G4cout << "or any confining volume may not overlap with" << G4endl;
        G4cout << "the source distribution or any confining volumes" << G4endl;
        G4cout << "may not exist"<< G4endl;
        G4cout << "If you have set confine then this will be ignored" <<G4endl;
        G4cout << "for this event." << G4endl;
        G4cout << "*************************************" << G4endl;
        srcconf = true; //Avoids an infinite loop
      }
  }

  // Angular stuff
  if(AngDistType == "iso")
    GenerateIsotropicFlux();
  else if(AngDistType == "direction")
    SetParticleMomentumDirection(particle_momentum_direction);
  else
    G4cout << " LegendParticleSource **** Error: AngDistType has unusual value" << G4endl;
  // Energy stuff
  if(EnergyDisType == "Mono")
    GenerateMonoEnergetic();
  else if(EnergyDisType == "Ar39") {
    GenAr39Energy();
    hAr39Gen->Fill( GetParticleEnergy());
  } else 
    G4cout << " LegendParticleSource **** Error: EnergyDisType has unknown value = " << EnergyDisType << G4endl;
  
  // create a new vertex
  G4PrimaryVertex* vertex = 
    new G4PrimaryVertex(particle_position,particle_time);

  if(verbosityLevel >= 2)
    G4cout << " LegendParticleSource **** Creating primaries and assigning to vertex" << G4endl;
  // create new primaries and set them to the vertex
  G4double mass =  particle_definition->GetPDGMass();
  G4double energy = particle_energy + mass;
  G4double pmom = std::sqrt(energy*energy-mass*mass);
  G4double px = pmom*particle_momentum_direction.x();
  G4double py = pmom*particle_momentum_direction.y();
  G4double pz = pmom*particle_momentum_direction.z();
  
  if(verbosityLevel >= 1){
    G4cout << " LegendParticleSource **** Particle name: " 
	   << particle_definition->GetParticleName() << G4endl; 
    G4cout << "       Energy: "<<particle_energy << G4endl;
    G4cout << "     Position: "<<particle_position<< G4endl; 
    G4cout << "    Direction: "<<particle_momentum_direction << G4endl;
    G4cout << " NumberOfParticlesToBeGenerated: "
	   << NumberOfParticlesToBeGenerated << G4endl;
  }


  for( G4int i=0; i<NumberOfParticlesToBeGenerated; i++ ) {
    G4PrimaryParticle* particle =
      new G4PrimaryParticle(particle_definition,px,py,pz);
    particle->SetMass( mass );
    particle->SetCharge( particle_charge );
    particle->SetPolarization(particle_polarization.x(),
			      particle_polarization.y(),
			      particle_polarization.z());
    vertex->SetPrimary( particle );
  }
  evt->AddPrimaryVertex( vertex );
  if(verbosityLevel > 1)
    G4cout << " LegendParticleSource ****  Primary Vetex generated "<< G4endl;   
}




