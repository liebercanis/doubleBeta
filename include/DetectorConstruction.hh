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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "LegendAnalysis.hh"
#include "PMTSD.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4PhysicalConstants.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"

class G4VPhysicalVolume;
class DetectorMessenger;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    void UpdateGeometry();

    void SetInnerShield (G4bool);
    void SetOuterShield (G4bool);
	  void SetOverlapsCheck(G4bool);
		void SetShieldStyle(G4String);
		void SetFillMaterial(G4String);
		void ConstructSDandField();
    void PlacePMT(G4ThreeVector r,double top_or_bot,int num);
    
    
    static const G4double LambdaE;
    G4double DetectorDimension[35][3] =
    {{32.11*mm,15.49*mm,692},{34.46*mm,25.42*mm,690},{35.35*mm,21.00*mm,688},{34.69*mm,24.88*mm,640},{0,0,0},
      {35.00*mm,15.00*mm,584},{32.94*mm,22.12*mm,674},{35.08*mm,16.99*mm,576},{30.10*mm,23.76*mm,680},{0,0,0},
      {35.00*mm,15.00*mm,676},{34.70*mm,23.90*mm,616},{33.66*mm,24.48*mm,614},{35.85*mm,23.80*mm,610},{0,0,0},
      {35.00*mm,15.00*mm,608},{35.00*mm,15.00*mm,598},{35.00*mm,15.00*mm,600},{35.00*mm,15.00*mm,594},{35.00*mm,15.00*mm,592},
      {35.00*mm,15.00*mm,664},{32.60*mm,32.30*mm,662},{34.70*mm,17.10*mm,656},{34.15*mm,26.25*mm,696},{0,0,0},
      {33.93*mm,20.42*mm,628},{35.85*mm,17.10*mm,626},{33.80*mm,19.64*mm,624},{30.58*mm,19.83*mm,632},{0,0,0},
      {35.00*mm,15.00*mm,646},{33.95*mm,19.70*mm,644},{30.45*mm,19.80*mm,642},{33.91*mm,26.27*mm,630},{0,0,0},
    };

 
    //LAr Construction Fucntions
    void ArgonOpticalProperties();
    void WLSOpticalProperties();
    G4double LArEpsilon(const G4double lambda);
    G4double LArRefIndex(const G4double lambda);
    G4double LArRayLength(const G4double lambda,const G4double temp);
    G4double ArScintillationSpectrum(const G4double kk);

    // TPB Spectrum from TGraph
    G4double TPBEmissionSpectrum(G4double energy) { 
      G4double wavelength = LambdaE/energy/nm;
      G4double eff=0;
      if(wavelength>350.0 && wavelength < 650.0) eff =fTPBspec->Eval(wavelength);
      return eff;
    }
    
  private:
    G4bool checkOverlaps;
    G4NistManager* nist; 
    DetectorMessenger* detectorMessenger;  // pointer to the Messenger
    TDirectory* fDir;
    TGraph* fTPBspec;
    TH1F* hWLSPhotonE;
    TH1F* hWLSPhotonWavelength;
    TH1F* hArPhotonE;
    TH1F* hArPhotonWavelength;
    
    G4String innerVessel_FillMaterial;

    G4Box* solidWorld;
    G4VPhysicalVolume* physicalWorld;
    G4LogicalVolume* logicalWorld;
    G4LogicalVolume* larSourceLogical;

    //PMT parameters
    G4double grouprmax, groupzmax, WLSHalfThickness; 
    G4LogicalVolume* logicalPmtHousing;  
    G4LogicalVolume* logicalPmtGlass;  
    G4LogicalVolume* logicalPMTWLS; 
    G4VPhysicalVolume* larPhysical;

    //Materials...most are found in materials files
    G4Material* mat_fill;
    G4Material*  fTPB;
    G4OpticalSurface* fPMTGlassOptSurface;
    G4MaterialPropertiesTable *tpbTable;
    G4MaterialPropertiesTable *fPMTGlassOptTable;
    G4Material* mat_ArLiq;
    G4Material* mat_ArCold;
    G4Material* mat_NCold;
    G4Material* mat_NLiq;
    G4Material* Det_mat;
 
    // needed for SDandField
    G4LogicalVolume* logicalWLSDisk;

     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

