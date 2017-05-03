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
/// \file optical/Legend/include/UserTrackInformation.hh
/// \brief Definition of the UserTrackInformation class
//
#include "G4VUserTrackInformation.hh"
#include "globals.hh"

#ifndef UserTrackInformation_h
#define UserTrackInformation_h 1
//IF you are going to add a status it must be a power of 2 because of the bitwise
//functionallity of the AddTrackStatus Function
//--Neil
/*  note these defined in enum G4TrackStatus
  fAlive,             // Continue the tracking
  fStopButAlive,      // Invoke active rest physics processes and
                      // and kill the current track afterward
  fStopAndKill,       // Kill the current track
  fKillTrackAndSecondaries,
                      // Kill the current track and also associated
                      // secondaries.
  fSuspend,           // Suspend the current track
  fPostponeToNextEvent
                      // Postpones the tracking of thecurrent track 
                      // to the next event.

 and these are in enum G4OpBoundaryProcessStatus {  Undefined,
                                  Transmission, FresnelRefraction,
                                  FresnelReflection, TotalInternalReflection,
                                  LambertianReflection, LobeReflection,
                                  SpikeReflection, BackScattering,
                                  Absorption, Detection, NotAtBoundary,
                                  SameMaterial, StepTooSmall, NoRINDEX,
                                  PolishedLumirrorAirReflection,
                                  PolishedLumirrorGlueReflection,
                                  PolishedAirReflection,
                                  PolishedTeflonAirReflection,
                                  PolishedTiOAirReflection,
                                  PolishedTyvekAirReflection,
                                  PolishedVM2000AirReflection,
                                  PolishedVM2000GlueReflection,
                                  EtchedLumirrorAirReflection,
                                  EtchedLumirrorGlueReflection,
                                  EtchedAirReflection,
                                  EtchedTeflonAirReflection,
                                  EtchedTiOAirReflection,
                                  EtchedTyvekAirReflection,
                                  EtchedVM2000AirReflection,
                                  EtchedVM2000GlueReflection,
                                  GroundLumirrorAirReflection,
                                  GroundLumirrorGlueReflection,
                                  GroundAirReflection,
                                  GroundTeflonAirReflection,
                                  GroundTiOAirReflection,
                                  GroundTyvekAirReflection,
                                  GroundVM2000AirReflection,
                                  GroundVM2000GlueReflection,
                                  Dichroic };
*/

/*TrackStatus:
  active: still being tracked
  hitPMT: stopped by being detected in a PMT
  absorbed: stopped by being absorbed with G4OpAbsorbtion
  boundaryAbsorbed: stopped by being aborbed with G4OpAbsorbtion
  hitSphere: track hit the sphere at some point
  inactive: track is stopped for some reason
   -This is the sum of all stopped flags so can be used to remove stopped flags
*/
enum TrackStatus { active=1, hitPMT=2, absorbed=4, boundaryAbsorbed=8,
                      absorbedLAr=16, inactive=32, hitWLS = 64, totalInternal=128, backScatter=256, notBoundary=512,
                      scint=2*notBoundary, 
                      eIoni=2*scint, 
                      hIoni=2*eIoni, 
                      ionIoni=2*hIoni, 
                      hitGe=2*ionIoni, isBad=2*hitGe};


enum TrackBit {MaxHistogramBit=17};

class UserTrackInformation : public G4VUserTrackInformation
{
  public:

    UserTrackInformation();
    virtual ~UserTrackInformation();

    //Sets the track status to s (does not check validity of flags)
    void SetTrackStatusFlags(int s){fStatus=s;}
    //Does a smart add of track status flags (disabling old flags that conflict)
    //If s conflicts with itself it will not be detected
    void AddTrackStatusFlag(int s);
 
    int GetTrackStatus()const {return fStatus;}

    // just return first bit of interest for histograming
    int GetTrackBit() const { 
      //photon 
      if(fStatus&hitPMT) return 1;
      else if(fStatus&absorbed) return 2; 
      else if(fStatus&boundaryAbsorbed) return 3;
      else if(fStatus&absorbedLAr) return 4;
      else if(fStatus&hitWLS) return 5;
      else if(fStatus&totalInternal) return 6;
      else if(fStatus&backScatter) return 7;
      // ionization
      else if(fStatus&hitGe) return 11;
      else if(fStatus&eIoni) return 12;
      else if(fStatus&hIoni) return 13;
      else if(fStatus&ionIoni) return 14;
      else if(fStatus&scint) return 15;
      else if(fStatus&notBoundary) return 16;
      else return MaxHistogramBit;
    }
 
    void IncReflections(){fReflections++;}
    G4int GetReflectionCount()const {return fReflections;}

    void SetForceDrawTrajectory(G4bool b){fForcedraw=b;}
    G4bool GetForceDrawTrajectory(){return fForcedraw;}

    inline virtual void Print() const{};

    void SetPrimary() { fPrimary=true; }
    G4bool IsPrimary() { return fPrimary;}
    void SetParentId( G4int id) { fParentId = id;}
    G4int GetParentId() { return fParentId;}

  private:

    G4int fStatus;
    G4int fParentId;
    G4bool fPrimary;
    G4int  fReflections;
    G4bool fForcedraw;
};

#endif
