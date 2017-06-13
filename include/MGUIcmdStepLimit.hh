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
// $Id: MGUIcmdStepLimit.hh,v 1.1 2009-01-17 15:26:59 jasondet Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef MGUIcmdStepLimit_H
#define MGUIcmdStepLimit_H 1

#include "G4UIcommand.hh"

// class description:

class MGUIcmdStepLimit : public G4UIcommand
{
  public: 
    MGUIcmdStepLimit(const char* commandPath, G4UImessenger* messenger);
    G4double GetStepSize(const char* paramString);
    G4String GetParticleName(const char* paramString);
    G4String GetVolumeName(const char* paramString);
};

#endif
