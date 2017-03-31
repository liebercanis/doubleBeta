//
// ********************************************************************
// * License and Disclaimer                                           *
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "UserEventInformation.hh"
#include "UserTrackInformation.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
:detector(det), eventaction(evt)
{ 

  // create directory 
  fDir = LegendAnalysis::Instance()->topDir()->mkdir("step");
  fDir->cd();
  G4cout<<" StepAction working root directory  is  " << G4endl;  
  gDirectory->pwd();
  G4cout << " ... " << G4endl;
  G4double LowE = 1.7712*eV;//700 nm
  G4double HighE = 12.3984*eV;//100 nm
  hWLSPhotonE = new TH1F("StepWLSPhotonE"," photon energy from WLS",1000,LowE,HighE);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

  // get volume of the current step
  G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4String volumename = 			volume->GetName();
  G4Track* aTrack = 					step->GetTrack();
  G4ThreeVector pos = 				aTrack->GetPosition();
  G4ThreeVector dir = 				aTrack->GetMomentumDirection();
  G4double length = 					step->GetStepLength();

  //G4cout << " stepping volume " << volumename << " length " << length << G4endl; 

  // check if we are in scoring volume
  if (volumename(0,4) == "source"){
    G4ThreeVector dir2 = eventaction->direction;
    if (dir2.x()*dir2.x()+dir2.y()*dir2.y()+dir2.z()*dir2.z() == 0){
      eventaction->direction = dir;
      eventaction->position = pos;
      //G4cout << " d " <<dir.cosTheta() << " " << dir.phi() << G4endl;
      //G4cout << " p " <<pos.cosTheta() << " " << pos.phi() << G4endl;
    }
  }


  if (volumename.contains("P") || volumename.contains("B")){
    G4int k = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
    //G4cout <<  k << " " <<length/cm << G4endl;
    eventaction->FillDetector(k,length);
  }

  /*************************************
  ** from LegendStepping Action 
  *************************************/

  if ( aTrack->GetCurrentStepNumber() == 1 ) fExpectedNextStatus = Undefined;
  UserTrackInformation* trackInformation = (UserTrackInformation*) aTrack->GetUserInformation();
  UserEventInformation* eventInformation = (UserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();
  

  G4StepPoint* thePrePoint = step->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

  G4StepPoint* thePostPoint = step->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
  
  if(!thePostPV){//out of the world
    G4cout<<"SteppingAction::Primary Vertex is out of this world \n\t Ending Stepping Action!"<<G4endl;
    fExpectedNextStatus=Undefined;
    return;
  }

  //This is a primary track 
  // did we miss any secondaries from the primary track?
  if(aTrack->GetParentID()==0){
    G4TrackVector* fSecondary = fpSteppingManager->GetfSecondary();
    G4int tN2ndariesTot = fpSteppingManager->GetfN2ndariesAtRestDoIt()
      + fpSteppingManager->GetfN2ndariesAlongStepDoIt()
      + fpSteppingManager->GetfN2ndariesPostStepDoIt();

    //If we havent already found the conversion position and there were
    //Loop over all 2ndaries that have not been found with N2ndariesTot
    if(!eventInformation->IsConvPosSet() && tN2ndariesTot>0 ){
      for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; lp1<(*fSecondary).size(); lp1++){
        const G4VProcess* creator=(*fSecondary)[lp1]->GetCreatorProcess();
        if(creator){
          G4String creatorName=creator->GetProcessName();
          //Added Scint to list -Neil
          if(creatorName=="phot"||creatorName=="compt"||creatorName=="conv"||creatorName=="Scintillation"){
            //since this is happening before the secondary is being tracked
            //the Vertex position has not been set yet(set in initial step)
            //so set Conversion Position
            eventInformation->SetConvPos((*fSecondary)[lp1]->GetPosition());
          } else if(!(creatorName=="eIoni"||creatorName=="eBrem")) G4cout << " SteppingAction unknown creatorName " << creatorName << G4endl;
        }
      }
    }
  }


  G4OpBoundaryProcessStatus boundaryStatus = Undefined;
  static G4ThreadLocal G4OpBoundaryProcess* boundary = NULL;

  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm = step->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i = 0; i < nprocesses; i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary" ){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }

  //Used to find othe non optical processes
  const G4VProcess * process = aTrack->GetCreatorProcess();
  
  G4String processName;
  if(process) processName = process->GetProcessName();

  G4ParticleDefinition* particleType = aTrack->GetDefinition();

  //Optical Photons
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){

    //Need local definition for ScintSDHit processing
    G4Step* step = const_cast<G4Step*>(step);

    //Kill photons exiting cryostat
    if(thePostPV->GetName()=="phy_World"){
      aTrack->SetTrackStatus(fStopAndKill);
      //eventInformation->IncPhotonCount_Escape();
      return;
    }

    //The photon was absorbed at another place other than a boundry
    if(thePostPoint->GetProcessDefinedStep()->GetProcessName()=="OpAbsorption"){
      //eventInformation->IncAbsorption();
      trackInformation->AddTrackStatusFlag(absorbed);
      //if the photon was absorbed in LAr ProcessHit
      /*
      ** if(thePrePV->GetName()=="phy_fillGas"){
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdName="ScintSD";
        LegendScintSD* ScintSD = (LegendScintSD*)SDman->FindSensitiveDetector(sdName);
        if(ScintSD){ 
          ScintSD->ProcessHits(step,NULL);
        }
        trackInformation->AddTrackStatusFlag(absorbedLAr);
      }
        */
    }
    
    boundaryStatus=boundary->GetStatus();
    
    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid
    //Prior to Geant4.6.0-p1 this would not have been enough to check
    /* enum G4OpBoundaryProcessStatus {  Undefined,
                                  Transmission, FresnelRefraction,
                                  FresnelReflection, TotalInternalReflection,
                                  LambertianReflection, LobeReflection,
                                  SpikeReflection, BackScattering,
                                  Absorption, Detection, NotAtBoundary,
                                  SameMaterial, StepTooSmall, NoRINDEX,
      .... and more in G4OpBoundaryProcess.hh
    */
    if(thePostPoint->GetStepStatus()==fGeomBoundary){
      if(fExpectedNextStatus==StepTooSmall){
        if(boundaryStatus!=StepTooSmall){
          G4cout<< "SteppingAction::UserSteppingAction(): No reallocation step after reflection!"<<G4endl;          
          G4cout<<"SteppinAction:: thePrePV of Process is :: "<< thePrePV->GetName()<<G4endl;
          G4cout<<"SteppinAction:: thePostPV of Process is :: "<< thePostPV->GetName()<<G4endl;
          G4cout<<"\t    >>>>>>>>>> Something is wrong with the surface normal or geometry....Track is killed"<<G4endl;

          aTrack->SetTrackStatus(fStopAndKill);
        }
      }
      fExpectedNextStatus=Undefined;
      switch(boundaryStatus){
        case Absorption:
          {
            //This all Transportation
            trackInformation->AddTrackStatusFlag(boundaryAbsorbed);
            eventInformation->IncBoundaryAbsorption();
            break;
          }
        case Detection:
          {
            //Note, this assumes that the volume causing detection
            //is the photocathode because it is the only one with
            //non-zero efficiency
            //Triger sensitive detector manually since photon is
            //absorbed but status was Detection
            /*G4SDManager* SDman = G4SDManager::GetSDMpointer();
              G4String sdName="PhotoCathode";//"/LegendDet/pmtSD";
              LegendPMTSD* pmtSD = (LegendPMTSD*)SDman->FindSensitiveDetector(sdName);
              if(pmtSD) pmtSD->ProcessHits_constStep(step,NULL);
               */

            trackInformation->AddTrackStatusFlag(hitPMT);
            break;
          }

        case FresnelReflection:
        case TotalInternalReflection:
        case LambertianReflection:
        case LobeReflection:
        case SpikeReflection:
        case BackScattering:
          trackInformation->IncReflections();
          fExpectedNextStatus=BackScattering;//: StepTooSmall;
          break;
          //added by Neil
        case NotAtBoundary:
        default:
          break;
      }
    }  //end of if(thePostPoint->GetStepStatus()==fGeomBoundary)
    //end of if OpticalPhoton
  } else if(processName == "phot" ){ 
  } else if(processName == "eIoni"){
  } else if(processName == "compt"){
  } else if(processName != ""){
    G4cout<<"LegendSteppingAction:: Process Name that Neil could not find is ... "<<processName<<" !!!"<<G4endl;
  }
}


