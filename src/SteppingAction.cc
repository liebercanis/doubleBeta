//
// ********************************************************************
// * License and Disclaimer                                           *
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "PMTSD.hh"
#include "GermaniumSD.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "UserEventInformation.hh"
#include "UserTrackInformation.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
:detector(det), eventaction(evt)
{ 

  // create directory 
  fDir = LegendAnalysis::Instance()->topHistDir()->mkdir("step");
  fDir->cd();
  G4cout<<" StepAction working root directory  is  " << G4endl;  
  gDirectory->pwd();
  G4cout << " ... " << G4endl;
  hBoundaryStatus = new TH1F("StepBoundaryStatus"," boundary status ",Dichroic,0,Dichroic); // last in enum G4OpBoundaryProcessStatus
  hParticleType = new TH1F("StepParticleType"," step particle type ",100,0,100);

  // must be in top directory for ChangeFile to work
  LegendAnalysis::Instance()->topTreeDir()->cd();
  ntStep = new TNtuple("ntStep"," step variables ","ev:parent:pdg:status:microsec:length:energy");
  //ntGeStep = new TNtuple("ntGeStep"," step variables ","num:pdg:length:energy");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

  LegendAnalysis::Instance()->getHistFile();
  // get volume of the current step
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4double length = 					step->GetStepLength();
  G4Track* aTrack = 					step->GetTrack();
  //G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* volume = theTouchable->GetVolume();
  G4String volumename = 			volume->GetName();
  G4ThreeVector pos = 				aTrack->GetPosition();
  G4ThreeVector dir = 				aTrack->GetMomentumDirection();
  G4ParticleDefinition* particleType = aTrack->GetDefinition();
  //Used to find other non optical processes
  const G4VProcess * process = aTrack->GetCreatorProcess();
  G4String processName;
  if(process) processName = process->GetProcessName();


  // check if we are in Ge volume
  G4bool  inGeDetector = false;
  G4int GeDetectorNumber=-1;
  G4int GePostNumber =0;
  if (  (volumename.find("B8") != string::npos) ||(volumename.find("P4") != string::npos ) ) {
    GePostNumber = step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber();
    if(GePostNumber!=0) inGeDetector = true;  // remains inside detetor, otherwise it is reflected
    /* no need to do this as 
    This method is invoked by G4SteppingManager when a step is composed in the G4LogicalVolume which has the pointer to this sensitive detector. 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    GermaniumSD* geSD = dynamic_cast<GermaniumSD*>(SDman->FindSensitiveDetector(G4String("GeDetector")));
    */
  }

  if(inGeDetector) {
    GeDetectorNumber = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
    G4cout <<  " stepping action in ge det " << volumename 
      << " copy " << GeDetectorNumber << " post " << GePostNumber << " pdg " << particleType->GetPDGEncoding() << "  " <<  processName <<G4endl;
    //eventaction->FillDetector(GeDetectorNumber,length);
    //ntGeStep->Fill(GeDetectorNumber,particleType->GetPDGEncoding(),length,aTrack->GetKineticEnergy());
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
    G4cout<<"SteppingAction:: WARNING Primary Vertex is out of this world \n\t Ending Stepping Action!"<<G4endl;
    fExpectedNextStatus=Undefined;
    return;
  }

  //This is a primary track 
  // did we miss any secondaries from the primary track?
  trackInformation->SetParentId(aTrack->GetParentID());
  if(aTrack->GetParentID()==0){
    //G4cout<<"SteppingAction::Primary Vertex found "<<G4endl;
    trackInformation->SetPrimary();
    eventInformation->SetPrimaryPhysicalVolume(thePostPV);
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
          } //else if(!(creatorName=="eIoni"||creatorName=="eBrem")) G4cout << " SteppingAction unknown creatorName " << creatorName << G4endl;
        }
      }
    }
  }

 
  // find the optical boundary process only once
  // this is a list of all available processes 
  G4OpBoundaryProcessStatus boundaryStatus = Undefined;
  static G4ThreadLocal G4OpBoundaryProcess* boundary = NULL;

  if(!boundary){
      G4ProcessManager* pm = step->GetTrack()->GetDefinition()->GetProcessManager();
      G4int nprocesses = pm->GetProcessListLength();
      G4ProcessVector* pv = pm->GetProcessList();
     G4cout << "  Stepping action looking for OpBoundary process " << G4endl;
     for(G4int i = 0; i < nprocesses; i++) G4cout << "\t" << i << " process  " << (*pv)[i]->GetProcessName()<< G4endl ;
     for(G4int i = 0; i < nprocesses; i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary" ){
        boundary = dynamic_cast<G4OpBoundaryProcess*>( (*pv)[i] );
        G4cout << "  Stepping action has found OpBoundary " << G4endl;
        break;
      }
    }
  }
  
 
  //Optical Photons
  /*
  G4cout<<"\t SteppingAction:: particleType " 
    << particleType->GetParticleName() 
    << "  type  " << particleType->GetParticleType() 
    << "  PDG " << particleType->GetPDGEncoding() << G4endl;
    */
  hParticleType->Fill(particleType->GetPDGEncoding());

  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){

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
      //if the photon was absorbed in LAr ProcessHit with name defined in DetectorConstruction 
      if(thePrePV->GetName()=="larPhysical"){
        trackInformation->AddTrackStatusFlag(absorbedLAr);
        /*G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdName="ScintSD";
        LegendScintSD* ScintSD = (LegendScintSD*)SDman->FindSensitiveDetector(sdName);
        if(ScintSD){ 
          ScintSD->ProcessHits(step,NULL);
        }
        */
      }
    }
    
    boundaryStatus=boundary->GetStatus();
    
    hBoundaryStatus->Fill(boundaryStatus);
    // G4cout << " Stepping geom boundary process " << boundaryStatus  << G4endl;
    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid
    //Prior to Geant4.6.0-p1 this would not have been enough to check
    /* enum G4OpBoundaryProcessStatus   Undefined,
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
            aTrack->SetTrackStatus(fStopAndKill);
            break;
          }
        case Detection:
          {
            aTrack->SetTrackStatus(fStopAndKill);
            trackInformation->AddTrackStatusFlag(hitPMT);
            //Note, this assumes that the volume causing detection
            // is the photocathode because it is the only one with non-zero efficiency
            //Triger sensitive detector manually since photon is absorbed but status was Detection
            G4SDManager* SDman = G4SDManager::GetSDMpointer();
            G4String sdName="PhotoCathode";//"/LegendDet/pmtSD";
            PMTSD* pmtSD = dynamic_cast<PMTSD*>(SDman->FindSensitiveDetector(sdName));
            if(pmtSD) pmtSD->ProcessHits_constStep(step,NULL);
            else G4cout << " SteppingAction ERROR!!!!!   cannot find PhotoCathode " << G4endl;
            break;
          }
        case FresnelReflection:
        case TotalInternalReflection:
            trackInformation->AddTrackStatusFlag(totalInternal);
        case LambertianReflection:
        case LobeReflection:
        case SpikeReflection:
        case BackScattering:
          trackInformation->IncReflections();
          fExpectedNextStatus=BackScattering;//: StepTooSmall;
          trackInformation->AddTrackStatusFlag(backScatter);
          break;
          //added by Neil
        case NotAtBoundary:
            trackInformation->AddTrackStatusFlag(notBoundary);
        default:
          break;
      }
      // WLS is optical but doesnt seem to correspond to above boundary case 
      if(processName == "OpWLS" ){ 
        trackInformation->AddTrackStatusFlag(hitWLS);
        aTrack->SetTrackStatus(fStopAndKill);
      } 
    }  //end of if(thePostPoint->GetStepStatus()==fGeomBoundary)
  } 
  /*else if(processName == "phot" ){ 
  } else if(processName == "compt"){
  } else if(processName == "eBrem"){
  } else if(processName == "conv"){
  } else if(processName == "Cerenkov"){
  } */



  // scint
  if(processName =="Scintillation") {
    trackInformation->AddTrackStatusFlag(scint);
  }
  // ionizing process
  if(processName == "eIoni" ) {   
    trackInformation->AddTrackStatusFlag(eIoni);
    //if(inGeDetector) G4cout<<"SteppingAction:: eIoni Process Name ... "<<processName<<" boundaryStatus " << boundaryStatus <<G4endl;
  }

  // ionizing process
  if(processName == "hIoni" ) {   
    trackInformation->AddTrackStatusFlag(hIoni);
    if(inGeDetector) G4cout<<"SteppingAction:: hIoni Process Name ... "<<processName<<" boundaryStatus " << boundaryStatus <<G4endl;
  };

  // ionizing process
  if(processName == "ionIoni" ) {   
    trackInformation->AddTrackStatusFlag(ionIoni);
    if(inGeDetector) G4cout<<"SteppingAction:: ionIoni Process Name ... "<<processName<<" boundaryStatus " << boundaryStatus <<G4endl;
  }

  // compt process
  if(processName == "compt" ) {   
    trackInformation->AddTrackStatusFlag(compton);
    if(inGeDetector) G4cout<<"SteppingAction:: compt Process Name ... "<<processName<<" boundaryStatus " << boundaryStatus <<G4endl;
  }

           
  if(inGeDetector) trackInformation->AddTrackStatusFlag(hitGe);
  //if(trackInformation->GetTrackStatus()&hitGe) G4cout << " SteppingAction hitGe  " << G4endl;
  G4int eventId = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  ntStep->Fill(eventId,trackInformation->GetParentId(),particleType->GetPDGEncoding(),
      trackInformation->GetTrackStatus(),aTrack->GetGlobalTime()/microsecond,length,aTrack->GetKineticEnergy());
}


