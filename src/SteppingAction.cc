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
   GeDebug = false;
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
  ntStep = new TNtuple("ntStep"," step variables ","ev:parent:pdg:status:tglobal:tlocal:length:energy");
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

  // determine if step is at Ge boundary
  G4String preName = thePrePV->GetName();
  G4String postName = thePostPV->GetName();
  G4bool isPreGroup1   =  preName.contains("group1Physical") && thePrePoint->GetStepStatus() == fGeomBoundary;
  G4bool isPostGroup1  =  postName.contains("group1Physical") && thePostPoint->GetStepStatus() == fGeomBoundary;
  G4bool isPreGe  = (preName.contains("B8")||preName.contains("P4")) && thePrePoint->GetStepStatus() == fGeomBoundary;
  G4bool isPostGe = (postName.contains("B8")||postName.contains("P4")) && thePostPoint->GetStepStatus() == fGeomBoundary;
  G4bool isInToGe = isPreGroup1&&isPostGe;
  G4bool isOutOfGe = isPreGe&&isPostGroup1;
  
  if(!thePostPV){//out of the world
    G4cout<<"SteppingAction:: WARNING Primary Vertex is out of this world \n\t Ending Stepping Action!"<<G4endl;
    fExpectedNextStatus=Undefined;
    return;
  }

  //This is a primary track 
  // did we miss any secondaries from the primary track?
  trackInformation->SetParentId(aTrack->GetParentID());
  trackInformation->SetProcessName(processName);
  trackInformation->SetPreName(thePrePV->GetName());
  trackInformation->SetPostName(thePostPV->GetName());
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
     //G4cout << "  Stepping action looking for OpBoundary process " << G4endl;
     //for(G4int i = 0; i < nprocesses; i++) G4cout << "\t" << i << " process  " << (*pv)[i]->GetProcessName()<< G4endl ;
     for(G4int i = 0; i < nprocesses; i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary" ){
        boundary = dynamic_cast<G4OpBoundaryProcess*>( (*pv)[i] );
        G4cout << "  Stepping action has found OpBoundary " << G4endl;
        break;
      }
    }
  }
  if(step == NULL){
    G4cout<<"NULL Step in steppingaction!"<<G4endl;
  }

  // this is private boundary->BoundaryProcessVerbose();
  // check if we are in Ge volume
  G4bool  inGeDetector = false;
  G4int GeDetectorNumber=-1;
  G4int GePostNumber =0;
  if (  (volumename.find("B8") != string::npos) ||(volumename.find("P4") != string::npos ) ) {
    GePostNumber = step->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber();
    if(GePostNumber!=0){
      inGeDetector = true;
      /*if(step->GetTotalEnergyDeposit() > 1.0*eV){
        G4cout <<  " SteppingAction:: Step inside GeDet...Energy Deposited = "<<step->GetTotalEnergyDeposit()<<
          ", by "<<particleType->GetPDGEncoding()<<", with process name "<<processName<<
          ", Total energy for Track "<<aTrack->GetTotalEnergy()<<G4endl;
      }*/
      //eventaction->FillDetector(GeDetectorNumber,length);
      //ntGeStep->Fill(GeDetectorNumber,particleType->GetPDGEncoding(),length,aTrack->GetKineticEnergy());
      
    } 
  }


  if(boundary) {
    boundaryStatus=boundary->GetStatus();
    hBoundaryStatus->Fill(boundaryStatus);
    trackInformation->AddBoundaryProcessStatus(boundaryStatus);
    trackInformation->AddBoundaryName(preName);
  }

  // add position and energy points of step to vectors.
  trackInformation->AddPositionHistory(aTrack->GetPosition());
  trackInformation->AddPositionEnergy(step->GetTotalEnergyDeposit());
  trackInformation->AddStepLength(step->GetStepLength());
  
 
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
      aTrack->SetTrackStatus(fStopAndKill); // DO I NEED TO DO THIS BY HAND????
      trackInformation->AddTrackStatusFlag(absorbed);
      eventInformation->IncAbsorption();
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
    
 
    //G4cout << " Stepping Action  " <<  boundaryStatus << " " << trackInformation->GetBoundaryProcessStatus() << G4endl;
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
      // case are evaluated if boundaryStatus == case 

      //if( isInToGe) G4cout<<"SteppingAction::  isInTo Ge geom boundary "<< thePrePV->GetName()<< "->" << thePostPV->GetName() << " status " << boundaryStatus << G4endl;
      //if( isOutOfGe) G4cout<<"SteppingAction::  isOutOf Ge geom boundary "<< thePrePV->GetName()<< "->"  << thePostPV->GetName() << " status "<< boundaryStatus << G4endl;

      if( isInToGe) trackInformation->IncInToGe();
      if( isOutOfGe) trackInformation->IncOutOfGe();
      // breaks not needed in case statement, but leaving in ... M.Gold
      switch(boundaryStatus){
        case Absorption: 
          {
            //This all Transportation
            trackInformation->AddTrackStatusFlag(boundaryAbsorbed);
            eventInformation->IncAbsorption();
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
          {
  //G4cout<<"SteppingAction::  FresnelReflection reflectivity = "<< boundary->GetReflectivity() << "  " << thePrePV->GetName()<< "->"  << thePostPV->GetName() << " status "<< boundaryStatus << G4endl;
            
            trackInformation->AddTrackStatusFlag(fresnelReflect);
            break;
          }
        case TotalInternalReflection:
          {
            trackInformation->AddTrackStatusFlag(totalInternal);
            break;
          }
        case LambertianReflection:
          break;
        case LobeReflection:
          break;
        case SpikeReflection:
          {
            trackInformation->IncSpikeReflection();
            /*
              G4double LambdaE = 2.0*TMath::Pi()*1.973269602e-10 * mm * MeV;// in default Mev-mm units
              G4cout<<"SteppingAction::  SpikeReflection wave " << LambdaE/aTrack->GetKineticEnergy()/nanometer <<  " (nm)  reflectivity = "
              << boundary->GetReflectivity() << " angle " <<  boundary->GetIncidentAngle()<< " @ "  
              << thePrePV->GetName()<< "->"  << thePostPV->GetName() << " status "<< boundaryStatus << G4endl;
              */
            break;
          }
        case StepTooSmall:  // what to do here?gg
          break;
        case BackScattering:
          {
            trackInformation->IncReflections();
            fExpectedNextStatus= BackScattering;//: StepTooSmall; MG problem with VM2000? if StepTooSmall ???
            trackInformation->AddTrackStatusFlag(backScatter);
            break;
          }
          //added by Neil
        case NotAtBoundary:
          {
            trackInformation->AddTrackStatusFlag(notBoundary);
            break;
          }
        default:
          break;
      }


      // WLS is optical but doesnt seem to correspond to above boundary case
      if(processName == "OpWLS" ){ 
        trackInformation->AddTrackStatusFlag(hitWLS);
        // G4cout << " SteppingAction OpWLS boundary status = %i " << boundaryStatus << G4endl;
        //aTrack->SetTrackStatus(fStopAndKill);  why stop and kill this?
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
    if(inGeDetector && GeDebug) G4cout<<"SteppingAction:: eIoni Process Name ... "<<processName<<" boundaryStatus " << boundaryStatus <<G4endl;
  }

  // ionizing process
  if(processName == "hIoni" ) {   
    trackInformation->AddTrackStatusFlag(hIoni);
    if(inGeDetector && GeDebug) G4cout<<"SteppingAction:: hIoni Process Name ... "<<processName<<" boundaryStatus " << boundaryStatus <<G4endl;
  }

  // ionizing process
  if(processName == "ionIoni" ) {   
    trackInformation->AddTrackStatusFlag(ionIoni);
    if(inGeDetector && GeDebug) G4cout<<"SteppingAction:: ionIoni Process Name ... "<<processName<<" boundaryStatus " << boundaryStatus <<G4endl;
  }

  // compt process
  if(processName == "compt" ) {   
    trackInformation->AddTrackStatusFlag(compton);
    if(inGeDetector && GeDebug) G4cout<<"SteppingAction:: compt Process Name ... "<<processName<<" boundaryStatus " << boundaryStatus <<G4endl;
  }

           
  if(inGeDetector) trackInformation->AddTrackStatusFlag(hitGe);
  if(trackInformation->GetTrackStatus()&hitGe && GeDebug) G4cout << " SteppingAction hitGe...TrackStatus "<<trackInformation->GetTrackStatus() << G4endl;
  
  
  G4int eventId = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  ntStep->Fill(eventId,trackInformation->GetParentId(),particleType->GetPDGEncoding(),
      trackInformation->GetTrackStatus(),aTrack->GetGlobalTime()/microsecond,aTrack->GetLocalTime()/microsecond,length,aTrack->GetKineticEnergy());
}


