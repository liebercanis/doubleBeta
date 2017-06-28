/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "LegendAnalysis.hh"
#include "LegendTrajectory.hh"
#include "UserEventInformation.hh"
#include "RunAction.hh"

#include "G4VisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"


#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventAction::EventAction(RunAction* run)
: runAct(run),fPrintModulo(100)
{
  fDir = LegendAnalysis::Instance()->topHistDir()->mkdir("event");
  fDir->cd();
  G4cout<<" EventAction working root directory  is  " << G4endl;  
    
  hEventPhotonLambda = new TH1F("EventPhotonLambda"," wls yield photons/kev ",700,50,750);
  hEventPhotonLambda->GetYaxis()->SetTitle("  photons ");
  hEventPhotonLambda->GetXaxis()->SetTitle("  wavelength (nm) ");
  hEventPhotonCount = new TH1F("EventPhotonCount"," photons detected per event  ",1000,0,1000);
  hEventPhotonCount->GetYaxis()->SetTitle("  events ");
  hEventPhotonCount->GetXaxis()->SetTitle("  photons detected  ");
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
   
  G4EventManager::GetEventManager()->SetUserInformation(new UserEventInformation);
  
  G4int eventNb = event->GetEventID();
  if (eventNb%fPrintModulo == 0) {
    G4cout << "\n************ Begin of event: " << eventNb << G4endl;
  }
  
  // clear the LegendAnalysis event branch 
  LegendAnalysis::Instance()->getEvent()->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* anEvent)
{
  //G4EventManager::GetEventManager()->KeepTheCurrentEvent();
  // fill the analysis tree
  LegendAnalysis::Instance()->anaEvent( anEvent );
  
  G4int nEntries = LegendAnalysis::Instance()->getTree()->GetEntries();
  G4int eventNb = anEvent->GetEventID();
  G4TrajectoryContainer* trajectoryContainer = anEvent->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  if (eventNb%fPrintModulo == 0) {
    G4cout << " **********  EndOfEventAction ********** event " << anEvent->GetEventID() << " **** size of tree *** " << nEntries << G4endl;
      G4cout << "\t number of primary verticies = "<< anEvent->GetNumberOfPrimaryVertex() 
      << " number of trajectories = "<< n_trajectories <<G4endl;
    LegendAnalysis::Instance()->getEvent()->print();
  }
   
  // extract the trajectories and draw them
  /*
  G4int hitCount=0;
  if(G4VVisManager::GetConcreteInstance()) {
    for (G4int i=0; i<n_trajectories; i++) {
      LegendTrajectory *trajectory = dynamic_cast<LegendTrajectory*>((*(anEvent->GetTrajectoryContainer()))[i]);
      const G4Track* aTrack= trajectory->GetTrack();
      if(!aTrack) {
        G4cout << "  EventAction:: no track for this trajectory  " << i << endl;
        continue;
      }
        G4cout << " EventAction dynamic cast  " << i << G4endl;
      UserTrackInformation* trackInformation= dynamic_cast<UserTrackInformation*>(aTrack->GetUserInformation());
      if(trackInformation->GetPostStepStatusLast() != TrackPostStepStatus::isGeomBoundary)  {
        G4cout << " EventAction skipping track with step status " << trackInformation->GetPostStepStatusLast() << G4endl;
        continue;
      }
      
      //if(trajectory->IsWLS()) trajectory->ShowTrajectory(); // print out to G4cout
      if(trajectory->IsPmtHit()) {
        G4double ke = trajectory->GetInitialMomentum().mag();
        hEventPhotonLambda->Fill(CLHEP::h_Planck*CLHEP::c_light/ke/nm); 
        ++hitCount;
      }
      //trajectory->DrawTrajectory();
    }
    hEventPhotonCount->Fill(double(hitCount));
  }
  */

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::FillDetector(G4int no, G4double l)
{
}
