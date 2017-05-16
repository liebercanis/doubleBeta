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
: runAct(run),
  fPrintModulo(100)
{
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

	direction.setX(0);
	direction.setY(0);
	direction.setZ(0);
	position.setX(0);
	position.setY(0);
	position.setZ(0);
	counter = 0;

  for (int k=0;k<300;k++) length[k] = 0.;

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
  if (eventNb%fPrintModulo == 0) {
    G4cout << " **********  EndOfEventAction ********** event " << anEvent->GetEventID() << " **** size of tree *** " << nEntries << G4endl;
    G4TrajectoryContainer* trajectoryContainer = anEvent->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if(trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    G4cout << "\t number of primary verticies = "<< anEvent->GetNumberOfPrimaryVertex() 
      << " number of trajectories = "<< n_trajectories <<G4endl;
    LegendAnalysis::Instance()->getEvent()->print();
  }
  // extract the trajectories and draw them
  /*if(G4VVisManager::GetConcreteInstance()) {
    for (G4int i=0; i<n_trajectories; i++) {
      LegendTrajectory *trajectory = dynamic_cast<LegendTrajectory*>((*(anEvent->GetTrajectoryContainer()))[i]);
      //if(trajectory->IsWLS()) trajectory->ShowTrajectory(); // print out to G4cout
      if(trajectory->GetParticleName()=="opticalphoton") {
        //trajectory->SetForceDrawTrajectory(true);
      }
      //trajectory->DrawTrajectory();
    }
  }*/

   if(LegendAnalysis::Instance()->getEvent()->nGeHits>0) LegendAnalysis::Instance()->getEvent()->print();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::FillDetector(G4int no, G4double l)
{
  length[no] = l;
	counter++;
}
