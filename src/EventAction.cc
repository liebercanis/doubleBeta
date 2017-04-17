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
  fPrintModulo(100000)
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4cout << " **********  end of event action ********** " << anEvent->GetEventID() << G4endl;
  G4cout << "\t number of primary verticies = "<< anEvent->GetNumberOfPrimaryVertex() <<G4endl;
  G4TrajectoryContainer* trajectoryContainer = anEvent->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  G4cout << "\t number of trajectories = "<< n_trajectories <<G4endl;
 // extract the trajectories and draw them
  G4int nOptical =0;
  if(G4VVisManager::GetConcreteInstance()) {
    for (G4int i=0; i<n_trajectories; i++) {
      LegendTrajectory *trj = dynamic_cast<LegendTrajectory*>((*(anEvent->GetTrajectoryContainer()))[i]);
      
      if(trj->IsWLS()) trj->ShowTrajectory(); // print out to G4cout
      
      if(trj->GetParticleName()=="opticalphoton") {
        //trj->SetForceDrawTrajectory(true);
        ++nOptical;
      }
      //trj->DrawTrajectory();
    }
  }
  G4cout << "\t number of optical  = "<< nOptical <<G4endl;
  
  // fill the analysis tree
  LegendAnalysis::Instance()->anaEvent( anEvent );
  //anEvent->Print();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::FillDetector(G4int no, G4double l)
{
  length[no] = l;
	counter++;
}
