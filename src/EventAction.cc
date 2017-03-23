/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

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
  G4int eventNb = event->GetEventID();
  if (eventNb%fPrintModulo == 0) {
    G4cout << "\n---> Begin of event: " << eventNb << G4endl;
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

void EventAction::EndOfEventAction(const G4Event* /*event*/)
{
	if (counter) runAct->Fill(direction,position,counter,length);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::FillDetector(G4int no, G4double l)
{
  length[no] = l;
	counter++;
}
