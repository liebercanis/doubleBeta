/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"


//root includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <sys/time.h>

#include "G4ThreeVector.hh"

class G4Run;
class DetectorConstruction;
class PrimaryGeneratorAction;
class RunMessenger;
using namespace std;

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

	void Fill(G4ThreeVector, G4ThreeVector, G4int, G4double [300]);

  G4String filename;
  private:


	TFile* fRootFout;
	TTree* fRootTree;

	DetectorConstruction   *detector;
	PrimaryGeneratorAction *primary;
   RunMessenger           *runMessenger;


    timeval start, end;


	double start_x;
	double start_y;
	double start_z;
	double start_theta;
	double start_phi;

	double start_px;
	double start_py;
	double start_pz;
	double start_ptheta;
	double start_pphi;
	vector<int> 		D;
	vector<double> 	L;
	int number;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

