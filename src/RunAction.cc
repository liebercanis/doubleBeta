#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "RunMessenger.hh"
//root includes

#include <iostream>
#include <time.h>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det,PrimaryGeneratorAction* prim)
:detector(det),primary(prim),G4UserRunAction()
{
 filename = "test";
 runMessenger = new RunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  G4String filename2 = filename + ".root";
  fRootFout = new TFile(filename2,"RECREATE",filename) ;
  fRootTree = new TTree("tree",filename2);
  fRootTree->SetMaxTreeSize(1073741824);


  fRootTree->Branch("number",&number,"number/I");

	fRootTree->Branch("start_x",&start_x,"start_x/D");
	fRootTree->Branch("start_y",&start_y,"start_y/D");
  fRootTree->Branch("start_z",&start_z,"start_z/D");
	fRootTree->Branch("start_phi",&start_phi,"start_phi/D");
  fRootTree->Branch("start_theta",&start_theta,"start_theta/D");

	fRootTree->Branch("start_px",&start_px,"start_px/D");
	fRootTree->Branch("start_py",&start_py,"start_py/D");
  fRootTree->Branch("start_pz",&start_pz,"start_pz/D");
	fRootTree->Branch("start_pphi",&start_pphi,"start_pphi/D");
  fRootTree->Branch("start_ptheta",&start_ptheta,"start_ptheta/D");

	fRootTree->Branch("detector",&D);
  fRootTree->Branch("length",&L);

  gettimeofday(&start, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{

  G4int nofEvents = aRun->GetNumberOfEvent();
  gettimeofday(&end, 0);

  if (nofEvents == 0) return;

  fRootFout->cd();
  //fRootTree->Print();
  fRootTree->Write();
  fRootFout->Close();

  G4cout <<nofEvents<<" Events in "<< end.tv_sec-start.tv_sec<<"s"<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::Fill(G4ThreeVector r_dir, G4ThreeVector r_pos, G4int r_counter, G4double r_length[80])
{

	D.clear();
	L.clear();

	for (int i = 0; i<300;i++){
		if (r_length[i]){
			L.push_back(r_length[i]);
			D.push_back(i);
		}
	}
	number = r_counter;
	start_x =  r_pos.x();
	start_y =  r_pos.y();
	start_z =  r_pos.z();
	start_theta =  r_pos.theta();
	start_phi =  r_pos.phi();

	start_px =  r_dir.x();
	start_py =  r_dir.y();
	start_pz =  r_dir.z();
	start_ptheta =  r_dir.theta();
	start_pphi =  r_dir.phi();

	fRootTree->Fill();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
