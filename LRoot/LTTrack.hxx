/**
** MG, May 2017 
**/
#ifndef LTTRACK_DEFINED
#define LTTRACK_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <TTree.h>
#include <TVector3.h>
#include <vector>
// LRoot classes
//#include "LTPVertex.hxx"
//#include "TBits.h"

using namespace std;


// class to store info for the event

class LTTrack: public TNamed {
	public:
		LTTrack();
		~LTTrack();
    void clear();
    void print();
    void addPositionHistory( Double_t px, Double_t py, Double_t pz); 
    // data elements

    Int_t evId;    // event id
    Int_t trkId;
    Int_t parentId;
    Int_t status;
    Int_t pdg;
    std::vector<Int_t> boundaryStatus;
    std::vector<std::string> boundaryName;
    Int_t nstep;
    Int_t copy;
    Double_t length;
    Double_t stepLength;
    Double_t ke;   // kinetic energy electronvolts
    Double_t edep; // energy deposited in step (electronvolts) 
    Double_t time;   //  event time, microseconds
    Double_t trkTime; // time since track began, microseconds
    TVector3 position; //end
    std::vector<TVector3> positionHistory;
    std::vector<Double_t> positionEnergy;
    TVector3 vertPosition;
    TString process;
    TString physVolName;
    TString particleName;
    TString preName;
    TString postName;
    Int_t nInToGe;
    Int_t nOutOfGe;
    Int_t nSpike;
    bool isLeaving;

ClassDef(LTTrack,11)
};
#endif

