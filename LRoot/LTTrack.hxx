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

using namespace std;


// class to store info for the event

class LTTrack: public TNamed {
	public:
		LTTrack();
		~LTTrack();

    void clear();
    void print();
    // data elements

    Int_t evId;    // event id
    Int_t trkId;
    Int_t parentId;
    Int_t status;
    Int_t nstep;
    Double_t length;
    Double_t stepLength;
    Double_t ke;   // kinetic energy electronvolts
    Double_t edep; // energy deposited in step (electronvolts) 
    Double_t time;   //  event time, microseconds
    Double_t trkTime; // time since track began, microseconds
    TVector3 position;
    TVector3 vertPosition;
    TString physVolName;
    TString particleName;
    TString copy;
    bool isLeaving;
    

ClassDef(LTTrack,1)
};
#endif

