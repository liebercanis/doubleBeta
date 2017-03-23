#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TProof.h"
#include "TVectorD.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

# define M_PI           3.14159265358979323846  /* pi */
using namespace std;

void Detectorposition(){
  ifstream in ("Detectorposition.txt");
  TString bla;
  TString name;
  vector<TString> names;
  while(in >> bla >> bla >> bla >> bla >> bla >> bla >> name ){
    //cout << name << endl;
    names.push_back(name);
  }
  in.close();


  int i = 0;
  int j = 0;
  for(j = 0; j<names.size();j++){
    ifstream in2 ("MJD.JSON");
    name = "\"" + names[j] + "\":{";

    while(in2 >> bla){
      if(bla.Contains(name)){
        cout << "---------------found " << names[j];  
        do{
          in2 >> bla;
          if (bla.Contains("r\":\"")) cout << " " << bla; 
          if (bla.Contains("z\":\"")) cout << " " << bla;
        }
        while(!(bla.Contains("]")));
        cout << endl;
        break;
      }
    }
    in2.close();
  }



}
