using namespace std;
vector<TString> names;

void detpos()
{
  string input;

  double mm = 1E+1;
  int Detector_number;
  TString Detector_name;
  TString Detector_name_sol;
  TString Detector_name_log;
  TVector3 Detector_position;
  int Detector_slices;
  
  vector<string> Detectors;
  vector<double> Detector_r;
  vector<double> Detector_z;

  TH1F* hRpos=new TH1F("Rpos"," r pos ",500,0,500);
  TH1F* hZpos=new TH1F("Zpos"," z pos ",2000,-1000,1000);
   
  ifstream infile("Detectorposition.txt",ios_base::in);
  if (! infile.is_open()) cout<< "DetectorConstruction unable to open Detectorposition.txt" << endl;
  else  cout<< "DetectorConstruction opened Detectorposition.txt" << endl;

  while(getline(infile,input)){
    Detectors.push_back(input);
  }
  infile.close();
  cout << "DetectorConstruction Detector size = " << Detectors.size() <<endl;

  int j =0;
  double sum_x = 0;
  double sum_y = 0;
  double sum_z = 0;
  
  
  for (int i =0;i< Detectors.size();i++){
    printf(" Detectors[%i]=%s\n",i,Detectors[i].c_str());
    Detector_number = 0;
    Detector_r.clear();
    Detector_z.clear();
    stringstream ss(Detectors[i]);
    string s;
    j = 0;
    while (ss >> s) {
     j++;
     if (j==1) Detector_number += atoi(s.c_str())*100;
     else if (j==2) Detector_number += atoi(s.c_str())*10;
     else if (j==3) Detector_number += atoi(s.c_str())*1;
     else if (j==4) Detector_position.SetX(atof(s.c_str())*mm);
     else if (j==5) Detector_position.SetY(atof(s.c_str())*mm);
     else if (j==6) Detector_position.SetZ(atof(s.c_str())*mm);
     else if (j==7) Detector_name = s;
     else if (j==8) Detector_slices = atoi(s.c_str());
     else if (j==9)  Detector_r.push_back((double) atof(s.c_str())*mm);
     else if (j==10) Detector_z.push_back((double) atof(s.c_str())*mm);
     else if (j==11) Detector_r.push_back((double) atof(s.c_str())*mm);
     else if (j==12) Detector_z.push_back((double) atof(s.c_str())*mm);
     else if (j==13) Detector_r.push_back((double) atof(s.c_str())*mm);
     else if (j==14) Detector_z.push_back((double) atof(s.c_str())*mm);
    }
    //const G4double r[] = {Detector_r[0],Detector_r[1],Detector_r[2]};
    //const G4double z[] = {Detector_z[0],Detector_z[1],Detector_z[2]};
    printf(" r %f %f %f \n",Detector_r[0],Detector_r[1],Detector_r[2]);
    printf(" z %f %f %f \n",Detector_z[0],Detector_z[1],Detector_z[2]);
    hRpos->Fill(Detector_r[0]);
    hRpos->Fill(Detector_r[1]);
    hRpos->Fill(Detector_r[2]);
    hZpos->Fill(Detector_z[0]);
    hZpos->Fill(Detector_z[1]);
    hZpos->Fill(Detector_z[2]);
    
    
    Detector_name_sol = Detector_name + "_sol";
    Detector_name_log = Detector_name + "_log";
    if(Detector_slices==2){
      Detector_r.push_back(Detector_r[1]+0.00001*mm);
      Detector_z.push_back(Detector_z[1]+0.00001*mm);
    }
    sum_x += Detector_position.x() /Detectors.size();
    sum_y += Detector_position.y() /Detectors.size();
    sum_z += Detector_position.z() /Detectors.size();
    // max extent of this detector is 
  }

  cout<< "\t *******************************************" <<endl;
  cout<< "\t DetectorConstruction -- Construct finished" <<endl;
  printf(" mean x %f mean y %f mean z %f \n",sum_x,sum_y,sum_z);
  cout<< "\t *******************************************" <<endl;

}
