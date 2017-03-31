using namespace std;
vector<TString> names;

void readJson() 
{
  int i = 0;
  int j = 0;
  TString bla;
  TString name;
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


void Detectorposition(){
  ifstream in ("Detectorposition.txt");
  TString bla;
  TString name;
  char line[180];
  while(in >> bla >> bla >> bla >> bla >> bla >> bla >> name ){
    names.push_back(name);
    in.getline(line,100);
  }
  in.close();

  printf(" number of detectors is %i\n",(int) names.size());
  for(unsigned in=0; in<names.size(); ++in) printf(" %s \n",names[in].Data());

  readJson();
}

