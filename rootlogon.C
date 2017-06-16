{
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl;
  gSystem->Load("/home/gold/doubleBeta/LRoot/libLRoot.so");
  gInterpreter->AddIncludePath("." );
  gInterpreter->AddIncludePath("./LRoot" );
// many are old and no longer work, but i use the simple ones often
  gROOT->LoadMacro("/home/gold/root_macros/UseRoot.C");
  cout << " loading: ";
  ("\n this is ROOT \n");
}
