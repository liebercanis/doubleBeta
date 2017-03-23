#include "LTHitSegment.hxx"
ClassImp(LTHitSegment);

LTHitSegment::LTHitSegment() 
{
  clear();
} 

LTHitSegment::~LTHitSegment()
{ 
}

void LTHitSegment::clear() 
{
  
  Contributors.clear();
  PrimaryId=0;
  EnergyDeposit=0;
  SecondaryDeposit=0;
  TrackLength=0;
  StartT=0;
  StopT=0;
  trajIndex=0;
  dedx=0;
  Start.Clear();
  Stop.Clear();
  pdg=0; 
}

void LTHitSegment::print(int iseg){
  printf(" %i) LTHitSegment:  primaryID %i pdg %i contributors %i EnergyDespoit %f  SecondaryDeposit %f TrackLength %f  StopT %f \n",
      iseg,PrimaryId,pdg,int(Contributors.size()),EnergyDeposit,SecondaryDeposit,TrackLength,StopT);
}
