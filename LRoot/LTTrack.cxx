#include "LTTrack.hxx"
ClassImp(LTTrack)

LTTrack::LTTrack(): TNamed("LTTrack","LTTrack")
{
  clear();
}

LTTrack::~LTTrack()
{
}

void LTTrack::clear()
{

  evId=0;    // event id
  trkId=0;    // track id
  parentId=0;    // parent id
  status=0;
  time=0;  
  ke=0;
  length=0;
  nstep=0;
  stepLength=0;
  position.Clear();
  vertPosition.Clear();
  traject.Clear(); 
}
void LTTrack::print(){
  printf(" ********************  LTTrack ********************* \n");
  printf(" \tevent  %i id %i parent %i step %i time %1.3E KE %1.3E \n",evId,trkId,parentId,nstep,time,ke);
  printf(" *************************************************** \n");
}