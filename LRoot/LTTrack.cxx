#include "LTTrack.hxx"
ClassImp(LTTrack)

LTTrack::LTTrack(): TNamed("LTTrack","LTTrack")
{
  clear();
}

LTTrack::~LTTrack()
{
}

void LTTrack::addPositionHistory( Double_t px, Double_t py, Double_t pz) 
{
  TVector3 p(px,py,pz);
  positionHistory.push_back(p);
}

void LTTrack::clear()
{
  evId=0;    // event id
  trkId=0;    // track id
  pdg=0;
  parentId=0;    // parent id
  status=0;
  time=0; 
  trkTime=0;
  ke=0;
  edep=0;
  length=0;
  nstep=0;
  trkStepLength=0;
  copy=-1;
  nInToGe=0;
  nOutOfGe=0;
  nSpike=0;
  
  // root Clears
  position.Clear();
  vertPosition.Clear();
  process.Clear();
  physVolName.Clear();
  particleName.Clear();
  preName.Clear();
  postName.Clear();

  // std vector clears
  boundaryStatus.clear();
  boundaryName.clear();
  positionHistory.clear();
  positionEnergy.clear();
  stepLength.clear();
  stepKE.clear();
  postStatus.clear();
  
}
void LTTrack::print(){
  printf(" ********************  LTTrack ********************* \n");
  printf(" \tevent  %i id %i parent %i step %i time %1.3E KE %1.3E \n",evId,trkId,parentId,nstep,time,ke);
  printf(" \tOptical nInToGe %i nOutOfGe %i nSpike %i \n",nInToGe,nOutOfGe,nSpike);
  printf(" *************************************************** \n");
}
