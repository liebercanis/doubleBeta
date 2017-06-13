#include "LTEvent.hxx"
ClassImp(LTEvent)

LTEvent::LTEvent(): TNamed("LTEvent","LTEvent")
{
  clear();
}

LTEvent::~LTEvent()
{
}

void LTEvent::clear()
{

  evId=0;    // event id
  nPVert=0; // number primary verticies
  pvertex.clear();
  nTraj=0;
  traject.clear(); // vector<LTTraject> trajec

  nPmtHits=0;
  nGeHits=0;
  nOptPhotons=0;
  nArScint=0;
  nWlsScint=0;
  nCherenkov=0;
  nAbsorbed=0;
  nAbsorbedBoundary=0;
  

  nTrajPmtHits=0;
  nTrajOptPhotons=0;
  nTrajArScint=0;
  nTrajWlsScint=0;
  nTrajGeHits=0;

  PDG=0;
  ePrimary=0;
  eCharged=0;
  eNeutral=0;
  eOptical=0;
  ePmt=0;
  eMaxDeposit=0;

  hasConversion;//true (initial) converstion position
  positionEWeight.Clear();
  positionConversion.Clear();
  positionMax.Clear();


}
void LTEvent::print(){
  printf(" ********************  LTevent ********************* \n");
  printf(" \tPrimary pdg %i E Primary %1.3E E charged %1.3E E neutral %1.3E E optical %1.3E E PMT %1.3E \n\tPrimaryVerticies %i Trajectories %i  \n",
      PDG,ePrimary,eCharged,eNeutral,eOptical,ePmt,nPVert,nTraj);
  if(nPVert>0) printf("\tPrimary vertex at position,t (%f,%f,%f,%f) num part %i \n",
      pvertex[0].Position.X(),pvertex[0].Position.Y(),pvertex[0].Position.Z(),pvertex[0].Position.T(),pvertex[0].nParticles);
  else  printf("\tNO primary vertex \n");
  printf(" \tOptical %i Ar scint %i WLS scint %i PMT hits %i Ge hits %i \n",nOptPhotons,nArScint,nWlsScint,nPmtHits,nGeHits);
  printf(" *************************************************** \n");
}
