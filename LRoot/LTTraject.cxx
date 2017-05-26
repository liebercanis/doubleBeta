#include "LTTraject.hxx"
ClassImp(LTTraject);

LTTraject::LTTraject() 
{
  clear();
} 

LTTraject::~LTTraject()
{ 
}

void LTTraject::clear() 
{
  evId=0;    // event id
  trkId=0;    // track id
  trkParentId=0;    // parent id
  trkStatus=0;
  trkStep=0;
  trkLength=0;
  trkEdep=0;
  time=0;
  trkTime=0;
  PDG=0;
  mass=0;
  charge=0;
  ke=0;
  physVolName.Clear();
  copy=0;
  type= LTTrajectType::UNK;
  name.Clear();
  // each element in vector corresponds to a point on the particle path
  position.clear();
  momentum.clear();
  region.clear();
  segments.clear();
}

void LTTraject::print(int itraj){
  printf(" %i) LTTraject: id %i  parent %i  primary %i  PDG %i  mass %.3f  charge %i positions %i momenta %i type %i \n ",itraj,
      trajId, parentId, primaryId, PDG, mass, charge, (int) position.size(), (int) momentum.size(), type  );
  if(itraj==0) return;
  // detailed info 
  for(unsigned int ip=0; ip< position.size(); ++ip)
    printf("     traj region %i  time/position (%.2f,%.2f,%.2f,%.2f) 3 momentum (%.2f %.2f %.2f) \n",region[ip],
        position[ip].T(),position[ip].X(),position[ip].Y(),position[ip].Z(),momentum[ip].X(),momentum[ip].Y(),momentum[ip].Z() );
  printf("      traj  %i with %i segments:   \n",itraj, (int) segments.size() );
  for(int iseg =0; iseg< int(segments.size()) ; ++ iseg ) segments[iseg].print(iseg);
  
}
