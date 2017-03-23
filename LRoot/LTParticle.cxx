#include "LTParticle.hxx"

ClassImp(LTParticle);

LTParticle::LTParticle() 
{
  clear();
} 

LTParticle::~LTParticle()
{ 
}

void LTParticle::print(int ip) {
  printf("Particle %i) Track Id %i VertexId %i PDG %i, Mass %f, Charge %i MomentumEP(%f,%f,%f,%f) \n",
      ip,TrackId,VertexId,PDG,Mass,Charge,Momentum.T(),Momentum.X(),Momentum.Y(),Momentum.Z());
}

void LTParticle::clear() 
{
  TrackId=0;
  VertexId=0;
  PDG=0;
  Mass=0;
  Charge=0;
  KEnergy=0;
  Momentum.SetPxPyPzE(0,0,0,0);
}
