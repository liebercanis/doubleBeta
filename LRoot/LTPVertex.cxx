#include "LTPVertex.hxx"

ClassImp(LTPVertex);

LTPVertex::LTPVertex() 
{
  clear();
} 

LTPVertex::~LTPVertex()
{ 
}

void LTPVertex::print(int iv){
  printf(" %i) vertex volume %s id %i particles %i t/position(%.2f,%.2f,%.2f,%.2f) \n",iv,physVolumeName.Data(),nParticles,
      VertexId,Position.T(),Position.X(),Position.Y(),Position.Z());

  printf("         with %i particles: ", (int) particle.size() );
  for(int ipp =0; ipp< int(particle.size()) ; ++ ipp ) particle[ipp].print(ipp);
  printf("\n");
}

void LTPVertex::clear() 
{
   VertexId=0;
   Position.SetXYZT(0,0,0,0);
   RelPosition.SetXYZT(0,0,0,0);
   nParticles=0;
   particle.clear();
   physVolumeName.Clear(); 
}
