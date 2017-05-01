#ifndef LTPVertex_hxx_seen
#define LTPVertex_hxx_seen

#include <string>
#include <vector>

#include <TObject.h>
#include <TLorentzVector.h>
#include <TString.h>
#include "LTParticle.hxx"

class LTPVertex: public TObject {
  public :
    LTPVertex();
    virtual ~LTPVertex();
    void clear();
    void print(int i);

    Int_t VertexId;
    TString physVolumeName;
    TLorentzVector Position;
    TLorentzVector RelPosition;
    Int_t nParticles;
    std::vector<LTParticle> particle;
    
ClassDef(LTPVertex,1)
};
#endif
