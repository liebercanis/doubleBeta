#ifndef LTParticle_hxx_seen
#define LTParticle_hxx_seen

#include <string>
#include <vector>

#include <TObject.h>
#include <TLorentzVector.h>


/// A class to save a G4 trajectory into a root output file without linking to
/// geant.  A trajectory is the truth information about the path of a particle
/// through the G4 simulation. It saves the parent trajectory that generated
/// this particle, the initial momentum of the particle, and the path followed
/// by the particle in the detector.  
class LTParticle: public TObject {
  public :
    LTParticle();
    virtual ~LTParticle();
    void clear();
    void print(int ip);

    /// The vertex id that this particle comes from.
    Int_t TrackId;
    Int_t VertexId;
    Int_t PDG;
    Float_t Mass;
    Int_t Charge;
    Float_t KEnergy;
    TLorentzVector Momentum;
ClassDef(LTParticle,1)    
 };
#endif
