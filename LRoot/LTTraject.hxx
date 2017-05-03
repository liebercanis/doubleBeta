#ifndef LTTraject_hxx_seen
#define LTTraject_hxx_seen

#include <string>
#include <vector>

#include <TObject.h>
#include <TString.h>
#include <TLorentzVector.h>
#include "LTHitSegment.hxx"


/// A class to save a G4 trajectory into a root output file without linking to
/// geant.  A trajectory is the truth information about the path of a particle
/// through the G4 simulation. It saves the parent trajectory that generated
/// this particle, the initial momentum of the particle, and the path followed
/// by the particle in the detector.  
//
enum LTTrajectType {UNK,PRI,SCI,WLS,PMTHIT,GEHIT};
class LTTraject: public TObject {
  public :
    LTTraject();
    virtual ~LTTraject();
    void clear();
    void print(int mode=0);
    
    Int_t           TrajId;   //[Trajectory_]
    Int_t           ParentId;   //[Trajectory_]
    Int_t           PrimaryId;   //[Trajectory_]
    Int_t           PDG;   //[Trajectory_]
    Float_t         Mass;   //[Trajectory_]
    Float_t         Charge;   //[Trajectory_]
    Float_t         KE;
    TString         Name;
    Int_t           Type;
    // each element in std::vector corresponds to a point on the particle path
    std::vector<TLorentzVector> Position;
    std::vector<TVector3> Momentum;
    std::vector<Int_t>   Region;
    std::vector<LTHitSegment> segments;
ClassDef(LTTraject,1)
};
#endif
