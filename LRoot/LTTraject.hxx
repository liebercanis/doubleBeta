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
    
    Int_t           evId;
    Int_t           trkId;
    Int_t           trkParentId;
    Int_t           trkStatus;
    Int_t           trkStep;
    Float_t         trkLength;
    Float_t         trkEdep;
    Float_t         time;
    Float_t         trkTime;
    Int_t           trajId;   //[Trajectory_]
    Int_t           parentId;   //[Trajectory_]
    Int_t           primaryId;   //[Trajectory_]
    Int_t           PDG;   //[Trajectory_]
    Float_t         mass;   //[Trajectory_]
    Float_t         charge;   //[Trajectory_]
    Float_t         ke;
    TString         name;
    TString         physVolName;
    Int_t           copy;
    Int_t           type;

    // each element in std::vector corresponds to a point on the particle path
    std::vector<TLorentzVector> position;
    std::vector<TVector3> momentum;
    std::vector<Int_t>   region;
    std::vector<LTHitSegment> segments;
ClassDef(LTTraject,2)
};
#endif
