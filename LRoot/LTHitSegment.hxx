#ifndef LTHitSegment_hxx_seen
#define LTHitSegment_hxx_seen

#include <string>
#include <vector>

#include <TNamed.h>
#include <TLorentzVector.h>


/// same content as TG4HitSegment.hxx
// **********************************
/// This is a base class used to save G4Hit objects into a root output file
/// when the detector simulation expects the electronics to be smulated off
/// line.  It contains the global position of the starting point and stopping
/// point of the track segment that created the hit, plus the fields provided
/// by TG4VHit.  The detector simulation expects the electronics response
/// simulation to assign the sensor identifiers to each hit.
///
/// In the simulation, this class is multiply inherited with a G4VHit derived
/// class that can be used in a G4HitCollection.  The behavior is designed so
/// that TG4HitSegment contains all of the payload for the hit, and detector
/// simulation class contains all of the interface related to the G4
/// infrastructure.

class LTHitSegment: public TNamed {
  public :
    LTHitSegment();
    virtual ~LTHitSegment();
    void clear();
    void print(int iseg=0);


   /// The TrackID for each trajectory that contributed to this hit.  This
    /// could contain the TrackID of the primary particle, but not
    /// necessarily.  
    std::vector<Int_t> Contributors;

    /// The track id of the primary particle.
    Int_t PrimaryId;

    /// The total energy deposit in this hit.  
    Float_t EnergyDeposit;

    /// The "secondary" energy deposit in this hit.  This is used to help
    /// simulate the recombination of electrons, and is part of the total
    /// energy deposit.
    ///
    /// DETSIM is expected to use this field to save the amount of energy
    /// deposited as opticalphotons.  The remaining energy will be deposited
    /// as ionization.  In this model (in argon), the mean number of quanta
    /// created will be <N_q> = (fEnergyDeposit)/(19.5*eV), N_q should be
    /// fluctuated around <N_q>, N_ph = N_q*fSecondaryDeposit/fEnergyDeposit,
    /// and N_e = N_q - N_ph.  Thd fSecondaryDeposit value already includes
    /// the binomial fluctuation, so don't fluctuate N_ph or N_e.
    Float_t SecondaryDeposit;

    Float_t dedx;
    
    /// The total charged track length in this hit.  This includes the
    /// contribution from all of the secondary particles (e.g. delta-rays)
    /// that are included in this hit.
    Float_t TrackLength;

    /// The position of the hit starting point.
    TVector3 Start;

    /// The time of the hit starting point.
    Float_t StartT;

    /// The position of the hit stopping point.
    TVector3 Stop;

    /// The time of the hit stopping point.
    Float_t StopT;

    // for Truth, matches TCapTrajectory TrajId; 
    Int_t   trajIndex; 

    Int_t pdg; 

    enum pdgcode {UNKNOWN=0,NEUTRON=2112,EPLUS=-11,EMINUS=11,GAMMA=22,PIZERO=111,PIPLUS=211,PIMINUS=-211}; 
      
  ClassDef(LTHitSegment,1)
};
#endif
