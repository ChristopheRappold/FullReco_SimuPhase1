#ifndef TDECAYVERTEX
#define TDECAYVERTEX

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "THyphiAttributes.h"
#include "TRandom3.h"

typedef TDataProcess<FullRecoEvent, MCAnaEventG4Sol> TDataProcessInterface;


class TDecayVertex final : public TDataProcessInterface
{
public:
  const THyphiAttributes& att;

  TDecayVertex(const THyphiAttributes& attr);
  ~TDecayVertex();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderDecayVertex(FullRecoEvent& RecoEvent);

  void FragmentTracksFinder(std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                            int& pdgParticle,
                            std::vector<DecayTrackInfo>& FragmentTracks);

  void PionTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results,
                        std::vector<DecayTrackInfo>& PionTracks);

  void CloseDist(DecayTrackInfo& FragmentTrack, DecayTrackInfo& PionTrack, double& distance, double& z);

  double f_function(DecayTrackInfo& DecayTrack, TVector3& PosXYZ);

  double V_function(std::vector<double>& f_vector);

  void SpaceDiscretization(double& Xi, double& Xf, size_t& NstepsX, double& Yi, double& Yf, size_t& NstepsY,
                                          double& Zi, double& Zf, size_t& NstepsZ, size_t& border, std::vector<TVector3>& PosXYZ);

  void TrackstoDecayVertex(std::vector<DecayTrackInfo>& FragmentTracks, std::vector<DecayTrackInfo>& PionTracks,
                                          TVector3& PrimVtxRecons, TVector3& DecayVertexRecons);





  double Zo_target = 24.5; // in cm
  double Zf_target = 25.5; // in cm
  double Xo_target = -1.;
  double Xf_target = 1.;
  double Yo_target = -1.;
  double Yf_target = +1.;

  double Zo_minifibers = 45;

  size_t NstepsdiscretXY      = 7;
  size_t NstepsdiscretZ       = 21;
  size_t Nstepsdiscretbox     = 5;
  size_t nTimesDiscretization = 7;
  double boxDistXY        = 3.;

  double k_factor       = 3.; // multiplies the fragment track f function value

  double min_f_value = 0.5;

  double EnergyThreshold       = 0.001; // in MeV
  double MaxEnergyDiffStrips   = 0.2;   // in MeV
  double MaxEnergyDiffSilicons = 0.1;   // in MeV

  double MaxDistTracks = 0.2;


  TRandom3* rand;
  PDG_fromName pid_fromName;

  struct LocalHists
  {
    TH1F* h_Pt_fragments;
    TH1F* h_Pz_fragments;

    TH1F* h_Pt_pions;
    TH1F* h_Pz_pions;
    TH1F* h_Chi2ndf_pions;

    TH1F* h_Closedist_Distance;
    TH1F* h_Closedist_PosZ;

    TH1F* h_DecayVertexDistance;
    TH1F* h_DecayVertexDistanceX;
    TH1F* h_DecayVertexDistanceY;
    TH1F* h_DecayVertexDistanceZ;

    TH1F* h_DecayVtxstats;
  };
  
  LocalHists LocalHisto;
};



#endif