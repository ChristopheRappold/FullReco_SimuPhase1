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

  void RealTracksFinder(std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                            int& pdgParticle, int cutConditions,
                            std::vector<DecayTrackInfo>& RealTracks);

  void PionTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results,
                        std::vector<DecayTrackInfo>& PionTracks);

  void CloseDist(DecayTrackInfo& FragmentTrack, DecayTrackInfo& PionTrack, double& distance, TVector3& centroid);

  double f_function(DecayTrackInfo& DecayTrack, TVector3& PosXYZ);

  double V_function(std::vector<double>& f_vector);

  void SpaceDiscretization(double& Xi, double& Xf, size_t& NstepsX, double& Yi, double& Yf, size_t& NstepsY,
                                          double& Zi, double& Zf, size_t& NstepsZ, size_t& border, std::vector<TVector3>& PosXYZ);

  void TrackstoDecayVertex(std::vector<DecayTrackInfo>& FragmentTracks, std::vector<DecayTrackInfo>& PionTracks,
                                          TVector3& PrimVtxRecons, TVector3& DecayVertexRecons);

  void Dist_DecayTrackPrimVtx(DecayTrackInfo& PionTrack, TVector3& PrimVtxRecons, double& distance);





  double Zo_target = 24.5; // in cm
  double Zo_minifibers = 45;

  double k_factor       = 3.; // multiplies the fragment track f function value

  size_t NstepsdiscretXY      = 7;
  size_t NstepsdiscretZ       = 21;
  size_t Nstepsdiscretbox     = 5;
  size_t nTimesDiscretization = 7;
  double boxDistXY        = 3.;

  int He3_pdg = 10003;
  int pi_pdg = -211;


  TRandom3* rand;
  PDG_fromName pid_fromName;

  struct LocalHists
  {
    TH1F* h_Pt_fragments;
    TH1F* h_Pz_fragments;

    TH1F* h_Pt_pions;
    TH1F* h_Pz_pions;
    TH1F* h_Chi2ndf_pions;

    TH1F* h_Pt_realpions;
    TH1F* h_Pz_realpions;

    TH1F* h_Pt_cutpions;
    TH1F* h_Pz_cutpions;

    TH1F* h_Nrealpions;
    TH1F* h_Ncutpions;
    TH1F* h_Npions;

    TH1F* h_Closedist_Distance;
    TH1F* h_Closedist_PosZ;
    TH1F* h_Dist_DecayTrackPrimVtx;

    TH1F* h_Closedist_realDistance;
    TH1F* h_Closedist_realPosZ;
    TH1F* h_Dist_realDecayTrackPrimVtx;

    TH1F* h_Closedist_cutDistance;
    TH1F* h_Closedist_cutPosZ;
    TH1F* h_Dist_cutDecayTrackPrimVtx;

    TH1F* h_DecayVertexDistance;
    TH1F* h_DecayVertexDistanceX;
    TH1F* h_DecayVertexDistanceY;
    TH1F* h_DecayVertexDistanceZ;

    TH1F* h_DecayVertexDistance_centroid;
    TH1F* h_DecayVertexDistanceX_centroid;
    TH1F* h_DecayVertexDistanceY_centroid;
    TH1F* h_DecayVertexDistanceZ_centroid;

    TH1F* h_DecayVertexDistance_2centroid_average;
    TH1F* h_DecayVertexDistanceX_2centroid_average;
    TH1F* h_DecayVertexDistanceY_2centroid_average;
    TH1F* h_DecayVertexDistanceZ_2centroid_average;

    TH1F* h_DecayVertexDistance_2centroid_closest;
    TH1F* h_DecayVertexDistanceX_2centroid_closest;
    TH1F* h_DecayVertexDistanceY_2centroid_closest;
    TH1F* h_DecayVertexDistanceZ_2centroid_closest;

    TH1F* h_DecayVertexrealDistance;
    TH1F* h_DecayVertexrealDistanceX;
    TH1F* h_DecayVertexrealDistanceY;
    TH1F* h_DecayVertexrealDistanceZ;

    TH1F* h_DecayVertexcutDistance;
    TH1F* h_DecayVertexcutDistanceX;
    TH1F* h_DecayVertexcutDistanceY;
    TH1F* h_DecayVertexcutDistanceZ;

    TH1F* h_DecayVertexPosZ_real;
    TH1F* h_DecayVertexPosZ_vfunction;
    TH1F* h_DecayVertexPosZ_centroid;

    TH1F* h_HypInvariantMass;

    TH1F* h_DecayVtxstats;
  };

  LocalHists LocalHisto;
};



#endif