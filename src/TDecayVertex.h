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

  void FragmentSelector(std::vector<DecayTrackInfo>& FragmentTracks_All, TVector3& PrimVtxRecons, std::vector<DecayTrackInfo>& FragmentTracks);

  void PionTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results,
                        std::vector<DecayTrackInfo>& PionTracks);

  void PionSelector(std::vector<DecayTrackInfo>& PionTracks_All, TVector3& PrimVtxRecons, std::vector<DecayTrackInfo>& PionTracks);

  void CloseDist(DecayTrackInfo& FragmentTrack, DecayTrackInfo& PionTrack, double& distance, TVector3& centroid);

  double f_function(DecayTrackInfo& DecayTrack, TVector3& PosXYZ);

  double V_function(std::vector<double>& f_vector);

  void SpaceDiscretization(double& Xi, double& Xf, size_t& NstepsX, double& Yi, double& Yf, size_t& NstepsY,
                                          double& Zi, double& Zf, size_t& NstepsZ, size_t& border, std::vector<TVector3>& PosXYZ);

  void TrackstoDecayVertex(std::vector<DecayTrackInfo>& FragmentTracks, std::vector<DecayTrackInfo>& PionTracks,
                                          TVector3& PrimVtxRecons, TVector3& DecayVertexRecons);

  void ThetaDist_TrackPrimVtx(DecayTrackInfo& Track, TVector3& PrimVtxRecons, double& theta, double& distance);

  void MotherTracksRecons(std::vector<DecayTrackInfo>& FragmentTracks, std::vector<DecayTrackInfo>& PionTracks,
                          TVector3& PrimVtxRecons, TVector3& DecayVtxRecons, std::vector<DecayTrackInfo>& MotherTracks,
                          std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks);

  void MotherTrack_SiHit(TVector3& PrimVtxRecons, TVector3& DecayVtxRecons, double& Z_plane,
                          std::vector<std::vector<double> >& Hits_Si, TVector3& Mother_Sihit);

  void MotherTrackSiliconHits(TVector3& PrimVtxRecons, TVector3& DecayVtxRecons, std::vector<std::vector<double> >& Hits_Si1,
                               std::vector<std::vector<double> >& Hits_Si2, DecayTrackInfo& Si_MotherTrack);

  void AllTrackstoDecayVertex_Vfunction(std::vector<DecayTrackInfo>& AllTracks, TVector3& Old_DecayVertexRecons, TVector3& DecayVertexRecons);

  void AllTrackstoDecayVertex_Centroids(std::vector<DecayTrackInfo>& AllTracks, TVector3& DecayVertexRecons);



  double Zo_target = 24.; // in cm
  double Zo_minifibers = 45;

  double Z_plane_Si1    = 27.;  // in cm
  double Z_plane_Si2    = 30.;  // in cm
  double Dist_to_Silicons = 0.5;

  double k_factor       = 3.; // multiplies the fragment track f function value

  size_t NstepsdiscretXY      = 7;
  size_t NstepsdiscretZ       = 21;
  size_t Nstepsdiscretbox     = 5;
  size_t nTimesDiscretization = 7;
  double boxDistXY            = 3.;

  size_t NstepsdiscretboxXYZ  = 9;
  double boxXYZ               = 2.;
  size_t nTimesBoxXYZ         = 5;


  int He3_pdg = 10003;
  int pi_pdg = -211;

  double pi_mass = 0.13957018; //in GeV

  double MinDist_FragmentTracksPrimVtx = 0.;
  double MinDist_PionTracksPrimVtx = 0.05;

  double Max_DaughtersTracks_closedist = 0.5;
  double Max_Dist_MotherTrackPrimVtx = 20.;
  double Max_Theta_MotherTrackPrimVtx = 30.;

  double Max_dist_Mother_SiHit = 0.5;
  double Min_EnergyDeposition_Si = 0.08;

  TRandom3* rand;
  PDG_fromName pid_fromName;

  struct LocalHists
  {
    TH1F* h_Pt_fragments;
    TH1F* h_Pz_fragments;
    TH1F* h_Dist_FragmentTrackPrimVtx;

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
    TH2F* h_Dist_DecayTrackPrimVtx;

    TH1F* h_Closedist_realDistance;
    TH1F* h_Closedist_realPosZ;
    TH2F* h_Dist_realDecayTrackPrimVtx;

    TH1F* h_Closedist_cutDistance;
    TH1F* h_Closedist_cutPosZ;
    TH2F* h_Dist_cutDecayTrackPrimVtx;

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

    TH1F* h_DecayVertexDistance_2centroid_IPCheck;
    TH1F* h_DecayVertexDistanceX_2centroid_IPCheck;
    TH1F* h_DecayVertexDistanceY_2centroid_IPCheck;
    TH1F* h_DecayVertexDistanceZ_2centroid_IPCheck;

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
    TH1F* h_DecayVertexPosZ_AllVfunc;
    TH1F* h_DecayVertexPosZ_AllCentroid;


    TH1F* h_Dist_MotherTrackPrimVtx;
    TH1F* h_Theta_MotherTrackPrimVtx;
    TH1F* h_HypInvariantMass;

    TH1F* h_N_Si_MotherTracks;

    TH1F* h_DecayVertexDistance_AllVfunc;
    TH1F* h_DecayVertexDistanceX_AllVfunc;
    TH1F* h_DecayVertexDistanceY_AllVfunc;
    TH1F* h_DecayVertexDistanceZ_AllVfunc;

    TH1F* h_DecayVertexDistance_AllCentroid;
    TH1F* h_DecayVertexDistanceX_AllCentroid;
    TH1F* h_DecayVertexDistanceY_AllCentroid;
    TH1F* h_DecayVertexDistanceZ_AllCentroid;

    TH1F* h_DecayVtxstats;
  };

  LocalHists LocalHisto;
};



#endif