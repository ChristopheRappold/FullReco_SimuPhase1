#ifndef TDECAYVERTEX
#define TDECAYVERTEX

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "THyphiAttributes.h"
#include "TRandom3.h"
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFParticleSIMD.h"
#include "KFParticleBaseSIMD.h"

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
                         int& pdgParticle, int& cutConditions,
                         std::vector<KFParticle>& RealTracks);

  void FragmentSelector(std::vector<KFParticle>& FragmentTracks_All, TVector3& PrimVtxRecons, std::vector<KFParticle>& FragmentTracks);

  void PionTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results,
                        std::vector<KFParticle>& PionTracks);

  void PionSelector(std::vector<KFParticle>& PionTracks_All, TVector3& PrimVtxRecons, std::vector<KFParticle>& PionTracks);

  void CloseDist(KFParticle& FragmentTrack, KFParticle& PionTrack, double& distance, TVector3& centroid);

  double f_function(KFParticle& DecayTrack, TVector3& PosXYZ);

  double V_function(std::vector<double>& f_vector);

  void SpaceDiscretization(double& Xi, double& Xf, size_t& NstepsX, double& Yi, double& Yf, size_t& NstepsY,
                                          double& Zi, double& Zf, size_t& NstepsZ, size_t& border, std::vector<TVector3>& PosXYZ);

  void TrackstoDecayVertex(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                                          TVector3& PrimVtxRecons, TVector3& DecayVertexRecons);

  void ThetaDist_TrackPrimVtx(KFParticle& Track, TVector3& PrimVtxRecons, double& theta, double& distance);

  void MotherTracksRecons(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                          TVector3& PrimVtxRecons, std::array<double,6> Cov_PrimVtx, std::vector<KFParticle>& MotherTracks,
                          std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks);

  void MotherTrack_SiHit(TVector3& PrimVtxRecons, TVector3& DecayVtxRecons, double& Z_plane,
                          std::vector<std::vector<double> >& Hits_Si, TVector3& Mother_Sihit);

  void MotherTrackSiliconHits(TVector3& PrimVtxRecons, TVector3& DecayVtxRecons, std::vector<std::vector<double> >& Hits_Si1,
                               std::vector<std::vector<double> >& Hits_Si2, KFParticle& Mother, KFParticle& Si_MotherTrack);

  void AllTrackstoDecayVertex_Vfunction(std::vector<KFParticle>& AllTracks, TVector3& Old_DecayVertexRecons, TVector3& DecayVertexRecons);

  void AllTrackstoDecayVertex_Centroids(std::vector<KFParticle>& AllTracks, TVector3& DecayVertexRecons);



  double Zo_target = 24.; // in cm
  double Zo_minifibers = 45;

  double Z_plane_Si1    = 27.;  // in cm
  double widthStrip_Si1 = 0.03; // in cm

  double Z_plane_Si2    = 30.;  // in cm
  double widthStrip_Si2 = 0.03; // in cm

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


  PDG_fromName pid_fromName;

  int H3L_pdg = pid_fromName("H3L");
  int H3L_charge = TDatabasePDG::Instance()->GetParticle(H3L_pdg)->Charge()/3.;
  double H3L_mass = TDatabasePDG::Instance()->GetParticle(H3L_pdg)->Mass();

  int He3_pdg = pid_fromName("He3");
  int He3_charge = TDatabasePDG::Instance()->GetParticle(He3_pdg)->Charge()/3.;
  double He3_mass = TDatabasePDG::Instance()->GetParticle(He3_pdg)->Mass();

  int pi_pdg = pid_fromName("pi-");
  int pi_charge = TDatabasePDG::Instance()->GetParticle(pi_pdg)->Charge()/3.;
  double pi_mass = TDatabasePDG::Instance()->GetParticle(pi_pdg)->Mass();

  int No_cutconditions = 0;
  int Yes_cutconditions = 1;

  double MinDist_FragmentTracksPrimVtx = 0.;
  double MinDist_PionTracksPrimVtx = 0.05;

  double Max_DaughtersTracks_closedist = 0.5;
  double Max_Dist_MotherTrackPrimVtx = 20.;
  double Max_Theta_MotherTrackPrimVtx = 30.;

  int KFPart_fConstructMethod = 2;
  double KFPart_fMassHypo = H3L_mass;

  double Max_dist_Mother_SiHit = 0.5;
  double Min_EnergyDeposition_Si = 0.08;

  TRandom3* rand;

  std::vector<KFParticle> KFPart_Pions;
  std::vector<KFParticle> KFPart_Fragments;
  std::vector<KFParticle> KFPart_Daughters;
  std::vector<KFParticle> KFPart_Mother;

  struct LocalHists
  {
    TH1F* h_Pt_fragments;
    TH1F* h_Pz_fragments;
    TH1F* h_Dist_FragmentTrackPrimVtx;

    TH1F* h_Pt_pions;
    TH1F* h_Pz_pions;
    TH1F* h_Chi2ndf_pions;

    TH1F* h_Pt_cutpions;
    TH1F* h_Pz_cutpions;

    TH1F* h_Nrealpions;
    TH1F* h_Ncutpions;
    TH1F* h_Npions;


    TH1F* h_Closedist_Distance;
    TH1F* h_Closedist_PosZ;
    TH2F* h_Dist_DecayTrackPrimVtx;

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

    TH1F* h_DecayVertexcutDistance;
    TH1F* h_DecayVertexcutDistanceX;
    TH1F* h_DecayVertexcutDistanceY;
    TH1F* h_DecayVertexcutDistanceZ;


    TH1F* h_DecayVertexPosZ_real;
    TH1F* h_DecayVertexPosZ_vfunction;
    TH1F* h_DecayVertexPosZ_centroid;
    TH1F* h_DecayVertexPosZ_AllVfunc;
    TH1F* h_DecayVertexPosZ_AllCentroid;
    TH1F* h_DecayVertexPosZ_AllKFPart;


    TH2F* h_N_MotherTracks;
    TH2F* h_Dist_MotherTrackPrimVtx;
    TH2F* h_Theta_MotherTrackPrimVtx;
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

    TH1F* h_DecayVertexDistance_AllKFPart;
    TH1F* h_DecayVertexDistanceX_AllKFPart;
    TH1F* h_DecayVertexDistanceY_AllKFPart;
    TH1F* h_DecayVertexDistanceZ_AllKFPart;

    TH1F* h_DecayVtxstats;
  };

  LocalHists LocalHisto;
};

class OutputContainer
{
  
 public:
  
  OutputContainer() = default;
  virtual ~OutputContainer() = default;  
  
  //  candidate parameters setters for two daugthers
  void SetChi2PrimPos(float value) {chi2_prim_pos_ = value;};
  void SetChi2PrimNeg(float value) {chi2_prim_neg_ = value;};
  void SetDistance(float value) {distance_ = value;};
  void SetCosineDaughterPos(float value) {cosine_daughter_pos_ = value;};
  void SetCosineDaughterNeg(float value) {cosine_daughter_neg_ = value;};
  void SetChi2Geo(float value) {chi2_geo_ = value;};
  void SetL(float value) {l_ = value;};
  void SetLdL(float value) {ldl_ = value;};
  void SetIsFromPV(int value) {is_from_pv_ = value;};
  void SetCosineTopo(float value) {cosine_topo_ = value;};
  void SetSigmaMassRatio(float value) {sigma_mass_ratio_ = value;};
  void SetChi2Topo(float value) {chi2_topo_ = value;};
  void SetNHitsPos(int value) {nhits_pos_ = value;};
  void SetNHitsNeg(int value) {nhits_neg_ = value;};

  //  candidate parameters setters for third daugther
  void SetChi2PrimThird(float value) {chi2_prim_third_ = value;};
  void SetDistanceThird(float value) {distance_third_ = value;};
  void SetCosineDaughterThird(float value) {cosine_daughter_third_ = value;};
  void SetChi2GeoThree(float value) {chi2_geo_three_ = value;};
  void SetCosineTopoThree(float value) {cosine_topo_three_ = value;};
  void SetChi2TopoThree(float value) {chi2_topo_three_ = value;};
  void SetNHitsThird(int value) {nhits_third_ = value;};
  
  void SetParticle(const KFParticle& particle) {particle_ = particle;};
  
  //  candidate parameters getters for two daugthers
  float GetChi2PrimPos() const {return chi2_prim_pos_;};
  float GetChi2PrimNeg() const {return chi2_prim_neg_;};
  float GetDistance() const {return distance_;};
  float GetCosineDaughterPos() const {return cosine_daughter_pos_;};
  float GetCosineDaughterNeg() const {return cosine_daughter_neg_;};
  float GetChi2Geo() const {return chi2_geo_;};
  float GetL() const {return l_;};
  float GetLdL() const {return ldl_;};
  int   GetIsFromPV() const {return is_from_pv_;};
  float GetCosineTopo() const {return cosine_topo_;};
  float GetSigmaMassRatio() const {return sigma_mass_ratio_;};
  float GetChi2Topo() const {return chi2_topo_;};
  int   GetNHitsPos() const {return nhits_pos_;};
  int   GetNHitsNeg() const {return nhits_neg_;};

  //  candidate parameters getters for third daugther
  float GetChi2PrimThird() const {return chi2_prim_third_;};
  float GetDistanceThird() const {return distance_third_;};
  float GetCosineDaughterThird() const {return cosine_daughter_third_;};
  float GetChi2GeoThree() const {return chi2_geo_three_;};
  float GetCosineTopoThree() const {return cosine_topo_three_;};
  float GetChi2TopoThree() const {return chi2_topo_three_;};
  int   GetNHitsThird() const {return nhits_third_;};

  const KFParticle& GetParticle() const {return particle_;};

 protected:
   
  // candidate selection parameters (to be cut) for two daughters
  float chi2_prim_pos_ {-1.};       ///< \f$\chi^2\f$ of the positive track to the primary vertex (PV)
  float chi2_prim_neg_ {-1.};       ///< \f$\chi^2\f$ of the negative track to the PV
  float distance_ {-1.};            ///< Distance between daughter tracks in their closest approach
  float cosine_daughter_pos_ {-1.}; ///< Cosine of the angle between positive daughter's and mother's momenta
  float cosine_daughter_neg_ {-1.}; ///< Cosine of the angle between negative daughter's and mother's momenta
  float chi2_geo_ {-1.};            ///< \f$\chi^2\f$ of daughters' tracks in their closest approach
  float l_ {-1.};                   ///< Lenght of interpolated track from secondary to primary vertex
  float ldl_ {-1.};                 ///< Distance between primary and secondary vertices divided by error 
  int   is_from_pv_ {-1};           ///< Flag variable whether mother particle comes from the PV (1-yes, 0-no)
  float cosine_topo_{-1.};          ///< Cosine of the angle between reconstructed mother's momentum and mother's radius vector beginning in the PV
  float sigma_mass_ratio_ {-1.};    ///< Difference between invariant and real mother's mass divided by the error (not used now)
  float chi2_topo_ {-1.};           ///< \f$\chi^2\f$ of the mother's track to the PV
  int   nhits_pos_{-1};
  int   nhits_neg_{-1};

  // candidate selection parameters (to be cut) for third daughter
  float chi2_prim_third_ {-1.};       ///< \f$\chi^2\f$ of the third track to the PV
  float distance_third_ {-1.};        ///< Distance between third daughter track and SV in their closest approach
  float cosine_daughter_third_ {-1.}; ///< Cosine of the angle between third daughter's and mother's momenta
  float chi2_geo_three_ {-1.};        ///< \f$\chi^2\f$ of all three daughters' tracks in their closest approach
  float cosine_topo_three_{-1.};      ///< Cosine of the angle between reconstructed mother's momentum of three particles and mother's radius vector beginning in the PV
  float chi2_topo_three_ {-1.};       ///< \f$\chi^2\f$ of the mother's track of three particles to the PV
  int   nhits_third_{-1};

  KFParticle particle_;
  
};


#endif