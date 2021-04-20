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

  void KFPart_PrimaryVertex(TVector3& PrimVtxRecons, std::array<double,6> Cov_PrimVtx, KFParticleSIMD& temp_PrimVtx);

  void MotherTracksRecons(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                          const KFParticleSIMD* pointer_PrimVtx, std::vector<KFParticle>& MotherTracks,
                          std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks);

  void MotherSelector(std::vector<KFParticle>& MotherTracks_All, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks_All,
                      std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks, TVector3& PrimVtxRecons,
                      std::vector<KFParticle>& MotherTracks, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks);

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
  //double widthStrip_Si2 = 0.03; // in cm

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

  double c_light_speed = 299792458.; //in m/s
  double c_light_speed_cmps = c_light_speed * 1.e-10; //in cm/ps


  //Cut conditions on real pions
  int No_cutconditions = 0;
  int Yes_cutconditions = 1;


  //Cut conditions on reconstructed fragments
  int ifCut_MaxChi2ndf_FragmentTracks = 1;
  double MaxChi2ndf_FragmentTracks = 3.; //Change !

  int ifCut_MinDist_FragmentTracksPrimVtx = 0;
  double MinDist_FragmentTracksPrimVtx = 0.; //Change !


  //Cut conditions on reconstructed pions
  int ifCut_MaxChi2ndf_PionTracks = 1;
  double MaxChi2ndf_PionTracks = 3.; //Change !

  int ifCut_MinDist_PionTracksPrimVtx = 0;
  double MinDist_PionTracksPrimVtx = 0.05; //Change !


  //Cut conditions on reconstructed hypernuclei
  int ifSet_ProductionVertex = 1;
  int ifSet_MassConstraint = 0;

  int KFPart_fConstructMethod = 2;
  double KFPart_fMassHypo = H3L_mass;

  int ifCut_MaxClosedist_DaughterTracks = 0; 
  double MaxClosedist_DaughterTracks = 10.; //Change !

  int ifCut_MaxAngle_MotherFragment = 0;
  double MaxAngle_MotherFragment = 90.; //In degrees (º) Change !

  int ifCut_MaxAngle_MotherPion = 0;
  double MaxAngle_MotherPion = 90.; //In degrees (º) Change !

  int ifCut_MaxChi2ndf = 0;
  double MaxChi2ndf = 7.; //Change !

  int ifCut_MaxDist_MotherTrackPrimVtx = 0;
  double MaxDist_MotherTrackPrimVtx = 10.; //Change !

  int ifCut_MaxAngle_MotherTrackPrimVtx = 0;
  double MaxAngle_MotherTrackPrimVtx = 60.; //Change !

  int ifCut_MaxPosZ_DecayVertex = 0;
  double MaxPosZ_DecayVertex = 80.;

  int ifCut_MinPosZ_DecayVertex = 0;
  double MinPosZ_DecayVertex = 10.;

  int ifCut_ArmenterosPodolanski = 0;

  double Max_dist_Mother_SiHit = 0.5;
  double Min_EnergyDeposition_Si = 0.08;


  //TRandom3* rand;

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

    TH1F* h_DecayVertexDistance_KFPart;
    TH1F* h_DecayVertexDistanceX_KFPart;
    TH1F* h_DecayVertexDistanceY_KFPart;
    TH1F* h_DecayVertexDistanceZ_KFPart;

    TH1F* h_DecayVertexDistance_KFPart_PrimVtx;
    TH1F* h_DecayVertexDistanceX_KFPart_PrimVtx;
    TH1F* h_DecayVertexDistanceY_KFPart_PrimVtx;
    TH1F* h_DecayVertexDistanceZ_KFPart_PrimVtx;

    TH1F* h_DecayVertexcutDistance;
    TH1F* h_DecayVertexcutDistanceX;
    TH1F* h_DecayVertexcutDistanceY;
    TH1F* h_DecayVertexcutDistanceZ;

    TH1F* h_DecayVertexcutDistance_KFPart;
    TH1F* h_DecayVertexcutDistanceX_KFPart;
    TH1F* h_DecayVertexcutDistanceY_KFPart;
    TH1F* h_DecayVertexcutDistanceZ_KFPart;

    TH1F* h_DecayVertexcutDistance_KFPart_PrimVtx;
    TH1F* h_DecayVertexcutDistanceX_KFPart_PrimVtx;
    TH1F* h_DecayVertexcutDistanceY_KFPart_PrimVtx;
    TH1F* h_DecayVertexcutDistanceZ_KFPart_PrimVtx;


    TH1F* h_DecayVertexPosZ_real;
    TH1F* h_DecayVertexPosZ_vfunction;
    TH1F* h_DecayVertexPosZ_centroid;
    TH1F* h_DecayVertexPosZ_KFPart;
    TH1F* h_DecayVertexPosZ_KFPart_PrimVtx;
    TH1F* h_DecayVertexPosZ_AllVfunc;
    TH1F* h_DecayVertexPosZ_AllCentroid;
    TH1F* h_DecayVertexPosZ_AllKFPart;


    TH2F* h_N_MotherTracks;
    TH2F* h_Dist_MotherTrackPrimVtx;
    TH2F* h_Theta_MotherTrackPrimVtx;
    
    TH1F* h_HypInvariantMass;
    TH1F* h_HypInvariantMassCheck;
    TH1F* h_HypErrorInvariantMass;
    TH1F* h_HypLifeTime_PrimVtx;
    TH1F* h_HypErrorLifeTime_PrimVtx;
    TH1F* h_HypcutLifeTime_PrimVtx;
    TH2F* h_Hyp_ArmenterosPodolanski;
    TH2F* h_Hyp_CutArmenterosPodolanski;

    TH1F* h_HypInvariantMass_LorentzVect;
    TH1F* h_HypInvariantMass_CutLorentzVect;

    
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

#endif