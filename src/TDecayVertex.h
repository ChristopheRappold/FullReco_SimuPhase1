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

  void StudyCaseSelector(std::string StudyCase, int& Hyp_pdg, int& Fragment_pdg);

  void RealTracksFinder(std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                         int& pdgParticle, int& cutConditions,
                         std::vector<KFParticle>& RealTracks);

  void FragmentMDCTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results, int& fragment_pdg,
                                std::vector<KFParticle>& FragmentMDCTracks);

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

  void SiHitsFinder(KFParticle& Track, std::vector<std::vector<double> >& Hits_Si, std::vector<std::vector<double> >& Track_Sihit);

  void MotherDaughtersTrack_SiHits(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                  std::vector<KFParticle>& MotherTracks, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks,
                  std::vector<std::vector<double> >& Hits_Si1, std::vector<std::vector<double> >& Hits_Si2,
                  std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>>& Fragment_SiHits,
                  std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>>& Pion_SiHits,
                  std::vector<std::tuple<int, std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>>>& Mother_SiHits);

  void SiHitsFinder2(KFParticle& Track, int idSilicon, int stripDirection, //Strip direction: 1 -> X, 2 -> Y 
                    std::vector<std::tuple<double, size_t> >& Hits_Si, std::vector<std::vector<double> >& Track_Sihit);

  void MotherDaughtersTrack_SiHits2(std::vector<KFParticle>& FragmentTracks, std::vector<KFParticle>& PionTracks,
                    std::vector<KFParticle>& MotherTracks, std::vector<std::tuple<size_t, size_t>>& RefDaughtersTracks,
                    std::vector<std::tuple<double, size_t> >& HitsX_Si1, std::vector<std::tuple<double, size_t> >& HitsY_Si1,
                    std::vector<std::tuple<double, size_t> >& HitsX_Si2, std::vector<std::tuple<double, size_t> >& HitsY_Si2,
                    std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>>& Fragment_SiHits,
                    std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>>& Pion_SiHits,
                    std::vector<std::tuple<int, std::vector<std::vector<std::vector<double>>>>>& Mother_SiHits);

  void AllTrackstoDecayVertex_Vfunction(std::vector<KFParticle>& AllTracks, TVector3& Old_DecayVertexRecons, TVector3& DecayVertexRecons);

  void AllTrackstoDecayVertex_Centroids(std::vector<KFParticle>& AllTracks, TVector3& DecayVertexRecons);



  double Zo_target = 24.5; // in cm
  double Zf_target = 25.5; // in cm
  double Zo_minifibers = 45; // in cm

  double Z_plane_Si1    = 27.;  // in cm
  //double widthStrip_Si1 = 0.019; // in cm

  double Z_plane_Si2    = 30.;  // in cm
  //double widthStrip_Si2 = 0.019; // in cm

  double Dist_to_Target = 0.5; // in cm
  double Dist_to_Silicons = 0.5; // in cm
  double Dist_to_Minifibers = 0.5; // in cm

  double k_factor       = 3.; // multiplies the fragment track f function value

  size_t NstepsdiscretXY      = 7;
  size_t NstepsdiscretZ       = 21;
  size_t Nstepsdiscretbox     = 5;
  size_t nTimesDiscretization = 7;
  double boxDistXY            = 3.;

  size_t NstepsdiscretboxXYZ  = 9;
  double boxXYZ               = 2.;
  size_t nTimesBoxXYZ         = 5;

  float fieldMDCParameters[10] = {0.f, 0.f, 0.f, // Bx(dz=0), Bx'(dz=0), Bx''(dz)
                                  0.f, 0.f, 0.f, // By(dz=0), By'(dz=0), By''(dz)
                                  0.f, 0.f, 0.f, // Bz(dz=0), Bz'(dz=0), Bz''(dz)
                                  0.f};          // z_0

  //KFParticleFieldRegion fieldMDC(fieldMDCParameters); // Correct field
  //fieldMDC.SetOneEntry(fieldMDCParameters);


  PDG_fromName pid_fromName;

  int Fragment_pdg;

  int He3_pdg = pid_fromName("He3");
  int He4_pdg = pid_fromName("alpha");
  int deuteron_pdg = pid_fromName("deuteron");
  int lambda_pdg = pid_fromName("lambda");
  int proton_pdg = pid_fromName("proton");
  int pi_pdg = pid_fromName("pi-");

  int H3L_pdg = 20001;
  int H4L_pdg = 20002;

  int Hyp_pdg;
  int Hyp_charge;
  double Hyp_mass;

  double c_light_speed = 299792458.; //in m/s
  double c_light_speed_cmps = c_light_speed * 1.e-10; //in cm/ps


  //Cut conditions on real pions
  int No_cutconditions = 0;
  int Yes_cutconditions = 1;


  //Methods for reconstructed fragments
  int recons_from_FRS_MDC = -1; // 1-> FRS, 2-> MDC

  int ref_RealFragment = -1;
  bool ifOnlyRealFragment = false;

  //Cut conditions on reconstructed fragments
  int ifCut_MaxChi2ndf_FragmentTracks = 1;
  double MaxChi2ndf_FragmentTracks = 3.; //Change !

  int ifCut_MinDist_FragmentTracksPrimVtx = 0;
  double MinDist_FragmentTracksPrimVtx = 0.; //Change !

  int ifCut_MinMomZ_FragmentTracks = 0;
  double MinMomZ_FragmentTracks = 0.; //Change !

  int ifCut_MaxTheta_FragmentTracks = 0;
  double MaxTheta_FragmentTracks = 5; //Change !
  double MaxTheta_FragmentMDCTracks = 45; //Change !


  //Cut conditions on reconstructed pions
  int ifCut_MaxChi2ndf_PionTracks = 1;
  double MaxChi2ndf_PionTracks = 3.; //Change !

  int ifCut_MinDist_PionTracksPrimVtx = 0;
  double MinDist_PionTracksPrimVtx = 0.05; //Change !

  int ifCut_MinMomZ_PionTracks = 0;
  double MinMomZ_PionTracks = 0.; //Change !


  //Cut conditions on reconstructed hypernuclei
  int ifSet_ProductionVertex;
  int ifSet_MassConstraint;

  int KFPart_fConstructMethod = 2;

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


  //Cut conditions on silicons hits from reconstructed tracks
  int ifCut_MaxDist_SiHit = 1;
  double MaxDist_SiHit = 0.05; //Change !

  int ifCut_MinEnergyDeposition_SiHit = 1;
  double MinEnergyDeposition_SiHit = 0.1; //Change !


  //Silicon details (copy from TPrimaryVertex.h)
  double Z_plane_Si1x    = 27.; // in cm
  double Z_plane_Si1y    = 27.05; // in cm
  double widthStrip_Si1 = 0.019; // in cm
  double lenghtSi_Si1   = 9.728; // in cm
  double thicknessSi_Si1 = 0.03; // in cm
  bool restrict_actlenght_Si1 = true;
  double actlenghtX_Si1 = 6.08; // in cm
  double actlenghtY_Si1 = 6.08; // in cm

  double Z_plane_Si2x    = 30.; // in cm
  double Z_plane_Si2y    = 30.05; // in cm
  double widthStrip_Si2 = 0.019; // in cm
  double lenghtSi_Si2   = 9.728; // in cm
  double thicknessSi_Si2 = 0.03; // in cm
  bool restrict_actlenght_Si2 = true;
  double actlenghtX_Si2 = 6.08; // in cm
  double actlenghtY_Si2 = 6.08; // in cm

  double sigma_Si1      = widthStrip_Si1 / std::sqrt(12.);
  double sigma_Si2      = widthStrip_Si2 / std::sqrt(12.);

  float Z_plane;
  double widthStrip;
  double lenghtSi;
  bool restrict_actlenght;
  double actlenghtX;
  double actlenghtY;
  double sigma;
  int nStrips;



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

    TH1F* h_DecayVertexDistance_KFPart_PrimVtx_Mass;
    TH1F* h_DecayVertexDistanceX_KFPart_PrimVtx_Mass;
    TH1F* h_DecayVertexDistanceY_KFPart_PrimVtx_Mass;
    TH1F* h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass;

    TH1F* h_DecayVertexcutDistance;
    TH1F* h_DecayVertexcutDistanceX;
    TH1F* h_DecayVertexcutDistanceY;
    TH1F* h_DecayVertexcutDistanceZ;

/*
    TH1F* h_DecayVertexcutDistance_KFPart;
    TH1F* h_DecayVertexcutDistanceX_KFPart;
    TH1F* h_DecayVertexcutDistanceY_KFPart;
    TH1F* h_DecayVertexcutDistanceZ_KFPart;
*/

    TH1F* h_DecayVertexcutDistance_KFPart_PrimVtx;
    TH1F* h_DecayVertexcutDistanceX_KFPart_PrimVtx;
    TH1F* h_DecayVertexcutDistanceY_KFPart_PrimVtx;
    TH1F* h_DecayVertexcutDistanceZ_KFPart_PrimVtx;


    TH1F* h_DecayVertexPosZ_real;
    TH1F* h_DecayVertexPosZ_vfunction;
    TH1F* h_DecayVertexPosZ_centroid;
    TH1F* h_DecayVertexPosZ_KFPart;
    TH1F* h_DecayVertexPosZ_AllVfunc;
    TH1F* h_DecayVertexPosZ_AllCentroid;
    TH1F* h_DecayVertexPosZ_AllKFPart;


    TH2F* h_N_MotherTracks;
    TH2F* h_Dist_DaughterTracks;
    TH2F* h_Angle_MotherFragment;
    TH2F* h_Angle_MotherPion;
    TH2F* h_Chi2ndf_MotherTracks;
    TH2F* h_Dist_MotherTrackPrimVtx;
    TH2F* h_Theta_MotherTrackPrimVtx;
    TH2F* h_DecayVertexPosZ_KFPart_PrimVtx;
    TH2F* h_DecayFragmentMomZ_KFPart_PrimVtx;
    TH2F* h_DecayPionMomZ_KFPart_PrimVtx;
    TH2F* h_Hyp_ArmenterosPodolanski;
    TH2F* h_Hyp_CutArmenterosPodolanski;


    TH1F* h_HypInvariantMass;
    TH1F* h_HypErrorInvariantMass;

    TH1F* h_Hyp_RealLifeTime;
    TH1F* h_HypLifeTime_PrimVtx;
    TH1F* h_HypErrorLifeTime_PrimVtx;
    TH1F* h_HypcutLifeTime_PrimVtx;

    TH2F* h_HypInvariantMassCheck;
    TH2F* h_HypInvariantErrorMassCheck;

    TH1F* h_HypInvariantMass_LorentzVect;
    TH1F* h_HypInvariantMass_CutLorentzVect;

    TH1F* h_EffPosZ_real;
    TH1F* h_EffPosZ_preKF;
    TH1F* h_EffPosZ_postKF;
    TH1F* h_EffPosZ_preKFPart;
    TH1F* h_EffPosZ_postKFPart;

    TH2F* h_N_SiHits_ReconsTracks;
/*
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
*/
    TH1F* h_DecayVtxstats;
  };

  LocalHists LocalHisto;
};

#endif