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
#include "KFParticleSIMD.h"

//typedef TDataProcess<FullRecoEvent, MCAnaEventG4Sol> TDataProcessInterface;
template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent, Out>;


template<class Out>
class TDecayVertex final : public TDataProcessInterface<Out>
{
public:
  const THyphiAttributes& att;

  TDecayVertex(const THyphiAttributes& attr, int pi_type=0);
  ~TDecayVertex();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, Out* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderDecayVertex(FullRecoEvent& RecoEvent);

  void StudyCaseSelector_Hyp(std::string StudyCase);

  void RealTracksFinder(std::unordered_map<int, std::vector<std::vector<SimHit> > >& TrackDAFSim,
                         int& pdgParticle, int& cutConditions,
                         std::vector<KFParticle>& RealTracks,
                         std::vector<KFFitInfo>& Vect_FitInfo);

  void FragmentMDCTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results, int& fragment_pdg,
                                std::vector<FragmentTrack>& FragmentMDCTracks);

//  void FragmentMDCTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results, int& fragment_pdg,
//                                std::vector<KFParticle>& FragmentMDCTracks);

  void FragmentSelector(std::vector<KFParticle>& FragmentTracks_All, TVector3& PrimVtxRecons, std::vector<KFParticle>& FragmentTracks);

  void PionTracksFinder(std::unordered_map<int, ResSolDAF>& DAF_results,
                        std::vector<KFParticle>& PionTracks,
                        std::vector<KFFitInfo>& Vect_FitInfo);

  void PionSelector(std::vector<KFParticle>& PionTracks_All, TVector3& PrimVtxRecons, std::vector<KFParticle>& PionTracks,
                      std::vector<KFFitInfo>& Vect_FitInfo);

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
/*
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
*/
  void AllTrackstoDecayVertex_Vfunction(std::vector<KFParticle>& AllTracks, TVector3& Old_DecayVertexRecons, TVector3& DecayVertexRecons);

  void AllTrackstoDecayVertex_Centroids(std::vector<KFParticle>& AllTracks, TVector3& DecayVertexRecons);

/*
  //Silicon details (copy from TPrimaryVertex.h)

  double Z_plane_Si1x   = 27.4; // in cm
  double Z_plane_Si1y   = 27.45; // in cm
  double Z_plane_Si1    = (Z_plane_Si1x + Z_plane_Si1y)/2.;
  size_t combineStrips_Si1 = 4; // power of 2
  double widthStrip_Si1 = 0.008 * combineStrips_Si1; // in cm
  double sigma_Si1      = widthStrip_Si1 / std::sqrt(12.);
  double lenghtSi_Si1   = 2.048; // in cm
  double thicknessSi_Si1 = 0.0285; // in cm
  int nStrips_Si1 = static_cast<int>(lenghtSi_Si1 / widthStrip_Si1);

  double Z_plane_Si2x    = 28.2; // in cm
  double Z_plane_Si2y    = 28.25; // in cm
  double Z_plane_Si2     = (Z_plane_Si2x + Z_plane_Si2y)/2.;  // in cm
  size_t combineStrips_Si2 = 4; // power of 2
  double widthStrip_Si2 = 0.008 * combineStrips_Si2; // in cm
  double sigma_Si2       = widthStrip_Si2 / std::sqrt(12.);
  double lenghtSi_Si2   = 2.048; // in cm
  double thicknessSi_Si2 = 0.0285; // in cm
  int nStrips_Si2 = static_cast<int>(lenghtSi_Si2 / widthStrip_Si2);

  float Z_plane;
  double widthStrip;
  double lenghtSi;
  bool restrict_actlenght;
  double actlenghtX;
  double actlenghtY;
  double sigma;
  int nStrips;
*/

  double Zo_target; // in cm
  double Zf_target; // in cm
  double Zo_minifibers; // in cm


//  double Dist_to_Target = 0.5; // in cm
//  double Dist_to_Silicons = 0.5; // in cm
//  double Dist_to_Minifibers = 0.5; // in cm

  double k_factor       = 3.; // multiplies the fragment track f function value

  size_t NstepsdiscretXY      = 7;
  size_t NstepsdiscretZ       = 21;
  size_t Nstepsdiscretbox     = 5;
  size_t nTimesDiscretization = 7;
  double boxDistXY            = 3.;

  size_t NstepsdiscretboxXYZ  = 9;
  double boxXYZ               = 2.;
  size_t nTimesBoxXYZ         = 5;

  KFParticleFieldRegionHypHI fieldWASA;
  
  PDG_fromName pid_fromName;

  int He3_pdg = pid_fromName("He3");
  int He4_pdg = pid_fromName("alpha");
  int deuteron_pdg = pid_fromName("deuteron");
  int lambda_pdg = pid_fromName("lambda");
  int proton_pdg = pid_fromName("proton");
  int pi_pdg;
  int piminus_pdg = pid_fromName("pi-");
  int piplus_pdg = pid_fromName("pi+");

  int pion_type = -1;

  int H3L_pdg = 20001;
  int H4L_pdg = 20002;
  int nnL_pdg = 0; //CHECK Change!

  int Hyp_pdg;
  int Hyp_charge;
  double Hyp_mass;

  double c_light_speed = 299792458.; //in m/s
  double c_light_speed_cmps = c_light_speed * 1.e-10; //in cm/ps


  //Cut conditions on real pions
  int No_cutconditions = 0;
  int Yes_cutconditions = 1;


  //Methods for reconstructed fragments
  int Fragment_pdg;
  int recons_from_FRS_MDC = -1; // 1-> FRS, 2-> MDC

  int ref_RealFragment = -1;
  bool ifOnlyRealFragment = false;

  //Cut conditions on reconstructed fragments
  bool ifCut_MaxChi2ndf_FragmentTracks = false;
  double MaxChi2ndf_FragmentTracks = 3.; //Change !

  bool ifCut_MinDist_FragmentTracksPrimVtx = false;
  double MinDist_FragmentTracksPrimVtx = 0.; //Change !

  bool ifCut_MinMomZ_FragmentTracks = false;
  double MinMomZ_FragmentTracks = 0.; //Change !

  bool ifCut_MaxTheta_FragmentTracks = false;
  double MaxTheta_FragmentTracks = 5; //Change !
  double MaxTheta_FragmentMDCTracks = 45; //Change !


  //Cut conditions on reconstructed pions
  bool ifCut_MaxChi2ndf_PionTracks = true;
  double MaxChi2ndf_PionTracks = 300.; //Change !

  bool ifCut_MinDist_PionTracksPrimVtx = false;
  double MinDist_PionTracksPrimVtx = 0.05; //Change !

  bool ifCut_MinMomZ_PionTracks = false;
  double MinMomZ_PionTracks = 0.; //Change !


  //Cut conditions on reconstructed hypernuclei
  bool ifSet_ProductionVertex;
  bool ifSet_MassConstraint;

  int KFPart_fConstructMethod = 2;

  bool ifRemoveNaNMother = true;

  bool ifCut_MaxClosedist_DaughterTracks = false; 
  double MaxClosedist_DaughterTracks = 10.; //Change !

  bool ifCut_MaxAngle_MotherFragment = false;
  double MaxAngle_MotherFragment = 90.; //In degrees (º) Change !

  bool ifCut_MaxAngle_MotherPion = false;
  double MaxAngle_MotherPion = 90.; //In degrees (º) Change !

  bool ifCut_MaxChi2ndf = false;
  double MaxChi2ndf = 7.; //Change !

  bool ifCut_MaxDist_MotherTrackPrimVtx = false;
  double MaxDist_MotherTrackPrimVtx = 10.; //Change !

  bool ifCut_MaxAngle_MotherTrackPrimVtx = false;
  double MaxAngle_MotherTrackPrimVtx = 60.; //Change !

  bool ifCut_MaxPosZ_DecayVertex = false;
  double MaxPosZ_DecayVertex = 80.;

  bool ifCut_MinPosZ_DecayVertex = false;
  double MinPosZ_DecayVertex = 10.;

  bool ifCut_ArmenterosPodolanski = false;

/*
  //Cut conditions on silicons hits from reconstructed tracks
  bool ifCut_MaxDist_SiHit = true;
  double MaxDist_SiHit = 0.05; //Change !

  bool ifCut_MinEnergyDeposition_SiHit = true;
  double MinEnergyDeposition_SiHit = 0.1; //Change !
*/


  struct LocalHists
  {
    TH1F* h_P_fragments[2];
    TH1F* h_Pt_fragments[2];
    TH1F* h_Pz_fragments[2];
    TH1F* h_Dist_FragmentTrackPrimVtx[2];

    TH1F* h_P_pions[2];
    TH1F* h_Pt_pions[2];
    TH1F* h_Pz_pions[2];
    TH1F* h_Chi2ndf_pions[2];

    TH1F* h_Pt_cutpions[2];
    TH1F* h_Pz_cutpions[2];

    TH1F* h_Nrealpions[2];
    TH1F* h_Ncutpions[2];
    TH1F* h_Npions[2];


    TH1F* h_Closedist_Distance[2];
    TH1F* h_Closedist_PosZ[2];
    TH2F* h_Dist_DecayTrackPrimVtx[2];

    TH1F* h_Closedist_cutDistance[2];
    TH1F* h_Closedist_cutPosZ[2];
    TH2F* h_Dist_cutDecayTrackPrimVtx[2];


    TH1F* h_DecayVertexDistance[2];
    TH1F* h_DecayVertexDistanceX[2];
    TH1F* h_DecayVertexDistanceY[2];
    TH1F* h_DecayVertexDistanceZ[2];

    TH1F* h_DecayVertexDistance_centroid[2];
    TH1F* h_DecayVertexDistanceX_centroid[2];
    TH1F* h_DecayVertexDistanceY_centroid[2];
    TH1F* h_DecayVertexDistanceZ_centroid[2];

    TH1F* h_DecayVertexDistance_KFPart[2];
    TH1F* h_DecayVertexDistanceX_KFPart[2];
    TH1F* h_DecayVertexDistanceY_KFPart[2];
    TH1F* h_DecayVertexDistanceZ_KFPart[2];

    TH1F* h_DecayVertexDistance_KFPart_PrimVtx[2];
    TH1F* h_DecayVertexDistanceX_KFPart_PrimVtx[2];
    TH1F* h_DecayVertexDistanceY_KFPart_PrimVtx[2];
    TH1F* h_DecayVertexDistanceZ_KFPart_PrimVtx[2];

    TH1F* h_DecayVertexDistance_KFPart_PrimVtx_Mass[2];
    TH1F* h_DecayVertexDistanceX_KFPart_PrimVtx_Mass[2];
    TH1F* h_DecayVertexDistanceY_KFPart_PrimVtx_Mass[2];
    TH1F* h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass[2];

    TH1F* h_DecayVertexcutDistance[2];
    TH1F* h_DecayVertexcutDistanceX[2];
    TH1F* h_DecayVertexcutDistanceY[2];
    TH1F* h_DecayVertexcutDistanceZ[2];

/*
    TH1F* h_DecayVertexcutDistance_KFPart[2];
    TH1F* h_DecayVertexcutDistanceX_KFPart[2];
    TH1F* h_DecayVertexcutDistanceY_KFPart[2];
    TH1F* h_DecayVertexcutDistanceZ_KFPart[2];
*/

    TH1F* h_DecayVertexcutDistance_KFPart_PrimVtx[2];
    TH1F* h_DecayVertexcutDistanceX_KFPart_PrimVtx[2];
    TH1F* h_DecayVertexcutDistanceY_KFPart_PrimVtx[2];
    TH1F* h_DecayVertexcutDistanceZ_KFPart_PrimVtx[2];


    TH1F* h_DecayVertexPosZ_real[2];
    TH1F* h_DecayVertexPosZ_vfunction[2];
    TH1F* h_DecayVertexPosZ_centroid[2];
    TH1F* h_DecayVertexPosZ_KFPart[2];
    TH1F* h_DecayVertexPosZ_AllVfunc[2];
    TH1F* h_DecayVertexPosZ_AllCentroid[2];
    TH1F* h_DecayVertexPosZ_AllKFPart[2];


    TH2F* h_N_MotherTracks[2];
    TH2F* h_Dist_DaughterTracks[2];
    TH2F* h_Angle_MotherFragment[2];
    TH2F* h_Angle_MotherPion[2];
    TH2F* h_Chi2ndf_MotherTracks[2];
    TH2F* h_Dist_MotherTrackPrimVtx[2];
    TH2F* h_Theta_MotherTrackPrimVtx[2];
    TH2F* h_DecayVertexPosZ_KFPart_PrimVtx[2];
    TH2F* h_DecayFragmentMomZ_KFPart_PrimVtx[2];
    TH2F* h_DecayPionMomZ_KFPart_PrimVtx[2];
    TH2F* h_Hyp_ArmenterosPodolanski[2];
    TH2F* h_Hyp_CutArmenterosPodolanski[2];


    TH1F* h_HypInvariantMass[2];
    TH1F* h_HypInvariantMass_Z05[2];
    TH1F* h_HypInvariantMass_Z10[2];
    TH1F* h_HypInvariantMass_Z15[2];
    TH1F* h_HypInvariantMass_Z20[2];
    TH1F* h_HypErrorInvariantMass[2];

    TH1F* h_Hyp_RealLifeTime[2];
    TH1F* h_HypLifeTime_PrimVtx[2];
    TH1F* h_HypErrorLifeTime_PrimVtx[2];
    TH1F* h_HypcutLifeTime_PrimVtx[2];

    TH2F* h_HypInvariantMassCheck[2];
    TH2F* h_HypInvariantErrorMassCheck[2];

    TH1F* h_HypInvariantMass_LorentzVect[2];
    TH1F* h_HypInvariantMass_CutLorentzVect[2];

    TH1F* h_EffPosZ_real[2];
    TH1F* h_EffPosZ_preKF[2];
    TH1F* h_EffPosZ_postKF[2];
    TH1F* h_EffPosZ_preKFPart[2];
    TH1F* h_EffPosZ_postKFPart[2];

    TH2F* h_EffPosZPosR_real[2];
    TH2F* h_EffPosZPosR_postKFPart[2];
/*
    TH2F* h_N_SiHits_ReconsTracks[2];

    TH1F* h_N_Si_MotherTracks[2];


    TH1F* h_DecayVertexDistance_AllVfunc[2];
    TH1F* h_DecayVertexDistanceX_AllVfunc[2];
    TH1F* h_DecayVertexDistanceY_AllVfunc[2];
    TH1F* h_DecayVertexDistanceZ_AllVfunc[2];

    TH1F* h_DecayVertexDistance_AllCentroid[2];
    TH1F* h_DecayVertexDistanceX_AllCentroid[2];
    TH1F* h_DecayVertexDistanceY_AllCentroid[2];
    TH1F* h_DecayVertexDistanceZ_AllCentroid[2];

    TH1F* h_DecayVertexDistance_AllKFPart[2];
    TH1F* h_DecayVertexDistanceX_AllKFPart[2];
    TH1F* h_DecayVertexDistanceY_AllKFPart[2];
    TH1F* h_DecayVertexDistanceZ_AllKFPart[2];
*/
    TH1F* h_DecayVtxstats[2];
  };

  LocalHists LocalHisto;
};

#endif
