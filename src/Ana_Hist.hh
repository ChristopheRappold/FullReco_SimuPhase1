#ifndef HYPHIRHISTO_H
#define HYPHIRHISTO_H

//#include "ShadowAna_Hist.hh"

//#include <iostream>

#include <string>
#include <unordered_map>
#include <vector>

#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TGraphErrors.h"
#include "TEllipse.h"

//#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMath.h"

#include "spdlog/logger.h"

enum StateHist : int
{
  DAF = 0,
  VERTEX,
  DCPROJ,
  FINDING,
  RIEMANN,
  HOUGH,
  SIMU,
  BUILDER,
  PRIMVTX,
  PRIMVTX_SI,
  DECAYVTX,
  DECAYVTX_PIPLUS,
  FRAGMENT,
  WASA,
  SIZEOF_STATEHIST
};

template <typename T>
struct Hist {
  T* h;
  std::vector<TH1*> store;
  void emplace_back(T* _h) {
    h = _h;
    h->SetDirectory(0);
    store.emplace_back(_h);
  }
};

class Ana_Hist
{
  std::shared_ptr<spdlog::logger> _logger;
  
  public:
  bool HaveBeenWritten;
  std::vector<bool> EnableState;

  Hist<TH1I> h_stats;
  Hist<TH2I> h_statsLess3Mes;
  Hist<TH2I> h_statsInvalid;
  Hist<TH1I> h_task_exit;


  //Field
  Hist<TH2F> FieldXY[3];
  Hist<TH2F> FieldXZ[3];
  Hist<TH2F> FieldYZ[3];
  Hist<TH2F> FieldXYmax[3];
  Hist<TH2F> FieldXZmax[3];
  Hist<TH2F> FieldYZmax[3];
  Hist<TH2F> FieldXY_n[3];
  Hist<TH2F> FieldXZ_n[3];
  Hist<TH2F> FieldYZ_n[3];


  // Finder
  Hist<TH2F> h_xy;
  Hist<TH2F> h_PxPy;
  Hist<TH2F> h_xy_extrap;
  Hist<TH2F> h_PxPy_extrap;
  Hist<TH2F> h_TrackFindingStat;
  Hist<TH2F> h_MDC_Dphi;
  Hist<TH2F> h_MDC_NextLayer;
  Hist<TH2F> h_MDC_DiffLayer;
  Hist<TH2F> h_MDC_DiffLayer2;
  Hist<TH2F> h_MDC_InLayer;
  
  Hist<TH2F> h_SolenoidGeo[3];
  std::vector<TEllipse*> geoSolenoid;

  Hist<TH1F> h_RZStats;
  Hist<TH2F> h_RZ;
  Hist<TH1F> h_RZfit_mom;
  Hist<TH2F> h_RZfit_Chi2;
  Hist<TH2F> h_XYfit_miniF;
  Hist<TH2F> h_MDC_Z_residu;
  Hist<TH2F> h_MDC_R_residu;
  Hist<TH2F> h_MDC_Z_pull;
  Hist<TH2F> h_MDC_R_pull;


  // Kalman:
  Hist<TH1F> h_pv;
  Hist<TH1F> h_chi2;
  Hist<TH1F> hd_pv[2];
  Hist<TH1F> hd_chi[2];
  
  Hist<TH1F> h_Path;
  Hist<TH1F> h_Path_Back;
  Hist<TH1F> h_MeanPath;
  Hist<TH1F> h_dpath;
  
  Hist<TH1F> h_beta;
  Hist<TH1F> h_beta2;
  Hist<TH1F> h_beta3;
  
  Hist<TH1F> h_Mass_All;
  Hist<TH1F> h_Mass_All2;
  Hist<TH1F> h_Mass_All3;
  
  Hist<TH2F> h_Mass_charge_All;
  Hist<TH2F> h_Mass_charge_All2;
  Hist<TH2F> h_Mass_charge_All3;
  
  Hist<TH2F> h_beta_mom;
  Hist<TH2F> h_beta_mom2;
  Hist<TH2F> h_beta_mom3;
  Hist<TH2F> h_beta_momcharge;
  Hist<TH2F> h_beta_momcharge2;
  Hist<TH2F> h_beta_momcharge3;

  Hist<TH2F> h_pv_mom;
  Hist<TH2F> h_pv_beta;
  Hist<TH2F> h_pv_mass;
  
  Hist<TH2F> h_path_tof;
  
  Hist<TH2F> h_mom_tof_cut;
  Hist<TH2F> h_path_mom_cut;
  Hist<TH2F> h_path_tof_cut;
  
  Hist<TH1F> h_Mass[4];
  Hist<TH1F> h_chi2_particle[4];
  Hist<TH1F> h_pv_particle[4];
  
  Hist<TH2F> h_mom_res[5];
  Hist<TH1D> h_ResPull[5][10];
  Hist<TH2F> h_ResPull_normal[5][10];

  Hist<TH2F> h_total_dE;


  // residual
  Hist<TH1F> h_ResFiber[9];
  Hist<TH1F> h_ResMiniFiber[6];
  Hist<TH1F> h_ResMDC[17][3];
  Hist<TH1F> h_ResPSCE[2];


  // Riemann Finder
  Hist<TH2F> h_RiemannChi2;
  Hist<TH2F> h_RiemannResidus;


  // PerfFinder
  Hist<TH2F> h_PerfFinder;
  Hist<TH2F> h_PerfFinderLevenshtein;


  // WASABuilder
  Hist<TH1I> h_Builderstats;

  Hist<TH2D> h10[7][3];
  Hist<TH1D> h11[7][3];
  Hist<TH1D> h12[7][3];
  Hist<TH1D> h13[7][3];
  Hist<TH1D> h14[7][3];
  Hist<TH1D> h15[7][3];
  Hist<TH2D> h75[7][3];
  Hist<TH2D> h16[7];
  Hist<TH1D> h17[7];
  Hist<TH1D> h17_2[7];
  Hist<TH1D> h18_3_1;
  Hist<TH1D> h18_3_2;
  Hist<TH2D> h18_3_3;
  Hist<TH2D> h18_3_4;
  Hist<TH1D> h18_3_5;
  Hist<TH1D> h18_3_6;
  Hist<TH1D> h18_3_7;
  Hist<TH1D> h18_3_8;
  Hist<TH1D> hfiber_13_0[7][3];
  Hist<TH1D> hfiber_13_1[7][3];
  Hist<TH1D> hfiber_13_2[7][3];
  Hist<TH2D> hfiber_13_3[7][3];
  Hist<TH2D> hfiber_13_4[7][3];
  Hist<TH2D> h51[3][3][2];

  //MFT12
  Hist<TH1D> hfiber_1_1;
  Hist<TH1D> hfiber_1_2;
  Hist<TH1D> hfiber_1_3;
  Hist<TH1D> hfiber_1_4;
  Hist<TH1D> hfiber_1_5;
  Hist<TH1D> hfiber_1_6;
  Hist<TH1D> hfiber_1_7;
  Hist<TH1D> hfiber_1_9;
  Hist<TH1D> hfiber_2_1_1;
  Hist<TH1D> hfiber_2_1_2;
  Hist<TH1D> hfiber_2_2_1;
  Hist<TH1D> hfiber_2_2_2;
  Hist<TH1D> hfiber_2_3;
  Hist<TH1D> hfiber_3_0;
  Hist<TH1D> hfiber_3_0_2;
  Hist<TH1D> hfiber_6_1;
  Hist<TH1D> hfiber_6_2;
  Hist<TH2D> hfiber_6_3;
  Hist<TH2D> hfiber_6_4;
  Hist<TH2D> hfiber_12_1_1;
  Hist<TH2D> hfiber_12_2_1;
  Hist<TH2D> hfiber_12_3_1;
  Hist<TH2D> hfiber_12_1_2;
  Hist<TH2D> hfiber_12_2_2;
  Hist<TH2D> hfiber_12_3_2;

  //DFT12
  Hist<TH1D> hfiber_4_1;
  Hist<TH2D> hfiber_4_2_1;
  Hist<TH2D> hfiber_4_3_1;
  Hist<TH2D> hfiber_4_4_1;
  Hist<TH1D> hfiber_4_5_1;
  Hist<TH2D> hfiber_4_2_2;
  Hist<TH2D> hfiber_4_3_2;
  Hist<TH2D> hfiber_4_4_2;
  Hist<TH1D> hfiber_4_5_2;
  Hist<TH1D> hfiber_4_1_3;
  Hist<TH2D> hfiber_4_2_3;
  Hist<TH2D> hfiber_4_3_3;
  Hist<TH2D> hfiber_4_4_3;
  Hist<TH1D> hfiber_4_5_3;
  Hist<TH1D> hfiber_5_1;
  Hist<TH1D> hfiber_5_2;
  Hist<TH1D> hfiber_5_3;
  Hist<TH1D> hfiber_5_4;
  Hist<TH1D> hfiber_5_5;
  Hist<TH1D> hfiber_5_6;
  Hist<TH1D> hfiber_5_7;

  Hist<TH1D> hpsb_0_1;
  Hist<TH1D> hpsb_0_2;
  Hist<TH1D> hpsb_0_3;
  Hist<TH2D> hpsb_0_4;
  Hist<TH1D> hpsb_1_1;
  Hist<TH1D> hpsb_2[46];
  Hist<TH2D> hpsb_3[46];
  Hist<TH1D> hpsb_4[46];
  Hist<TH2D> h76;

  Hist<TH1D> hpsfe_0_1;
  Hist<TH1D> hpsfe_0_2;
  Hist<TH1D> hpsfe_0_3;
  Hist<TH2D> hpsfe_0_4;

  Hist<TH1D> hpsbe_0_1;
  Hist<TH1D> hpsbe_0_2;
  Hist<TH1D> hpsbe_0_3;
  Hist<TH2D> hpsbe_0_4;
  Hist<TH2D> hpsbe_1_0;

  Hist<TH1D> ht0_0_1;
  Hist<TH1D> ht0_0_2;
  Hist<TH1D> ht0_0_3;
  Hist<TH1D> ht0_0_4;
  Hist<TH1D> ht0_1[28];

  Hist<TH1D> hmdc_0_1;
  Hist<TH2D> hmdc_0_2;
  Hist<TH1D> hmdc_0_3;
  Hist<TH1D> hmdc_0_4;
  Hist<TH1D> hmdc_0_5;
  Hist<TH2D> hmdc_0_6;
  Hist<TH2D> hmdc_0_9;
  Hist<TH1D> hmdc_1[17];
  Hist<TH1D> hmdc_2[17];
  Hist<TH1D> hmdc_3[17];
  Hist<TH1D> hmdc_2_3[17];
  Hist<TH1D> hmdc_3_3[17];

  Hist<TH1D> hmwdc_1_1;
  Hist<TH1D> hmwdc_1_2;
  Hist<TH1D> hmwdc_1_3;
  Hist<TH1D> hmwdc_1_4;
  Hist<TH1D> hmwdc_1_5;
  Hist<TH1D> hmwdc_1_6;

  Hist<TH1D> hs4sci_1_1;
  Hist<TH1D> hs4sci_1_2;
  Hist<TH1D> hs4sci_1_3;
  Hist<TH1D> hs4sci_1_4;
  Hist<TH2D> hs4sci_2_1;
  Hist<TH2D> hs4sci_2_2;
  Hist<TH2D> hs4sci_2_3;
  Hist<TH2D> hs4sci_2_4;

  Hist<TH1D> htrig_0;
  Hist<TH2D> htrig_1;
  Hist<TH2D> htrig_2;
  Hist<TH1D> htrig_3;
  Hist<TH2D> htrig_4;


  // FragmentFinder
  Hist<TH1D> hopt_1_1;
  Hist<TH1D> hopt_1_2;
  Hist<TH1D> hopt_1_3;
  Hist<TH1D> hopt_1_4;
  Hist<TH2D> hopt_2_1;
  Hist<TH2D> hopt_2_2;
  Hist<TH1D> hopt_2_3;
  Hist<TH1D> hopt_2_4;


  // WASAFinder
  Hist<TH2D> h23_1;
  Hist<TH2D> h23_2;
  Hist<TH2D> h24_1;
  Hist<TH2D> h24_2;
  Hist<TH2D> h24_3[17];
  Hist<TH2D> h24_4[17];
  Hist<TH1D> h24_2_1;
  Hist<TH1D> h24_2_2;
  Hist<TH1D> h24_2_3[17];
  Hist<TH1D> h24_2_4[17];
  Hist<TH1D> hmdc_2_2[17];
  Hist<TH1D> hmdc_3_2[17];

  // Primary Vertex
  Hist<TH1F> h_HitMultiplicity_Si1;
  Hist<TH1F> h_HitMultiplicityRecons_Si1;
  Hist<TH1F> h_HitMultiplicityDiff_Si1;
  Hist<TH2F> h_HitMultiplicityDiffNHits_Si1;

  Hist<TH1F> h_EnergyDiffStrips_Si1;
  Hist<TH1F> h_nEventsGoodrecons_Si1;
  Hist<TH1F> h_nEventsGhost_Si1 ;
  Hist<TH2F> h_nEventsGoodreconsGhost_Si1 ;
  Hist<TH2F> h_nEventsRealGoodrecons_Si1 ;
  Hist<TH2F> h_nEventsRealRejectPad_Si1 ;

  Hist<TH1F> h_HitMultiplicity_Si2 ;
  Hist<TH1F> h_HitMultiplicityRecons_Si2 ;
  Hist<TH1F> h_HitMultiplicityDiff_Si2 ;
  Hist<TH2F> h_HitMultiplicityDiffNHits_Si2 ;

  Hist<TH1F> h_EnergyDiffStrips_Si2 ;
  Hist<TH1F> h_nEventsGoodrecons_Si2 ;
  Hist<TH1F> h_nEventsGhost_Si2 ;
  Hist<TH2F> h_nEventsGoodreconsGhost_Si2 ;
  Hist<TH2F> h_nEventsRealGoodrecons_Si2 ;
  Hist<TH2F> h_nEventsRealRejectPad_Si2 ;

  Hist<TH1F> h_MFCheck_Theta_MomSi1MomSi2 ;
  Hist<TH1F> h_MFCheck_Dist_MomSi1HitSi2 ;
  Hist<TH1F> h_MFCheck_Dist_MomSi2HitSi1 ;
  Hist<TH1F> h_MFCheck_Dist_MomSi1HitIP ;
  Hist<TH1F> h_MFCheck_Dist_MomSi2HitIP ;

  Hist<TH2F> h_EnergyStripEnergyTotalReal ;
  Hist<TH2F> h_EnergyStripEnergyTotal ;
  Hist<TH2F> h_EnergyDiffSilicons ;

  Hist<TH2F> h_EnergyDepositionMother ;
  Hist<TH2F> h_EnergyDepositionDaughters ;


  Hist<TH2F> h_NFiberMult;
  Hist<TH2F> h_NCombsXUV_UFT12;

  Hist<TH2F> h_NCombsXUV[5];
  Hist<TH2F> h_CombsXUV_dvalue_theta[5];
  Hist<TH2F> h_CombsXUV_dvalue_phi[5];
  Hist<TH2F> h_CombsXUV_dvalue_phi_theta5[5];
  Hist<TH2F> h_CombsXUV_dvalue_phi_theta10[5];

  Hist<TH2F> h_NHits_PrimaryTracks ;

  Hist<TH2F> h_nTrackCandidates ;
  Hist<TH2F> h_DistanceBeamTracks ;
  Hist<TH2F> h_PosZBeamTracks ;
  Hist<TH2F> h_thetaTracks ;
  Hist<TH2F> h_chi2ndfTracks ;
  Hist<TH1F> h_thetaResol ;

  Hist<TH1F> h_Acc_ThetaCandidates ;
  Hist<TH1F> h_Acc_ThetaAllReal ;

  Hist<TH2F> h_nCandidatesRealTracks ;
  Hist<TH2F> h_nCandidatesRealTracks_IfRecons ;

  Hist<TH1F> h_nHypernucleiTrack ;
  Hist<TH1F> h_fvalues ;

  Hist<TH1F> h_InteractionPointPosX ;
  Hist<TH1F> h_InteractionPointPosY ;
  Hist<TH1F> h_InteractionPointPosZ ;

  Hist<TH2F> h_InteractionPointDistance_V_value ;

  Hist<TH1F> h_InteractionPointDistance ;
  Hist<TH1F> h_InteractionPointDistanceX ;
  Hist<TH1F> h_InteractionPointDistanceY ;
  Hist<TH1F> h_InteractionPointDistanceZ ;

  Hist<TH1F> h_InteractionPointDistanceX_pull ;
  Hist<TH1F> h_InteractionPointDistanceY_pull ;
  Hist<TH1F> h_InteractionPointDistanceZ_pull ;

  Hist<TH1F> h_CovarianceSigmaX ;
  Hist<TH1F> h_CovarianceSigmaY ;
  Hist<TH1F> h_CovarianceSigmaZ ;

  Hist<TH1F> h_IP_DecayDistance ;
  Hist<TH1F> h_IP_DecayDistanceX ;
  Hist<TH1F> h_IP_DecayDistanceY ;
  Hist<TH1F> h_IP_DecayDistanceZ ;

  Hist<TH1F> h_DecayPositionDistance ;
  Hist<TH1F> h_DecayPositionDistanceX ;
  Hist<TH1F> h_DecayPositionDistanceY ;
  Hist<TH1F> h_DecayPositionDistanceZ ;
  
  Hist<TH2F> h_PrimStatus;
  Hist<TH1F> h_PrimVtxstats ;


  // Decay Vertex
  Hist<TH1F> h_P_fragments[2] ;
  Hist<TH1F> h_Pt_fragments[2] ;
  Hist<TH1F> h_Pz_fragments[2] ;
  Hist<TH1F> h_Dist_FragmentTrackPrimVtx[2] ;

  Hist<TH1F> h_P_pions[2] ;
  Hist<TH1F> h_Pt_pions[2] ;
  Hist<TH1F> h_Pz_pions[2] ;
  Hist<TH1F> h_Chi2ndf_pions[2] ;

  Hist<TH1F> h_Pt_cutpions[2] ;
  Hist<TH1F> h_Pz_cutpions[2] ;

  Hist<TH1F> h_Nrealpions[2] ;
  Hist<TH1F> h_Ncutpions[2] ;
  Hist<TH1F> h_Npions[2] ;

  Hist<TH1F> h_Closedist_Distance[2] ;
  Hist<TH1F> h_Closedist_PosZ[2] ;
  Hist<TH2F> h_Dist_DecayTrackPrimVtx[2] ;

  Hist<TH1F> h_Closedist_cutDistance[2] ;
  Hist<TH1F> h_Closedist_cutPosZ[2] ;
  Hist<TH2F> h_Dist_cutDecayTrackPrimVtx[2] ;

  Hist<TH1F> h_DecayVertexDistance[2] ;
  Hist<TH1F> h_DecayVertexDistanceX[2] ;
  Hist<TH1F> h_DecayVertexDistanceY[2] ;
  Hist<TH1F> h_DecayVertexDistanceZ[2] ;

  Hist<TH1F> h_DecayVertexDistance_centroid[2] ;
  Hist<TH1F> h_DecayVertexDistanceX_centroid[2] ;
  Hist<TH1F> h_DecayVertexDistanceY_centroid[2] ;
  Hist<TH1F> h_DecayVertexDistanceZ_centroid[2] ;

  Hist<TH1F> h_DecayVertexDistance_KFPart[2] ;
  Hist<TH1F> h_DecayVertexDistanceX_KFPart[2] ;
  Hist<TH1F> h_DecayVertexDistanceY_KFPart[2] ;
  Hist<TH1F> h_DecayVertexDistanceZ_KFPart[2] ;

  Hist<TH1F> h_DecayVertexDistance_KFPart_PrimVtx[2] ;
  Hist<TH1F> h_DecayVertexDistanceX_KFPart_PrimVtx[2] ;
  Hist<TH1F> h_DecayVertexDistanceY_KFPart_PrimVtx[2] ;
  Hist<TH1F> h_DecayVertexDistanceZ_KFPart_PrimVtx[2] ;

  Hist<TH1F> h_DecayVertexDistance_KFPart_PrimVtx_Mass[2] ;
  Hist<TH1F> h_DecayVertexDistanceX_KFPart_PrimVtx_Mass[2] ;
  Hist<TH1F> h_DecayVertexDistanceY_KFPart_PrimVtx_Mass[2] ;
  Hist<TH1F> h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass[2] ;

  Hist<TH1F> h_DecayVertexcutDistance[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceX[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceY[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceZ[2] ;

/*
  Hist<TH1F> h_DecayVertexcutDistance_KFPart[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceX_KFPart[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceY_KFPart[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceZ_KFPart[2] ;
*/

  Hist<TH1F> h_DecayVertexcutDistance_KFPart_PrimVtx[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceX_KFPart_PrimVtx[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceY_KFPart_PrimVtx[2] ;
  Hist<TH1F> h_DecayVertexcutDistanceZ_KFPart_PrimVtx[2] ;

  Hist<TH1F> h_DecayVertexPosZ_real[2] ;
  Hist<TH1F> h_DecayVertexPosZ_vfunction[2] ;
  Hist<TH1F> h_DecayVertexPosZ_centroid[2] ;
  Hist<TH1F> h_DecayVertexPosZ_KFPart[2] ;
  Hist<TH1F> h_DecayVertexPosZ_AllVfunc[2] ;
  Hist<TH1F> h_DecayVertexPosZ_AllCentroid[2] ;
  Hist<TH1F> h_DecayVertexPosZ_AllKFPart[2] ;

  Hist<TH2F> h_N_MotherTracks[2] ;
  Hist<TH2F> h_Dist_DaughterTracks[2] ;
  Hist<TH2F> h_Angle_MotherFragment[2] ;
  Hist<TH2F> h_Angle_MotherPion[2] ;
  Hist<TH2F> h_Chi2ndf_MotherTracks[2] ;
  Hist<TH2F> h_Dist_MotherTrackPrimVtx[2] ;
  Hist<TH2F> h_Theta_MotherTrackPrimVtx[2] ;
  Hist<TH2F> h_DecayVertexPosZ_KFPart_PrimVtx[2] ;
  Hist<TH2F> h_DecayFragmentMomZ_KFPart_PrimVtx[2] ;
  Hist<TH2F> h_DecayPionMomZ_KFPart_PrimVtx[2] ;
  Hist<TH2F> h_Hyp_ArmenterosPodolanski[2] ;
  Hist<TH2F> h_Hyp_CutArmenterosPodolanski[2] ;
  
  Hist<TH1F> h_HypInvariantMass[2] ;
  Hist<TH1F> h_HypInvariantMass_Z05[2] ;
  Hist<TH1F> h_HypInvariantMass_Z10[2] ;
  Hist<TH1F> h_HypInvariantMass_Z15[2] ;
  Hist<TH1F> h_HypInvariantMass_Z20[2] ;
  Hist<TH1F> h_HypErrorInvariantMass[2] ;

  Hist<TH1F> h_Hyp_RealLifeTime[2] ;
  Hist<TH1F> h_HypLifeTime_PrimVtx[2] ;
  Hist<TH1F> h_HypErrorLifeTime_PrimVtx[2] ;
  Hist<TH1F> h_HypcutLifeTime_PrimVtx[2] ;

  Hist<TH2F> h_HypInvariantMassCheck[2] ;
  Hist<TH2F> h_HypInvariantErrorMassCheck[2] ;

  Hist<TH1F> h_HypInvariantMass_LorentzVect[2] ;
  Hist<TH1F> h_HypInvariantMass_CutLorentzVect[2] ;

  Hist<TH1F> h_EffPosZ_real[2] ;
  Hist<TH1F> h_EffPosZ_preKF[2] ;
  Hist<TH1F> h_EffPosZ_postKF[2] ;
  Hist<TH1F> h_EffPosZ_preKFPart[2] ;
  Hist<TH1F> h_EffPosZ_postKFPart[2] ;

  Hist<TH2F> h_EffPosZPosR_real[2] ;
  Hist<TH2F> h_EffPosZPosR_postKFPart[2] ;

/*
  Hist<TH1F> h_N_Si_MotherTracks[2] ;

  Hist<TH1F> h_DecayVertexDistance_AllVfunc[2] ;
  Hist<TH1F> h_DecayVertexDistanceX_AllVfunc[2] ;
  Hist<TH1F> h_DecayVertexDistanceY_AllVfunc[2] ;
  Hist<TH1F> h_DecayVertexDistanceZ_AllVfunc[2] ;

  Hist<TH1F> h_DecayVertexDistance_AllCentroid[2] ;
  Hist<TH1F> h_DecayVertexDistanceX_AllCentroid[2] ;
  Hist<TH1F> h_DecayVertexDistanceY_AllCentroid[2] ;
  Hist<TH1F> h_DecayVertexDistanceZ_AllCentroid[2] ;

  Hist<TH1F> h_DecayVertexDistance_AllKFPart[2] ;
  Hist<TH1F> h_DecayVertexDistanceX_AllKFPart[2] ;
  Hist<TH1F> h_DecayVertexDistanceY_AllKFPart[2] ;
  Hist<TH1F> h_DecayVertexDistanceZ_AllKFPart[2] ;
*/

  Hist<TH1F> h_DecayVtxstats[2] ;


  // CheckFiberXUV
  Hist<TH1F> h_ResidualFiberHitX[7];
  Hist<TH1F> h_ResidualFiberHitY[7];
  Hist<TH1F> h_ResidualFiberHitR[7];
  Hist<TH2F> h_ResidualFiberHitXY[7];
  Hist<TH2F> h_ResidualFiberHitX_Theta[7];
  Hist<TH2F> h_ResidualFiberHitY_Theta[7];
  Hist<TH2F> h_ResidualFiberHitX_HitX[7];
  Hist<TH2F> h_ResidualFiberHitY_HitY[7];
  Hist<TH2F> h_ResidualFiberHitR_Theta[7];
  Hist<TH1F> h_ResidualSingleFiberHitX[7];
  Hist<TH1F> h_ResidualSingleFiberHitY[7];
  Hist<TH1F> h_ResidualSingleFiberHitR[7];
  Hist<TH2F> h_ResidualSingleFiberHitX_Theta[7];
  Hist<TH2F> h_ResidualSingleFiberHitY_Theta[7];
  Hist<TH2F> h_ResidualSingleFiberHitR_Theta[7];
  Hist<TH2F> h_EfficiencyFiberHit;
  Hist<TH2F> h_EfficiencySingleFiberHit;
  Hist<TH2F> h_EfficiencyFiberHit_Theta[7];
  Hist<TH2F> h_EfficiencyFiberHit_dvalue[7];
  Hist<TH2F> h_EfficiencyFiberHit_mult[7];
  Hist<TH2F> h_EfficiencySingleFiberHit_Theta[7];
  Hist<TH2F> h_EfficiencySingleFiberHit_dvalue[7];
  Hist<TH2F> h_NumFiberHit_GoodReco[7];
  Hist<TH2F> h_NumFiberHit_Ghost[7];
  Hist<TH1F> h_FiberHit_dvalue[7];
  Hist<TH1F> h_FiberHitSingle_dvalue[7];
  Hist<TH1F> h_FiberHitReal_dvalue[7];
  Hist<TH2F> h_FiberHitReal_dvalue_Theta[7];
  Hist<TH2F> h_FiberHitReal_dvalue_Phi[7];
  Hist<TH2F> h_FiberHitReal_dvalue_Theta03_Phi[7];
  Hist<TH2F> h_FiberHitReal_dvalue_Theta310_Phi[7];
  Hist<TH2F> h_FiberHitReal_dvalue_Theta1020_Phi[7];
  Hist<TH2F> h_FiberHitReal_dvalue_HitX[7];
  Hist<TH2F> h_FiberHitReal_dvalue_HitY[7];
  Hist<TH2F> h_FiberHitReal_dvalue_PosX[7];
  Hist<TH2F> h_FiberHitReal_dvalue_tanThetacosPhi[7];
  Hist<TH2F> h_FiberHitReal_dvalue_dfunction[7];
  Hist<TH1F> h_FiberHit_Residualdvalue[7];
  Hist<TH2F> h_FiberHit_Residualdvalue_Realdvalue[7];


  // CheckFiberTrack
  Hist<TH2F> h_ResidualFiberX;
  Hist<TH2F> h_ResidualFiberY;

  Hist<TH2F> h_ResidualFiberX_Angle[11][2];
  Hist<TH2F> h_ResidualFiberY_Angle[11][2];

  Hist<TH2F> h_ResidualFiberDzDphi;
  Hist<TH2F> h_ResidualFiberDzDtheta;

  Hist<TH2F> h_ResidualMDC_dZ1;
  Hist<TH2F> h_ResidualMDC_dZ2;

  Hist<TH2F> h_RPhiZMDC_Sigma;


  Hist<TH2F> h_ResidualMDC_dZ_PSB;
  Hist<TH2F> h_ResidualMDC_dZ_PSBE;
  Hist<TH2F> h_ResidualMDC_dZ_PSFE;
  Hist<TH2F> h_ResidualMDC_dZ_More6;

  std::unordered_map<std::string, std::tuple<std::vector<std::vector<TH1*>*>, int> > HistRegisteredByDir;

  Ana_Hist(bool Daf = true, bool Vertex = true, bool DCproject = true, bool Finding = true, bool Riemann = true, bool Hough = true, bool Simu = false,
            bool Builder = false, bool PrimVtx = true, bool PrimVtx_Si = false, bool DecayVtx = true,  bool DecayVtx_piplus = true,
            bool FragmentFinder = true, bool WASAFinder = true);
  ~Ana_Hist();

  int Write(TFile*);
  int WriteTemp(char* tempfile);

  void DebugHists();
  
  template <typename T>
  T* CloneAndRegister(Hist<T>& hist) {
    T* newHist = dynamic_cast<T*>(hist.h->Clone());
    hist.store.emplace_back(newHist);
    newHist->SetDirectory(0);
    return newHist;
  }

  TDirectory* GetDir(TFile* out_file, const std::string& name_dir)
  {
    if(!out_file->GetDirectory(name_dir.c_str()))
      return out_file->mkdir(name_dir.c_str());
    else
      return out_file->GetDirectory(name_dir.c_str());
  }
  
};

#endif
