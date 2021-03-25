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
  HOUGH,
  SIMU,
  PRIMVTX,
  DECAYVTX,
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
  // K>lman:
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

  //residual
  Hist<TH1F> h_ResFiber[9];
  Hist<TH1F> h_ResMiniFiber[6];
  Hist<TH1F> h_ResMDC[17][3];
  Hist<TH1F> h_ResPSCE[2];

  //Primary Vertex
  Hist<TH1F> h_HitMultiplicity_Si1;
  Hist<TH1F> h_HitMultiplicityRecons_Si1;
  Hist<TH1F> h_HitMultiplicityDiff_Si1;
  Hist<TH2F> h_HitMultiplicityDiffNHits_Si1;

  Hist<TH1F> h_EnergyDiffStrips_Si1;

  Hist<TH1F> h_nEventsGoodrecons_Si1;
  Hist<TH1F> h_nEventsGhost_Si1 ;
  Hist<TH2F> h_nEventsGoodreconsGhost_Si1 ;
  Hist<TH2F> h_nEventsRealGoodrecons_Si1 ;
  


  Hist<TH1F> h_HitMultiplicity_Si2 ;
  Hist<TH1F> h_HitMultiplicityRecons_Si2 ;
  Hist<TH1F> h_HitMultiplicityDiff_Si2 ;
  Hist<TH2F> h_HitMultiplicityDiffNHits_Si2 ;

  Hist<TH1F> h_EnergyDiffStrips_Si2 ;

  Hist<TH1F> h_nEventsGoodrecons_Si2 ;
  Hist<TH1F> h_nEventsGhost_Si2 ;
  Hist<TH2F> h_nEventsGoodreconsGhost_Si2 ;
  Hist<TH2F> h_nEventsRealGoodrecons_Si2 ;



  Hist<TH2F> h_EnergyStripEnergyTotalReal ;
  Hist<TH2F> h_EnergyStripEnergyTotal ;
  Hist<TH2F> h_EnergyDiffSilicons ;

  Hist<TH2F> h_EnergyDepositionMother ;
  Hist<TH2F> h_EnergyDepositionDaughters ;

  Hist<TH2F> h_nTrackCandidates ;
  Hist<TH2F> h_DistanceBeamTracks ;
  Hist<TH2F> h_PosZBeamTracks ;
  Hist<TH2F> h_thetaTracks ;

  //Hist<TH1F> h_nHypTrackReal;
  Hist<TH1F> h_nHypernucleiTrack ;
  Hist<TH1F> h_fvalues ;

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

  //Decay Vertex
  Hist<TH1F> h_Pt_fragments ;
  Hist<TH1F> h_Pz_fragments ;
  Hist<TH1F> h_Dist_FragmentTrackPrimVtx ;

  Hist<TH1F> h_Pt_pions ;
  Hist<TH1F> h_Pz_pions ;
  Hist<TH1F> h_Chi2ndf_pions ;

  Hist<TH1F> h_Pt_realpions ;
  Hist<TH1F> h_Pz_realpions ;

  Hist<TH1F> h_Pt_cutpions ;
  Hist<TH1F> h_Pz_cutpions ;

  Hist<TH1F> h_Nrealpions ;
  Hist<TH1F> h_Ncutpions ;
  Hist<TH1F> h_Npions ;

  Hist<TH1F> h_Closedist_Distance ;
  Hist<TH1F> h_Closedist_PosZ ;
  Hist<TH2F> h_Dist_DecayTrackPrimVtx ;

  Hist<TH1F> h_Closedist_realDistance ;
  Hist<TH1F> h_Closedist_realPosZ ;
  Hist<TH2F> h_Dist_realDecayTrackPrimVtx ;

  Hist<TH1F> h_Closedist_cutDistance ;
  Hist<TH1F> h_Closedist_cutPosZ ;
  Hist<TH2F> h_Dist_cutDecayTrackPrimVtx ;

  Hist<TH1F> h_DecayVertexDistance ;
  Hist<TH1F> h_DecayVertexDistanceX ;
  Hist<TH1F> h_DecayVertexDistanceY ;
  Hist<TH1F> h_DecayVertexDistanceZ ;

  Hist<TH1F> h_DecayVertexDistance_centroid ;
  Hist<TH1F> h_DecayVertexDistanceX_centroid ;
  Hist<TH1F> h_DecayVertexDistanceY_centroid ;
  Hist<TH1F> h_DecayVertexDistanceZ_centroid ;

  Hist<TH1F> h_DecayVertexDistance_2centroid_average ;
  Hist<TH1F> h_DecayVertexDistanceX_2centroid_average ;
  Hist<TH1F> h_DecayVertexDistanceY_2centroid_average ;
  Hist<TH1F> h_DecayVertexDistanceZ_2centroid_average ;

  Hist<TH1F> h_DecayVertexDistance_2centroid_closest ;
  Hist<TH1F> h_DecayVertexDistanceX_2centroid_closest ;
  Hist<TH1F> h_DecayVertexDistanceY_2centroid_closest ;
  Hist<TH1F> h_DecayVertexDistanceZ_2centroid_closest ;

  Hist<TH1F> h_DecayVertexDistance_2centroid_IPCheck ;
  Hist<TH1F> h_DecayVertexDistanceX_2centroid_IPCheck ;
  Hist<TH1F> h_DecayVertexDistanceY_2centroid_IPCheck ;
  Hist<TH1F> h_DecayVertexDistanceZ_2centroid_IPCheck ;

  Hist<TH1F> h_DecayVertexrealDistance ;
  Hist<TH1F> h_DecayVertexrealDistanceX ;
  Hist<TH1F> h_DecayVertexrealDistanceY ;
  Hist<TH1F> h_DecayVertexrealDistanceZ ;

  Hist<TH1F> h_DecayVertexcutDistance ;
  Hist<TH1F> h_DecayVertexcutDistanceX ;
  Hist<TH1F> h_DecayVertexcutDistanceY ;
  Hist<TH1F> h_DecayVertexcutDistanceZ ;

  Hist<TH1F> h_DecayVertexPosZ_real ;
  Hist<TH1F> h_DecayVertexPosZ_vfunction ;
  Hist<TH1F> h_DecayVertexPosZ_centroid ;

  Hist<TH1F> h_Dist_MotherTrackPrimVtx ;
  Hist<TH1F> h_Theta_MotherTrackPrimVtx ;
  Hist<TH1F> h_HypInvariantMass ;
  Hist<TH1F> h_N_Si_MotherTracks ;

  Hist<TH1F> h_DecayVtxstats ;

  std::unordered_map<std::string, std::tuple<std::vector<std::vector<TH1*>*>, int> > HistRegisteredByDir;

  Ana_Hist(bool Daf = true, bool Vertex = true, bool DCproject = true, bool Finding = true, bool Hough = true, bool Simu = false, bool PrimVtx = true, bool DecayVtx = true);
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
