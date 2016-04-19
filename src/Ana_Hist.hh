#ifndef HYPHIRHISTO_H
#define HYPHIRHISTO_H

//#include "ShadowAna_Hist.hh"

#include <iostream>

#include <vector>
#include <map>
#include <string>

#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH1I.h"

//#include "TCanvas.h"
#include "TMath.h"
#include "TDirectory.h"

enum StateHist {DAF,VERTEX,DCPROJ,FINDING,HOUGH,SIMU,SIZEOF_STATEHIST};

class Ana_Hist {

public:
  bool HaveBeenWritten;
  std::vector<bool> EnableState;
  TH1I* h_task_exit;

  std::vector<TH2F*> h_2utrPos_XY; // Utrackers positions
  std::vector<TH2F*> h_2utrPos_XZ; // Utrackers positions
  std::vector<TH2F*> h_2utrPos_YZ; // Utrackers positions
  std::vector<TH1F*> h_utrEnergy; // Utrackers energy deposit
  std::vector<TH1F*> h_utrTime;
  std::vector<TH1F*> h_pull_X_det;
  std::vector<TH1F*> h_pull_Y_det;

  std::vector<TH1F*> h_layer_x;
  //std::vector<TH1F*> h_layer_y;
  //std::vector<TH1F*> h_layer_z;
  //std::vector<TH1F*> h_layer_t;
  //std::vector<TH1F*> h_layer_E;
  //std::vector<TH1F*> h_layer_PID;

  TH2F* h_2utrRotated_XZ;
  TH2F* h_2utrNotRotated_XZ;
  TH2F* h_2utrFrameMag_XZ;

  TH2F* h_TOFP_res;
  TH2F* h_TOFP_diff_ResInter;
  TH2F* h_TOFP_Rdiff_ResInter;
  TH2F* h_TOFP_res_Path[96];

  //std::vector<TH2F*> fHough_Z400;
  //std::vector<TH2F*> fHough_Z750;


  TH1F *pattern;
  // TH1F *alpha_Ek;
  // TH1F *alpha_Theta;
  // TH1F *alpha_Phi;
  // TH1F *decay_pos_z;
  //
  //TH2F *h_2TR1; // Position for fiber tracker 1
  //TH2F *h_2TR2; // Position for fiber tracker 1
  //TH2F *h_2AladinTOF; // Position ALADIN TOF wall

  TH2I* h_alpha_rejection;
  TH1I* h_presort; 
  TH1I* h_presort_comb; 
  TH2F* h_chi2_prefitY;
  TH1F* h_rchi2_prefitY;
  // TH1F *h_dmom1[4];
  // TH1F *h_dmom2[4];
  // TH2F *h_dmom3[4];
  // TH2F *h_dmom4[4];

  // TH1F *h_dmom_dmass[4];
  // TH1F *h_dmom_dmass_best[4];
  // TH1F *h_dmom_vertex[4];
  // TH1F *h_dmom_vertex_best[4];

  TH1F *h_Mass[4];

  // TH1F* h_slope_yz[4];
  // TH1F* h_dy;

  TH1I *h_stats;
  // TH1I *h_stats_Mom;
  TH1I *h_detection_eff_pion;
  TH1I *h_detection_eff_alpha;
  TH1I *h_detection_eff_proton;

  TH1F *h_beta;
  TH1F *h_Mass_All;

  // TH1I *h_hough;
  // TH1I *h_hough_pi;
  // TH2I *h_hough_pi_dist_tof;
  // TH2I *h_hough_pi_dist_dc;
  // TH2I *h_hough_pi_dc;
  // TH2I *h_hough_pi_tof;

  // TH2I *h_hough_pi_dc_true;
  // TH2I *h_hough_pi_dc_false;
  // TH2I *h_hough_pi_dc_true_tof_true;
  // TH2I *h_hough_pi_tof_true;
  // TH2I *h_hough_pi_tof_false;
  // TH2I *h_hough_pi_tof_true_dc_true;


  // TH2I *h_hough_2d;
  // TH1I *h_hough_optimal;
  // TH2I *h_hough_2d_optimal;
  TH2F* h_hough_eff;
  TH1F* h_hough_diff;
  TH1F* h_hough_TOF_dist;
  TH1F* h_hough_DC_dist;
  TH1F* h_hough_TOF_dist_2;
  TH1F* h_hough_DC_dist_2;

  TH1F* h_hough_tot_dist;
  TH1F* h_hough_tot_dist_2;

  TH1F* h_hough_tot_pv;

  //TH2F* h_hough_each;
  TH2F* h_patern_decay;

  TH1F* h_finding_chi2_tot_pid[7][3];
  TH1F* h_finding_chi2_x_pid[7][3];
  TH1F* h_finding_chi2_y_pid[7][3];
  TH2F* h_finding_chi2_2D[7];
  TH2F* h_finding_chi2_2D_ghost[7];

  TH2F* h_finding_eff;
  TH2F* h_finding_multi; 
  // TH1F* h_mass_invariant_cut_mass;
  // TH1F* h_mass_invariant_best_cut_mass;
  // TH1F* h_mass_invariant_cut_vertex;
  // TH1F* h_mass_invariant_best_cut_vertex;
  TH2F* h_beta_mom;
  TH2F* h_beta_mom2;
  TH2F* h_beta_mom3;

  TH2F* h_mom_check;
  TH2F* h_mom_check2;
  TH2F* h_tofdE;

  TH2F* h_path_tof;
  TH2F* h_pv_mom;
  TH2F* h_pv_beta;
  TH2F* h_pv_mass;

  // TH1F* h_vtx_x;
  // TH1F* h_vtx_y;
  // TH1F* h_vtx_z;

  // TH1F* h_vtx_x_cut;
  // TH1F* h_vtx_y_cut;
  // TH1F* h_vtx_z_cut;

  // TH1F* h_vtx_x_cut_best;
  // TH1F* h_vtx_y_cut_best;
  // TH1F* h_vtx_z_cut_best;

  //kalman

  TH1F* h_pv;
  TH1F* h_chi2;

  TH1F* h_ResidualX[7][4];
  TH1F* h_ResidualY[7][4];

  TH1F* h_weigth_max[4][19];
  TH1F* h_weigth_min[4][19];
  TH1F* h_weigth_sum[4][19];
  TH1F* h_weigth_one[4][19];
  TH1F* h_chi2_det[4][19];
  TH1F* h_chi2_det_one[4][19];

  TH1F* h_dB;
  TH1F* h_Bchi2;
  TH2F* h_Bz;

  //TH1F* h_chi2_smooth;
  //TH1F* h_pull[5];
   
  TH1F* h_chi2_particle[4];
  TH1F* h_pv_particle[4];

  TH1F* h_Path;
  TH1F* h_Path_Back;
  TH1F* h_MeanPath;

  TH1F* h_Path_newL;
  TH1F* h_Path_LTr1toTOF;
  TH1F* h_Path_LTr0toTr1;
  TH1F* h_Path_check1;
  TH1F* h_Path_check2;

  TH2F* h_dpath_new;
  TH2F* h_Time_Estimation;
  TH1F* h_Time;
  

  TH1F* h_beta2;
  TH1F* h_Mass_All2;
  TH1F* h_Mass2[3][4];

  TH1F* h_beta3;
  TH1F* h_Mass_All3;

  TH2F* h_Mass_charge_All;
  TH2F* h_Mass_charge_All2;
  TH2F* h_Mass_charge_All3;

  TH2F* h_mom_tof_cut;
  TH2F* h_path_tof_cut;
  TH2F* h_path_mom_cut;
  TH1F* h_dpath;
  TH2F* h_dpath2_chi2;
  TH2F* h_dpath2_pv;

  TH2F* h_mass_xTof;
  TH2F* h_mom_xTof ;
  TH2F* h_beta_xTof;
  TH2F* h_chi2_xTof;
  TH2F* h_pv_xTof;

  // TH1F* h_dmom_kalman[3][4];
  // TH1F* h_dmom_normal_kalman[3][4];
  // TH2F* h_dmom_mom_kalman[3][4];

  // TH1F *h_dmom_dmass2[3][4];
  // TH1F *h_dmom_dmass_best2[3][4];
  // TH1F *h_dmom_vertex2[3][4];
  // TH1F *h_dmom_vertex_best2[3][4];

  TH1F* h_mass_invariant_cut_mass2[3];
  TH1F* h_mass_invariant_best_cut_mass2[3];
  TH1F* h_mass_invariant_cut_vertex2[3];
  TH1F* h_mass_invariant_best_cut_vertex2[3];

  TH2F* h_momD1_momD2_cut_mass[3][4];
  TH2F* h_momD1_momD2_best_cut_mass[3][4];
  TH2F* h_momD1_momD2_cut_vertex[3][4];
  TH2F* h_momD1_momD2_best_cut_vertex[3][4];

  // TH1F* h_vtx_x2[3];
  // TH1F* h_vtx_y2[3];
  // TH1F* h_vtx_z2[3];

  // TH1F* h_vtx_x_cut2[3];
  // TH1F* h_vtx_y_cut2[3];
  // TH1F* h_vtx_z_cut2[3];

  // TH1F* h_vtx_x_cut_best2[3];
  // TH1F* h_vtx_y_cut_best2[3];
  // TH1F* h_vtx_z_cut_best2[3];

  TH1F* h_chi2_cut[3];
  TH1F* h_chi2_cut_best_cut_mass[3];
  TH1F* h_chi2_cut_cut_vtx[3];
  TH1F* h_chi2_cut_best_cut_vtx[3];
  TH1F* h_chi2_cut_cut_mass[3];
  //TH1F* h_chi2_cut[3];
  
  //Chamber Projection
  TH2I* h_dc1_project;
  TH2I* h_dc1x_project;
  TH2I* h_dc1u_project;
  TH2I* h_dc1v_project;

  TH2I* h_dc2_project;
  TH2I* h_dc2x_project;
  TH2I* h_dc2y_project;
  TH2I* h_dc2u_project;

  //Daisuke
  //#ifdef DAISUKE
  TH1F *h_status;
  TH2F *h_MassZ1BarBest;
  TH2F *h_MassZ1BarBestTofs;
  TH2F *h_MassZ1BarBest_ini;
  TH2F *h_MassZ1BarBest_iniTofs;

//   TH2F *h_MassZ1BarBest3;
//   TH2F *h_MassZ1BarBestTofs3;
//   TH2F *h_MassZ1BarBest_ini3;
//   TH2F *h_MassZ1BarBest_iniTofs3;

  TH2F *h_MassZ2BarBest;
  TH2F *h_MassZ2BarBestTofs;
  TH2F *h_MassZ2BarBest_ini;
  TH2F *h_MassZ2BarBest_iniTofs;

//   TH2F *h_MassZ2BarBest3;
//   TH2F *h_MassZ2BarBestTofs3;
//   TH2F *h_MassZ2BarBest_ini3;
//   TH2F *h_MassZ2BarBest_iniTofs3;

  TH2F *h_MassZ1BarBest_pi;
  TH2F *h_MassZ1BarBest_piTofs;
  TH2F *h_MassZ1BarBest_ini_pi;
  TH2F *h_MassZ1BarBest_ini_piTofs;

//   TH2F *h_MassZ1BarBest_pi3;
//   TH2F *h_MassZ1BarBest_piTofs3;
//   TH2F *h_MassZ1BarBest_ini_pi3;
//   TH2F *h_MassZ1BarBest_ini_piTofs3;

  TH2F *hd_MassPv[2];
  TH2F *h_MomCor[3];
  //two body
//   TH2F *hd_secz_all[3][3];
//   TH2F *hd_secz[3][3];
//   TH2F *hd_secz_cut[3][3];
//   TH2F *hd_secz_cut2[3][3];
//   TH2F *hd_secz_cut3[3][3];
//   TH2F *hd_secz_mompi[3][3];
//   TH1F *hd_invmass[4][3][3];
//   TH2F *hd_invmass_secz[4][3][3];
  
  TH2F *hd_theWorld;
  TH1F *hd_chi[2];
  TH1F *hd_pv[2];
//   TH2F* hd_BetaMom_Tp[32];
//   TH2F* hd_BetaMom_Bg[18];
  TH2F* hd_BetaMom_Tp_all;
  TH2F* hd_BetaMom_Tp_allTofs;
  TH2F* hd_BetaMom_Tp_allTofs_Z2;
  TH2F* hd_BetaMom_Tp_ini_all;
  TH2F* hd_BetaMom_Tp_ini_allTofs;
  TH2F* hd_BetaMom_Tp_kal_allTofs;
  TH2F* hd_BetaMom_Tp_ini_allTofs_Z2;
  TH2F* hd_BetaMom_Bg_all;
  TH2F* hd_BetaMom_Bg_allTofs;
  TH2F* hd_BetaMom_Bg_ini_all;
  TH2F* hd_BetaMom_Bg_ini_allTofs;
  TH1F* hd_TimeCalTp[32];
  TH2F *hd_TrackFrom[2];
  TH2F *hd_Momini_rchiy_Tofs;

  TH2F *hd_PoQ_E;
  TH2F *hd_PoQ_E_ini;
  TH2F *hd_PoQ_ETofs;
  TH2F *hd_PoQ_E_iniTofs;
  
  TH2F *hd_Mom_ini_vs_rchi_ini;

  //   TH2F *hd_PoQ_E3;
  //   TH2F *hd_PoQ_E_ini3;
  //   TH2F *hd_PoQ_ETofs3;
  //   TH2F *hd_PoQ_E_iniTofs3;
  
  TH2F *hd_IniMom_DiffMom_allTofs;
  TH2F *hd_NewMom_DiffMom_allTofs;
  TH2F *hd_IniMom_NewMom_all;
  TH2F *hd_IniMom_NewMom_all_tofs;
  TH2F *hd_IniMom_NewMom_all_tofs_ini;
  TH2F *hd_IniMom_NewMom_all_tofsZ1;
  TH2F *hd_IniMom_NewMom_all_tofsZ1_ini;
  TH2F *hd_Bar_DiffMom_allTofs;

  TH2F *hd_IniNewMomdiff_vs_IniNewMomplus_allTofs;
  TH2F *hd_charge_ini_vs_IniNewDiv;
  TH2F *hd_charge_vs_IniNewDiv;
  TH2F *hd_IniNewMomdiff_vs_IniNewMomplus_allTofs_Z1;
  TH2F *hd_IniNewMomdiff_vs_IniNewMomplus_allTofs_ini_Z1;
  TH2F *hd_best_pos_hitvstrack_tofpZ1;
  TH1F *hd_best_pos_hitdifftrack_tofpZ1;


  //TH2F *hd_IniMom_DiffMom_all;
  //TH2F *hd_NewMom_DiffMom_all;
  //TH2F *hd_IniMom_NewMom_all;
  //TH2F *hd_Bar_DiffMom_all;

  //#endif


  TH2F* h_MC_PxPy;
  TH1F* h_MC_Pt;
  TH2F* h_MC_Pt_y;
  TH2F* h_MC_TargetXZ;
  TH2F* h_MC_TargetYZ;
  TH1F* h_MC_stats;
  
  TH2F* h_MC_dE_mom;
  TH2F* h_MC_time_mom;
  TH2F* h_MC_beta_mom;
  TH2F* h_MC_length;
  TH2F* h_MC_length_Z;

  std::vector<TH2F*> Material_XX0_y;
  std::vector<TH2F*> Material_dE_y;
  

  Ana_Hist(bool Daf=true, bool Vertex=true, bool DCproject=true, bool Finding=true, bool Hough=true,bool Simu = false);
  ~Ana_Hist();

  int Write(TFile*);
  int WriteTemp(char *tempfile);
};

#endif
