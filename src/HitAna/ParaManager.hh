#ifndef PARA_MANAGER_HH
#define PARA_MANAGER_HH

#include <string>
#include <cstring>
#include <array>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"


struct InfoTDC
{
  int cmp_id = -1;
  int cmp_block = -1;
  int mdc_layer = -1;
  int mdc_sector = -1;
  int wire_id_begin = -1;
  int wire_id_end = -1;
  int wire_nb = -1;
  int wire_real_start = -1;
  int wire_real_end = -1;
};

struct InfoMDC
{
  float Dpsi1 = 0.;
  float Dpsi0b = 0.;
  float Dpsi0f = 0.;
  float dpsi0 = 0.;
  float rmax = -1.;
  int nb_wire = -1;
  float wire_size = -1.;
  float z_back = 0.;
  float z_front = 0.;
  float l_back = 0.;
  float l_front = 0.;
  float drift_par0 = 0.;
  float drift_par1 = 0.;
  float drift_par2 = 0.;
  float drift_par3 = 0.;
  float drift_par4 = 0.;
  float drift_par5 = 0.;
  float drift_par6 = 0.;
};

class ParaManager
{

  public:
  ParaManager(const std::unordered_map<std::string,std::string>& ParamFiles);
  ~ParaManager();
  //static ParaManager* Instance();

  //  Flag  //////
  bool write_hist;
  bool write_hist_all;
  bool write_canvas;
  bool flag_mft12;
  bool flag_debug;
  bool flag_dft12_cut;
  bool flag_dft12_combi;
  bool flag_dup_xuv_mft12;
  bool flag_dup_mft12_combi;
  bool flag_mft12_inclusive;
  bool flag_dup_xuv_uft3;
  bool flag_mft12_allcombi;
  bool flag_mft12_posang;
  bool flag_mft12_combi;
  bool flag_mft12_pair;
  bool flag_mft12_xuv_psb;
  bool flag_dup_mft12_xuv;

  //  Trig  //////
  bool trig_main;
  bool trig_clock;
  bool trig_t0;
  bool trig_sc41;

  //  Cut  //////
  double cut_chi2_gf;
  double cut_chi2_mft12;
  double cut_chi2_dft12;
  double cut_chi2_uft12;
  int    cut_num_mdc;
  double cut_psfe_phi;
  double cut_psbe_phi;
  double cut_t0_time;
  double cut_mft12_combi;
  double cut_dft12_combi;
  double cut_mft12_hit;
  double cut_dft12_tot_mean;
  double cut_dft12_tot_sig;
  double cut_dft12_tot_max;
  double cut_dft12_time_mean;
  double cut_dft12_time_sig;

  //  Field  //////
  bool   field_flag;
  double field;

  //  WASA  //////
  double wasa_tof_offset;

  //  Fiber  //////
  double fiber_mft_cut_min  ;  double fiber_mft_cut_max  ; // SFP0
  double fiber_uft3_cut_min ;  double fiber_uft3_cut_max ; // SFP1
  double fiber_uft12_cut_min;  double fiber_uft12_cut_max; // SFP2
  double fiber_dft12_cut_min;  double fiber_dft12_cut_max; // SFP3

  double fiber_tot_mft_cut_min  ;  double fiber_tot_mft_cut_max  ;
  double fiber_tot_uft3_cut_min ;  double fiber_tot_uft3_cut_max ;
  double fiber_tot_uft12_cut_min;  double fiber_tot_uft12_cut_max;
  double fiber_tot_dft12_cut_min;  double fiber_tot_dft12_cut_max;

  double fiber_time_mft_cut_min  ;  double fiber_time_mft_cut_max  ;
  double fiber_time_uft3_cut_min ;  double fiber_time_uft3_cut_max ;
  double fiber_time_uft12_cut_min;  double fiber_time_uft12_cut_max;
  double fiber_time_dft12_cut_min;  double fiber_time_dft12_cut_max;
  
  double fiber_mft_cut_d ;
  double fiber_uft1_cut_d;
  double fiber_uft2_cut_d;
  double fiber_uft3_cut_d;
  double fiber_dft1_cut_d;
  double fiber_dft2_cut_d;

  double fiber_ch2ns;
  double fiber_res;

  double fiber_tgt_pos_x;   double fiber_tgt_pos_y;   double fiber_tgt_pos_z; //CHECK possibly duplicated
  double fiber_tgt_size_x;  double fiber_tgt_size_y;  double fiber_tgt_size_z;
  double target_pos_x;      double target_pos_y;      double target_pos_z;

  double fiber_mft1_pos_x;  double fiber_mft1_pos_y;  double fiber_mft1_pos_z;
  double fiber_mft2_pos_x;  double fiber_mft2_pos_y;  double fiber_mft2_pos_z;
  double fiber_uft1_pos_x;  double fiber_uft1_pos_y;  double fiber_uft1_pos_z;
  double fiber_uft2_pos_x;  double fiber_uft2_pos_y;  double fiber_uft2_pos_z;
  double fiber_uft3_pos_x;  double fiber_uft3_pos_y;  double fiber_uft3_pos_z;
  double fiber_dft1_pos_x;  double fiber_dft1_pos_y;  double fiber_dft1_pos_z;
  double fiber_dft2_pos_x;  double fiber_dft2_pos_y;  double fiber_dft2_pos_z;

  double fiber_mft1_off_x1; double fiber_mft1_off_u1; double fiber_mft1_off_v1;
  double fiber_mft1_off_x2; double fiber_mft1_off_u2; double fiber_mft1_off_v2;
  double fiber_mft2_off_x1; double fiber_mft2_off_u1; double fiber_mft2_off_v1;
  double fiber_mft2_off_x2; double fiber_mft2_off_u2; double fiber_mft2_off_v2;
  double fiber_uft1_off_x;  double fiber_uft1_off_u;  double fiber_uft1_off_v;
  double fiber_uft2_off_x;  double fiber_uft2_off_u;  double fiber_uft2_off_v;
  double fiber_uft3_off_x;  double fiber_uft3_off_u;  double fiber_uft3_off_v;
  double fiber_dft1_off_x;  double fiber_dft1_off_u;  double fiber_dft1_off_v;
  double fiber_dft2_off_x;  double fiber_dft2_off_u;  double fiber_dft2_off_v;

  double fiber_mft1_step_x1; double fiber_mft1_step_u1; double fiber_mft1_step_v1;
  double fiber_mft1_step_x2; double fiber_mft1_step_u2; double fiber_mft1_step_v2;
  double fiber_mft2_step_x1; double fiber_mft2_step_u1; double fiber_mft2_step_v1;
  double fiber_mft2_step_x2; double fiber_mft2_step_u2; double fiber_mft2_step_v2;


  std::string fiber_name_offset;
  std::array<std::array<std::array<double, 384>, 3>, 7> fiber_offset;
  std::string fiber_name_timeoffset;
  std::array<std::array<std::array<double, 384>, 3>, 7> fiber_time_offset;
  std::string fiber_name_angleoffset;
  std::array<std::array<std::array<double, 2>, 3>, 7> fiber_angle_offset;

  double fiber_mft1_off_ang_x1; double fiber_mft1_off_ang_u1; double fiber_mft1_off_ang_v1;
  double fiber_mft1_off_ang_x2; double fiber_mft1_off_ang_u2; double fiber_mft1_off_ang_v2;
  double fiber_mft2_off_ang_x1; double fiber_mft2_off_ang_u1; double fiber_mft2_off_ang_v1;
  double fiber_mft2_off_ang_x2; double fiber_mft2_off_ang_u2; double fiber_mft2_off_ang_v2;

  std::string fiber_name_mftcor;
  std::array<std::array<std::array<std::array<double, 3>, 2>, 3>, 2> fiber_mft_cor_par;


  //  PSB  //////
  double psb_tcut_min; double psb_tcut_max;
  double psb_pos_x;  double psb_pos_y;  double psb_pos_z;
  double psb_res_phi;
  double psb_res_z;
  std::string psb_name_time;
  double psb_off_time[46][2];
  double psb_ch2ns[46][2];
  double psb_zpar[46];

  //  PSFE  //////
  double psfe_tcut_min; double psfe_tcut_max;
  double psfe_rmin;     double psfe_rmax;
  double psfe_pos_z;
  double psfe_res_phi;
  double psfe_ch2ns;

  //  PSBE  //////
  double psbe_tcut_min; double psbe_tcut_max;
  double psbe_rmin;     double psbe_rmax;
  double psbe_pos_z;

  //  T0  //////
  double t0_tcut_min; double t0_tcut_max;
  double t0_pos_z;
  std::string t0_name_time;
  double t0_off_time[28][2];
  double t0_ch2ns[28][2];


  //  MDC  //////
  double mdc_cut_min;
  double mdc_cut_max1, mdc_cut_max2, mdc_cut_max3;
  double mdc_cut_tot;

  double mdc_pos_x;
  double mdc_pos_y;
  double mdc_pos_z;

  double mdc_rot_x;
  double mdc_rot_y;
  double mdc_rot_z;

  double mdc_ch2ns;
  std::array<std::array<InfoTDC, 8>, 16> mdc_MappingMDC;
  std::array<InfoMDC, 17> mdc_PhysMDC;
  std::string mdc_name_map;
  std::string mdc_name_phys;

  std::string mdc_name_drift;

  double mdc_t0_off[17];
  double mdc_t0_off_wir[17][148];
  std::string mdc_name_t0;
  std::string mdc_name_t0_wir;

  double mdc_res;

  bool InitMDCParameter(const std::unordered_map<std::string,std::string>& ParamFilesMDC);
  bool mdc_init_done = false;


  //  MWDC  //////
  std::string mwdc_name_dtdx;
  char mwdc_name_dtdxtable[256];

  void SetMWDCdtdxFromTH1D(char*, float**, int*, int*);
  float dist_focS4;     /* All distances from    */
  float dist_MWDC_zref; /*  MWDC reference point   */
  float dist_s4test;    /*  for testing something. */
  float dist_SC41;      /*        in mm          */
  float dist_SC42;      /*        in mm          */
  float dist_SC43;      /*        in mm          */

  int mwdc_lt_valid_min;
  int mwdc_lt_valid_max;
  int id_plane[8][8];     //[id_board][id_group(16ch)]  to  id_plane
  int id_wiregroup[8][8]; //[id_board][id_group(16ch)]  to  id_plane
  char plane_name[16][32];
  int participate_tracking_plane[16];
  float mwdc_zpos[16], mwdc_center_id[16], mwdc_plane_angle[16], mwdc_plane_sign[16];
  float mwdc_assumed_plane_resolution;
  float mwdc_tracking_chi2cut_max;
  int mwdc_max_hit_combination_for_tracking;
  int mwdc_min_plane_with_hit_for_tracking;
  float mwdc_dtdxtable[16][1000]; // for our MWDC41/42, ~300 is enough. 1000 is just in case.
  int mwdc_dtdxtable_dtmin[16];
  int mwdc_dtdxtable_dtmax[16];
  float mwdc_shift_x_alignment; // shift due to detector misalignment at S4
  float mwdc_shift_a_alignment; // shift due to detector misalignment at S4
  float mwdc_shift_y_alignment; // shift due to detector misalignment at S4
  float mwdc_shift_b_alignment; // shift due to detector misalignment at S4
  float cut_a_s4;
  float mwdc_dtdxconversion(int, int); // i_plane, tdc_val(int)
  int mwdc_tdccut_timing_threshold[2]; // used in TFRSUserProc.cxx [0]:begin, [1]: end

  //// mwdc drift time to drift length convert parameter, later fit for each run with Fermi function
  double mwdc_dtdx_par[16][8];
  double mwdc_t0_off[16];
  double mwdc_tmax[16];
  double mwdc_ch2ns;

  /// --- S3, S4 Scintillators -- //
  double mtdc_ch2ns;
  int tcut_sc31_min;
  int tcut_sc31_max;
  int tcut_sc41_1_min;
  int tcut_sc41_1_max;
  int tcut_sc41_2_min;
  int tcut_sc41_2_max;
  int tcut_sc41_3_min;
  int tcut_sc41_3_max;
  int tcut_sc42_min;
  int tcut_sc42_max;
  int tcut_sc43_1_min;
  int tcut_sc43_1_max;
  int tcut_sc43_2_min;
  int tcut_sc43_2_max;
  int tcut_sc43_3_min;
  int tcut_sc43_3_max;
  int tcut_sc43_4_min;
  int tcut_sc43_4_max;
  double offset_tof_sc3141;
  double offset_tof_sc4143;


  std::unordered_map<std::string,std::string>::const_iterator itr_ParamFiles;

};


#endif
