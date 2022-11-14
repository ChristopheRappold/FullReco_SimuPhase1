#ifndef PARA_MANAGER_HH
#define PARA_MANAGER_HH

#include <string>
#include <array>
#include <unordered_map>

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

  //  Cut  //////
  double cut_chi2_gf;
  double cut_phi_fm;
  int    cut_num_mdc;
  double cut_psb_phi;
  double cut_psb_z;
  double cut_psfe_phi;
  double cut_psbe_phi;

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

  double fiber_mft_cut_d ;
  double fiber_uft1_cut_d;
  double fiber_uft2_cut_d;
  double fiber_uft3_cut_d;
  double fiber_dft1_cut_d;
  double fiber_dft2_cut_d;

  double fiber_tgt_pos_x;   double fiber_tgt_pos_y;   double fiber_tgt_pos_z;

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

  double fiber_mft1_off_ang_x1; double fiber_mft1_off_ang_u1; double fiber_mft1_off_ang_v1;
  double fiber_mft1_off_ang_x2; double fiber_mft1_off_ang_u2; double fiber_mft1_off_ang_v2;
  double fiber_mft2_off_ang_x1; double fiber_mft2_off_ang_u1; double fiber_mft2_off_ang_v1;
  double fiber_mft2_off_ang_x2; double fiber_mft2_off_ang_u2; double fiber_mft2_off_ang_v2;

  double fiber_res;


  //  PSB  //////
  double psb_tcut_min; double psb_tcut_max;
  double psb_pos_z;
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
  double mdc_t0_off;
  std::array<std::array<InfoTDC, 8>, 16> mdc_MappingMDC;
  std::array<InfoMDC, 17> mdc_PhysMDC;
  std::string mdc_name_map;
  std::string mdc_name_phys;
  std::string mdc_name_par;

  double mdc_res;

  bool InitMDCParameter(const std::unordered_map<std::string,std::string>& ParamFilesMDC);
  bool mdc_init_done = false;

  std::unordered_map<std::string,std::string>::const_iterator itr_ParamFiles;

};


#endif
