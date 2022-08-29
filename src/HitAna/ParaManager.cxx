#include "ParaManager.hh"

ParaManager *fixInstance=0;

ParaManager::ParaManager(){

  //  Flag  //////
  write_hist = true;
  write_hist_all = true;
  write_canvas = true;
  flag_mft12 = true;

  //  Cut  //////
  cut_chi2_gf = 50;
  cut_phi_fm  = 0.3;
  cut_num_mdc = 3;
  cut_psb_phi = 0.6;
  cut_psb_z   = 120;
  cut_psfe_phi = 0.6;
  cut_psbe_phi = 0.6;

  //  Field  //////
  field_flag = true; // true : measuremnt field , false : constant field
  field = 1; // Tesla

  //  WASA  //////
  wasa_tof_offset = 0;

  //  Fiber  //////
  double s2fiber_sfp0_mean = -1510; // mft12
  double s2fiber_sfp1_mean = -1590; // uft3
  double s2fiber_sfp2_mean = -1590; // uft12
  double s2fiber_sfp3_mean = -1610; // dft12
  double s2fiber_cut = 50;

  fiber_mft_cut_min   = s2fiber_sfp0_mean - s2fiber_cut;  fiber_mft_cut_max   = s2fiber_sfp0_mean + s2fiber_cut;; // SFP0
  fiber_uft3_cut_min  = s2fiber_sfp1_mean - s2fiber_cut;  fiber_uft3_cut_max  = s2fiber_sfp1_mean + s2fiber_cut;; // SFP1
  fiber_uft12_cut_min = s2fiber_sfp2_mean - s2fiber_cut;  fiber_uft12_cut_max = s2fiber_sfp2_mean + s2fiber_cut;; // SFP2
  fiber_dft12_cut_min = s2fiber_sfp3_mean - s2fiber_cut;  fiber_dft12_cut_max = s2fiber_sfp3_mean + s2fiber_cut;; // SFP3

  fiber_mft_cut_d  = 4.;
  fiber_uft1_cut_d = 4.;
  fiber_uft2_cut_d = 4.;
  fiber_uft3_cut_d = 20.;
  fiber_dft1_cut_d = 6.;
  fiber_dft2_cut_d = 6.;

  fiber_tgt_pos_x = 0.;  fiber_tgt_pos_y  = 0.;   fiber_tgt_pos_z  = 1961.2;

  fiber_mft1_pos_x = 0;  fiber_mft1_pos_y =   0;  fiber_mft1_pos_z = 2269.3;
  fiber_mft2_pos_x = 0;  fiber_mft2_pos_y =   0;  fiber_mft2_pos_z = 2309.3;
  fiber_uft1_pos_x = 0;  fiber_uft1_pos_y =   0;  fiber_uft1_pos_z = 1547.7;
  fiber_uft2_pos_x = 0;  fiber_uft2_pos_y =   0;  fiber_uft2_pos_z = 1711.7;
  fiber_uft3_pos_x = 0;  fiber_uft3_pos_y =   0;  fiber_uft3_pos_z = 1996.7;
  fiber_dft1_pos_x = 0;  fiber_dft1_pos_y =  12;  fiber_dft1_pos_z = 3960.0;
  fiber_dft2_pos_x = 0;  fiber_dft2_pos_y = -12;  fiber_dft2_pos_z = 4056.3;

  fiber_mft1_off_x1 = 0; fiber_mft1_off_u1 = 0; fiber_mft1_off_v1 = 0;
  fiber_mft1_off_x2 = 0; fiber_mft1_off_u2 = 0; fiber_mft1_off_v2 = 0;
  fiber_mft2_off_x1 = 0; fiber_mft2_off_u1 = 0; fiber_mft2_off_v1 = 0;
  fiber_mft2_off_x2 = 0; fiber_mft2_off_u2 = 0; fiber_mft2_off_v2 = 0;
  fiber_uft1_off_x  = 0; fiber_uft1_off_u  = 0; fiber_uft1_off_v  = 0;
  fiber_uft2_off_x  = 0; fiber_uft2_off_u  = 0; fiber_uft2_off_v  = 0;
  fiber_uft3_off_x  = 0; fiber_uft3_off_u  = 0; fiber_uft3_off_v  = 0;
  fiber_dft1_off_x  = 0; fiber_dft1_off_u  = 0; fiber_dft1_off_v  = 0;
  fiber_dft2_off_x  = 0; fiber_dft2_off_u  = 0; fiber_dft2_off_v  = 0;

  fiber_mft1_step_x1 = 1.1; fiber_mft1_step_u1 = 1.1; fiber_mft1_step_v1 = 1.1;
  fiber_mft1_step_x2 = 1.1; fiber_mft1_step_u2 = 1.1; fiber_mft1_step_v2 = 1.1;
  fiber_mft2_step_x1 = 1.1; fiber_mft2_step_u1 = 1.1; fiber_mft2_step_v1 = 1.1;
  fiber_mft2_step_x2 = 1.1; fiber_mft2_step_u2 = 1.1; fiber_mft2_step_v2 = 1.1;

  fiber_name_offset  = "setup/fiber_offset/fiber_offset.csv";
  for(int i=0; i<7; ++i){
    for(int j=0; j<3; ++j){
      for(int k=0; k<384; ++k){
        fiber_offset[i][j][k] = 0;
      }
    }
  }

  fiber_mft1_off_ang_x1 =  0.0; fiber_mft1_off_ang_u1 =  0.0; fiber_mft1_off_ang_v1 =  0.0;
  fiber_mft1_off_ang_x2 =  0.0; fiber_mft1_off_ang_u2 =  0.0; fiber_mft1_off_ang_v2 =  0.0;
  fiber_mft2_off_ang_x1 =  0.0; fiber_mft2_off_ang_u1 =  0.0; fiber_mft2_off_ang_v1 =  0.0;
  fiber_mft2_off_ang_x2 =  0.0; fiber_mft2_off_ang_u2 =  0.0; fiber_mft2_off_ang_v2 =  0.0;

  fiber_res = 0.15;

  //  PSB  //////
  psb_tcut_min = -20200; psb_tcut_max = -19400;
  psb_pos_z = 2760.;
  psb_res_phi = 11; //mm
  psb_res_z   = 10; //mm

  psb_name_time  = "setup/psb_param/psb_time.csv";
  for(int i=0; i<46; ++i){
    psb_zpar[i] = 70;
    for(int j=0; j<2; ++j){
      psb_off_time[i][j] = 0.;
      psb_ch2ns[i][j]    = 0.025;
    }
  }

  //  PSFE  //////
  psfe_tcut_min = -202000; psfe_tcut_max = 194000;
  psfe_rmin = 85.;         psfe_rmax = 187.5;
  psfe_pos_z = 3040.;
  psfe_res_phi  = 10; //mm
  psfe_ch2ns  = 0.025;

  //  PSBE  //////
  psbe_tcut_min = -202000; psbe_tcut_max = 194000;
  psbe_rmin = 117.;         psbe_rmax = 187.5;
  psbe_pos_z = 2480.;

  //  T0  //////
  t0_tcut_min = -200000 ; t0_tcut_max = 20000;
  t0_pos_z = 100; //tmp
  psb_name_time  = "setup/t0_param/t0_time.csv";
  for(int i=0; i<28; ++i){
    for(int j=0; j<2; ++j){
      t0_off_time[i][j] = 0.;
      t0_ch2ns[i][j]    = 0.025;
    }
  }


  //  MDC  //////
  mdc_cut_min  = -1350;
  mdc_cut_max1 = -1150; mdc_cut_max2  = -1080; mdc_cut_max3  = -1000;
  mdc_cut_tot =   180;

  mdc_ch2ns  = 0.4;
  mdc_t0_off = -530;

  mdc_pos_x = 0.;
  mdc_pos_y = 0.;
  mdc_pos_z = 2760.;

  mdc_rot_x = 0.;
  mdc_rot_y = 0.;
  mdc_rot_z = 0.;

  mdc_res = 0.2;

  InitMDCParameter();


}

ParaManager::~ParaManager() {}

ParaManager* ParaManager::Instance(){

  if(fixInstance==0) {
    fixInstance = new ParaManager();
  }
  return fixInstance;

};


bool ParaManager::InitMDCParameter()
{
  if(mdc_init_done)
    return true;

  mdc_name_map  = "./mapping/MDC_channelmap.csv";
  mdc_name_phys = "./mapping/MDC_PhysicalMap.csv";
  mdc_name_par  = "./mdc_driftparam/MDC_DriftParam.txt";

  for(int ctdc_id=0; ctdc_id<16; ++ctdc_id){
    for(int ctdc_block=0; ctdc_block<8; ++ctdc_block){
      mdc_MappingMDC[ctdc_id][ctdc_block].cmp_id = -1;
      mdc_MappingMDC[ctdc_id][ctdc_block].cmp_block = -1;
      mdc_MappingMDC[ctdc_id][ctdc_block].mdc_layer = -1;
      mdc_MappingMDC[ctdc_id][ctdc_block].mdc_sector = -1;
      mdc_MappingMDC[ctdc_id][ctdc_block].wire_id_begin = -1;
      mdc_MappingMDC[ctdc_id][ctdc_block].wire_id_end = -1;
      mdc_MappingMDC[ctdc_id][ctdc_block].wire_nb = -1;
    }
  }

  for(int layerID=0; layerID<17; ++layerID){

    mdc_PhysMDC[layerID].Dpsi1 = 0;
    mdc_PhysMDC[layerID].Dpsi0b = 0;
    mdc_PhysMDC[layerID].Dpsi0f = 0;
    mdc_PhysMDC[layerID].dpsi0 = 0;
    mdc_PhysMDC[layerID].rmax = 0;
    mdc_PhysMDC[layerID].wire_size = 0;
    mdc_PhysMDC[layerID].z_back = 0;
    mdc_PhysMDC[layerID].z_front = 0;
    mdc_PhysMDC[layerID].l_back = 0;
    mdc_PhysMDC[layerID].l_front = 0;
    mdc_PhysMDC[layerID].drift_par0 = 0;
    mdc_PhysMDC[layerID].drift_par1 = 0;
    mdc_PhysMDC[layerID].drift_par2 = 0;
    mdc_PhysMDC[layerID].drift_par3 = 0;
    mdc_PhysMDC[layerID].drift_par4 = 0;
    mdc_PhysMDC[layerID].drift_par5 = 0;
    mdc_PhysMDC[layerID].drift_par6 = 0;

  }
  mdc_PhysMDC[0].nb_wire = 52;
  mdc_PhysMDC[1].nb_wire = 64;
  mdc_PhysMDC[2].nb_wire = 76;
  mdc_PhysMDC[3].nb_wire = 88;
  mdc_PhysMDC[4].nb_wire = 100;
  mdc_PhysMDC[5].nb_wire = 76;
  mdc_PhysMDC[6].nb_wire = 86;
  mdc_PhysMDC[7].nb_wire = 96;
  mdc_PhysMDC[8].nb_wire = 106;
  mdc_PhysMDC[9].nb_wire = 116;
  mdc_PhysMDC[10].nb_wire = 126;
  mdc_PhysMDC[11].nb_wire = 104;
  mdc_PhysMDC[12].nb_wire = 112;
  mdc_PhysMDC[13].nb_wire = 120;
  mdc_PhysMDC[14].nb_wire = 130;
  mdc_PhysMDC[15].nb_wire = 138;
  mdc_PhysMDC[16].nb_wire = 148;

  for(size_t i = 0; i< mdc_MappingMDC.size() ; ++i){
    for(size_t j = 0 ;j < mdc_MappingMDC[i].size(); ++j){
      mdc_MappingMDC[i][j].wire_real_start = -1;
      mdc_MappingMDC[i][j].wire_real_end = -1;
    }
  }

  return true;
}
