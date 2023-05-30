#include "ParaManager.hh"

//ParaManager *fixInstance=0;

ParaManager::ParaManager(const std::unordered_map<std::string,std::string>& ParamFiles){

  //  Flag  //////
  write_hist = true;
  write_hist_all = true;
  write_canvas = true;
  flag_mft12 = true;
  flag_debug = false;
  flag_dft12_cut = true;
  flag_dft12_combi = true;
  flag_dup_xuv_mft12    = false;
  flag_dup_mft12_combi = false;
  flag_mft12_inclusive = true;
  flag_dup_xuv_uft3     = true;
  flag_mft12_allcombi = false;
  flag_mft12_posang = true;
  flag_mft12_combi = true;
  flag_mft12_pair = true;
  flag_mft12_xuv_psb = false;
  flag_dup_mft12_xuv = false;
  flag_uft12_combi = false;

  //  Trig  //////
  trig_main  = false;
  trig_clock = false;
  trig_t0    = false;
  trig_sc41  = false;

  //  Cut  //////
  cut_chi2_gf = 50;
  cut_chi2_dft12 = 10.;
  cut_chi2_mft12 = 10.;
  cut_chi2_uft12 = 10.;
  cut_num_mdc = 3;
  cut_psfe_phi = 0.6;
  cut_psbe_phi = 0.6;
  cut_mft12_combi = 1e5;
  cut_dft12_combi = 1e5;
  cut_dft12_tot_mean  = 88.8;
  cut_dft12_tot_sig   =  3.2;
  cut_dft12_tot_max   = 75.;
  cut_dft12_time_mean = 0.04;
  cut_dft12_time_sig  = 1.09;

  //  Field  //////
  field_flag = true; // true : measurement field , false : constant field
  field = 1; // Tesla

  //  WASA  //////
  wasa_tof_offset = 3.402;

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

  fiber_tot_mft_cut_min   = -1;  fiber_tot_mft_cut_max   = 1e+10;
  fiber_tot_uft3_cut_min  = -1;  fiber_tot_uft3_cut_max  = 80.;
  fiber_tot_uft12_cut_min = -1;  fiber_tot_uft12_cut_max = 1e+10;
  fiber_tot_dft12_cut_min = -1;  fiber_tot_dft12_cut_max = 1e+10;

  fiber_time_uft12_cut_min = -5.;  fiber_time_uft12_cut_max = 8.;
  fiber_time_uft3_cut_min  = -5.;  fiber_time_uft3_cut_max  = 10.;
  fiber_time_mft_cut_min   = -5.;  fiber_time_mft_cut_max   = 10.;
  fiber_time_dft12_cut_min = -5.;  fiber_time_dft12_cut_max = 7.;

  fiber_mft_cut_d  = 4.;
  fiber_uft1_cut_d = 4.;
  fiber_uft2_cut_d = 4.;
  fiber_uft3_cut_d = 8.;
  fiber_dft1_cut_d = 6.;
  fiber_dft2_cut_d = 6.;

  fiber_ch2ns = 0.416666;
  fiber_res = 0.2;

  double buf_mft1 = -6.0;
  double buf_mft2 = -6.0;
  double buf_mft_dist =  0.4;

  fiber_tgt_pos_x = 0.;   fiber_tgt_pos_y  = 0.;   fiber_tgt_pos_z  = 1961.2;
  fiber_tgt_size_x = 30.; fiber_tgt_size_y  = 30.; fiber_tgt_size_z  = 30.;

  fiber_mft1_pos_x = 2.1;fiber_mft1_pos_y = 2.5;  fiber_mft1_pos_z = 2269.3 + buf_mft1 - buf_mft_dist/2.;
  fiber_mft2_pos_x = 2.1;fiber_mft2_pos_y = 2.5;  fiber_mft2_pos_z = 2309.3 + buf_mft2 + buf_mft_dist/2.;
  fiber_uft1_pos_x = 0;  fiber_uft1_pos_y =   0;  fiber_uft1_pos_z = 1547.7;
  fiber_uft2_pos_x = 0;  fiber_uft2_pos_y =   0;  fiber_uft2_pos_z = 1711.7;
  fiber_uft3_pos_x = 0;  fiber_uft3_pos_y =   0;  fiber_uft3_pos_z = 1996.7;
  fiber_dft1_pos_x = 0;  fiber_dft1_pos_y =  12;  fiber_dft1_pos_z = 3960.0;
  fiber_dft2_pos_x = 0;  fiber_dft2_pos_y = -12;  fiber_dft2_pos_z = 4056.3;

  fiber_mft1_off_x1 =  0.93 + 0.7 ;  fiber_mft1_off_u1 = 1.49 - 0.45;         fiber_mft1_off_v1 =  0.76 + 0.05 + 0.66;
  fiber_mft1_off_x2 = -0.22 - 0.05;  fiber_mft1_off_u2 = 0.25 - 0.05 - 0.17;  fiber_mft1_off_v2 = -0.52 - 0.29;
  fiber_mft2_off_x1 =  0.70 + 0.93;  fiber_mft2_off_v1 = 0.61 + 0.8;          fiber_mft2_off_u1 =  1.62 - 0.55;
  fiber_mft2_off_x2 = -0.13 - 0.15;  fiber_mft2_off_v2 = 0.10 - 0.50;         fiber_mft2_off_u2 = -0.08 - 0.20 - 0.08;

  fiber_uft1_off_x  =  0.294; fiber_uft1_off_u  =  0.102; fiber_uft1_off_v  =  0.273;
  fiber_uft2_off_x  = -0.258; fiber_uft2_off_u  = -0.100; fiber_uft2_off_v  = -0.247;
  fiber_uft3_off_x  = -1.194; fiber_uft3_off_u  = -0.422; fiber_uft3_off_v  = -0.887;
  fiber_dft1_off_x  =  0.040; fiber_dft1_off_u  = -0.287; fiber_dft1_off_v  =  0.892;
  fiber_dft2_off_x  =  0.183; fiber_dft2_off_u  =  0.436; fiber_dft2_off_v  = -0.377;

  fiber_mft1_step_x1 = 1.1; fiber_mft1_step_u1 = 1.1; fiber_mft1_step_v1 = 1.1;
  fiber_mft1_step_x2 = 1.1; fiber_mft1_step_u2 = 1.1; fiber_mft1_step_v2 = 1.1;
  fiber_mft2_step_x1 = 1.1; fiber_mft2_step_u1 = 1.1; fiber_mft2_step_v1 = 1.1;
  fiber_mft2_step_x2 = 1.1; fiber_mft2_step_u2 = 1.1; fiber_mft2_step_v2 = 1.1;

  itr_ParamFiles = ParamFiles.find("fiber_offset_file");
  if(itr_ParamFiles != ParamFiles.end())
    fiber_name_offset  = itr_ParamFiles->second;

  for(int i=0; i<7; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<384; ++k)
        fiber_offset[i][j][k] = 0.;


  std::ifstream ifs_fiber ( fiber_name_offset );
  if(ifs_fiber.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs_fiber,temp_line))
        {
          std::stringstream stream(temp_line);
          std::string testComment(stream.str());
          std::size_t it_comment = testComment.find(CommentSymbol);
          if(it_comment!=std::string::npos)
            {
              //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
              continue;
            }

          int det_id, lay_id, fib_id;
          double offset;

          stream >> det_id >> lay_id >> fib_id >> offset;
          //printf("%d %d %d : %.2f\n", det_id, lay_id, fib_id, offset );

          fiber_offset[det_id][lay_id][fib_id] = offset;
        }
      //std::cout << "done " << mdc_name_map << std::endl;
      printf("fiber offset loaded : %s\n", fiber_name_offset.c_str());
    }
  else
    {
      //std::cout << " ! fail to open " << mdc_name_map << std::endl;
      printf(" ! fail to open  : %s\n", fiber_name_offset.c_str());
      exit(-1);
    }


  itr_ParamFiles = ParamFiles.find("fiber_timeoffset_file");
  if(itr_ParamFiles != ParamFiles.end())
    fiber_name_timeoffset  = itr_ParamFiles->second;

  for(int i=0; i<7; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<384; ++k)
        fiber_time_offset[i][j][k] = 0.;


  std::ifstream ifs_timefiber ( fiber_name_timeoffset );
  if(ifs_timefiber.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs_timefiber,temp_line))
        {
          std::stringstream stream(temp_line);
          std::string testComment(stream.str());
          std::size_t it_comment = testComment.find(CommentSymbol);
          if(it_comment!=std::string::npos)
            {
              //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
              continue;
            }

          int det_id, lay_id, fib_id;
          double offset;

          stream >> det_id >> lay_id >> fib_id >> offset;
          //printf("%d %d %d : %.2f\n", det_id, lay_id, fib_id, offset );

          fiber_time_offset[det_id][lay_id][fib_id] = offset;
        }
      //std::cout << "done " << mdc_name_map << std::endl;
      printf("fiber time offset loaded : %s\n", fiber_name_timeoffset.c_str());
    }
  else
    {
      //std::cout << " ! fail to open " << mdc_name_map << std::endl;
      printf(" ! fail to open  : %s\n", fiber_name_timeoffset.c_str());
      exit(-1);
    }



  itr_ParamFiles = ParamFiles.find("fiber_angleoffset_file");
  if(itr_ParamFiles != ParamFiles.end())
    fiber_name_angleoffset  = itr_ParamFiles->second;

  for(int i=0; i<7; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<2; ++k)
        fiber_angle_offset[i][j][k] = 0.;


  std::ifstream ifs_anglefiber ( fiber_name_angleoffset );
  if(ifs_anglefiber.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs_anglefiber,temp_line))
        {
          std::stringstream stream(temp_line);
          std::string testComment(stream.str());
          std::size_t it_comment = testComment.find(CommentSymbol);
          if(it_comment!=std::string::npos)
            {
              //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
              continue;
            }

          int det_id, lay_id, seg_id;
          double angle;

          stream >> det_id >> lay_id >> seg_id >> angle;
          //printf("%d %d %d : %.2f\n", det_id, lay_id, fib_id, offset );

          fiber_angle_offset[det_id][lay_id][seg_id] = angle;
        }
      //std::cout << "done " << mdc_name_map << std::endl;
      printf("fiber angle offset loaded : %s\n", fiber_name_angleoffset.c_str());
    }
  else
    {
      //std::cout << " ! fail to open " << mdc_name_map << std::endl;
      printf(" ! fail to open  : %s\n", fiber_name_angleoffset.c_str());
      exit(-1);
    }

  fiber_mft1_off_ang_x1 = fiber_angle_offset[3][0][0];  fiber_mft1_off_ang_x2 = fiber_angle_offset[3][0][1];
  fiber_mft1_off_ang_u1 = fiber_angle_offset[3][1][0];  fiber_mft1_off_ang_u2 = fiber_angle_offset[3][1][1];
  fiber_mft1_off_ang_v1 = fiber_angle_offset[3][2][0];  fiber_mft1_off_ang_v2 = fiber_angle_offset[3][2][1];
  fiber_mft2_off_ang_x1 = fiber_angle_offset[4][0][0];  fiber_mft2_off_ang_x2 = fiber_angle_offset[4][0][1];
  fiber_mft2_off_ang_v1 = fiber_angle_offset[4][1][0];  fiber_mft2_off_ang_v2 = fiber_angle_offset[4][1][1];
  fiber_mft2_off_ang_u1 = fiber_angle_offset[4][2][0];  fiber_mft2_off_ang_u2 = fiber_angle_offset[4][2][1];



  itr_ParamFiles = ParamFiles.find("fiber_mftcor_file");
  if(itr_ParamFiles != ParamFiles.end())
    fiber_name_mftcor  = itr_ParamFiles->second;

  for(int i=0; i<2; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<2; ++k)
        for(int l=0; l<3; ++l)
          fiber_mft_cor_par[i][j][k][l] = 0.;


  std::ifstream ifs_mftcorfiber ( fiber_name_mftcor );
  if(ifs_mftcorfiber.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs_mftcorfiber,temp_line))
        {
          std::stringstream stream(temp_line);
          std::string testComment(stream.str());
          std::size_t it_comment = testComment.find(CommentSymbol);
          if(it_comment!=std::string::npos)
            {
              //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
              continue;
            }

          int det, lay, seg;
          double p0, p1, p2;

          stream >> det >> lay >> seg >> p0 >> p1 >> p2;
          //printf("%d %d %d : %.2f\n", det_id, lay_id, fib_id, offset );

          fiber_mft_cor_par[det][lay][seg][0] = p0;
          fiber_mft_cor_par[det][lay][seg][1] = p1;
          fiber_mft_cor_par[det][lay][seg][2] = p2;
        }
      //std::cout << "done " << mdc_name_map << std::endl;
      printf("fiber mft corrections loaded : %s\n", fiber_name_mftcor.c_str());
    }
  else
    {
      //std::cout << " ! fail to open " << mdc_name_map << std::endl;
      printf(" ! fail to open  : %s\n", fiber_name_mftcor.c_str());
      exit(-1);
    }



  //  PSB  //////
  psb_tcut_min = -20200; psb_tcut_max = -19400;
  psb_pos_x = 0.5; psb_pos_y = 5.5; psb_pos_z = 2760.;
  psb_rot_z = -0.4;
  cut_psb_phi = 0.4;
  cut_psb_z   = 150;
  psb_res_phi = 1000.; //mm
  psb_res_z   = 10.; //mm

  itr_ParamFiles = ParamFiles.find("psb_time_file");
  if(itr_ParamFiles != ParamFiles.end())
    psb_name_time  = itr_ParamFiles->second;
  
  for(int i=0; i<46; ++i){
    psb_zpar[i] = 70;
    for(int j=0; j<2; ++j){
      psb_off_time[i][j] = 0.;
      psb_ch2ns[i][j]    = 0.025;
    }
  }


  std::ifstream ifs_psb_time ( psb_name_time );
  if(ifs_psb_time.is_open())
  {
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs_psb_time,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos)
      {
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }

      int seg;
      double ns_u, ns_d, off_u, off_d, zpar;

      stream >> seg >> ns_u >> ns_d >> off_u >> off_d >> zpar;

      psb_ch2ns[seg][0] = ns_u;
      psb_ch2ns[seg][1] = ns_d;
      psb_off_time[seg][0] = off_u;
      psb_off_time[seg][1] = off_d;
      psb_zpar[seg] = zpar;
    }
    printf("PSB time file done : %s\n", psb_name_time.c_str());
  }else
  {
    printf(" ! fail to open  : %s\n", psb_name_time.c_str());
    exit(-1);
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
  t0_tcut_min = -26200 ; t0_tcut_max = -25800;
  t0_pos_z = 100; //tmp
  cut_t0_time = 2.;

  itr_ParamFiles = ParamFiles.find("t0_time_file");
  if(itr_ParamFiles != ParamFiles.end())
    t0_name_time  = itr_ParamFiles->second;

  for(int i=0; i<28; ++i){
    for(int j=0; j<2; ++j){
      t0_off_time[i][j] = 0.;
      t0_ch2ns[i][j]    = 0.025;
    }
  }

  std::ifstream ifs_t0 ( t0_name_time );
  if(ifs_t0.is_open())
  {
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs_t0,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos)
      {
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }

      int seg;
      double ns_u, ns_d, off_u, off_d;

      stream >> seg >> ns_u >> ns_d >> off_u >> off_d;

      t0_ch2ns[seg][0] = ns_u;
      t0_ch2ns[seg][1] = ns_d;
      t0_off_time[seg][0] = off_u;
      t0_off_time[seg][1] = off_d;

    }
    //std::cout << "done " << mdc_name_map << std::endl;
    printf("t0 offset loaded : %s\n", t0_name_time.c_str());
  }
  else
  {
    //std::cout << " ! fail to open " << mdc_name_map << std::endl;
    printf(" ! fail to open  : %s\n", t0_name_time.c_str());
    exit(-1);
  }


  //  MDC  //////
  mdc_cut_min  = -1350;
  mdc_cut_max1 = -1150; mdc_cut_max2  = -1080; mdc_cut_max3  = -1000;
  mdc_cut_tot =   180;

  mdc_ch2ns  = 0.416666;

  mdc_pos_x = 1.5;
  mdc_pos_y = 5.5;
  mdc_pos_z = 2760. - 4.;

  mdc_rot_x = -0.32;
  mdc_rot_y = -0.1;
  mdc_rot_z = -0.4;

  mdc_res = 0.3;

  InitMDCParameter(ParamFiles);


  /// ---- S4 Scintillators ---///
  mtdc_ch2ns = 0.015625; // exactly 1/64.
  tcut_sc31_min = 21800; // run802, He3 setting, maybe not needed?
  tcut_sc31_max = 22500;
  tcut_sc41_1_min = 11200;
  tcut_sc41_1_max = 11700;
  tcut_sc41_2_min = 10700;
  tcut_sc41_2_max = 11500;
  tcut_sc41_3_min = 10800;
  tcut_sc41_3_max = 11500;
  tcut_sc42_min = 12000;
  tcut_sc42_max = 12700;
  tcut_sc43_1_min = 12900;
  tcut_sc43_1_max = 13500;
  tcut_sc43_2_min = 12900;
  tcut_sc43_2_max = 13500;
  tcut_sc43_3_min = 12900;
  tcut_sc43_3_max = 13500;
  tcut_sc43_4_min = 12900;
  tcut_sc43_4_max = 13500;

  // 3He setting (Brho=13.1308Tm, beta=0.941871, velocity=282.366mm/ns)
  // offset_tof_sc3141 = 61.231 - (-168.5); // [ns]
  // offset_tof_sc4143 = 8.6352 - 28.8; // [ns]   --> ~run523, see setup.C of wasa_unpack
  offset_tof_sc3141 = 62.287 - (-169.8); // [ns]
  offset_tof_sc4143 = 8.7842 - 28.99; // [ns]

  /// --- mwdc --- ///
  dist_focS4 = 3200.0; //etaprime, hyphi
  dist_MWDC_zref = 992.0; //fisrt one plate. second one is set as relative +1000 by MWDCParameter
  dist_SC41 = 2190.0;
  dist_SC42 = 2321.0;
  dist_SC43 = 4668.0;
  dist_s4test = 3200.0; // y waist for etaprime, hyphi, and near chiara's detector test table
  mwdc_shift_x_alignment = 0.0;
  mwdc_shift_a_alignment = 0.0;
  mwdc_shift_y_alignment = -5.5;
  mwdc_shift_b_alignment = -7.0 + 5.5;
  mwdc_lt_valid_min = -1000.; 
  mwdc_lt_valid_max = -500.;//mm This is still initial value. you can change it in setup.C
  mwdc_assumed_plane_resolution = 0.5;  //mm This is still initial value. you can change it in setup.C
  mwdc_tracking_chi2cut_max     = 20.0; //   This is still initial value. you can change it in setup.C
  participate_tracking_plane[0] = 1; //for fast analysis in online Go4, maybe we should remove some planes.
  participate_tracking_plane[1] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[2] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[3] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[4] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[5] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[6] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[7] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[8] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[9] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[10] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[11] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[12] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[13] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[14] = 1;//just initial value. Change it in setup.C
  participate_tracking_plane[15] = 1;//just initial value. Change it in setup.C
  mwdc_max_hit_combination_for_tracking = 10;//just initial value. Change it in setup.C
  mwdc_min_plane_with_hit_for_tracking = 12;//just initial value. Change it in setup.C

  mwdc_ch2ns = 1.0; //modify later

  itr_ParamFiles = ParamFiles.find("mwdc_name_dtdx");
  if(itr_ParamFiles != ParamFiles.end())
    mwdc_name_dtdx  = itr_ParamFiles->second;

  for(int i=0; i<16; ++i)
  {
    mwdc_t0_off[i] = -9999.;
    mwdc_tmax[i] = -9999.;
    for(int j=0; j<8; ++j) mwdc_dtdx_par[i][j] = 0.;
  }
  

  std::ifstream ifs_mwdc_dtdx ( mwdc_name_dtdx );
  if(ifs_mwdc_dtdx.is_open())
  {
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs_mwdc_dtdx,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos)
      {
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }

      int plane_tmp = -99;
      double mwdc_t0_off_ = -9999.;
      double mwdc_tmax_ = -9999.;
      double par0_tmp, par1_tmp, par2_tmp, par3_tmp, par4_tmp, par5_tmp, par6_tmp, par7_tmp;

      stream >> plane_tmp >> mwdc_t0_off_ >> mwdc_tmax_ >> par0_tmp >> par1_tmp >> par2_tmp >> par3_tmp >> par4_tmp >> par5_tmp >> par6_tmp >> par7_tmp;

      mwdc_t0_off[plane_tmp] = mwdc_t0_off_;
      mwdc_tmax[plane_tmp] = mwdc_tmax_;

      mwdc_dtdx_par[plane_tmp][0] = par0_tmp;
      mwdc_dtdx_par[plane_tmp][1] = par1_tmp;
      mwdc_dtdx_par[plane_tmp][2] = par2_tmp;
      mwdc_dtdx_par[plane_tmp][3] = par3_tmp;
      mwdc_dtdx_par[plane_tmp][4] = par4_tmp;
      mwdc_dtdx_par[plane_tmp][5] = par5_tmp;
      mwdc_dtdx_par[plane_tmp][6] = par6_tmp;
      mwdc_dtdx_par[plane_tmp][7] = par7_tmp;
    }
    //std::cout << "done " << mdc_name_map << std::endl;
    printf("mwdc dtdx loaded : %s\n", mwdc_name_dtdx.c_str());
  }
  else
  {
    //std::cout << " ! fail to open " << mdc_name_map << std::endl;
    printf(" ! fail to open  : %s\n", mwdc_name_dtdx.c_str());
    exit(-1);
  }






  // clocktdcc to MWDC channel assignment. This will most probably not changed. So I fix this here.
  // [id_board_of_clocktdc][id_group_clock_tdc(16ch as a unit)]  to id_plane(mwdc) and id_wiregroup(ASD group:0,1,2)
  id_plane[0][0] = 0;   id_wiregroup[0][0] = 2; // MWDC41-X1-gr3
  id_plane[0][1] = 2;   id_wiregroup[0][1] = 2; // MWDC41-X2gr3
  id_plane[0][2] = 6;   id_wiregroup[0][2] = 0; // MWDC41-V-gr1
  id_plane[0][3] = 4;   id_wiregroup[0][3] = 0; // MWDC41-U-gr1
  id_plane[0][4] = 2;   id_wiregroup[0][4] = 0; // MWDC41-X2-gr1
  id_plane[0][5] = 0;   id_wiregroup[0][5] = 0; // MWDC41-X1-gr1
  id_plane[0][6] = 4;   id_wiregroup[0][6] = 1; // MWDC41-U-gr2
  id_plane[0][7] = 6;   id_wiregroup[0][7] = 1; // MWDC41-V-gr2
  //
  id_plane[1][0] = 5;   id_wiregroup[1][0] = 1;  // MWDC41-U'-gr2
  id_plane[1][1] = 7;   id_wiregroup[1][1] = 1;  // MWDC41-V'-gr2
  id_plane[1][2] = 3;   id_wiregroup[1][2] = 0;  // MWDC41-X2'-gr1
  id_plane[1][3] = 1;   id_wiregroup[1][3] = 0;  // MWDC41-X1'-gr1
  id_plane[1][4] = 2;   id_wiregroup[1][4] = 1;  // MWDC41-X2-gr2
  id_plane[1][5] = 0;   id_wiregroup[1][5] = 1;  // MWDC41-X1-gr2
  id_plane[1][6] = 4;   id_wiregroup[1][6] = 2;  // MWDC41-U-gr3
  id_plane[1][7] = 6;   id_wiregroup[1][7] = 2;  // MWDC41-V-gr3
  //
  id_plane[2][0] = 5;   id_wiregroup[2][0] = 2;  // MWDC41-U'-gr3
  id_plane[2][1] = 7;   id_wiregroup[2][1] = 2;  // MWDC41-V'-gr3
  id_plane[2][2] = 3;   id_wiregroup[2][2] = 1;  // MWDC41-X2'-gr2
  id_plane[2][3] = 1;   id_wiregroup[2][3] = 1;  // MWDC41-X1'-gr2
  id_plane[2][4] = 7;   id_wiregroup[2][4] = 0;  // MWDC41-V'-gr1
  id_plane[2][5] = 5;   id_wiregroup[2][5] = 0;  // MWDC41-U'-gr1
  id_plane[2][6] = 1;   id_wiregroup[2][6] = 2;  // MWDC41-X1'-gr3
  id_plane[2][7] = 3;   id_wiregroup[2][7] = 2;  // MWDC41-X2'-gr3
  //
  id_plane[3][0] = 8;   id_wiregroup[3][0] = 2; // MWDC42-X1-gr3
  id_plane[3][1] =10;   id_wiregroup[3][1] = 2; // MWDC42-X2gr3
  id_plane[3][2] =14;   id_wiregroup[3][2] = 0; // MWDC42-V-gr1
  id_plane[3][3] =12;   id_wiregroup[3][3] = 0; // MWDC42-U-gr1
  id_plane[3][4] =10;   id_wiregroup[3][4] = 0; // MWDC42-X2-gr1
  id_plane[3][5] = 8;   id_wiregroup[3][5] = 0; // MWDC42-X1-gr1
  id_plane[3][6] =12;   id_wiregroup[3][6] = 1; // MWDC42-U-gr2
  id_plane[3][7] =14;   id_wiregroup[3][7] = 1; // MWDC42-V-gr2
  //
  id_plane[4][0] = 9;   id_wiregroup[4][0] = 0;  // MWDC42-X1'-gr1
  id_plane[4][1] =11;   id_wiregroup[4][1] = 0;  // MWDC42-X2'-gr1
  id_plane[4][2] =15;   id_wiregroup[4][2] = 1;  // MWDC42-V'-gr2
  id_plane[4][3] =13;   id_wiregroup[4][3] = 1;  // MWDC42-U'-gr2
  id_plane[4][4] =10;   id_wiregroup[4][4] = 1;  // MWDC42-X2-gr2
  id_plane[4][5] = 8;   id_wiregroup[4][5] = 1;  // MWDC42-X1-gr2
  id_plane[4][6] =12;   id_wiregroup[4][6] = 2;  // MWDC42-U-gr3
  id_plane[4][7] =14;   id_wiregroup[4][7] = 2;  // MWDC42-V-gr3
  //
  id_plane[5][0] = 9;   id_wiregroup[5][0] = 1;  // MWDC42-X1'-gr2
  id_plane[5][1] =11;   id_wiregroup[5][1] = 1;  // MWDC42-X2'-gr2
  id_plane[5][2] =15;   id_wiregroup[5][2] = 2;  // MWDC42-V'-gr3
  id_plane[5][3] =13;   id_wiregroup[5][3] = 2;  // MWDC42-U'-gr3
  id_plane[5][4] =11;   id_wiregroup[5][4] = 2;  // MWDC42-X2'-gr3
  id_plane[5][5] = 9;   id_wiregroup[5][5] = 2;  // MWDC42-X1'-gr3
  id_plane[5][6] =13;   id_wiregroup[5][6] = 0;  // MWDC42-U'-gr1
  id_plane[5][7] =15;   id_wiregroup[5][7] = 0;  // MWDC42-V'-gr1
  //
  id_plane[6][0] = -1;   id_wiregroup[6][0] = -1;  //
  id_plane[6][1] = -1;   id_wiregroup[6][1] = -1;  //
  id_plane[6][2] = -1;   id_wiregroup[6][2] = -1;  //
  id_plane[6][3] = -1;   id_wiregroup[6][3] = -1;  //
  id_plane[6][4] = -1;   id_wiregroup[6][4] = -1;  //
  id_plane[6][5] = -1;   id_wiregroup[6][5] = -1;  //
  id_plane[6][6] = -1;   id_wiregroup[6][6] = -1;  //
  id_plane[6][7] = -1;   id_wiregroup[6][7] = -1;  //
  //
  id_plane[7][0] = -1;   id_wiregroup[7][0] = -1;  //
  id_plane[7][1] = -1;   id_wiregroup[7][1] = -1;  //
  id_plane[7][2] = -1;   id_wiregroup[7][2] = -1;  //
  id_plane[7][3] = -1;   id_wiregroup[7][3] = -1;  //
  id_plane[7][4] = -1;   id_wiregroup[7][4] = -1;  //
  id_plane[7][5] = -1;   id_wiregroup[7][5] = -1;  //
  id_plane[7][6] = -1;   id_wiregroup[7][6] = -1;  //
  id_plane[7][7] = -1;   id_wiregroup[7][7] = -1;  //
  //

  //Configuration for WASA-2022
  //This is identical to S436/S437 July 2014 (MWDC41="B" and 42="A")
  //Only COSY-2014Jan experiment, everything was opposite
  // plane information: angle, sign ([i_order] from upstream) // !!! COSY 2014/Jan. with xy plane seeing from upstream side (Left-handed-xyz) !!!
  float angle_X =   0.0;  // unit is rad.
  float angle_V = -15.0 * 3.14159265/180.0; // sign of angle is depending on the MWDC direction. (+ for COSY, - for FRS)
  float angle_U = +15.0 * 3.14159265/180.0; //  (- for COSY, + for FRS)
  mwdc_plane_angle[15]  = angle_V ;  mwdc_plane_sign[15] = +1.0; // - for cosy
  mwdc_plane_angle[14]  = angle_V ;  mwdc_plane_sign[14] = -1.0;
  mwdc_plane_angle[13]  = angle_U ;  mwdc_plane_sign[13] = +1.0;
  mwdc_plane_angle[12]  = angle_U ;  mwdc_plane_sign[12] = -1.0;
  mwdc_plane_angle[11]  = angle_X ;  mwdc_plane_sign[11] = +1.0;
  mwdc_plane_angle[10]  = angle_X ;  mwdc_plane_sign[10] = -1.0;
  mwdc_plane_angle[9]   = angle_X ;  mwdc_plane_sign[9]  = +1.0;
  mwdc_plane_angle[8]   = angle_X ;  mwdc_plane_sign[8]  = -1.0;
  mwdc_plane_angle[7]   = angle_V ;  mwdc_plane_sign[7]  = +1.0;
  mwdc_plane_angle[6]   = angle_V ;  mwdc_plane_sign[6]  = -1.0;
  mwdc_plane_angle[5]   = angle_U ;  mwdc_plane_sign[5]  = +1.0;
  mwdc_plane_angle[4]   = angle_U ;  mwdc_plane_sign[4]  = -1.0;
  mwdc_plane_angle[3]   = angle_X ;  mwdc_plane_sign[3]  = +1.0;
  mwdc_plane_angle[2]   = angle_X ;  mwdc_plane_sign[2]  = -1.0;
  mwdc_plane_angle[1]   = angle_X ;  mwdc_plane_sign[1]  = +1.0;
  mwdc_plane_angle[0]   = angle_X ;  mwdc_plane_sign[0]  = -1.0;

  //Configuration for WASA-2022
  //This is identical to S436/S437 July 2014 (MWDC41="B" and 42="A")
  //Only COSY-2014Jan experiment, everything was opposite

  float z_pos_B = 0.0; // MWDC-41! in mm (see definition in )
  float z_pos_A = 1000.0;// MWDC-42! in mm (see definition in )
  float z_shift_VP = 0.0;
  float z_shift_V  = -4.8;
  float z_shift_UP = -24.4;
  float z_shift_U  = -29.2;
  float z_shift_XP = -48.8;
  float z_shift_X  = -53.6;
  float z_shift_XPN= -73.2;
  float z_shift_XN = -78.0;

  sprintf(plane_name[15],"42_Vp") ; mwdc_zpos[15] = z_pos_A + z_shift_VP ; mwdc_center_id[15] = 23.5-0.015-0.003-0.010903+0.001762-0.001015;//+
  sprintf(plane_name[14],"42_V")  ; mwdc_zpos[14] = z_pos_A + z_shift_V  ; mwdc_center_id[14] = 24.0+0.004-0.003+0.010925-0.001998+0.001349;//-
  sprintf(plane_name[13],"42_Up") ; mwdc_zpos[13] = z_pos_A + z_shift_UP ; mwdc_center_id[13] = 23.5-0.015-0.014808+0.004182-0.001338;//+
  sprintf(plane_name[12],"42_U")  ; mwdc_zpos[12] = z_pos_A + z_shift_U  ; mwdc_center_id[12] = 23.0-0.005+0.007541-0.000596-0.000224;//- //center_id[3] = 24.0;
  sprintf(plane_name[11],"42_Xp") ; mwdc_zpos[11] = z_pos_A + z_shift_XP ; mwdc_center_id[11] = 23.5+0.018-0.004561-0.001762-0.000085;//+
  sprintf(plane_name[10],"42_X")  ; mwdc_zpos[10] = z_pos_A + z_shift_X  ; mwdc_center_id[10] = 24.0-0.003+0.002888+0.001891+0.000045;//-
  sprintf(plane_name[9] ,"42_Xpn"); mwdc_zpos[9]  = z_pos_A + z_shift_XPN; mwdc_center_id[9]  = 23.5+0.025177-0.003251+0.000202;
  sprintf(plane_name[8] ,"42_Xn") ; mwdc_zpos[8]  = z_pos_A + z_shift_XN ; mwdc_center_id[8]  = 24.0-0.014913+0.002477+0.000136;
  sprintf(plane_name[7] ,"41_Vp") ; mwdc_zpos[7]  = z_pos_B + z_shift_VP ; mwdc_center_id[7]  = 23.5-0.001-0.008790+0.001490-0.000479;//+
  sprintf(plane_name[6] ,"41_V")  ; mwdc_zpos[6]  = z_pos_B + z_shift_V  ; mwdc_center_id[6]  = 24.0+0.031+0.012508-0.003461+0.001316;//-
  sprintf(plane_name[5] ,"41_Up") ; mwdc_zpos[5]  = z_pos_B + z_shift_UP ; mwdc_center_id[5]  = 23.5-0.021-0.013236+0.003383-0.001517;//+
  sprintf(plane_name[4] ,"41_U")  ; mwdc_zpos[4]  = z_pos_B + z_shift_U  ; mwdc_center_id[4]  = 23.0+0.005+0.007646-0.001212+0.000387;//-// center_id[11]= 24.0;
  sprintf(plane_name[3] ,"41_Xp") ; mwdc_zpos[3]  = z_pos_B + z_shift_XP ; mwdc_center_id[3]  = 23.5+0.031-0.011035-0.001411+0.000153;//+
  sprintf(plane_name[2] ,"41_X")  ; mwdc_zpos[2]  = z_pos_B + z_shift_X  ; mwdc_center_id[2]  = 24.0-0.026+0.005554+0.000983-0.000103;//-
  sprintf(plane_name[1] ,"41_Xpn"); mwdc_zpos[1]  = z_pos_B + z_shift_XPN; mwdc_center_id[1]  = 23.5+0.025290-0.002695+0.000274;
  sprintf(plane_name[0] ,"41_Xn") ; mwdc_zpos[0]  = z_pos_B + z_shift_XN ; mwdc_center_id[0]  = 24.0-0.015982+0.002609-0.000248;

  float *ptable[16];
  for (int ii = 0; ii < 16; ii++)
  {
    ptable[ii] = (mwdc_dtdxtable[ii]);
  }
  // dtdx table
  
  //char filename_dtdxtable[256];
  //sprintf(filename_dtdxtable, "setup/mwdc_dtdxtable/mwdcdtdx_run0735.root");

  itr_ParamFiles = ParamFiles.find("mwdc_name_dtdxtable");
  if(itr_ParamFiles != ParamFiles.end())
    strcpy(mwdc_name_dtdxtable, (itr_ParamFiles->second).c_str());

  SetMWDCdtdxFromTH1D(mwdc_name_dtdxtable, ptable, (mwdc_dtdxtable_dtmin), (mwdc_dtdxtable_dtmax));

}


ParaManager::~ParaManager() {}


bool ParaManager::InitMDCParameter(const std::unordered_map<std::string,std::string>& ParamFilesMDC)
{
  if(mdc_init_done)
    return true;

  itr_ParamFiles = ParamFilesMDC.find("mdc_map_file");
  if(itr_ParamFiles != ParamFilesMDC.end())
    mdc_name_map  = itr_ParamFiles->second;

  itr_ParamFiles = ParamFilesMDC.find("mdc_phys_file");
  if(itr_ParamFiles != ParamFilesMDC.end())
    mdc_name_phys  = itr_ParamFiles->second;

  itr_ParamFiles = ParamFilesMDC.find("mdc_drift_file");
  if(itr_ParamFiles != ParamFilesMDC.end())
    mdc_name_drift  = itr_ParamFiles->second;

  itr_ParamFiles = ParamFilesMDC.find("mdc_parT0_file");
  if(itr_ParamFiles != ParamFilesMDC.end())
    mdc_name_t0  = itr_ParamFiles->second;

  itr_ParamFiles = ParamFilesMDC.find("mdc_parT0wir_file");
  if(itr_ParamFiles != ParamFilesMDC.end())
    mdc_name_t0_wir  = itr_ParamFiles->second;

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


  for(int i=0; i<17; ++i){
    mdc_t0_off[i] = 530.;
    for(int j=0; j<148; ++j){
      mdc_t0_off_wir[i][j] = 0.;
    }
  }

  std::ifstream ifs_mdc ( mdc_name_map );
  if(ifs_mdc.is_open())
  {
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs_mdc,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos)
      {
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }

      int cmp_id, cmp_block, mdc_layer, mdc_sector, wire_id_begin, wire_id_end, wire_nb, ctdc_id, ctdc_block;

      stream >> cmp_id >> cmp_block >> mdc_layer >> mdc_sector >> wire_id_begin >> wire_id_end >> wire_nb >> ctdc_id >> ctdc_block;
      //std::cout<<"cmp ["<<cmp_id<<" / "<<cmp_block<<"] MDC layer"<<mdc_layer<<" "<<mdc_sector<<" wire ["<<wire_id_begin<<", "<< wire_id_end<<"] #"<<wire_nb<<" : cTDC = ["<<ctdc_id<<" "<<ctdc_block<<"] \n";

      /// TRIG_REF: cmp_id = 27, cmp_block = 3 -> ctdc_id =  6, ctdc_block = 6 -> MDC layer, etc = -10
      /// Empty:    cmp_id = 11, cmp_block = 0 -> ctdc_id = 10, ctdc_block = 5 -> MDC layer, etc = -1
      /// Empty:    cmp_id = 11, cmp_block = 3 -> ctdc_id =  2, ctdc_block = 6 -> MDC layer, etc = -1
      /// Empty:    cmp_id = 31, cmp_block = 2 -> ctdc_id =  7, ctdc_block = 5 -> MDC layer, etc = -1

      if(mdc_layer < 0)
        continue;

      mdc_MappingMDC[ctdc_id][ctdc_block].cmp_id = cmp_id;
      mdc_MappingMDC[ctdc_id][ctdc_block].cmp_block = cmp_block;
      mdc_MappingMDC[ctdc_id][ctdc_block].mdc_layer = mdc_layer;
      mdc_MappingMDC[ctdc_id][ctdc_block].mdc_sector = mdc_sector;
      mdc_MappingMDC[ctdc_id][ctdc_block].wire_id_begin = wire_id_begin;
      mdc_MappingMDC[ctdc_id][ctdc_block].wire_id_end = wire_id_end;
      mdc_MappingMDC[ctdc_id][ctdc_block].wire_nb = wire_nb;

    }
    //std::cout << "done " << mdc_name_map << std::endl;
    printf("done MDC Channel Map: %s\n", mdc_name_map.c_str());
  }else
  {
    //std::cout << " ! fail to open " << mdc_name_map << std::endl;
    printf(" ! fail to open  : %s\n", mdc_name_map.c_str());
    exit(-1);
  }

  std::ifstream ifs2_mdc ( mdc_name_phys );
  if(ifs2_mdc.is_open())
  {
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs2_mdc,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos)
      {
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }
      int layerID;
      double Dpsi1, Dpsi0b, Dpsi0f, dpsi0, rmax;
      int nb_wire;
      double wire_size, z_back, z_front, l_back, l_front;
      stream >> layerID;
      stream >> Dpsi1 >> Dpsi0b >> Dpsi0f >> dpsi0 >> rmax;
      stream >> nb_wire;
      stream >> wire_size >> z_back >> z_front >> l_back >> l_front;

      //std::cout<<"layer ID "<<layerID<<" psi ["<<Dpsi1<<" "<<Dpsi0b<<" "<<Dpsi0f<<" "<<dpsi0<<"] Nb wire:"<<nb_wire<<" "<<wire_size<<" :Z ["<<z_front<<" "<<z_back<<"] L ["<<l_front<<" "<<l_back<<"] \n";

      //double zoff = atof(getenv("ZOFF"));
      double zoff = 0;
      //double zoff = -6;
      //double zoff = -10.;
      //double zoff = -12.;

      mdc_PhysMDC[layerID-1].Dpsi1 = Dpsi1;
      mdc_PhysMDC[layerID-1].Dpsi0b = Dpsi0b;
      mdc_PhysMDC[layerID-1].Dpsi0f = Dpsi0f;
      mdc_PhysMDC[layerID-1].dpsi0 = dpsi0;
      mdc_PhysMDC[layerID-1].rmax = rmax;
      mdc_PhysMDC[layerID-1].nb_wire = nb_wire;
      mdc_PhysMDC[layerID-1].wire_size = wire_size;
      mdc_PhysMDC[layerID-1].z_back = z_back - zoff;
      mdc_PhysMDC[layerID-1].z_front = z_front + zoff;
      mdc_PhysMDC[layerID-1].l_back = l_back;
      mdc_PhysMDC[layerID-1].l_front = l_front;

    }

    //std::cout << "done " << mdc->name_phys << std::endl;
    printf("done MDC Physical Map: %s\n", mdc_name_phys.c_str());
  }else
  {
    //std::cout << " ! fail to open " << mdc->name_phys << std::endl;
    printf(" ! fail to open  : %s\n", mdc_name_phys.c_str());
    exit(-1);
  }


  for(size_t i = 0; i< mdc_MappingMDC.size() ; ++i)
    for(size_t j = 0 ;j < mdc_MappingMDC[i].size(); ++j)
    {
      if(mdc_MappingMDC[i][j].mdc_layer>=0)
      {
        int reduc = mdc_PhysMDC[mdc_MappingMDC[i][j].mdc_layer-1].nb_wire / 2 ;

        mdc_MappingMDC[i][j].wire_real_start = mdc_MappingMDC[i][j].wire_id_begin + reduc*(mdc_MappingMDC[i][j].mdc_sector-1);
        mdc_MappingMDC[i][j].wire_real_end   = mdc_MappingMDC[i][j].wire_id_end   + reduc*(mdc_MappingMDC[i][j].mdc_sector-1);

        //std::cout<<"block :TDC["<<i<<", "<<j<<"] : cmp ["<<mdc->MappingMDC[i][j].cmp_id<<" / "<<mdc->MappingMDC[i][j].cmp_block<<"]MDClayer : "<<mdc->MappingMDC[i][j].mdc_layer<<" "<<mdc->MappingMDC[i][j].mdc_sector<<" wire ["<<mdc->MappingMDC[i][j].wire_id_begin<<", "<<MappingMDC [i][j].wire_id_end<<"] \n"
      }
    }


  std::ifstream ifs3_mdc ( mdc_name_drift );
  if(ifs3_mdc.is_open()){
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs3_mdc,temp_line)){
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos){
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }
      int layerID;
      float par0, par1, par2, par3, par4, par5, par6;
      stream >> layerID;
      stream >> par0 >> par1 >> par2 >> par3 >> par4 >> par5 >> par6;

      mdc_PhysMDC[layerID].drift_par0 = par0;
      mdc_PhysMDC[layerID].drift_par1 = par1;
      mdc_PhysMDC[layerID].drift_par2 = par2;
      mdc_PhysMDC[layerID].drift_par3 = par3;
      mdc_PhysMDC[layerID].drift_par4 = par4;
      mdc_PhysMDC[layerID].drift_par5 = par5;
      mdc_PhysMDC[layerID].drift_par6 = par6;
    }

    printf("done MDC DriftParam File: %s\n", mdc_name_drift.c_str());
  }else
  {
    printf(" ! fail to open  : %s\n", mdc_name_drift.c_str());
    exit(-1);
  }

  std::ifstream ifs4_mdc ( mdc_name_t0 );
  if(ifs4_mdc.is_open()){
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs4_mdc,temp_line)){
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos){
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }
      int layerID;
      float par0;
      stream >> layerID;
      stream >> par0;

      mdc_t0_off[layerID] = par0;
    }

    printf("done MDC T0Param File: %s\n", mdc_name_t0.c_str());
  }else
  {
    printf(" ! fail to open  : %s\n", mdc_name_t0.c_str());
    exit(-1);
  }

  std::ifstream ifs5_mdc ( mdc_name_t0_wir );
  if(ifs5_mdc.is_open()){
    const std::string CommentSymbol("#");

    std::string temp_line;
    while(std::getline(ifs5_mdc,temp_line)){
      std::stringstream stream(temp_line);
      std::string testComment(stream.str());
      std::size_t it_comment = testComment.find(CommentSymbol);
      if(it_comment!=std::string::npos){
        //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
        continue;
      }
      int layerID;
      int wireID;
      float par0;
      stream >> layerID >> wireID;
      stream >> par0;

      mdc_t0_off_wir[layerID][wireID] = par0;
      //printf("%d %d %.2f\n", layerID, wireID, par0);
    }

    printf("done MDC T0wirParam File: %s\n", mdc_name_t0_wir.c_str());
  }else
  {
    printf(" ! fail to open  : %s\n", mdc_name_t0_wir.c_str());
    exit(-1);
  }

  return true;
}


float ParaManager::mwdc_dtdxconversion(int i_plane, int timeval){
  if(timeval < mwdc_dtdxtable_dtmin[i_plane] || mwdc_dtdxtable_dtmax[i_plane] < timeval){ return -999.; }
  int i_array = timeval - mwdc_dtdxtable_dtmin[i_plane];
  return mwdc_dtdxtable[i_plane][i_array];
}

void ParaManager::SetMWDCdtdxFromTH1D(char* filenametmp, float** ptable, int* dtmin, int* dtmax){
  printf("\n-----------------------------------------\n");
  printf("MWDC DTDX conversion table (%s)\n",filenametmp);
  TFile *_file0 = TFile::Open(filenametmp);
  TH1D *hist[16]; //32 for signal-IN + 4 for TRn
  for(int i=0; i<16; i++){
    hist[i]  = (TH1D*)((_file0->Get(Form("dtdx_graph_plane%d",i)))->Clone());
    dtmin[i] = (int)((hist[i] ->GetXaxis()->GetXmin())) ;
    dtmax[i] = (int)((hist[i] ->GetXaxis()->GetXmax())) ;
    TGraph tempgr(hist[i]);
    for(int j=0; j<(dtmax[i]-dtmin[i]+1); j++){
      double dt_here = (double)(dtmin[i] + j);
      ptable[i][j] = tempgr.Eval(dt_here);
      //        printf("iplane=%d : ipoint = %d : >> dl = %f\n",i,j,ptable[i][j]);
    }
    hist[i]->Delete();
    //    printf("\n\n");
  }
  _file0->Close();
  _file0->Delete();
  printf("... done\n\n");fflush(stdout);
}