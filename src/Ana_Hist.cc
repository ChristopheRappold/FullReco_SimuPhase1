#include "Ana_Hist.hh"
#include "TFile.h"
#include "Riostream.h"

#include "TGeoManager.h"
using namespace std;

Ana_Hist::~Ana_Hist()
{
  //cout<<"~Ana_Hist()"<<endl;
  if(HaveBeenWritten==false)
    {
      TFile* f= new TFile("Ana_Hist_Default.root","RECREATE");
      Write(f);
      f->Close();
    }
}

/********************************************************************/
Ana_Hist::Ana_Hist(bool Daf, bool Vertex, bool DCproject, bool Finding, bool Hough, bool Simu)
{
  EnableState.resize(SIZEOF_STATEHIST);
  EnableState[DAF] = Daf;
  EnableState[VERTEX] = Vertex;
  EnableState[DCPROJ] = DCproject;
  EnableState[FINDING] = Finding;
  EnableState[HOUGH] = Hough;
  EnableState[SIMU] = Simu;

  cout<<"Creating Histograms ";
  for(int i=0;i<SIZEOF_STATEHIST;++i)
    if(EnableState[i]==true)
      cout<<"Hist ["<<i<<"] ";

  HaveBeenWritten=false;

#ifdef ADD_HISTO

  h_task_exit = new TH1I("h_task_exit","h_task_exit",20,0,20);
  
  h_2utrPos_XZ.resize(21);
  h_2utrPos_YZ.resize(21);
  h_2utrPos_XY.resize(21);
  h_utrEnergy.resize(21);
  h_utrTime.resize(21);
  h_pull_X_det.resize(21);
  h_pull_Y_det.resize(21);
  
  h_layer_x.resize(21);

  TString tmpstr;
  //  TString nameDetIn[8] = {"Tr0","Tr1","Tr2","CbM","CaM","TofAladin","TofPlus","TofPlusPlus"};
 
 TString nameDetIn[21] = {"Tr0x","Tr0y","Tr1x","Tr1y","Tr2x","Tr2y",
                          "CbMu","CbMup","CbMv","CbMvp","CbMx","CbMxp",
                          "CaMu","CaMx","CaMxp","CaMy","CaMyp",
                          "TofAladin","TofPlus","TofPlusBar","TofPlusPlus"};  


  TString nameLayIn[21] = {"Sebxxxxxx","Tr0y","Tr1x","TR1y","TR2x","TR2y",
                          "CbMu","CbMup","CbMv","CbMvp","CbMx","CbMxp",
                          "CaMu","CaMx","CaMxp","CaMy","CaMyp",
                          "TofAladin","TofPlus","TofPlusBar","TofPlusPlus"};

  if(Simu)
    {
      nameDetIn[17] = "STOP"; 

    }
  //double sizeX[21]={40.,40.,140.,140.,250.,250.,250.,250.,250.,250.,140.,140.,1200.,2500.,1000.,1500.};
  //double sizeY[21]={40.,80.,120.,80.,1000.,1100.,1000.,1300.};
  //double PosX[21]={0.,25.,50.,25.,-200.,2000.,-100, 2000.};
 
 
 //  TString nameDetIn[8] = {"Tr0","Tr1","Tr2","CbM","CaM","TofAladin","TofPlus","TofPlusPlus"};
 //TString nameDetIn[18] = {"Tr0","Bfm1p","Bfm2p","Bfm1n","Bfm2n",
 //			  "CbMx","CbMxp","CbMu","CbMup","CbMv","CbMvp",
 //			  "CaMx","CaMxp","CaMy","CaMyp","CaMu",
 //			  "Afm2p","Afm2n"};

  double sizeX[21]={25.,25.,70.,70.,130.,130.,
		    300.,300.,300.,300.,300.,300.,
		    700.,700.,700.,700.,700.,
		    2000.,500.,500.,1000.};

  double sizeY[21]={25.,25.,40.,40.,60.,
		    200.,200.,200.,200.,200.,200.,
		    600.,600.,600.,600.,600.,
		    80.,70.,70.,80.};
  
  double PosX[21]={23.,23.,25.45,25.45,-0.9,-0.9,
		   -0.3,-0.3,-0.3,-0.3,-0.3,-0.3,
		   338.,338.,338.,338.,338.,
		   1709,-269.29,-269.29,1709.7};
  
  //std::vector<double> sizeX(21,0.);
  //std::vector<double> sizeY(21,0.);
  //std::vector<double> PosX(21,0.);

  for (Int_t n=0;n<21;n++)
    {
      tmpstr = nameDetIn[n];
      tmpstr+="Position_XY";
      h_2utrPos_XY[n] = new TH2F(tmpstr,tmpstr,200,PosX[n]-sizeX[n],PosX[n]+sizeX[n],
				 200,-sizeY[n],sizeY[n]);
      tmpstr = nameDetIn[n];
      tmpstr+="Position_YZ";
      h_2utrPos_YZ[n] = new TH2F(tmpstr,tmpstr,200,-2000.,2000.,200,-3000.,3000);
      
      tmpstr = nameDetIn[n];
      tmpstr+="Position_XZ";
      h_2utrPos_XZ[n] = new TH2F(tmpstr,tmpstr,200,-4000.,4000.,200,-3000.,3000.);
      tmpstr = nameDetIn[n];
      tmpstr+="Energy";
      h_utrEnergy[n] = new TH1F(tmpstr,tmpstr,200,0,200);
      tmpstr = nameDetIn[n];
      tmpstr+="Time";
      h_utrTime[n] = new TH1F(tmpstr,tmpstr,200,0,20);
      tmpstr = nameDetIn[n];
      tmpstr+="pull_X";
      h_pull_X_det[n] = new TH1F(tmpstr,tmpstr,1000,-10,10);
      tmpstr = nameDetIn[n];
      tmpstr+="pull_Y";
      h_pull_Y_det[n] = new TH1F(tmpstr,tmpstr,1000,-10,10);
      
    }

  for (Int_t n=0;n<21;n++)
    {
      tmpstr = nameLayIn[n];
      tmpstr+="testx";
      h_layer_x[n] = new TH1F(tmpstr,tmpstr,200,-4000.,4000.);
    }

  h_2utrRotated_XZ = new TH2F("Rotated_XZ","Rotated_XZ",200,-4000.,4000.,300,-3000.,6000.);
  h_2utrNotRotated_XZ = new TH2F("NotRotated_XZ","NotRotated_XZ",200,-4000.,4000.,300,-3000.,6000.);
  h_2utrFrameMag_XZ = new TH2F("Rotated_XZ_FrameMag","Rotated_XZ_FrameMag",200,-4000.,4000.,300,-3000.,6000.);

  double Aladin_Width = 1560.,Aladin_Length = 1760.;
  for(int i = 0; i < 10;i++)
    for(int j = 0 ;j <= 10;j++)
      {
	double x_temp = Aladin_Width/2.-j*Aladin_Width/10.;
	double z_temp = Aladin_Length/2.-j*Aladin_Length/10.;
	h_2utrFrameMag_XZ->Fill(Aladin_Width/2.,z_temp);
	h_2utrFrameMag_XZ->Fill(-Aladin_Width/2.,z_temp);
	h_2utrFrameMag_XZ->Fill(x_temp,Aladin_Length/2.);
	h_2utrFrameMag_XZ->Fill(x_temp,-Aladin_Length/2.);
      }
  
  h_TOFP_res = new TH2F("h_TOFP_res","h_TOFP_res",100,0.,100.,1000.,0.,4.5);
  h_TOFP_diff_ResInter = new TH2F("h_TOFP_diff_ResInter","h_TOFP_diff_ResInter",100,0.,100.,1000.,-2.,2.);
  h_TOFP_Rdiff_ResInter = new TH2F("h_TOFP_Rdiff_ResInter","h_TOFP_Rdiff_ResInter",100,0.,100.,1000.,-2.,2.);
  for(int j= 0;j<96;++j)
    {
      TString name_R("h_TOFP_res_Path_");
      name_R += j;
      h_TOFP_res_Path[j] = new TH2F(name_R,name_R,500,0.,0.5,1000,470.,520.);
    }
#endif
  // fHough_Z400.resize(50);
  // fHough_Z750.resize(50);

  // for(unsigned int i=0;i<fHough_Z750.size();++i)
  //   {
  //     TString name_hist1("AVHough_Z400_Event_");
  //     name_hist1+=i;

  //     TString name_hist2("AVHough_Z750_Event_");
  //     name_hist2+=i;


  //     fHough_Z400[i]=new TH2F(name_hist1,name_hist1,500,-165.,85.,240,-60.,60.);
  //     fHough_Z750[i]=new TH2F(name_hist2,name_hist2,500,-165.,85.,240,-60.,60.);
      
      //cout<<"Arrr done !"<<i<<endl;
  //    }


  if(DCproject == true)
    {
      h_dc1_project = new TH2I("DC1_Project","DC1_Project",96,-167.5,72.5,140,-70,70);
      h_dc1x_project = new TH2I("DC1x_Project","DC1x_Project",96,-167.5,72.5,140,-70,70);
      h_dc1u_project = new TH2I("DC1u_Project","DC1u_Project",96,-167.5,72.5,140,-70,70);
      h_dc1v_project = new TH2I("DC1v_Project","DC1v_Project",96,-167.5,72.5,140,-70,70);
      
      h_dc2_project = new TH2I("DC2_Project","DC2_Project",256,-258.75,893.25,864,-451.22,412.78);
      h_dc2x_project = new TH2I("DC2x_Project","DC2x_Project",256,-258.75,893.25,864,-451.22,412.78);
      h_dc2y_project = new TH2I("DC2y_Project","DC2y_Project",256,-258.75,893.25,864,-451.22,412.78);
      h_dc2u_project = new TH2I("DC2u_Project","DC2u_Project",256,-258.75,893.25,864,-451.22,412.78);
    }


#ifdef ADD_HISTO

  std::string name[4] = {"_proton","_alpha","_pi-","_he3"};
  std::string name_hyp[4] = {"_L","_H3L","_H4L","_He5L"};
  //double range[4] ={6.,18.,2.,18};
  std::string pv[3]={"_pv_0.01","_pv_0.05","_pv_0.1"};

  for (int i=0;i<4;i++)
    {
      std::string temp;
      //temp = "dmom_rho";
      //temp+=name[i];
      //h_dmom1[i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

      //std::string temp2 = "dmom_angle";
      //temp2+=name[i];
      //h_dmom2[i] = new TH1F(temp2.c_str(),temp2.c_str(),10000,-2,2);

      // temp = "mom_vs_mom_rho";
      // temp+=name[i];
      // temp2 = "mom_vs_mom_angle";
      // temp2+=name[i];

      // h_dmom3[i] = new TH2F(temp.c_str(),temp.c_str(),500, 0.,range[i],500,0.,range[i]);
      // h_dmom4[i] = new TH2F(temp2.c_str(),temp2.c_str(),500,0.,range[i],500,0.,range[i]);

      // temp = "dmom_rho_cut_mass";
      // temp+=name[i];
      // h_dmom_dmass[i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

      // temp = "dmom_rho_cut_mass_best";
      // temp+=name[i];
      // h_dmom_dmass_best[i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

      // temp = "dmom_rho_cut_vertex";
      // temp+=name[i];
      // h_dmom_vertex[i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

      // temp= "dmom_rho_cut_vertex_best";
      // temp+=name[i];
      // h_dmom_vertex_best[i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);


      temp="mass";
      temp+=name[i];
      h_Mass[i]= new TH1F(temp.c_str(),temp.c_str(),10000,0,10.);

      temp="chi2_";
      temp+=name[i];
      h_chi2_particle[i] = new TH1F(temp.c_str(),temp.c_str(),1000,0,10);

      temp="pv_";
      temp+=name[i];
      h_pv_particle[i] = new TH1F(temp.c_str(),temp.c_str(),100,0,1);

      for(int j=0;j<3;j++)
	{
	  // temp = "dmom_rho_cut_mass";
	  // temp+=name[i];
	  // temp+=pv[j];
	  // h_dmom_dmass2[j][i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

	  // temp = "dmom_rho_cut_mass_best";
	  // temp+=name[i];
	  // temp+=pv[j];
	  // h_dmom_dmass_best2[j][i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

	  // temp = "dmom_rho_cut_vertex";
	  // temp+=name[i];
	  // temp+=pv[j];
	  // h_dmom_vertex2[j][i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

	  // temp= "dmom_rho_cut_vertex_best";
	  // temp+=name[i];
	  // temp+=pv[j];
	  // h_dmom_vertex_best2[j][i] = new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

	  // temp= "dmom_after_mass_selection";
	  // temp+=name[i];
	  // temp+=pv[j];
	  // h_dmom_kalman[j][i]=new TH1F(temp.c_str(),temp.c_str(),10000,-2,2);

	  // temp= "dmom_normal_after_mass_selection";
	  // temp+=name[i];
	  // temp+=pv[j];
	  // h_dmom_normal_kalman[j][i]=new TH1F(temp.c_str(),temp.c_str(),10000,-1.,1.);
		
	  // temp= "dmom_mom_after_mass_selection";
	  // temp+=name[i];
	  // temp+=pv[j];
	  // h_dmom_mom_kalman[j][i]=new TH2F(temp.c_str(),temp.c_str(),1000,0.,20.,1000,0.,20.);

	  temp="mass";
	  temp+=name[i];
	  temp+=pv[j];
	  h_Mass2[j][i]= new TH1F(temp.c_str(),temp.c_str(),10000,0,10.);

	}

    }

  h_alpha_rejection = new TH2I("alpha_rejection","alpha_rejection",10,0,10,10,0,10);
  //h_presort = new TH1I("PreSort_Combinasion","PreSort_Combinasion",100,0,10000);
  h_presort_comb = new TH1I("PreSort_Combinasion","PreSort_Combinasion",100,0,10000);
  h_chi2_prefitY = new TH2F("chi2_prefitY","chi2_prefitY",100,0,1000.,5,0,5);
  //h_hough_each = new TH2F("chi2_each","chi2_each",10,0,10,1000,0,5000.);
  h_rchi2_prefitY = new TH1F("Rchi2_prefitY","Rchi2_prefitY",100,0,1000.);

  // h_hough_2d = new TH2I("hough_2d","hough_2d",10,0,10,10,0,10);
  // h_hough = new TH1I("hough","hough",2,0,2);
  // h_hough_pi = new TH1I("hough_pi","hough_pi",2,0,2);
  // h_hough_pi_dist_tof = new TH2I("hough_pi_tof_dist","hough_pi_tof_dist",10000,0,2,20,0,20);
  // h_hough_pi_dist_dc = new TH2I("hough_pi_dc_dist","hough_pi_dc_dist",10000,0,2,20,0,20);
  // h_hough_pi_tof = new TH2I("hough_pi_tof","hough_pi_tof",20,0,20,20,0,20);
  // h_hough_pi_dc = new TH2I("hough_pi_dc","hough_pi_dc",20,0,20,20,0,20);

  // h_hough_pi_dc_true = new TH2I("hough_pi_dc_true","hough_pi_dc_true",20,0,20,20,0,20);
  // h_hough_pi_dc_false = new TH2I("hough_pi_dc_false","hough_pi_dc_false",20,0,20,20,0,20);
  // h_hough_pi_dc_true_tof_true = new TH2I("hough_pi_dc_true_tof_true","hough_pi_dc_true_tof_true",20,0,20,20,0,20);
  // h_hough_pi_tof_true = new TH2I("hough_pi_tof_true","hough_pi_tof_true",20,0,20,20,0,20);
  // h_hough_pi_tof_false = new TH2I("hough_pi_tof_false","hough_pi_tof_false",20,0,20,20,0,20);
  // h_hough_pi_tof_true_dc_true = new TH2I("hough_pi_tof_true_dc_true","hough_pi_tof_true_dc_true",20,0,20,20,0,20);


  // h_hough_2d_optimal = new TH2I("hough_2d_optimal","hough_2d_optimal",10,0,10,10,0,10);
  // h_hough_optimal = new TH1I("hough_optimal","hough_optimal",2,0,2);
  if(Hough == true)
    {
      h_hough_eff = new TH2F("h_hough_eff","h_hough_eff",10,0,10.,10,0,10.);
      h_hough_diff = new TH1F("h_hough_diff","h_hough_diff",20,0,20);
      h_hough_TOF_dist = new TH1F("h_hough_TOF_dist","h_hough_TOF_dist",2000,0,20);
      h_hough_DC_dist = new TH1F("h_hough_DC_dist","h_hough_DC_dist",2000,0,20);
      h_hough_TOF_dist_2 = new TH1F("h_hough_TOF_dist_2","h_hough_TOF_dist_2",10000,0,0.01);
      h_hough_DC_dist_2 = new TH1F("h_hough_DC_dist_2","h_hough_DC_dist_2",10000,0,0.01);
  
      h_hough_tot_dist = new TH1F("h_hough_tot_dist","h_hough_tot_dist",2000,0,20);
      h_hough_tot_dist_2 = new TH1F("h_hough_tot_dist_2","h_hough_tot_dist_2",10000,0,0.01);
      h_hough_tot_pv = new TH1F("h_hough_tot_pv","h_hough_tot_pv",200,0.,1.);
    }

  h_patern_decay= new TH2F("patern_decay","patern_decay",10,0,10,6,0,6);
#endif

  h_stats = new TH1I("stats","stats",10,0,10);

#ifdef ADD_HISTO
  h_presort = new TH1I("presort","presort",100,0,10000);
  // h_stats_Mom = new TH1I("stats_Mom","stats_Mom",4,0,4);
  h_detection_eff_pion = new TH1I("Detection_Eff_pion","Detection_Eff_pion",10,0,10);
  h_detection_eff_proton = new TH1I("Detection_Eff_proton","Detection_Eff_proton",10,0,10);
  h_detection_eff_alpha = new TH1I("Detection_Eff_alpha","Detection_Eff_alpha",10,0,10);

  h_mom_check = new TH2F("mom_check_momVStof","mom_check_momVStof",200,0.,20.,90,10.,40.);
  h_mom_check2 = new TH2F("mom_check_rhoVStof","mom_check_rhoVStof",200,0.,10000.,90,10.,40.);
  h_tofdE = new TH2F("tofdE","tofdE",120,10.,40.,250,0.,1500.);

  h_beta = new TH1F("beta","beta",300,0.,3.);
  h_beta2 = new TH1F("beta_pvcut_0.75","beta_pvcut_0.75",300,0.,3.);
  h_beta3 = new TH1F("beta_pvcut_0.5","beta_pvcut_0.5",300,0.,3.);

  h_beta_mom = new TH2F("beta_mom","beta_mom",500,0.,15.,200,0.,2.5);
  h_beta_mom2 = new TH2F("beta_mom_pvcut_0.75","beta_mom_pvcut_0.75",500,0.,15.,100,0.,2.5);
  h_beta_mom3 = new TH2F("beta_mom_pvcut_0.5","beta_mom_pvcut_0.5",500,0.,15.,100,0.,2.5);

  h_Mass_All= new TH1F("mass_all","mass_all",500,0,10.);
  h_Mass_All2= new TH1F("mass_all_pvcut_0.75","mass_all_pvcut_0.75",500,0,10.);
  h_Mass_All3= new TH1F("mass_all_pvcut_0.5","mass_all_pvcut_0.5",500,0,10.);

  h_Mass_charge_All  = new TH2F("mass_charge_all","mass_charge_all",200,0.,10.,10,-5,5);
  h_Mass_charge_All2 = new TH2F("mass_charge_all_pvcut_0.75","mass_charge_all_pvcut_0.75",200,0.,10.,10,-5,5);
  h_Mass_charge_All3 = new TH2F("mass_charge_all_pvcut_0.5","mass_charge_all_pvcut_0.5",200,0.,10.,10,-5,5);

  h_Path =  new TH1F("path_length","path_length",200,400.,800.);
  h_Path_Back =  new TH1F("path_length_back","path_length_back",200,400.,800.);
  h_MeanPath =  new TH1F("Mean_path_length","Mean_path_length",200,400.,800.);
  h_dpath = new TH1F("Dpath","Dpath",200,-10.,10.);
  h_dpath2_chi2 = new TH2F("Dpath2_chi2","Dpath2_chi2",200,-10.,10.,200,0,100);
  h_dpath2_pv = new TH2F("Dpath2_pv","Dpath2_pv",200,-10.,10.,100,0.,1.);

  h_Path_newL = new TH1F("path_length_NEW","path_length_NEW",400,400.,600.);
  h_Path_LTr1toTOF = new TH1F("path_length_NEW_Tr1TOF","path_length_NEW_Tr1TOF",400,300.,500.);
  h_Path_LTr0toTr1 = new TH1F("path_length_NEW_Tr0Tr1","path_length_NEW_Tr0Tr1",300,-10.,50.);
  h_Path_check1 = new TH1F("path_length_NEW_CheckTr1TOF","path_length_NEW_CheckTr1TOF",200,-600.,600.);
  h_Path_check2 = new TH1F("path_length_NEW_CheckTr0Tr1","path_length_NEW_CheckTr0Tr1",200,-600.,600.);

  h_dpath_new = new TH2F("Dpath_NEW_2d","Dpath_NEW_2d",800,-20.,20.,10,0,10);
  h_Time_Estimation = new TH2F("Time_NEW_2d","Time_NEW_2d",800,-20.,20.,10,0,10);
  h_Time = new TH1F("Time","Time",1000,0.,50.);


  h_path_tof = new TH2F("Mean_path_TOF","Mean_path_TOF",500,15.,20.,90,10,40);
  h_pv_mom = new TH2F("mom_pv","mom_pv",200,0.,20.,250,0,1);
  h_pv_beta = new TH2F("beta_pv","beta_pv",250,0.,2.5,250,0,1);
  h_pv_mass = new TH2F("mass_pv","mass_pv",200,0.,10.,250,0,1);
  h_mom_tof_cut = new TH2F("mom_tof_pvcut_0.5","mom_tof_pvcut_0.5",200,0.,20.,500,0.,50.);
  h_path_tof_cut = new TH2F("path_tof_pvcut_0.5","path_tof_pvcut_0.5",500,15.,20.,500,0.,50.);
  h_path_mom_cut = new TH2F("path_mom_pvcut","path_mom_pvcut",500,15.,20.,200,0,20.);

  // h_slope_yz[0] = new TH1F("slope_yz_no_fit","slope_yz_no_fit",10000,-1.,1.);
  // h_slope_yz[1] = new TH1F("slope_yz_fit_1","slope_yz_fit_1",10000,-1.,1.);
  // h_slope_yz[2] = new TH1F("slope_yz_fit_2","slope_yz_fit_2",10000,-1.,1.);
  // h_slope_yz[3] = new TH1F("slope_yz_fit_3","slope_yz_fit_3",10000,-1.,1.);

  //h_dy = new TH1F("delta_y","delta y on tof",10000,-1.,1.);

  // h_mass_invariant_cut_mass = new TH1F("mass_invariant","mass_invariant",14000,1.,5.3);
  // h_mass_invariant_best_cut_mass = new TH1F("mass_invariant_best_cut_mass","mass_invariant_best_cut_mass",14000,1.,5.3);
  // h_mass_invariant_cut_vertex = new TH1F("mass_invariant_cut_vertex","mass_invariant_cut_vertex",14000,1.,5.3);
  // h_mass_invariant_best_cut_vertex = new TH1F("mass_invariant_best_cut_vertex","mass_invariant_best_cut_vertex",14000,1.,5.3);

  // h_vtx_x = new TH1F("vtx_reso_x","vtx_reso_x",500,-25.,25.);
  // h_vtx_y = new TH1F("vtx_reso_y","vtx_reso_y",500,-25.,25.);
  // h_vtx_z = new TH1F("vtx_reso_z","vtx_reso_z",500,-150.,150.);

  // h_vtx_x_cut = new TH1F("vtx_reso_x_cut","vtx_reso_x_cut",500,-25.,25.);
  // h_vtx_y_cut = new TH1F("vtx_reso_y_cut","vtx_reso_y_cut",500,-25.,25.);
  // h_vtx_z_cut = new TH1F("vtx_reso_z_cut","vtx_reso_z_cut",500,-150.,150.);

  //kalman
  h_pv = new TH1F("total_pvalue","total_pvalue",200,0,1);
  h_chi2 = new TH1F("total chi2","total chi2",400,0,40);

  TString name_residu[5]={"TR1","TR2","DC2","TOF+","TFW"};
  TString name_residu_Q[4]={"_pi-","_q1","_q2","_q3"};
  
  for(int ii=0;ii<5;++ii)
    for(int cc=0;cc<4;++cc)
      {
	TString name_temp_resX = "h_residualX_";
	TString name_temp_resY = "h_residualY_";
	name_temp_resX += name_residu[ii];
	name_temp_resY += name_residu[ii];
	name_temp_resX += name_residu_Q[cc];
	name_temp_resY += name_residu_Q[cc];
	
	h_ResidualX[ii][cc] = new TH1F(name_temp_resX,name_temp_resX,2000,-20.,20.);
	h_ResidualY[ii][cc] = new TH1F(name_temp_resY,name_temp_resY,2000,-20.,20.);
      }

  //if(DAF == true)
  //{
      /* TString name_residu[7]={"TR1x","TR1y","TR2x","TR2y","DC1","DC2","TOF"};

      for(int ii=0;ii<7;++ii)
	for(int cc=0;cc<4;++cc)
	  {
	    TString name_temp_resX = "h_residualX_";
	    TString name_temp_resY = "h_residualY_";
	    name_temp_resX += name_residu[ii];
	    name_temp_resY += name_residu[ii];
	    name_temp_resX += name_residu_Q[cc];
	    name_temp_resY += name_residu_Q[cc];
	    
	    h_ResidualX[ii][cc] = new TH1F(name_temp_resX,name_temp_resX,2000,-20.,20.);
	    h_ResidualY[ii][cc] = new TH1F(name_temp_resY,name_temp_resY,2000,-20.,20.);
	  }*/
  {
      TString DetectorName[] = {"TR0X","TR0Y","TR1X","TR1Y","DC1XP","DC1X","DC1UP","DC1U","DC1V","DC1VP","TR2X","TR2Y","DC2U","DC2YP","DC2Y","DC2XP","DC2X","TOFP","TFWN"};
      TString name_residu_Q[4]={"_pi-","_q1","_q2","_q3"};
      
      for(int cc=0;cc<4;++cc)
	for(int ii = 0;ii<19;++ii)
	  {
	    TString name_temp_weight1("h_weight_max_");
	    name_temp_weight1+=DetectorName[ii];
	    name_temp_weight1+=name_residu_Q[cc];
	    TString name_temp_weight2("h_weight_min_");
	    name_temp_weight2+=DetectorName[ii];
	    name_temp_weight2+=name_residu_Q[cc];
	    TString name_temp_weight3("h_weight_sum_");
	    name_temp_weight3+=DetectorName[ii];
	    name_temp_weight3+=name_residu_Q[cc];
	    TString name_temp_weight4("h_weight_unique_");
	    name_temp_weight4+=DetectorName[ii];
	    name_temp_weight4+=name_residu_Q[cc];
	    TString name_temp_chi2("h_chi2_det_");
	    name_temp_chi2+=DetectorName[ii];
	    name_temp_chi2+=name_residu_Q[cc];
	    TString name_temp_chi2_1("h_chi2_det_unique_");
	    name_temp_chi2_1+=DetectorName[ii];
	    name_temp_chi2_1+=name_residu_Q[cc];
	    h_weigth_max[cc][ii] = new TH1F(name_temp_weight1,name_temp_weight1,200,0.,1.);
	    h_weigth_min[cc][ii] = new TH1F(name_temp_weight2,name_temp_weight2,200,0.,1.);
	    h_weigth_sum[cc][ii] = new TH1F(name_temp_weight3,name_temp_weight3,200,0.,1.);
	    h_weigth_one[cc][ii] = new TH1F(name_temp_weight4,name_temp_weight4,200,0.,1.);
	    h_chi2_det[cc][ii] = new TH1F(name_temp_chi2,name_temp_chi2,200.,0.,10.);
	    h_chi2_det_one[cc][ii] = new TH1F(name_temp_chi2_1,name_temp_chi2_1,200.,0.,10.);
	  }
  }
  //h_chi2_smooth = new TH1F("total chi2 smoothed","total chi2 smoothed",10000,0,100);
			     

  h_dB = new TH1F("h_dB","h_dB",100,-0.05,0.05);
  h_Bchi2 = new TH1F("h_B_chi2","h_B_chi2",100,0,0.05);
  h_Bz = new TH2F("B_z","B_z",320,-160.,160.,100,0.,1.);

  // std::string temp_str[5]={"x","y","tx","ty","q/p"};
  // for(int i=0;i<5;i++)
  // {
  //     std::string temp ="pull_";
  //     temp+=temp_str[i];
  //     h_pull[i]= new TH1F(temp.c_str(),temp.c_str(),1000,-10,10);
  // }
  if(Vertex == true)
    {
      for(int i=0;i<3;i++)
	{
	  std::string temp = "mass_invariant";
	  temp +=pv[i];
	  h_mass_invariant_cut_mass2[i] = new TH1F(temp.c_str(),temp.c_str(),14000,1.,5.3);
	  
	  temp="mass_invariant_best_cut_mass";
	  temp +=pv[i];
	  h_mass_invariant_best_cut_mass2[i] = new TH1F(temp.c_str(),temp.c_str(),14000,1.,5.3);
	  
	  temp="mass_invariant_cut_vertex";
	  temp +=pv[i];
	  h_mass_invariant_cut_vertex2[i] = new TH1F(temp.c_str(),temp.c_str(),14000,1.,5.3);
	  
	  temp="mass_invariant_best_cut_vertex";
	  temp +=pv[i];
	  h_mass_invariant_best_cut_vertex2[i] = new TH1F(temp.c_str(),temp.c_str(),14000,1.,5.3);
	  
	  for(int j=0;j<4;j++)
	    {
	      temp="mom_D1_mom_D2_cut_mass";
	      temp +=name_hyp[j];
	      temp +=pv[i];
	      h_momD1_momD2_cut_mass[i][j] = new TH2F(temp.c_str(),temp.c_str(),200,0.,1.2,200,0.,13.);
	      
	      temp="mom_D1_mom_D2_best_cut_mass";
	      temp +=name_hyp[j];
	      temp +=pv[i];
	      h_momD1_momD2_best_cut_mass[i][j] = new TH2F(temp.c_str(),temp.c_str(),200,0.,1.2,200,0.,13.);
	      
	      temp="mom_D1_mom_D2_cut_vertex";
	      temp +=name_hyp[j];
	      temp +=pv[i];
	      h_momD1_momD2_cut_vertex[i][j] = new TH2F(temp.c_str(),temp.c_str(),200,0.,1.2,200,0.,13.);
	      
	      temp="mom_D1_mom_D2_best_cut_vertex";
	      temp +=name_hyp[j];
	      temp +=pv[i];
	      h_momD1_momD2_best_cut_vertex[i][j] = new TH2F(temp.c_str(),temp.c_str(),200,0.,1.2,200,0.,13.);
	    }
	  
	  // temp="vtx_reso_x";
	  // temp +=pv[i];
	  // h_vtx_x2[i] = new TH1F(temp.c_str(),temp.c_str(),500,-2.5,2.5);
	  
	  // temp="vtx_reso_y";
	  // temp +=pv[i];
	  // h_vtx_y2[i] = new TH1F(temp.c_str(),temp.c_str(),500,-2.5,2.5);
	  
	  // temp="vtx_reso_z";
	  // temp +=pv[i];
	  // h_vtx_z2[i] = new TH1F(temp.c_str(),temp.c_str(),500,-15.,15.);

	  // temp="vtx_reso_x_cut";
	  // temp +=pv[i];
	  // h_vtx_x_cut2[i] = new TH1F(temp.c_str(),temp.c_str(),500,-2.5,2.5);
	  
	  // temp="vtx_reso_y_cut";
	  // temp +=pv[i];
	  // h_vtx_y_cut2[i] = new TH1F(temp.c_str(),temp.c_str(),500,-2.5,2.5);
	  
	  // temp="vtx_reso_z_cut";
	  // temp +=pv[i];
	  // h_vtx_z_cut2[i] = new TH1F(temp.c_str(),temp.c_str(),500,-15.,15.);
      

	  temp="chi2_cut";
	  temp +=pv[i];
	  h_chi2_cut[i] = new TH1F(temp.c_str(),temp.c_str(),200,0.,20);
	  
	  temp="chi2_cut_best_cut_mass";
	  temp +=pv[i];
	  h_chi2_cut_best_cut_mass[i] = new TH1F(temp.c_str(),temp.c_str(),200,0.,20);
	  
	  temp="chi2_cut_cut_mass";
	  temp +=pv[i];
	  h_chi2_cut_cut_mass[i] = new TH1F(temp.c_str(),temp.c_str(),200,0.,20);
	  
	  temp="chi2_cut_cut_vtx";
	  temp +=pv[i];
	  h_chi2_cut_cut_vtx[i] = new TH1F(temp.c_str(),temp.c_str(),200,0.,20);

	  temp="chi2_cut_best_cut_vtx";
	  temp +=pv[i];
	  h_chi2_cut_best_cut_vtx[i] = new TH1F(temp.c_str(),temp.c_str(),200,0.,20);
	  
	}
    }

  h_mass_xTof = new TH2F("mass_xTof","mass_xTof",500,0.,25.,200,-800.,200.);
  h_mom_xTof  = new TH2F("mom_xTof","mom_xTof",500,0.,15.,200,-800.,200.);
  h_beta_xTof = new TH2F("beta_xTof","beta_xTof",400,0.,2.,200,-800.,200.);
  h_chi2_xTof = new TH2F("chi2_xTof","chi2_xTof",500,0.,15.,200,-800.,200.);
  h_pv_xTof = new TH2F("pv_xTof","pv_xTof",200,0.,1.,200,-800.,200.);

  if(Finding == true)
    {
      TString part[7]={"pi-","pi+","K+","alpha","proton","He3","deuteron"};
      TString ext[3]={"_","_noFound_","_Ghost_"};
      for(int i=0;i<7;++i)
	{
	  for(int j=0;j<3;++j)
	    {
	      TString name_temp("find_chi2_tot_pid");
	      name_temp+=ext[j];
	      name_temp+=part[i];
	      h_finding_chi2_tot_pid[i][j] = new TH1F(name_temp,name_temp,10000,0.,1000.);
	      
	      name_temp="find_chi2_x_pid";
	      name_temp+=ext[j];
	      name_temp+=part[i];
	      h_finding_chi2_x_pid[i][j] = new TH1F(name_temp,name_temp,10000,0.,1000.);
	      
	      name_temp="find_chi2_y_pid";
	      name_temp+=ext[j];
	      name_temp+=part[i];
	      h_finding_chi2_y_pid[i][j] = new TH1F(name_temp,name_temp,10000,0.,1000.);
	    }
	  TString name_temp("find_chi2_2D_");
	  name_temp+=part[i];
	  h_finding_chi2_2D[i] = new TH2F(name_temp,name_temp,1000,-1.,99.,1000,-1,99.);
	  name_temp="find_chi2_2D_Ghost_";
	  name_temp+=part[i];
	  h_finding_chi2_2D_ghost[i] = new TH2F(name_temp,name_temp,1000,-1.,199.,1000,-1,199.);
	  
	}
      
      h_finding_eff = new TH2F("find_eff_stat","find_eff_stat",10,0.,10.,10,0,10);
      h_finding_multi = new TH2F("find_multi_stat","find_multi_stat",20,0.,20.,100,0,100);
    }
#endif

  //Daisuke
  //#ifdef DAISUKE
  h_status = new TH1F("AnalysisStatus","AnalysiStatus",100,0,100);

  h_MassZ1BarBest = new TH2F("AAD_MassZ1Bar","AAD_MassZ1Bar",1000,0,10,35,0,35);//2000
  h_MassZ1BarBestTofs = new TH2F("AAD_MassZ1Bar_Tofs","AAD_MassZ1Bar_Tofs",1000,0,10,35,0,35);//2000
  //h_MassZ1BarBest_ini = new TH2F("AAD_MassZ1Bar_ini","AAD_MassZ1Bar_ini",1000,0,10,35,0,35);//2000
  //h_MassZ1BarBest_iniTofs = new TH2F("AAD_MassZ1Bar_iniTofs","AAD_MassZ1Bar_iniTofs",1000,0,10,35,0,35);//2000

//   h_MassZ1BarBest3 = new TH2F("AAD_MassZ1Bar3","AAD_MassZ1Bar3",2000,0,10,35,0,35);
//   h_MassZ1BarBestTofs3 = new TH2F("AAD_MassZ1Bar_Tofs3","AAD_MassZ1Bar_Tofs3",2000,0,10,35,0,35);
//   h_MassZ1BarBest_ini3 = new TH2F("AAD_MassZ1Bar_ini3","AAD_MassZ1Bar_ini3",2000,0,10,35,0,35);
//   h_MassZ1BarBest_iniTofs3 = new TH2F("AAD_MassZ1Bar_iniTofs3","AAD_MassZ1Bar_iniTofs3",2000,0,10,35,0,35);

  //h_MassZ2BarBest = new TH2F("AAD_MassZ2Bar","AAD_MassZ2Bar",1000,0,10,35,0,35);//2000
  //h_MassZ2BarBestTofs = new TH2F("AAD_MassZ2Bar_Tofs","AAD_MassZ2Bar_Tofs",1000,0,10,35,0,35);//2000
  //h_MassZ2BarBest_ini = new TH2F("AAD_MassZ2Bar_ini","AAD_MassZ2Bar_ini",1000,0,10,35,0,35);//2000
  //h_MassZ2BarBest_iniTofs = new TH2F("AAD_MassZ2Bar_iniTofs","AAD_MassZ2Bar_iniTofs",1000,0,10,35,0,35);//2000

//   h_MassZ2BarBest3 = new TH2F("AAD_MassZ2Bar3","AAD_MassZ2Bar3",2000,0,10,35,0,35);
//   h_MassZ2BarBestTofs3 = new TH2F("AAD_MassZ2Bar_Tofs3","AAD_MassZ2Bar_Tofs3",2000,0,10,35,0,35);
//   h_MassZ2BarBest_ini3 = new TH2F("AAD_MassZ2Bar_ini3","AAD_MassZ2Bar_ini3",2000,0,10,35,0,35);
//   h_MassZ2BarBest_iniTofs3 = new TH2F("AAD_MassZ2Bar_iniTofs3","AAD_MassZ2Bar_iniTofs3",2000,0,10,35,0,35);

  //h_MassZ1BarBest_pi = new TH2F("AAD_MassZ1Bar_pi","AAD_MassZ1Bar_pi",1000,0,5,20,0,20);
  //h_MassZ1BarBest_piTofs = new TH2F("AAD_MassZ1Bar_pi_Tofs","AAD_MassZ1Bar_pi_Tofs",1000,0,5,20,0,20);
  //h_MassZ1BarBest_ini_pi = new TH2F("AAD_MassZ1Bar_ini_pi","AAD_MassZ1Bar_ini_pi",1000,0,5,20,0,20);
  //h_MassZ1BarBest_ini_piTofs = new TH2F("AAD_MassZ1Bar_ini_piTofs","AAD_MassZ1Bar_ini_piTofs",1000,0,5,20,0,20);

//   h_MassZ1BarBest_pi3 = new TH2F("AAD_MassZ1Bar_pi3","AAD_MassZ1Bar_pi3",1000,0,5,20,0,20);
//   h_MassZ1BarBest_piTofs3 = new TH2F("AAD_MassZ1Bar_pi_Tofs3","AAD_MassZ1Bar_pi_Tofs3",1000,0,5,20,0,20);
//   h_MassZ1BarBest_ini_pi3 = new TH2F("AAD_MassZ1Bar_ini_pi3","AAD_MassZ1Bar_ini_pi3",1000,0,5,20,0,20);
//   h_MassZ1BarBest_ini_piTofs3 = new TH2F("AAD_MassZ1Bar_ini_piTofs3","AAD_MassZ1Bar_ini_piTofs3",1000,0,5,20,0,20);

  hd_TrackFrom[0] = new TH2F("ADD_Track_vs_Tofs","ADD_Track_vs_Tofs",15,0,15,400,-200,200);
  hd_TrackFrom[1] = new TH2F("ADD_Track_vs_Tofs_pi","ADD_Track_vs_Tofs_pi",15,0,15,400,-200,200);
  h_MomCor[0]=new TH2F("AAD_MomCor_X","AAD_MomCor_X",1000,0,5,1000,0,5);
  h_MomCor[1]=new TH2F("AAD_MomCor_Y","AAD_MomCor_Y",1000,0,5,1000,0,5);
  h_MomCor[2]=new TH2F("AAD_MomCor_Z","AAD_MomCor_Z",1000,0,5,1000,0,5);
  
  hd_theWorld = new TH2F("AAD_theWorld","AAD_theWorld",300,-1000,5000,300,-3000,3000);//3000,3000
  hd_chi[0] = new TH1F("AAD_Chi2_TOFp","AAD_Chi2_TOFp",500,0,100);
  hd_chi[1] = new TH1F("AAD_Chi2_BigTof","AAD_Chi2_BigTof",500,0,100);
  hd_pv[0] = new TH1F("AAD_Pval_TOFp","AAD_Pval_TOFp",500,0,1);
  hd_pv[1] = new TH1F("AAD_Pval_BigTof","AAD_Pval_BigTof",500,0,1);

  hd_MassPv[0]=new TH2F("AAD_Mass_Pv","AAD_Mass_Pv",1000,0,5,1000,0,1);
  hd_MassPv[1]=new TH2F("AAD_Mass_Pv_pi","AAD_Mass_Pv_pi",1000,0,5,1000,0,1);



  char tCstr[300];
  


  sprintf(tCstr,"AAD_BetaMom_TofP_bar_all");
  hd_BetaMom_Tp_all = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000
  sprintf(tCstr,"AAD_BetaMom_TofP_bar_allTofs");
  hd_BetaMom_Tp_allTofs = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000

  sprintf(tCstr,"AAD_BetaMom_TofP_bar_allTofs_Z2");
  hd_BetaMom_Tp_allTofs_Z2 = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000

  sprintf(tCstr,"AAD_BetaMom_TofP_bar_ini_all");
  hd_BetaMom_Tp_ini_all = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000
  sprintf(tCstr,"AAD_BetaMom_TofP_bar_ini_allTofs");
  hd_BetaMom_Tp_ini_allTofs = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000
  sprintf(tCstr,"AAD_BetaMom_TofP_bar_kal_allTofs");
  hd_BetaMom_Tp_kal_allTofs = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000

  sprintf(tCstr,"AAD_BetaMom_TofP_bar_ini_allTofs_Z2");
  hd_BetaMom_Tp_ini_allTofs_Z2 = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000

  //hd_Mom_ini_vs_rchi_ini= new TH2F("AAD_Mom_ini_vs_rchi_ini","AAD_Mom_ini_vs_rchi_ini",1000,0,20,1000,0,20);
  //hd_Momini_rchiy_Tofs= new TH2F("AAD_Momini_rchiy_tofs","AAD_Momini_rchiy_tofs",1000,0,20,1000,0,20);

  sprintf(tCstr,"AAD_BetaMom_BigTof_bar_all");
  hd_BetaMom_Bg_all = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000
  sprintf(tCstr,"AAD_BetaMom_BigTof_bar_allTofs");
  hd_BetaMom_Bg_allTofs = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000

  sprintf(tCstr,"AAD_BetaMom_BigTof_bar_ini_all");
  hd_BetaMom_Bg_ini_all = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000
  sprintf(tCstr,"AAD_BetaMom_BigTof_bar_ini_allTofs");
  hd_BetaMom_Bg_ini_allTofs = new TH2F(tCstr,tCstr,1000,0,20,400,0,2);//2000

  for(int bar=0;bar<32;bar++)
    {
//       sprintf(tCstr,"AAD_BetaMom_TofP_bar%i",bar);
//       hd_BetaMom_Tp[bar] = new TH2F(tCstr,tCstr,2000,0,20,400,0,2);
      sprintf(tCstr,"AAD_TimeCalcDiff_TofP_bar%i",bar);
      hd_TimeCalTp[bar] = new TH1F(tCstr,tCstr,200,-10,10);
    }
//   for(int bar=0;bar<18;bar++)
//     {
//       sprintf(tCstr,"AAD_BetaMom_BigTof_bar%i",bar);
//       hd_BetaMom_Bg[bar] = new TH2F(tCstr,tCstr,1000,0,10,400,0,2);
//     }
  //#endif

//     for(int mm=0;mm<3;mm++)
//       for(int dd=0;dd<3;dd++)
// 	{
// 	  sprintf(tCstr,"AAD_SecZX_all_pv%i_mother%i",dd,mm);
// 	  hd_secz_all[mm][dd] = new TH2F(tCstr,tCstr,1000,-300,700,200,-100,100);
// 	  sprintf(tCstr,"AAD_SecZX_pv%i_mother%i",dd,mm);
// 	  hd_secz[mm][dd] = new TH2F(tCstr,tCstr,1000,-300,700,200,-100,100);
// 	  sprintf(tCstr,"AAD_SecZX_cut_pv%i_mother%i",dd,mm);
// 	  hd_secz_cut[mm][dd] = new TH2F(tCstr,tCstr,1000,-300,700,200,-100,100);
// 	  sprintf(tCstr,"AAD_SecZX_cut2_pv%i_motehr%i",dd,mm);
// 	  hd_secz_cut2[mm][dd] = new TH2F(tCstr,tCstr,1000,-300,700,200,-100,100);
// 	  sprintf(tCstr,"AAD_SecZX_cut3_pv%i_mother%i",dd,mm);
// 	  hd_secz_cut3[mm][dd] = new TH2F(tCstr,tCstr,1000,-300,700,200,-100,100);
// 	  sprintf(tCstr,"AAD_SecZX_mompi_pv%i_mother%i",dd,mm);
// 	  hd_secz_mompi[mm][dd] = new TH2F(tCstr,tCstr,1000,-300,700,200,0,5);
// 	}

//   for(int bar=0;bar<4;bar++)
//     for(int jj=0;jj<3;jj++)
//       {
// 	sprintf(tCstr,"AAD_InvMassLambda%i_pv%i",bar,jj);
// 	hd_invmass[0][jj][bar] = new TH1F(tCstr,tCstr,2000,1.0,1.5);
	
// 	sprintf(tCstr,"AAD_InvMassLambdaSecZ%i_pv%i",bar,jj);
// 	hd_invmass_secz[0][jj][bar] = new TH2F(tCstr,tCstr,2000,1.0,1.5,1000,-300,700);

// 	sprintf(tCstr,"AAD_InvMassH3L%i_pv%i",bar,jj);
// 	hd_invmass[1][jj][bar] = new TH1F(tCstr,tCstr,2000,2.8,3.3);
	
// 	sprintf(tCstr,"AAD_InvMassH3LSecZ%i_pv%i",bar,jj);
// 	hd_invmass_secz[1][jj][bar] = new TH2F(tCstr,tCstr,2000,2.8,3.3,1000,-300,700);

// 	sprintf(tCstr,"AAD_InvMassH4L%i_pv%i",bar,jj);
// 	hd_invmass[2][jj][bar] = new TH1F(tCstr,tCstr,2000,3.8,4.3);
	
// 	sprintf(tCstr,"AAD_InvMassH4LSecZ%i_pv%i",bar,jj);
// 	hd_invmass_secz[2][jj][bar] = new TH2F(tCstr,tCstr,2000,3.8,4.3,1000,-300,700);

//       }
  


  


  //hd_PoQ_E = new TH2F("AAD_PoQ_E","AAD_PoQ_E",1000,0,20,600,0,300);//2000,600
  //hd_PoQ_E_ini = new TH2F("AAD_PoQ_E_ini","AAD_PoQ_E_ini",1000,0,20,600,0,300);//2000,600

  //hd_PoQ_ETofs = new TH2F("AAD_PoQ_E_tofs","AAD_PoQ_E_tofs",1000,0,20,600,0,300);//2000,600
  //hd_PoQ_E_iniTofs = new TH2F("AAD_PoQ_E_ini_tofs","AAD_PoQ_E_ini_tofs",1000,0,20,600,0,300);//2000,600

//   hd_PoQ_E3 = new TH2F("AAD_PoQ_E3","AAD_PoQ_E3",2000,0,20,1200,0,300);
//   hd_PoQ_E_ini3 = new TH2F("AAD_PoQ_E_ini3","AAD_PoQ_E_ini3",2000,0,20,1200,0,300);

//   hd_PoQ_ETofs3 = new TH2F("AAD_PoQ_E_tofs3","AAD_PoQ_E_tofs3",2000,0,20,1200,0,300);
//   hd_PoQ_E_iniTofs3 = new TH2F("AAD_PoQ_E_ini_tofs3","AAD_PoQ_E_ini_tofs3",2000,0,20,1200,0,300);

  //hd_IniMom_NewMom_all = new TH2F("AAD_X_IniMom_Y_NewMom_all","AAD_X_IniMom_Y_NewMom_all",1000,0,100,1000,0,100);
  //hd_IniMom_NewMom_all_tofs = new TH2F("AAD_X_IniMom_Y_NewMom_all_tofs","AAD_X_IniMom_Y_NewMom_all_tofs",1000,0,100,1000,0,100);
  //hd_IniMom_NewMom_all_tofs_ini = new TH2F("AAD_X_IniMom_Y_NewMom_all_tofs_ini","AAD_X_IniMom_Y_NewMom_all_tofs_ini",1000,0,100,1000,0,100);
  hd_IniMom_NewMom_all_tofsZ1 = new TH2F("AAD_X_IniMom_Y_NewMom_all_tofsZ1","AAD_X_IniMom_Y_NewMom_all_tofsZ1",1000,0,100,1000,0,100);
  hd_IniMom_NewMom_all_tofsZ1_ini = new TH2F("AAD_X_IniMom_Y_NewMom_all_tofsZ1_ini","AAD_X_IniMom_Y_NewMom_all_tofsZ1_ini",1000,0,100,1000,0,100);
  //hd_IniMom_DiffMom_allTofs = new TH2F("AAD_X_IniMom_Y_DiffMom_allTofs","AAD_X_IniMom_Y_DiffMom_allTofs",1000,0,100,1000,-100,100);
  //hd_NewMom_DiffMom_allTofs = new TH2F("AAD_X_NewMom_Y_DiffMom_allTofs","AAD_X_NewMom_Y_DiffMom_allTofs",1000,0,100,1000,-100,100);
  //hd_Bar_DiffMom_allTofs = new TH2F("AAD_X_Bar_Y_DiffMom_allTofs","AAD_X_Bar_Y_DiffMom_allTofs",32,0,32,1000,-100,100);

  hd_IniNewMomdiff_vs_IniNewMomplus_allTofs = new TH2F("AAD_IniNewMomdiff_vs_IniNewMomplus_allTofs","AAD_IniNewMomdiff_vs_IniNewMomplus_allTofs",1000,-50,50,1000,0,100);
  
  hd_charge_ini_vs_IniNewDiv = new TH2F("AAD_charge_ini_vs_IniNewDiv","AAD_charge_ini_vs_IniNewDiv",10,0,10,1000,0,50);
  hd_charge_vs_IniNewDiv = new TH2F("AAD_charge_vs_IniNewDiv","AAD_charge_vs_IniNewDiv",10,0,10,1000,0,50);
  
  hd_IniNewMomdiff_vs_IniNewMomplus_allTofs_Z1 = new TH2F("IniNewMomdiff_vs_IniNewMomplus_allTofs_Z1","IniNewMomdiff_vs_IniNewMomplus_allTofs_Z1",1000,-20,20,1000,0,50);
  hd_IniNewMomdiff_vs_IniNewMomplus_allTofs_ini_Z1 = new TH2F("IniNewMomdiff_vs_IniNewMomplus_allTofs_ini_Z1","IniNewMomdiff_vs_IniNewMomplus_allTofs_ini_Z1",1000,-20,20,1000,0,50);

  hd_best_pos_hitvstrack_tofpZ1 = new TH2F("best_pos_hitvstrack_tofpZ1","best_pos_hitvstrack_tofpZ1",1000,-1000,500,1000,-1000,500);
  hd_best_pos_hitdifftrack_tofpZ1 = new TH1F("best_pos_hitdifftrack_tofpZ1","best_pos_hitdifftrack_tofpZ1",1000,-500,500);

  //hd_IniMom_NewMom_all = new TH2F("AAD_X_IniMom_Y_NewMom_all","AAD_X_IniMom_Y_NewMom_all",1000,0,100,1000,0,100);
  //hd_IniMom_DiffMom_all = new TH2F("AAD_X_IniMom_Y_DiffMom_all","AAD_X_IniMom_Y_DiffMom_all",1000,0,100,1000,-100,100);
  //hd_NewMom_DiffMom_all = new TH2F("AAD_X_NewMom_Y_DiffMom_all","AAD_X_NewMom_Y_DiffMom_all",1000,0,100,1000,-100,100);
  //hd_Bar_DiffMom_all = new TH2F("AAD_X_Bar_Y_DiffMom_all","AAD_X_Bar_Y_DiffMom_all",32,0,32,1000,-100,100);



  h_MC_PxPy = new TH2F("MC_PxPy","MC_PxPy",200,-1.,1.,200,-1.,1.); 
  h_MC_Pt = new TH1F("MC_Pt","MC_Pt",200,0,0.4); 
  h_MC_Pt_y= new TH2F("MC_Pt_y","MC_Pt_y",400,-2,6,200,0.,0.4); 
  h_MC_TargetXZ = new TH2F("MC_TargetXZ","MC_TargetXZ",200,-5.,45.,200,-6.,6.); 
  h_MC_TargetYZ = new TH2F("MC_TargetYZ","MC_TargetYZ",200,-5.,45.,200,-6.,6.); 
  h_MC_stats = new TH1F("MC_stats","MC_stats",10,0,100); 
  
  h_MC_dE_mom   = new TH2F("MC_dE_mom","MC_dE_mom",1000,-5.,25.,200,0.,20.);
  h_MC_time_mom = new TH2F("MC_time_mom","MC_time_mom",1000,-5.,25.,200,0.,20.);
  h_MC_beta_mom = new TH2F("MC_beta_mom","MC_beta_mom",1000,-5.,25.,200,0.,5.);
  h_MC_length = new TH2F("MC_length","MC_length",6,-2,4,700,0.,700.);
  h_MC_length_Z = new TH2F("MC_length_Z","MC_length_Z",200,-5.,45.,14000,0.,700.);
  

  TObjArray* L_vol = gGeoManager->GetListOfVolumes();
  int n_volume = L_vol->GetEntries();
  //double max_bin_XX0[] = {10.,1000.,100.,1000.,1000.,1000.,1000.,0.1,10.,10.,10.,50.,10.,1.,10.,10.,1000.,50.,10.,1000.,1000.,1000.,2000.,100.,100.,2000.};
  for(int n_v = 0;n_v<n_volume;++n_v)
    {
      TString name_vol_temp(L_vol->At(n_v)->GetName());

      TString name_temp_vX("Mat_XX0_");
      name_temp_vX+=name_vol_temp;
      Material_XX0_y.push_back(new TH2F(name_temp_vX,name_temp_vX,200,-1,4,1000,0,1));

      TString name_temp_vE("Mat_DE_");
      name_temp_vE+=name_vol_temp;
      Material_dE_y.push_back(new TH2F(name_temp_vE,name_temp_vE,200,-1,4,1000,0,1.));
    }

  Material_XX0_y.push_back(new TH2F("Mat_XX0_Total","Mat_XX0_Total",200,-1,4,1000,0,1));
  
  Material_dE_y.push_back(new TH2F("Mat_DE_Total","Mat_DE_Total",200,-1,4,1000,0,1));


  std::cout<<" : done !"<<std::endl;

}
int Ana_Hist::Write(TFile* out_file)
{
  HaveBeenWritten=true;
  out_file->cd();
  //h_2TR1->Write();
  //h_2TR2->Write();
  //h_2AladinTOF->Write();
  TDirectory* det;
  TDirectory* Chamb;
  TDirectory* det_pull;
  TDirectory* mom1;
  TDirectory* mass;
  TDirectory* hough;
  TDirectory* find;
  TDirectory* find_part[7];
  TDirectory* kalman;
  TDirectory* particle;
  TDirectory* pv[3];
  TDirectory* momD1D2[3][4];
  TDirectory* Field;
  TDirectory* Daisuke;
  TDirectory* Daf;
  TDirectory* MC;
  TDirectory* Mat;

  TString part[7]={"pi-","pi+","K+","alpha","proton","He3","deuteron"};  
  if(!out_file->GetDirectory("Detector"))
    {
      std::cout<<"making directory "<<std::endl;
      det = out_file->mkdir("Detector");
      // AVhough = out_file->mkdir("AVHough");
      // AVHough_sub[10];

      Chamb = out_file->mkdir("ChamberProject");

      
  
      // for(int i=0;i<10;++i)
      //   {
      //     TString name("Event_");
      //     name+=i;
      //     AVHough_sub[i] = AVhough->mkdir(name);
      //   }

       det_pull = det->mkdir("pull");
       mom1 = out_file->mkdir("Momentum_Check");
      // mom2 = out_file->mkdir("Momentum_after_mass");
      // mom3 = out_file->mkdir("Momentum_after_vertex");
       mass = out_file->mkdir("Mass");
      // vertex = out_file->mkdir("Vertex");
       hough = out_file->mkdir("Hough");
       find = out_file->mkdir("Finding");

      for(int i=0;i<7;++i)
	find_part[i]=(TDirectory*) find->mkdir(part[i]);

      kalman = out_file->mkdir("Kalman");

      Daf = out_file->mkdir("DAF");

      particle= out_file->mkdir("Particle");
      pv[0] = (TDirectory*)out_file->mkdir("pv_0.01");
      pv[1] = (TDirectory*)out_file->mkdir("pv_0.05");
      pv[2] = (TDirectory*)out_file->mkdir("pv_0.1");
  
      for(int i=0;i<3;i++)
	{
	  momD1D2[i][0]= (TDirectory*) pv[i]->mkdir("L");
	  momD1D2[i][1]= (TDirectory*) pv[i]->mkdir("H3L");
	  momD1D2[i][2]= (TDirectory*) pv[i]->mkdir("H4L");
	  momD1D2[i][3]= (TDirectory*) pv[i]->mkdir("He5L");
	}
      Field = out_file->mkdir("Field");

      Daisuke = out_file->mkdir("Daisuke");
      MC = out_file->mkdir("MC");
      Mat = out_file->mkdir("Material");
    }
  else
    {
      std::cout<<"get directory "<<std::endl;
//       out_file->Write("",TObject::kOverwrite);
//       return 0;
       det = out_file->GetDirectory("Detector");
       if(det->cd())
	 {
	   std::cout<<"det "<<endl;
	 }
       else
	 {
	   std::cout<<"else det "<<endl;
	 }
       out_file->cd();
      // AVhough = out_file->GetDirectory("AVHough");
      // AVHough_sub[10];

       Chamb = out_file->GetDirectory("ChamberProject");
       
       
  
      // for(int i=0;i<10;++i)
      //   {
      //     TString name("Event_");
      //     name+=i;
      //     AVHough_sub[i] = AVhough->GetDirectory(name);
      //   }

       det_pull = det->GetDirectory("pull");
       mom1 = out_file->GetDirectory("Momentum_Check");
      // mom2 = out_file->GetDirectory("Momentum_after_mass");
      // mom3 = out_file->GetDirectory("Momentum_after_vertex");
       mass = out_file->GetDirectory("Mass");
      // vertex = out_file->GetDirectory("Vertex");
       hough = out_file->GetDirectory("Hough");
       find = out_file->GetDirectory("Finding");
       //      TString part[7]={"pi-","pi+","K+","alpha","proton","He3","deuteron"};
      for(int i=0;i<7;++i)
	find_part[i]=(TDirectory*) find->GetDirectory(part[i]);

       kalman = out_file->GetDirectory("Kalman");
       Daf = out_file->GetDirectory("DAF");
       particle= out_file->GetDirectory("Particle");
       pv[0] = out_file->GetDirectory("pv_0.01");
       pv[1] = out_file->GetDirectory("pv_0.05");
       pv[2] = out_file->GetDirectory("pv_0.1");
       for(int i=0;i<3;i++)
	 {
	   momD1D2[i][0]=  pv[i]->GetDirectory("L");
	   momD1D2[i][1]=  pv[i]->GetDirectory("H3L");
	   momD1D2[i][2]=  pv[i]->GetDirectory("H4L");
	   momD1D2[i][3]=  pv[i]->GetDirectory("He5L");
	 }
       Field = out_file->GetDirectory("Field");
       Daisuke = out_file->GetDirectory("Daisuke");       
       MC = out_file->GetDirectory("MC");
       Mat = out_file->GetDirectory("Material");
    }

  out_file->cd();

#ifdef ADD_HISTO

  h_task_exit->Write("",TObject::kOverwrite);
  
  h_alpha_rejection->Write("",TObject::kOverwrite);
  det->cd();
  for (Int_t n=0;n<21;n++)
    {
      h_2utrPos_XZ[n]->Write("",TObject::kOverwrite);
      h_2utrPos_XZ[n]->SetDirectory(det);
      h_2utrPos_YZ[n]->Write("",TObject::kOverwrite);
      h_2utrPos_YZ[n]->SetDirectory(det);
      h_2utrPos_XY[n]->Write("",TObject::kOverwrite);
      h_2utrPos_XY[n]->SetDirectory(det);
      h_utrEnergy[n]->Write("",TObject::kOverwrite);
      h_utrEnergy[n]->SetDirectory(det);
      h_utrTime[n]->Write("",TObject::kOverwrite);
      h_utrTime[n]->SetDirectory(det);
    }
  h_patern_decay->SetOption("text");
  h_patern_decay->SetMarkerSize(1.6);
  h_patern_decay->Write("",TObject::kOverwrite);
  h_2utrRotated_XZ->Write("",TObject::kOverwrite);
  h_2utrFrameMag_XZ->Write("",TObject::kOverwrite);
  h_2utrNotRotated_XZ->Write("",TObject::kOverwrite);

  h_TOFP_res->Write();
  h_TOFP_diff_ResInter->Write();
  h_TOFP_Rdiff_ResInter->Write();
  for(int j= 0;j<96;++j)
    {
      h_TOFP_res_Path[j]->Write();
    }

  det_pull->cd();
  for (Int_t n=0;n<18;n++)
    {
      h_pull_X_det[n]->Write("",TObject::kOverwrite);
      h_pull_X_det[n]->SetDirectory(det_pull);
      h_pull_Y_det[n]->Write("",TObject::kOverwrite);
      h_pull_Y_det[n]->SetDirectory(det_pull);
    }
  
  for (int i=0;i<4;i++)
    {
      out_file->cd();
      h_chi2_particle[i]->Write("",TObject::kOverwrite);
      h_pv_particle[i]->Write("",TObject::kOverwrite);

      // mom1->cd();
      // h_dmom1[i]->SetDirectory(mom1);
      // h_dmom1[i]->Write("",TObject::kOverwrite);
      // h_dmom2[i]->SetDirectory(mom1);
      // h_dmom2[i]->Write("",TObject::kOverwrite);
      // h_dmom3[i]->SetDirectory(mom1);
      // h_dmom3[i]->Write("",TObject::kOverwrite);
      // h_dmom4[i]->SetDirectory(mom1);
      // h_dmom4[i]->Write("",TObject::kOverwrite);
      // mom2->cd();
      // h_dmom_dmass[i]->SetDirectory(mom2);
      // h_dmom_dmass[i]->Write("",TObject::kOverwrite);
      // h_dmom_dmass_best[i]->SetDirectory(mom2);
      // h_dmom_dmass_best[i]->Write("",TObject::kOverwrite);
      // mom3->cd();
      // h_dmom_vertex[i]->SetDirectory(mom3);
      // h_dmom_vertex[i]->Write("",TObject::kOverwrite);
      // h_dmom_vertex_best[i]->SetDirectory(mom3);
      // h_dmom_vertex_best[i]->Write("",TObject::kOverwrite);

      mass->cd();
      h_Mass[i]->SetDirectory(mass);
      h_Mass[i]->Write("",TObject::kOverwrite);
	
    }

  h_Mass_All->Write("",TObject::kOverwrite);
  h_Mass_All2->Write("",TObject::kOverwrite);
  h_Mass_All3->Write("",TObject::kOverwrite);

  h_Mass_charge_All->Write("",TObject::kOverwrite);
  h_Mass_charge_All2->Write("",TObject::kOverwrite);
  h_Mass_charge_All3->Write("",TObject::kOverwrite);

  mom1->cd();
  h_mom_check->Write("",TObject::kOverwrite);
  h_mom_check->SetDirectory(mom1);
  h_mom_check2->Write("",TObject::kOverwrite);
  h_mom_check2->SetDirectory(mom1);
  h_tofdE->Write("",TObject::kOverwrite);
  h_tofdE->SetDirectory(mom1);

  out_file->cd();
  // for (int i=0;i<4;i++)
  //   h_slope_yz[i]->Write("",TObject::kOverwrite);


  h_beta->Write("",TObject::kOverwrite);
  h_beta2->Write("",TObject::kOverwrite);
  h_beta3->Write("",TObject::kOverwrite);

  h_beta_mom->Write("",TObject::kOverwrite);
  h_beta_mom2->Write("",TObject::kOverwrite);
  h_beta_mom3->Write("",TObject::kOverwrite);

  kalman->cd();
  h_pv->Write("",TObject::kOverwrite);
  h_chi2->Write("",TObject::kOverwrite);
  h_Path->Write("",TObject::kOverwrite);
  h_Path_Back->Write("",TObject::kOverwrite);
  h_MeanPath->Write("",TObject::kOverwrite);
  h_path_tof->Write("",TObject::kOverwrite);
  h_pv_mom->Write("",TObject::kOverwrite);
  h_pv_beta->Write("",TObject::kOverwrite);
  h_pv_mass->Write("",TObject::kOverwrite);
  h_mom_tof_cut->Write("",TObject::kOverwrite);
  h_path_tof_cut->Write("",TObject::kOverwrite);
  h_dpath->Write("",TObject::kOverwrite);
  h_dpath2_chi2->Write("",TObject::kOverwrite);
  h_dpath2_pv->Write("",TObject::kOverwrite);
  h_path_mom_cut->Write("",TObject::kOverwrite);


  h_Path_newL->Write("",TObject::kOverwrite);
  h_Path_LTr1toTOF->Write("",TObject::kOverwrite);
  h_Path_LTr0toTr1->Write("",TObject::kOverwrite);
  h_Path_check1->Write("",TObject::kOverwrite);
  h_Path_check2->Write("",TObject::kOverwrite);

  h_dpath_new->Write("",TObject::kOverwrite);

  h_Time_Estimation->Write("",TObject::kOverwrite);
  h_Time->Write("",TObject::kOverwrite);


  if(EnableState[DAF] == true)
    {
      //kalman->cd();
      //for(int ii=0;ii<7;++ii)
      //for(int cc=0;cc<4;++cc)
      //{
      //h_ResidualX[ii][cc]->Write("",TObject::kOverwrite);
      //h_ResidualY[ii][cc]->Write("",TObject::kOverwrite);
      //}
      Daf->cd();
      for(int cc=0;cc<4;++cc)
	for(int ii = 0;ii<19;++ii)
	  {
	    h_weigth_max[cc][ii]->Write("",TObject::kOverwrite);
	    h_weigth_min[cc][ii]->Write("",TObject::kOverwrite);
	    h_weigth_sum[cc][ii]->Write("",TObject::kOverwrite);
	    h_weigth_one[cc][ii]->Write("",TObject::kOverwrite);
	    h_chi2_det[cc][ii]->Write("",TObject::kOverwrite);
	    h_chi2_det_one[cc][ii]->Write("",TObject::kOverwrite);
	  }
    }
  
  //h_chi2_smooth->Write("",TObject::kOverwrite);
  //for(int i=0;i<5;i++)
  //   h_pull[i]->Write("",TObject::kOverwrite);

  particle->cd();
  
  h_mass_xTof->Write("",TObject::kOverwrite);
  h_mom_xTof->Write("",TObject::kOverwrite);
  h_beta_xTof->Write("",TObject::kOverwrite);
  h_chi2_xTof->Write("",TObject::kOverwrite);
  h_pv_xTof->Write("",TObject::kOverwrite);

  out_file->cd();
  // h_patern_decay->Write("",TObject::kOverwrite);
  h_stats->Write("",TObject::kOverwrite);
  h_presort->Write("",TObject::kOverwrite);
  // h_stats_Mom->Write("",TObject::kOverwrite);
  h_detection_eff_alpha->Write("",TObject::kOverwrite);
  h_detection_eff_pion->Write("",TObject::kOverwrite);
  h_detection_eff_proton->Write("",TObject::kOverwrite);
  
  if(EnableState[HOUGH] == true)
    {
      hough->cd();
      h_chi2_prefitY->Write("",TObject::kOverwrite);
      h_rchi2_prefitY->Write("",TObject::kOverwrite);
      h_hough_eff->Write("",TObject::kOverwrite);
      h_hough_diff->Write("",TObject::kOverwrite);
      h_hough_TOF_dist->Write("",TObject::kOverwrite);
      h_hough_DC_dist->Write("",TObject::kOverwrite);
      h_hough_TOF_dist_2->Write("",TObject::kOverwrite);
      h_hough_DC_dist_2->Write("",TObject::kOverwrite);
      h_hough_tot_dist->Write("",TObject::kOverwrite);
      h_hough_tot_dist_2->Write("",TObject::kOverwrite);
      h_hough_tot_pv->Write("",TObject::kOverwrite);
      out_file->cd();
    }
  //h_hough_each->Write("",TObject::kOverwrite);

  if(EnableState[FINDING] == true)
    {
      find->cd();
      h_finding_eff->SetOption("text");
      h_finding_eff->Write("",TObject::kOverwrite);
      h_finding_multi->Write("",TObject::kOverwrite);
      for(int i=0;i<7;++i)
	{
	  find_part[i]->cd();
	  for(int j=0;j<3;++j)
	    {
	      h_finding_chi2_tot_pid[i][j]->Write("",TObject::kOverwrite);
	      h_finding_chi2_x_pid[i][j]->Write("",TObject::kOverwrite);
	      h_finding_chi2_y_pid[i][j]->Write("",TObject::kOverwrite);
	    }
	  h_finding_chi2_2D[i]->Write("",TObject::kOverwrite);
	  h_finding_chi2_2D_ghost[i]->Write("",TObject::kOverwrite);
	}
      out_file->cd();
    }

  // AVhough->cd();
  // for(unsigned int i=0;i<10;++i)
  //   {
  //     AVHough_sub[i]->cd();
  //     for(unsigned int j=0;j<5;++j)
  // 	fHough_Z750[j+5*i]->Write("",TObject::kOverwrite);
      
  //for(unsigned int j=0;j<5;++j)
  //	fHough_Z400[j+5*i]->Write("",TObject::kOverwrite);
  //}




  if(EnableState[DCPROJ]==true)
    {
      Chamb->cd();
  
      h_dc1_project->Write("",TObject::kOverwrite);
      h_dc1x_project->Write("",TObject::kOverwrite);
      h_dc1u_project->Write("",TObject::kOverwrite);
      h_dc1v_project->Write("",TObject::kOverwrite);
      
      h_dc2_project->Write("",TObject::kOverwrite);
      h_dc2x_project->Write("",TObject::kOverwrite);
      h_dc2y_project->Write("",TObject::kOverwrite);
      h_dc2u_project->Write("",TObject::kOverwrite);
      out_file->cd();
	
    }




  Field->cd();
  h_Bz->Write("",TObject::kOverwrite);
  h_dB->Write("",TObject::kOverwrite);
  h_Bchi2->Write("",TObject::kOverwrite);

  /*hough->cd();
      h_hough_2d->Write("",TObject::kOverwrite);
      h_hough->Write("",TObject::kOverwrite);
      h_hough_pi->Write("",TObject::kOverwrite);
      h_hough_pi_dist_tof->Write("",TObject::kOverwrite);
      h_hough_pi_dist_dc->Write("",TObject::kOverwrite);
      
      h_hough_pi_tof->SetOption("text");
      h_hough_pi_tof->SetMarkerSize(1.6);
      h_hough_pi_tof->Write("",TObject::kOverwrite);
      
      h_hough_pi_dc->SetOption("text");
      h_hough_pi_dc->SetMarkerSize(1.6);
      h_hough_pi_dc->Write("",TObject::kOverwrite);
      
      h_hough_pi_dc_true->SetOption("text");
      h_hough_pi_dc_true->SetMarkerSize(1.6);
      h_hough_pi_dc_true->Write("",TObject::kOverwrite);
      
      h_hough_pi_dc_false->SetOption("text");
      h_hough_pi_dc_false->SetMarkerSize(1.6);
      h_hough_pi_dc_false->Write("",TObject::kOverwrite);
      
      h_hough_pi_dc_true_tof_true->SetOption("text");
      h_hough_pi_dc_true_tof_true->SetMarkerSize(1.6);
      h_hough_pi_dc_true_tof_true->Write("",TObject::kOverwrite);
      
      h_hough_pi_tof_true->SetOption("text");
      h_hough_pi_tof_true->SetMarkerSize(1.6);
      h_hough_pi_tof_true->Write("",TObject::kOverwrite);
       
      h_hough_pi_tof_false->SetOption("text");
      h_hough_pi_tof_false->SetMarkerSize(1.6);
      h_hough_pi_tof_false->Write("",TObject::kOverwrite);
      
      h_hough_pi_tof_true_dc_true->SetOption("text");
      h_hough_pi_tof_true_dc_true->SetMarkerSize(1.6);
      h_hough_pi_tof_true_dc_true->Write("",TObject::kOverwrite);
      
      
      h_hough_2d_optimal->Write("",TObject::kOverwrite);
      h_hough_optimal->Write("",TObject::kOverwrite);
  */

  out_file->cd();

//   h_mass_invariant_cut_mass->Write("",TObject::kOverwrite);
//   h_mass_invariant_best_cut_mass->Write("",TObject::kOverwrite);

//   h_mass_invariant_cut_vertex->Write("",TObject::kOverwrite);
//   h_mass_invariant_best_cut_vertex->Write("",TObject::kOverwrite);

//  h_dy->Write("",TObject::kOverwrite);

//  vertex->cd();
// h_vtx_x->Write("",TObject::kOverwrite);
// h_vtx_y->Write("",TObject::kOverwrite);
// h_vtx_z->Write("",TObject::kOverwrite);
// h_vtx_x_cut->Write("",TObject::kOverwrite);
// h_vtx_y_cut->Write("",TObject::kOverwrite);
// h_vtx_z_cut->Write("",TObject::kOverwrite);


  if(EnableState[VERTEX]== true)
    {
      for(int j=0;j<3;j++)
	{
	  pv[j]->cd();
	  
	  h_mass_invariant_cut_mass2[j]->SetDirectory(pv[j]);
	  h_mass_invariant_cut_mass2[j]->Write("",TObject::kOverwrite);
	  h_mass_invariant_best_cut_mass2[j]->SetDirectory(pv[j]);
	  h_mass_invariant_best_cut_mass2[j]->Write("",TObject::kOverwrite);
	  
	  h_mass_invariant_cut_vertex2[j]->SetDirectory(pv[j]);
	  h_mass_invariant_cut_vertex2[j]->Write("",TObject::kOverwrite);
	  h_mass_invariant_best_cut_vertex2[j]->SetDirectory(pv[j]);
	  h_mass_invariant_best_cut_vertex2[j]->Write("",TObject::kOverwrite);
	  
	  
	  for(int k=0;k<4;k++)
	    {
	      momD1D2[j][k]->cd();
	      h_momD1_momD2_cut_mass[j][k]->SetMarkerColor(k+1);
	      h_momD1_momD2_cut_mass[j][k]->SetFillColor(k+1);
	      h_momD1_momD2_cut_mass[j][k]->Write("",TObject::kOverwrite);
	      h_momD1_momD2_best_cut_mass[j][k]->SetMarkerColor(k+1);
	      h_momD1_momD2_best_cut_mass[j][k]->SetFillColor(k+1);
	      h_momD1_momD2_best_cut_mass[j][k]->Write("",TObject::kOverwrite);
	      h_momD1_momD2_cut_vertex[j][k]->SetMarkerColor(k+1);
	      h_momD1_momD2_cut_vertex[j][k]->SetFillColor(k+1);
	      h_momD1_momD2_cut_vertex[j][k]->Write("",TObject::kOverwrite);
	      h_momD1_momD2_best_cut_vertex[j][k]->SetMarkerColor(k+1);
	      h_momD1_momD2_best_cut_vertex[j][k]->SetFillColor(k+1);
	      //std::cout<<j<<" "<<k<<" "<<h_momD1_momD2_best_cut_vertex[j][k]<<std::endl;
	      h_momD1_momD2_best_cut_vertex[j][k]->Write("",TObject::kOverwrite);
	    }
	  //	 h_momD1_momD2_best_cut_vertex[0][0]->Write("",TObject::kOverwrite);
	  
	  pv[j]->cd();
	  
	  // h_vtx_x2[j]->SetDirectory(pv[j]);
	  // h_vtx_x2[j]->Write("",TObject::kOverwrite);
	  // h_vtx_y2[j]->SetDirectory(pv[j]);
	  // h_vtx_y2[j]->Write("",TObject::kOverwrite);
	  // h_vtx_z2[j]->SetDirectory(pv[j]);
	  // h_vtx_z2[j]->Write("",TObject::kOverwrite);
	  // h_vtx_x_cut2[j]->SetDirectory(pv[j]);
	  // h_vtx_x_cut2[j]->Write("",TObject::kOverwrite);
	  // h_vtx_y_cut2[j]->SetDirectory(pv[j]);
	  // h_vtx_y_cut2[j]->Write("",TObject::kOverwrite);
	  // h_vtx_z_cut2[j]->SetDirectory(pv[j]);
	  // h_vtx_z_cut2[j]->Write("",TObject::kOverwrite);
	  
	  for(int k=0;k<4;k++)
	    {
	      // h_dmom_kalman[j][k]->SetDirectory(pv[j]);
	      // h_dmom_kalman[j][k]->Write("",TObject::kOverwrite);
	      // h_dmom_normal_kalman[j][k]->SetDirectory(pv[j]);
	      // h_dmom_normal_kalman[j][k]->Write("",TObject::kOverwrite);
	      // h_dmom_mom_kalman[j][k]->SetDirectory(pv[j]);
	      // h_dmom_mom_kalman[j][k]->Write("",TObject::kOverwrite);
	      
	      h_Mass2[j][k]->SetDirectory(pv[j]);
	      h_Mass2[j][k]->Write("",TObject::kOverwrite);
	    }
	}
      out_file->cd();
    }


  
  out_file->cd();
  MC->cd();
  h_MC_PxPy->Write("",TObject::kOverwrite);
  h_MC_Pt->Write("",TObject::kOverwrite);
  h_MC_Pt_y->Write("",TObject::kOverwrite);
  h_MC_TargetXZ->Write("",TObject::kOverwrite);
  h_MC_TargetYZ->Write("",TObject::kOverwrite);
  h_MC_stats->Write("",TObject::kOverwrite);

  h_MC_dE_mom->Write("",TObject::kOverwrite);
  h_MC_time_mom->Write("",TObject::kOverwrite);
  h_MC_beta_mom->Write("",TObject::kOverwrite);
  h_MC_length->Write("",TObject::kOverwrite);
  h_MC_length_Z->Write("",TObject::kOverwrite);

  out_file->cd();
  Mat->cd();
  
  for(unsigned int n_v = 0;n_v<Material_XX0_y.size();++n_v)
    Material_XX0_y[n_v]->Write("",TObject::kOverwrite);

  for(unsigned int n_v = 0;n_v<Material_dE_y.size();++n_v)
    Material_dE_y[n_v]->Write("",TObject::kOverwrite);


 //Daisuke




 out_file->cd();
 Daisuke->cd();
 h_MassZ1BarBest->Write("",TObject::kOverwrite);
 h_MassZ1BarBestTofs->Write("",TObject::kOverwrite);
 //h_MassZ1BarBest_ini->Write("",TObject::kOverwrite);
 //h_MassZ1BarBest_iniTofs->Write("",TObject::kOverwrite);

//  h_MassZ1BarBest3->Write("",TObject::kOverwrite);
//  h_MassZ1BarBestTofs3->Write("",TObject::kOverwrite);
//  h_MassZ1BarBest_ini3->Write("",TObject::kOverwrite);
//  h_MassZ1BarBest_iniTofs3->Write("",TObject::kOverwrite);

 //h_MassZ2BarBest->Write("",TObject::kOverwrite);
 //h_MassZ2BarBestTofs->Write("",TObject::kOverwrite);
 //h_MassZ2BarBest_ini->Write("",TObject::kOverwrite);
 //h_MassZ2BarBest_iniTofs->Write("",TObject::kOverwrite);

//  h_MassZ2BarBest3->Write("",TObject::kOverwrite);
//  h_MassZ2BarBestTofs3->Write("",TObject::kOverwrite);
//  h_MassZ2BarBest_ini3->Write("",TObject::kOverwrite);
//  h_MassZ2BarBest_iniTofs3->Write("",TObject::kOverwrite);

 //h_MassZ1BarBest_pi->Write("",TObject::kOverwrite);
 //h_MassZ1BarBest_piTofs->Write("",TObject::kOverwrite);
 //h_MassZ1BarBest_ini_pi->Write("",TObject::kOverwrite);
 //h_MassZ1BarBest_ini_piTofs->Write("",TObject::kOverwrite);

//  h_MassZ1BarBest_pi3->Write("",TObject::kOverwrite);
//  h_MassZ1BarBest_piTofs3->Write("",TObject::kOverwrite);
//  h_MassZ1BarBest_ini_pi3->Write("",TObject::kOverwrite);
//  h_MassZ1BarBest_ini_piTofs3->Write("",TObject::kOverwrite);
 
 hd_theWorld->Write("",TObject::kOverwrite);
 hd_chi[0] ->Write("",TObject::kOverwrite);
 hd_chi[1] ->Write("",TObject::kOverwrite);
 hd_pv[0] ->Write("",TObject::kOverwrite);
 hd_pv[1] ->Write("",TObject::kOverwrite);
 hd_TrackFrom[0]->Write("",TObject::kOverwrite);
 hd_TrackFrom[1]->Write("",TObject::kOverwrite);
 hd_MassPv[0]->Write("",TObject::kOverwrite);
 hd_MassPv[1]->Write("",TObject::kOverwrite);
 for(int dd=0;dd<3;dd++)
   {
     h_MomCor[dd]->Write("",TObject::kOverwrite);
//      for(int mm=0;mm<3;mm++)
//        {
// 	 hd_secz[mm][dd]->Write("",TObject::kOverwrite);
// 	 hd_secz_cut[mm][dd]->Write("",TObject::kOverwrite);
// 	 hd_secz_cut2[mm][dd]->Write("",TObject::kOverwrite);
// 	 hd_secz_cut3[mm][dd]->Write("",TObject::kOverwrite);
// 	 hd_secz_all[mm][dd]->Write("",TObject::kOverwrite);
// 	 hd_secz_mompi[mm][dd]->Write("",TObject::kOverwrite);
//        }
   }
 for(int dd=0;dd<32;dd++)
   {
     //     hd_BetaMom_Tp[dd]->Write("",TObject::kOverwrite);
     hd_TimeCalTp[dd]->Write("",TObject::kOverwrite);
   }

 hd_BetaMom_Tp_all->Write("",TObject::kOverwrite);
 hd_BetaMom_Tp_allTofs->Write("",TObject::kOverwrite);
 hd_BetaMom_Tp_allTofs_Z2->Write("",TObject::kOverwrite);

 hd_BetaMom_Tp_ini_all->Write("",TObject::kOverwrite);
 hd_BetaMom_Tp_ini_allTofs->Write("",TObject::kOverwrite);
 hd_BetaMom_Tp_kal_allTofs->Write("",TObject::kOverwrite);
 hd_BetaMom_Tp_ini_allTofs_Z2->Write("",TObject::kOverwrite);
 //hd_Mom_ini_vs_rchi_ini->Write("",TObject::kOverwrite);
 //hd_Momini_rchiy_Tofs->Write("",TObject::kOverwrite);
//  for(int dd=0;dd<18;dd++)
//    hd_BetaMom_Bg[dd]->Write("",TObject::kOverwrite);

 hd_BetaMom_Bg_all->Write("",TObject::kOverwrite);
 hd_BetaMom_Bg_allTofs->Write("",TObject::kOverwrite);
 hd_BetaMom_Bg_ini_all->Write("",TObject::kOverwrite);
 hd_BetaMom_Bg_ini_allTofs->Write("",TObject::kOverwrite);
//  for(int dd=0;dd<4;dd++)
//    for(int jj=0;jj<3;jj++)
//      for(int mm=0;mm<3;mm++)
//        {
// 	 hd_invmass[mm][jj][dd]->Write("",TObject::kOverwrite);
// 	 hd_invmass_secz[mm][jj][dd]->Write("",TObject::kOverwrite);
//    }
 // out_file->Write("",TObject::kOverwrite);
 //hd_PoQ_E->Write("",TObject::kOverwrite);
  //hd_PoQ_E_ini->Write("",TObject::kOverwrite);
  //hd_PoQ_ETofs->Write("",TObject::kOverwrite);
  //hd_PoQ_E_iniTofs->Write("",TObject::kOverwrite);


//   hd_PoQ_E3->Write("",TObject::kOverwrite);
//   hd_PoQ_E_ini3->Write("",TObject::kOverwrite);
//   hd_PoQ_ETofs3->Write("",TObject::kOverwrite);
//   hd_PoQ_E_iniTofs3->Write("",TObject::kOverwrite);

#endif

 return 0;
}


int Ana_Hist::WriteTemp(char *tempfile)
{
  HaveBeenWritten=true;
//   char tempfile[300];
  
//   sprintf(tempfile,"Temp_OutFile_%i_%i_%s_%s",Now,HMS,HostName.c_str(),nameoutfile2.Data());
  
  TFile *ff = new TFile(tempfile,"RECREATE");
  cout<<"File = "<<tempfile<<endl;

  ff->cd();


  TDirectory *bmtp = ff->mkdir("BetaMom_Tof+"); 
  TDirectory *bmbg = ff->mkdir("BetaMom_BigTof+"); 
  TDirectory *Mother[3];//
  Mother[0] =(TDirectory *) ff->mkdir("Lambda");
  Mother[1] =(TDirectory *) ff->mkdir("H3L");
  Mother[2] =(TDirectory *) ff->mkdir("H4L");
  TDirectory *invpv[3][3];
  for(int mm=0;mm<3;mm++)
    {
      invpv[mm][0] = (TDirectory *)Mother[mm]->mkdir("pv_tight");
      invpv[mm][1] = (TDirectory *)Mother[mm]->mkdir("pv_mid");
      invpv[mm][2] = (TDirectory *)Mother[mm]->mkdir("pv_loose");
    }

  //#ifdef ADD_HISTO

  h_MassZ1BarBest->Write("",TObject::kOverwrite);

  //h_MassZ1BarBest_ini->Write("",TObject::kOverwrite);

  h_MassZ1BarBestTofs->Write("",TObject::kOverwrite);

  //h_MassZ1BarBest_iniTofs->Write("",TObject::kOverwrite);

//   h_MassZ1BarBest3->Write("",TObject::kOverwrite);
//   h_MassZ1BarBestTofs3->Write("",TObject::kOverwrite);
//   h_MassZ1BarBest_ini3->Write("",TObject::kOverwrite);
//   h_MassZ1BarBest_iniTofs3->Write("",TObject::kOverwrite);

  //h_MassZ2BarBest->Write("",TObject::kOverwrite);
  //h_MassZ2BarBestTofs->Write("",TObject::kOverwrite);
  //h_MassZ2BarBest_ini->Write("",TObject::kOverwrite);
  //h_MassZ2BarBest_iniTofs->Write("",TObject::kOverwrite);

//   h_MassZ2BarBest3->Write("",TObject::kOverwrite);
//   h_MassZ2BarBestTofs3->Write("",TObject::kOverwrite);
//   h_MassZ2BarBest_ini3->Write("",TObject::kOverwrite);
//   h_MassZ2BarBest_iniTofs3->Write("",TObject::kOverwrite);

  //h_MassZ1BarBest_pi->Write("",TObject::kOverwrite);
  //h_MassZ1BarBest_piTofs->Write("",TObject::kOverwrite);
  //h_MassZ1BarBest_ini_pi->Write("",TObject::kOverwrite);
  //h_MassZ1BarBest_ini_piTofs->Write("",TObject::kOverwrite);

//   h_MassZ1BarBest_pi3->Write("",TObject::kOverwrite);
//   h_MassZ1BarBest_piTofs3->Write("",TObject::kOverwrite);
//   h_MassZ1BarBest_ini_pi3->Write("",TObject::kOverwrite);
//   h_MassZ1BarBest_ini_piTofs3->Write("",TObject::kOverwrite);

  hd_theWorld->Write("",TObject::kOverwrite);
  hd_chi[0] ->Write("",TObject::kOverwrite);
  hd_chi[1] ->Write("",TObject::kOverwrite);
  hd_pv[0] ->Write("",TObject::kOverwrite);
  hd_pv[1] ->Write("",TObject::kOverwrite);
  hd_TrackFrom[0]->Write("",TObject::kOverwrite);
  hd_TrackFrom[1]->Write("",TObject::kOverwrite);
  hd_MassPv[0]->Write("",TObject::kOverwrite);
  hd_MassPv[1]->Write("",TObject::kOverwrite);
  for(auto & elem : h_MomCor)
    {
      elem->Write("",TObject::kOverwrite);
//       for(int mm=0;mm<3;mm++)
// 	{
// 	  Mother[mm]->cd();
// 	  hd_secz[mm][dd]->Write("",TObject::kOverwrite);
// 	  hd_secz_cut[mm][dd]->Write("",TObject::kOverwrite);
// 	  hd_secz_cut2[mm][dd]->Write("",TObject::kOverwrite);
// 	  hd_secz_cut3[mm][dd]->Write("",TObject::kOverwrite);
// 	  hd_secz_all[mm][dd]->Write("",TObject::kOverwrite);
// 	  hd_secz_mompi[mm][dd]->Write("",TObject::kOverwrite);
// 	}
 
    }
  for(auto & elem : hd_TimeCalTp)
    {
      bmtp->cd();
      //hd_BetaMom_Tp[dd]->SetDirectory(bmtp);
      //hd_BetaMom_Tp[dd]->Write("",TObject::kOverwrite);
      //      hd_TimeCalTp[dd]->SetDirectory(bmtp);
      elem->Write("",TObject::kOverwrite);


    }
  ff->cd();
  hd_BetaMom_Tp_all->Write("",TObject::kOverwrite);
  hd_BetaMom_Tp_allTofs->Write("",TObject::kOverwrite);
  hd_BetaMom_Tp_allTofs_Z2->Write("",TObject::kOverwrite);
  
  hd_BetaMom_Tp_ini_all->Write("",TObject::kOverwrite);
  hd_BetaMom_Tp_ini_allTofs->Write("",TObject::kOverwrite);
  hd_BetaMom_Tp_kal_allTofs->Write("",TObject::kOverwrite);
  hd_BetaMom_Tp_ini_allTofs_Z2->Write("",TObject::kOverwrite);
  //hd_Mom_ini_vs_rchi_ini->Write("",TObject::kOverwrite);
  //hd_Momini_rchiy_Tofs->Write("",TObject::kOverwrite);

  for(int dd=0;dd<18;dd++)
    {
      bmbg->cd();
      //hd_BetaMom_Bg[dd]->SetDirectory(bmbg);
      //hd_BetaMom_Bg[dd]->Write("",TObject::kOverwrite);

    }
  ff->cd();
  hd_BetaMom_Bg_all->Write("",TObject::kOverwrite);

  hd_BetaMom_Bg_allTofs->Write("",TObject::kOverwrite);
  hd_BetaMom_Bg_ini_all->Write("",TObject::kOverwrite);
  hd_BetaMom_Bg_ini_allTofs->Write("",TObject::kOverwrite);
//   for(int dd=0;dd<4;dd++)
//     for(int jj=0;jj<3;jj++)
//       for(int mm=0;mm<3;mm++)
// 	{
	  
// 	  invpv[mm][jj]->cd();
// 	  //	hd_invmass[jj][dd]->SetDirectory(invpv[jj]);
// 	  hd_invmass[mm][jj][dd]->Write("",TObject::kOverwrite);
	  
// 	  //hd_invmass_secz[jj][dd]->SetDirectory(invpv[jj]);
// 	  hd_invmass_secz[mm][jj][dd]->Write("",TObject::kOverwrite);
	  
//       }


  ff->cd();
  //hd_PoQ_E->Write("",TObject::kOverwrite);
  //hd_PoQ_E_ini->Write("",TObject::kOverwrite);
  //hd_PoQ_ETofs->Write("",TObject::kOverwrite);
  //hd_PoQ_E_iniTofs->Write("",TObject::kOverwrite);


//   hd_PoQ_E3->Write("",TObject::kOverwrite);
//   hd_PoQ_E_ini3->Write("",TObject::kOverwrite);
//   hd_PoQ_ETofs3->Write("",TObject::kOverwrite);
//   hd_PoQ_E_iniTofs3->Write("",TObject::kOverwrite);
//#endif

  //hd_IniMom_DiffMom_allTofs->Write("",TObject::kOverwrite);
  //hd_NewMom_DiffMom_allTofs->Write("",TObject::kOverwrite);
  //hd_IniMom_NewMom_all->Write("",TObject::kOverwrite);
  //hd_IniMom_NewMom_all_tofs->Write("",TObject::kOverwrite);
  //hd_IniMom_NewMom_all_tofs_ini->Write("",TObject::kOverwrite);
  hd_IniMom_NewMom_all_tofsZ1->Write("",TObject::kOverwrite);
  hd_IniMom_NewMom_all_tofsZ1_ini->Write("",TObject::kOverwrite);
  //hd_Bar_DiffMom_allTofs->Write("",TObject::kOverwrite);
  hd_IniNewMomdiff_vs_IniNewMomplus_allTofs->Write("",TObject::kOverwrite);
  hd_charge_ini_vs_IniNewDiv->Write("",TObject::kOverwrite);
  hd_charge_vs_IniNewDiv->Write("",TObject::kOverwrite);
  hd_IniNewMomdiff_vs_IniNewMomplus_allTofs_Z1->Write("",TObject::kOverwrite);
  hd_IniNewMomdiff_vs_IniNewMomplus_allTofs_ini_Z1->Write("",TObject::kOverwrite);
  hd_best_pos_hitvstrack_tofpZ1->Write("",TObject::kOverwrite);
  hd_best_pos_hitdifftrack_tofpZ1->Write("",TObject::kOverwrite);


  //hd_IniMom_DiffMom_all->Write("",TObject::kOverwrite);
  //hd_NewMom_DiffMom_all->Write("",TObject::kOverwrite);
  //hd_IniMom_NewMom_all->Write("",TObject::kOverwrite);
  //hd_Bar_DiffMom_all->Write("",TObject::kOverwrite);

  ff->Close();
  delete ff;
  //h_2TR1->Write();
  //h_2TR2->Write();
  //h_2AladinTOF->Write();

//   TDirectory* det = ff->mkdir("Detector");
//   //TDirectory* AVhough = ff->mkdir("AVHough");
//   //TDirectory* AVHough_sub[10];


//   TDirectory* Chamb = ff->mkdir("ChamberProject");


  
//   // for(int i=0;i<10;++i)
//   //   {
//   //     TString name("Event_");
//   //     name+=i;
//   //     AVHough_sub[i] = AVhough->mkdir(name);
//   //   }

//   TDirectory* det_pull = det->mkdir("pull");
//   TDirectory* mom1 = ff->mkdir("Momentum_Check");
//   //TDirectory* mom2 = ff->mkdir("Momentum_after_mass");
//   //TDirectory* mom3 = ff->mkdir("Momentum_after_vertex");
//   TDirectory* mass = ff->mkdir("Mass");
//   //TDirectory* vertex = ff->mkdir("Vertex");
//   TDirectory* hough = ff->mkdir("Hough");
//   TDirectory* find = ff->mkdir("Finding");
//   TDirectory* find_part[7];
//   TString part[7]={"pi-","pi+","K+","alpha","proton","He3","deuteron"};
//   for(int i=0;i<7;++i)
//     find_part[i]=(TDirectory*) find->mkdir(part[i]);

//   TDirectory* kalman = ff->mkdir("Kalman");
//   TDirectory* particle= ff->mkdir("Particle");
//   TDirectory* pv[3];
//   pv[0] = (TDirectory*)ff->mkdir("pv_0.01");
//   pv[1] = (TDirectory*)ff->mkdir("pv_0.05");
//   pv[2] = (TDirectory*)ff->mkdir("pv_0.1");
  
//   TDirectory* momD1D2[3][4];

//   for(int i=0;i<3;i++)
//     {
//       momD1D2[i][0]= (TDirectory*) pv[i]->mkdir("L");
//       momD1D2[i][1]= (TDirectory*) pv[i]->mkdir("H3L");
//       momD1D2[i][2]= (TDirectory*) pv[i]->mkdir("H4L");
//       momD1D2[i][3]= (TDirectory*) pv[i]->mkdir("He5L");
//     }
//   TDirectory* Field = ff->mkdir("Field");

//   ff->cd();
//   h_task_exit->Write();
  
//   h_alpha_rejection->Write();
//   det->cd();
//   for (Int_t n=0;n<18;n++)
//     {
//       h_2utrPos_XZ[n]->Write();
//       h_2utrPos_XZ[n]->SetDirectory(det);
//       h_2utrPos_YZ[n]->Write();
//       h_2utrPos_YZ[n]->SetDirectory(det);
//       h_2utrPos_XY[n]->Write();
//       h_2utrPos_XY[n]->SetDirectory(det);
//       h_utrEnergy[n]->Write();
//       h_utrEnergy[n]->SetDirectory(det);
//       h_utrTime[n]->Write();
//       h_utrTime[n]->SetDirectory(det);
//     }
//   h_patern_decay->SetOption("text");
//   h_patern_decay->SetMarkerSize(1.6);
//   h_patern_decay->Write();
//   h_2utrRotated_XZ->Write();
//   h_2utrFrameMag_XZ->Write();
//   h_2utrNotRotated_XZ->Write();

//   det_pull->cd();
//   for (Int_t n=0;n<18;n++)
//     {
//       h_pull_X_det[n]->Write();
//       h_pull_X_det[n]->SetDirectory(det_pull);
//       h_pull_Y_det[n]->Write();
//       h_pull_Y_det[n]->SetDirectory(det_pull);
//     }
  
//   for (int i=0;i<4;i++)
//     {
//       ff->cd();
//       h_chi2_particle[i]->Write();
//       h_pv_particle[i]->Write();

//       // mom1->cd();
//       // h_dmom1[i]->SetDirectory(mom1);
//       // h_dmom1[i]->Write();
//       // h_dmom2[i]->SetDirectory(mom1);
//       // h_dmom2[i]->Write();
//       // h_dmom3[i]->SetDirectory(mom1);
//       // h_dmom3[i]->Write();
//       // h_dmom4[i]->SetDirectory(mom1);
//       // h_dmom4[i]->Write();
//       // mom2->cd();
//       // h_dmom_dmass[i]->SetDirectory(mom2);
//       // h_dmom_dmass[i]->Write();
//       // h_dmom_dmass_best[i]->SetDirectory(mom2);
//       // h_dmom_dmass_best[i]->Write();
//       // mom3->cd();
//       // h_dmom_vertex[i]->SetDirectory(mom3);
//       // h_dmom_vertex[i]->Write();
//       // h_dmom_vertex_best[i]->SetDirectory(mom3);
//       // h_dmom_vertex_best[i]->Write();

//       mass->cd();
//       h_Mass[i]->SetDirectory(mass);
//       h_Mass[i]->Write();
	
//     }
//   h_Mass_All->Write();
//   h_Mass_All2->Write();
//   h_Mass_All3->Write();

//   h_Mass_charge_All->Write();
//   h_Mass_charge_All2->Write();
//   h_Mass_charge_All3->Write();

//   mom1->cd();
//   h_mom_check->Write();
//   h_mom_check->SetDirectory(mom1);
//   h_mom_check2->Write();
//   h_mom_check2->SetDirectory(mom1);
//   h_tofdE->Write();
//   h_tofdE->SetDirectory(mom1);

//   ff->cd();
//   // for (int i=0;i<4;i++)
//   //   h_slope_yz[i]->Write();


//   h_beta->Write();
//   h_beta2->Write();
//   h_beta3->Write();

//   h_beta_mom->Write();
//   h_beta_mom2->Write();
//   h_beta_mom3->Write();

//   kalman->cd();
//   h_pv->Write();
//   h_chi2->Write();
//   h_Path->Write();
//   h_Path_Back->Write();
//   h_MeanPath->Write();
//   h_path_tof->Write();
//   h_pv_mom->Write();
//   h_pv_beta->Write();
//   h_pv_mass->Write();
//   h_mom_tof_cut->Write();
//   h_path_tof_cut->Write();
//   h_dpath->Write();
//   h_path_mom_cut->Write();
//   //h_chi2_smooth->Write();
//   //for(int i=0;i<5;i++)
//   //   h_pull[i]->Write();

//   particle->cd();
  
//   h_mass_xTof->Write();
//   h_mom_xTof->Write();
//   h_beta_xTof->Write();
//   h_chi2_xTof->Write();
//   h_pv_xTof->Write();

//   ff->cd();
//   // h_patern_decay->Write();
//   h_stats->Write();
//   h_presort->Write();
//   // h_stats_Mom->Write();
//   h_detection_eff_alpha->Write();
//   h_detection_eff_pion->Write();
//   h_detection_eff_proton->Write();
//   hough->cd();
//   h_chi2_prefitY->Write();
//   h_rchi2_prefitY->Write();
//   h_hough_eff->Write();
//   h_hough_diff->Write();
//   h_hough_TOF_dist->Write();
//   h_hough_DC_dist->Write();
//   h_hough_TOF_dist_2->Write();
//   h_hough_DC_dist_2->Write();
//   h_hough_tot_dist->Write();
//   h_hough_tot_dist_2->Write();
//   h_hough_tot_pv->Write();
//   //h_hough_each->Write();

//   ff->cd();
//   find->cd();
//   h_finding_eff->SetOption("text");
//   h_finding_eff->Write();
//   h_finding_multi->Write();
//   for(int i=0;i<7;++i)
//     {
//       find_part[i]->cd();
//       for(int j=0;j<3;++j)
// 	{
// 	  h_finding_chi2_tot_pid[i][j]->Write();
// 	  h_finding_chi2_x_pid[i][j]->Write();
// 	  h_finding_chi2_y_pid[i][j]->Write();
// 	}
//       h_finding_chi2_2D[i]->Write();
//       h_finding_chi2_2D_ghost[i]->Write();
//     }


//   // AVhough->cd();
//   // for(unsigned int i=0;i<10;++i)
//   //   {
//   //     AVHough_sub[i]->cd();
//   //     for(unsigned int j=0;j<5;++j)
//   // 	fHough_Z750[j+5*i]->Write();
      
//   //for(unsigned int j=0;j<5;++j)
//   //	fHough_Z400[j+5*i]->Write();
//   //}





//   Chamb->cd();

//   h_dc1_project->Write();
//   h_dc1x_project->Write();
//   h_dc1u_project->Write();
//   h_dc1v_project->Write();

//   h_dc2_project->Write();
//   h_dc2x_project->Write();
//   h_dc2y_project->Write();
//   h_dc2u_project->Write();





//   Field->cd();
//   h_Bz->Write();
//   h_dB->Write();
//   h_Bchi2->Write();

//   // h_hough_2d->Write();
//   // h_hough->Write();
//   // h_hough_pi->Write();
//   // h_hough_pi_dist_tof->Write();
//   // h_hough_pi_dist_dc->Write();

//   // h_hough_pi_tof->SetOption("text");
//   // h_hough_pi_tof->SetMarkerSize(1.6);
//   // h_hough_pi_tof->Write();

//   // h_hough_pi_dc->SetOption("text");
//   // h_hough_pi_dc->SetMarkerSize(1.6);
//   // h_hough_pi_dc->Write();

//   // h_hough_pi_dc_true->SetOption("text");
//   // h_hough_pi_dc_true->SetMarkerSize(1.6);
//   // h_hough_pi_dc_true->Write();

//   // h_hough_pi_dc_false->SetOption("text");
//   // h_hough_pi_dc_false->SetMarkerSize(1.6);
//   // h_hough_pi_dc_false->Write();

//   // h_hough_pi_dc_true_tof_true->SetOption("text");
//   // h_hough_pi_dc_true_tof_true->SetMarkerSize(1.6);
//   // h_hough_pi_dc_true_tof_true->Write();

//   // h_hough_pi_tof_true->SetOption("text");
//   // h_hough_pi_tof_true->SetMarkerSize(1.6);
//   // h_hough_pi_tof_true->Write();

//   // h_hough_pi_tof_false->SetOption("text");
//   // h_hough_pi_tof_false->SetMarkerSize(1.6);
//   // h_hough_pi_tof_false->Write();

//   // h_hough_pi_tof_true_dc_true->SetOption("text");
//   // h_hough_pi_tof_true_dc_true->SetMarkerSize(1.6);
//   // h_hough_pi_tof_true_dc_true->Write();


//   // h_hough_2d_optimal->Write();
//   // h_hough_optimal->Write();
//   ff->cd();

// //   h_mass_invariant_cut_mass->Write();
// //   h_mass_invariant_best_cut_mass->Write();

// //   h_mass_invariant_cut_vertex->Write();
// //   h_mass_invariant_best_cut_vertex->Write();

// //  h_dy->Write();

// //  vertex->cd();
// // h_vtx_x->Write();
// // h_vtx_y->Write();
// // h_vtx_z->Write();
// // h_vtx_x_cut->Write();
// // h_vtx_y_cut->Write();
// // h_vtx_z_cut->Write();



//  for(int j=0;j<3;j++)
//    {
//      pv[j]->cd();
     
//      h_mass_invariant_cut_mass2[j]->SetDirectory(pv[j]);
//      h_mass_invariant_cut_mass2[j]->Write();
//      h_mass_invariant_best_cut_mass2[j]->SetDirectory(pv[j]);
//      h_mass_invariant_best_cut_mass2[j]->Write();
     
//      h_mass_invariant_cut_vertex2[j]->SetDirectory(pv[j]);
//      h_mass_invariant_cut_vertex2[j]->Write();
//      h_mass_invariant_best_cut_vertex2[j]->SetDirectory(pv[j]);
//      h_mass_invariant_best_cut_vertex2[j]->Write();
     

//      for(int k=0;k<4;k++)
//        {
// 	 momD1D2[j][k]->cd();
// 	 h_momD1_momD2_cut_mass[j][k]->SetMarkerColor(k+1);
// 	 h_momD1_momD2_cut_mass[j][k]->SetFillColor(k+1);
// 	 h_momD1_momD2_cut_mass[j][k]->Write();
// 	 h_momD1_momD2_best_cut_mass[j][k]->SetMarkerColor(k+1);
// 	 h_momD1_momD2_best_cut_mass[j][k]->SetFillColor(k+1);
// 	 h_momD1_momD2_best_cut_mass[j][k]->Write();
// 	 h_momD1_momD2_cut_vertex[j][k]->SetMarkerColor(k+1);
// 	 h_momD1_momD2_cut_vertex[j][k]->SetFillColor(k+1);
// 	 h_momD1_momD2_cut_vertex[j][k]->Write();
// 	 h_momD1_momD2_best_cut_vertex[j][k]->SetMarkerColor(k+1);
// 	 h_momD1_momD2_best_cut_vertex[j][k]->SetFillColor(k+1);
// 	 //std::cout<<j<<" "<<k<<" "<<h_momD1_momD2_best_cut_vertex[j][k]<<std::endl;
// 	 h_momD1_momD2_best_cut_vertex[j][k]->Write();
//        }

// //	 h_momD1_momD2_best_cut_vertex[0][0]->Write();

//      pv[j]->cd();

//      // h_vtx_x2[j]->SetDirectory(pv[j]);
//      // h_vtx_x2[j]->Write();
//      // h_vtx_y2[j]->SetDirectory(pv[j]);
//      // h_vtx_y2[j]->Write();
//      // h_vtx_z2[j]->SetDirectory(pv[j]);
//      // h_vtx_z2[j]->Write();
//      // h_vtx_x_cut2[j]->SetDirectory(pv[j]);
//      // h_vtx_x_cut2[j]->Write();
//      // h_vtx_y_cut2[j]->SetDirectory(pv[j]);
//      // h_vtx_y_cut2[j]->Write();
//      // h_vtx_z_cut2[j]->SetDirectory(pv[j]);
//      // h_vtx_z_cut2[j]->Write();
     
//      for(int k=0;k<4;k++)
//        {
// 	 // h_dmom_kalman[j][k]->SetDirectory(pv[j]);
// 	 // h_dmom_kalman[j][k]->Write();
// 	 // h_dmom_normal_kalman[j][k]->SetDirectory(pv[j]);
// 	 // h_dmom_normal_kalman[j][k]->Write();
// 	 // h_dmom_mom_kalman[j][k]->SetDirectory(pv[j]);
// 	 // h_dmom_mom_kalman[j][k]->Write();
	 
// 	 h_Mass2[j][k]->SetDirectory(pv[j]);
// 	 h_Mass2[j][k]->Write();
//        }

//    }

//  delete ff;
 return 0;
}






//

// int Ana_Hist::Write(TFile* out_file)
// {
//   HaveBeenWritten=true;
//   out_file->cd();
//   //h_2TR1->Write();
//   //h_2TR2->Write();
//   //h_2AladinTOF->Write();

//   TDirectory* det = out_file->mkdir("Detector");
//   //TDirectory* AVhough = out_file->mkdir("AVHough");
//   //TDirectory* AVHough_sub[10];


//   TDirectory* Chamb = out_file->mkdir("ChamberProject");


  
//   // for(int i=0;i<10;++i)
//   //   {
//   //     TString name("Event_");
//   //     name+=i;
//   //     AVHough_sub[i] = AVhough->mkdir(name);
//   //   }

//   TDirectory* det_pull = det->mkdir("pull");
//   TDirectory* mom1 = out_file->mkdir("Momentum_Check");
//   //TDirectory* mom2 = out_file->mkdir("Momentum_after_mass");
//   //TDirectory* mom3 = out_file->mkdir("Momentum_after_vertex");
//   TDirectory* mass = out_file->mkdir("Mass");
//   //TDirectory* vertex = out_file->mkdir("Vertex");
//   TDirectory* hough = out_file->mkdir("Hough");
//   TDirectory* find = out_file->mkdir("Finding");
//   TDirectory* find_part[7];
//   TString part[7]={"pi-","pi+","K+","alpha","proton","He3","deuteron"};
//   for(int i=0;i<7;++i)
//     find_part[i]=(TDirectory*) find->mkdir(part[i]);

//   TDirectory* kalman = out_file->mkdir("Kalman");
//   TDirectory* particle= out_file->mkdir("Particle");
//   TDirectory* pv[3];
//   pv[0] = (TDirectory*)out_file->mkdir("pv_0.01");
//   pv[1] = (TDirectory*)out_file->mkdir("pv_0.05");
//   pv[2] = (TDirectory*)out_file->mkdir("pv_0.1");
  
//   TDirectory* momD1D2[3][4];

//   for(int i=0;i<3;i++)
//     {
//       momD1D2[i][0]= (TDirectory*) pv[i]->mkdir("L");
//       momD1D2[i][1]= (TDirectory*) pv[i]->mkdir("H3L");
//       momD1D2[i][2]= (TDirectory*) pv[i]->mkdir("H4L");
//       momD1D2[i][3]= (TDirectory*) pv[i]->mkdir("He5L");
//     }
//   TDirectory* Field = out_file->mkdir("Field");

//   out_file->cd();
//   h_task_exit->Write();
  
//   h_alpha_rejection->Write();
//   det->cd();
//   for (Int_t n=0;n<18;n++)
//     {
//       h_2utrPos_XZ[n]->Write();
//       h_2utrPos_XZ[n]->SetDirectory(det);
//       h_2utrPos_YZ[n]->Write();
//       h_2utrPos_YZ[n]->SetDirectory(det);
//       h_2utrPos_XY[n]->Write();
//       h_2utrPos_XY[n]->SetDirectory(det);
//       h_utrEnergy[n]->Write();
//       h_utrEnergy[n]->SetDirectory(det);
//       h_utrTime[n]->Write();
//       h_utrTime[n]->SetDirectory(det);
//     }
//   h_patern_decay->SetOption("text");
//   h_patern_decay->SetMarkerSize(1.6);
//   h_patern_decay->Write();
//   h_2utrRotated_XZ->Write();
//   h_2utrFrameMag_XZ->Write();
//   h_2utrNotRotated_XZ->Write();

//   det_pull->cd();
//   for (Int_t n=0;n<18;n++)
//     {
//       h_pull_X_det[n]->Write();
//       h_pull_X_det[n]->SetDirectory(det_pull);
//       h_pull_Y_det[n]->Write();
//       h_pull_Y_det[n]->SetDirectory(det_pull);
//     }
  
//   for (int i=0;i<4;i++)
//     {
//       out_file->cd();
//       h_chi2_particle[i]->Write();
//       h_pv_particle[i]->Write();

//       // mom1->cd();
//       // h_dmom1[i]->SetDirectory(mom1);
//       // h_dmom1[i]->Write();
//       // h_dmom2[i]->SetDirectory(mom1);
//       // h_dmom2[i]->Write();
//       // h_dmom3[i]->SetDirectory(mom1);
//       // h_dmom3[i]->Write();
//       // h_dmom4[i]->SetDirectory(mom1);
//       // h_dmom4[i]->Write();
//       // mom2->cd();
//       // h_dmom_dmass[i]->SetDirectory(mom2);
//       // h_dmom_dmass[i]->Write();
//       // h_dmom_dmass_best[i]->SetDirectory(mom2);
//       // h_dmom_dmass_best[i]->Write();
//       // mom3->cd();
//       // h_dmom_vertex[i]->SetDirectory(mom3);
//       // h_dmom_vertex[i]->Write();
//       // h_dmom_vertex_best[i]->SetDirectory(mom3);
//       // h_dmom_vertex_best[i]->Write();

//       mass->cd();
//       h_Mass[i]->SetDirectory(mass);
//       h_Mass[i]->Write();
	
//     }
//   h_Mass_All->Write();
//   h_Mass_All2->Write();
//   h_Mass_All3->Write();

//   h_Mass_charge_All->Write();
//   h_Mass_charge_All2->Write();
//   h_Mass_charge_All3->Write();

//   mom1->cd();
//   h_mom_check->Write();
//   h_mom_check->SetDirectory(mom1);
//   h_mom_check2->Write();
//   h_mom_check2->SetDirectory(mom1);
//   h_tofdE->Write();
//   h_tofdE->SetDirectory(mom1);

//   out_file->cd();
//   // for (int i=0;i<4;i++)
//   //   h_slope_yz[i]->Write();


//   h_beta->Write();
//   h_beta2->Write();
//   h_beta3->Write();

//   h_beta_mom->Write();
//   h_beta_mom2->Write();
//   h_beta_mom3->Write();

//   kalman->cd();
//   h_pv->Write();
//   h_chi2->Write();
//   h_Path->Write();
//   h_Path_Back->Write();
//   h_MeanPath->Write();
//   h_path_tof->Write();
//   h_pv_mom->Write();
//   h_pv_beta->Write();
//   h_pv_mass->Write();
//   h_mom_tof_cut->Write();
//   h_path_tof_cut->Write();
//   h_dpath->Write();
//   h_path_mom_cut->Write();
//   //h_chi2_smooth->Write();
//   //for(int i=0;i<5;i++)
//   //   h_pull[i]->Write();

//   particle->cd();
  
//   h_mass_xTof->Write();
//   h_mom_xTof->Write();
//   h_beta_xTof->Write();
//   h_chi2_xTof->Write();
//   h_pv_xTof->Write();

//   out_file->cd();
//   // h_patern_decay->Write();
//   h_stats->Write();
//   h_presort->Write();
//   // h_stats_Mom->Write();
//   h_detection_eff_alpha->Write();
//   h_detection_eff_pion->Write();
//   h_detection_eff_proton->Write();
//   hough->cd();
//   h_chi2_prefitY->Write();
//   h_rchi2_prefitY->Write();
//   h_hough_eff->Write();
//   h_hough_diff->Write();
//   h_hough_TOF_dist->Write();
//   h_hough_DC_dist->Write();
//   h_hough_TOF_dist_2->Write();
//   h_hough_DC_dist_2->Write();
//   h_hough_tot_dist->Write();
//   h_hough_tot_dist_2->Write();
//   h_hough_tot_pv->Write();
//   //h_hough_each->Write();

//   out_file->cd();
//   find->cd();
//   h_finding_eff->SetOption("text");
//   h_finding_eff->Write();
//   h_finding_multi->Write();
//   for(int i=0;i<7;++i)
//     {
//       find_part[i]->cd();
//       for(int j=0;j<3;++j)
// 	{
// 	  h_finding_chi2_tot_pid[i][j]->Write();
// 	  h_finding_chi2_x_pid[i][j]->Write();
// 	  h_finding_chi2_y_pid[i][j]->Write();
// 	}
//       h_finding_chi2_2D[i]->Write();
//       h_finding_chi2_2D_ghost[i]->Write();
//     }


//   // AVhough->cd();
//   // for(unsigned int i=0;i<10;++i)
//   //   {
//   //     AVHough_sub[i]->cd();
//   //     for(unsigned int j=0;j<5;++j)
//   // 	fHough_Z750[j+5*i]->Write();
      
//   //for(unsigned int j=0;j<5;++j)
//   //	fHough_Z400[j+5*i]->Write();
//   //}





//   Chamb->cd();

//   h_dc1_project->Write();
//   h_dc1x_project->Write();
//   h_dc1u_project->Write();
//   h_dc1v_project->Write();

//   h_dc2_project->Write();
//   h_dc2x_project->Write();
//   h_dc2y_project->Write();
//   h_dc2u_project->Write();





//   Field->cd();
//   h_Bz->Write();
//   h_dB->Write();
//   h_Bchi2->Write();

//   // h_hough_2d->Write();
//   // h_hough->Write();
//   // h_hough_pi->Write();
//   // h_hough_pi_dist_tof->Write();
//   // h_hough_pi_dist_dc->Write();

//   // h_hough_pi_tof->SetOption("text");
//   // h_hough_pi_tof->SetMarkerSize(1.6);
//   // h_hough_pi_tof->Write();

//   // h_hough_pi_dc->SetOption("text");
//   // h_hough_pi_dc->SetMarkerSize(1.6);
//   // h_hough_pi_dc->Write();

//   // h_hough_pi_dc_true->SetOption("text");
//   // h_hough_pi_dc_true->SetMarkerSize(1.6);
//   // h_hough_pi_dc_true->Write();

//   // h_hough_pi_dc_false->SetOption("text");
//   // h_hough_pi_dc_false->SetMarkerSize(1.6);
//   // h_hough_pi_dc_false->Write();

//   // h_hough_pi_dc_true_tof_true->SetOption("text");
//   // h_hough_pi_dc_true_tof_true->SetMarkerSize(1.6);
//   // h_hough_pi_dc_true_tof_true->Write();

//   // h_hough_pi_tof_true->SetOption("text");
//   // h_hough_pi_tof_true->SetMarkerSize(1.6);
//   // h_hough_pi_tof_true->Write();

//   // h_hough_pi_tof_false->SetOption("text");
//   // h_hough_pi_tof_false->SetMarkerSize(1.6);
//   // h_hough_pi_tof_false->Write();

//   // h_hough_pi_tof_true_dc_true->SetOption("text");
//   // h_hough_pi_tof_true_dc_true->SetMarkerSize(1.6);
//   // h_hough_pi_tof_true_dc_true->Write();


//   // h_hough_2d_optimal->Write();
//   // h_hough_optimal->Write();
//   out_file->cd();

// //   h_mass_invariant_cut_mass->Write();
// //   h_mass_invariant_best_cut_mass->Write();

// //   h_mass_invariant_cut_vertex->Write();
// //   h_mass_invariant_best_cut_vertex->Write();

// //  h_dy->Write();

// //  vertex->cd();
// // h_vtx_x->Write();
// // h_vtx_y->Write();
// // h_vtx_z->Write();
// // h_vtx_x_cut->Write();
// // h_vtx_y_cut->Write();
// // h_vtx_z_cut->Write();



//  for(int j=0;j<3;j++)
//    {
//      pv[j]->cd();
     
//      h_mass_invariant_cut_mass2[j]->SetDirectory(pv[j]);
//      h_mass_invariant_cut_mass2[j]->Write();
//      h_mass_invariant_best_cut_mass2[j]->SetDirectory(pv[j]);
//      h_mass_invariant_best_cut_mass2[j]->Write();
     
//      h_mass_invariant_cut_vertex2[j]->SetDirectory(pv[j]);
//      h_mass_invariant_cut_vertex2[j]->Write();
//      h_mass_invariant_best_cut_vertex2[j]->SetDirectory(pv[j]);
//      h_mass_invariant_best_cut_vertex2[j]->Write();
     

//      for(int k=0;k<4;k++)
//        {
// 	 momD1D2[j][k]->cd();
// 	 h_momD1_momD2_cut_mass[j][k]->SetMarkerColor(k+1);
// 	 h_momD1_momD2_cut_mass[j][k]->SetFillColor(k+1);
// 	 h_momD1_momD2_cut_mass[j][k]->Write();
// 	 h_momD1_momD2_best_cut_mass[j][k]->SetMarkerColor(k+1);
// 	 h_momD1_momD2_best_cut_mass[j][k]->SetFillColor(k+1);
// 	 h_momD1_momD2_best_cut_mass[j][k]->Write();
// 	 h_momD1_momD2_cut_vertex[j][k]->SetMarkerColor(k+1);
// 	 h_momD1_momD2_cut_vertex[j][k]->SetFillColor(k+1);
// 	 h_momD1_momD2_cut_vertex[j][k]->Write();
// 	 h_momD1_momD2_best_cut_vertex[j][k]->SetMarkerColor(k+1);
// 	 h_momD1_momD2_best_cut_vertex[j][k]->SetFillColor(k+1);
// 	 //std::cout<<j<<" "<<k<<" "<<h_momD1_momD2_best_cut_vertex[j][k]<<std::endl;
// 	 h_momD1_momD2_best_cut_vertex[j][k]->Write();
//        }

// //	 h_momD1_momD2_best_cut_vertex[0][0]->Write();

//      pv[j]->cd();

//      // h_vtx_x2[j]->SetDirectory(pv[j]);
//      // h_vtx_x2[j]->Write();
//      // h_vtx_y2[j]->SetDirectory(pv[j]);
//      // h_vtx_y2[j]->Write();
//      // h_vtx_z2[j]->SetDirectory(pv[j]);
//      // h_vtx_z2[j]->Write();
//      // h_vtx_x_cut2[j]->SetDirectory(pv[j]);
//      // h_vtx_x_cut2[j]->Write();
//      // h_vtx_y_cut2[j]->SetDirectory(pv[j]);
//      // h_vtx_y_cut2[j]->Write();
//      // h_vtx_z_cut2[j]->SetDirectory(pv[j]);
//      // h_vtx_z_cut2[j]->Write();
     
//      for(int k=0;k<4;k++)
//        {
// 	 // h_dmom_kalman[j][k]->SetDirectory(pv[j]);
// 	 // h_dmom_kalman[j][k]->Write();
// 	 // h_dmom_normal_kalman[j][k]->SetDirectory(pv[j]);
// 	 // h_dmom_normal_kalman[j][k]->Write();
// 	 // h_dmom_mom_kalman[j][k]->SetDirectory(pv[j]);
// 	 // h_dmom_mom_kalman[j][k]->Write();
	 
// 	 h_Mass2[j][k]->SetDirectory(pv[j]);
// 	 h_Mass2[j][k]->Write();
//        }

//    }

//  return 0;
// }
