#include "Ana_Hist.hh"
#include "Riostream.h"
#include "TFile.h"

#include "TGeoManager.h"
using namespace std;

Ana_Hist::~Ana_Hist()
{
  // cout<<"~Ana_Hist()"<<endl;
  if(HaveBeenWritten == false)
    {
      TFile* f = new TFile("Ana_Hist_Default.root", "RECREATE");
      Write(f);
      f->Close();
      f->Delete();
      f = nullptr;
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

  cout << "Creating Histograms ";
  for(int i = 0; i < SIZEOF_STATEHIST; ++i)
    if(EnableState[i] == true)
      cout << "Hist [" << i << "] ";

  HaveBeenWritten = false;

  h_stats = new TH1I("stats", "stats", 10, 0, 10);
  h_task_exit = new TH1I("h_task_exit", "h_task_exit", 20, 0, 20);

  HistRegisteredByDir.insert(std::make_pair("stat", std::vector<TH1*>({h_stats, h_task_exit})));

  if(EnableState[DAF])
    {
      std::vector<TH1*> HistReg;
      h_pv = new TH1F("total_pvalue", "total_pvalue", 200, 0, 1);
      HistReg.emplace_back(h_pv);
      h_chi2 = new TH1F("total chi2", "total chi2", 400, 0, 40);
      HistReg.emplace_back(h_chi2);

      hd_pv[0] = new TH1F("AAD_Pval_TOFp", "AAD_Pval_TOFp", 500, 0, 1);
      HistReg.emplace_back(hd_pv[0]);
      hd_pv[1] = new TH1F("AAD_Pval_BigTof", "AAD_Pval_BigTof", 500, 0, 1);
      HistReg.emplace_back(hd_pv[1]);
      hd_chi[0] = new TH1F("AAD_Chi2_TOFp", "AAD_Chi2_TOFp", 500, 0, 100);
      HistReg.emplace_back(hd_chi[0]);
      hd_chi[1] = new TH1F("AAD_Chi2_BigTof", "AAD_Chi2_BigTof", 500, 0, 100);
      HistReg.emplace_back(hd_chi[0]);

      h_Path = new TH1F("path_length", "path_length", 200, 400., 800.);
      HistReg.emplace_back(h_Path);
      h_Path_Back = new TH1F("path_length_back", "path_length_back", 200, 400., 800.);
      HistReg.emplace_back(h_Path_Back);
      h_MeanPath = new TH1F("Mean_path_length", "Mean_path_length", 200, 400., 800.);
      HistReg.emplace_back(h_MeanPath);
      h_dpath = new TH1F("Dpath", "Dpath", 200, -10., 10.);
      HistReg.emplace_back(h_dpath);

      h_beta = new TH1F("beta", "beta", 300, 0., 3.);
      HistReg.emplace_back(h_beta);
      h_beta2 = new TH1F("beta_pvcut_0.75", "beta_pvcut_0.75", 300, 0., 3.);
      HistReg.emplace_back(h_beta2);
      h_beta3 = new TH1F("beta_pvcut_0.5", "beta_pvcut_0.5", 300, 0., 3.);
      HistReg.emplace_back(h_beta3);

      h_Mass_All = new TH1F("mass_all", "mass_all", 500, 0, 10.);
      HistReg.emplace_back(h_Mass_All);
      h_Mass_All2 = new TH1F("mass_all_pvcut_0.75", "mass_all_pvcut_0.75", 500, 0, 10.);
      HistReg.emplace_back(h_Mass_All2);
      h_Mass_All3 = new TH1F("mass_all_pvcut_0.5", "mass_all_pvcut_0.5", 500, 0, 10.);
      HistReg.emplace_back(h_Mass_All3);

      h_Mass_charge_All = new TH2F("mass_charge_all", "mass_charge_all", 200, 0., 10., 10, -5, 5);
      HistReg.emplace_back(h_Mass_charge_All);
      h_Mass_charge_All2 = new TH2F("mass_charge_all_pvcut_0.75", "mass_charge_all_pvcut_0.75", 200, 0., 10., 10, -5, 5);
      HistReg.emplace_back(h_Mass_charge_All2);
      h_Mass_charge_All3 = new TH2F("mass_charge_all_pvcut_0.5", "mass_charge_all_pvcut_0.5", 200, 0., 10., 10, -5, 5);
      HistReg.emplace_back(h_Mass_charge_All);

      h_beta_mom = new TH2F("beta_mom", "beta_mom", 500, 0., 15., 200, 0., 2.5);
      HistReg.emplace_back(h_beta_mom);
      h_beta_mom2 = new TH2F("beta_mom_pvcut_0.75", "beta_mom_pvcut_0.75", 500, 0., 15., 100, 0., 2.5);
      HistReg.emplace_back(h_beta_mom2);
      h_beta_mom3 = new TH2F("beta_mom_pvcut_0.5", "beta_mom_pvcut_0.5", 500, 0., 15., 100, 0., 2.5);
      HistReg.emplace_back(h_beta_mom3);

      h_pv_mom = new TH2F("mom_pv", "mom_pv", 200, 0., 20., 250, 0, 1);
      HistReg.emplace_back(h_pv_mom);
      h_pv_beta = new TH2F("beta_pv", "beta_pv", 250, 0., 2.5, 250, 0, 1);
      HistReg.emplace_back(h_pv_beta);
      h_pv_mass = new TH2F("mass_pv", "mass_pv", 200, 0., 10., 250, 0, 1);
      HistReg.emplace_back(h_pv_mass);

      h_path_tof = new TH2F("Mean_path_TOF", "Mean_path_TOF", 500, 15., 20., 90, 10, 40);
      HistReg.emplace_back(h_path_tof);

      h_mom_tof_cut = new TH2F("mom_tof_pvcut_0.5", "mom_tof_pvcut_0.5", 200, 0., 20., 500, 0., 50.);
      HistReg.emplace_back(h_mom_tof_cut);
      h_path_tof_cut = new TH2F("path_tof_pvcut_0.5", "path_tof_pvcut_0.5", 500, 15., 20., 500, 0., 50.);
      HistReg.emplace_back(h_path_tof_cut);
      h_path_mom_cut = new TH2F("path_mom_pvcut", "path_mom_pvcut", 500, 15., 20., 200, 0, 20.);
      HistReg.emplace_back(h_path_mom_cut);

      HistRegisteredByDir.insert(std::make_pair("Kalman", HistReg));

      std::string name[4] = {"_proton", "_alpha", "_pi-", "_he3"};
      for(int i = 0; i < 4; i++)
        {
          std::vector<TH1*> HistRegPID;
          std::string temp;
          temp = "mass";
          temp += name[i];
          h_Mass[i] = new TH1F(temp.c_str(), temp.c_str(), 10000, 0, 10.);
          HistRegPID.emplace_back(h_Mass[i]);

          temp = "chi2_";
          temp += name[i];
          h_chi2_particle[i] = new TH1F(temp.c_str(), temp.c_str(), 1000, 0, 10);
          HistRegPID.emplace_back(h_chi2_particle[i]);

          temp = "pv_";
          temp += name[i];
          h_pv_particle[i] = new TH1F(temp.c_str(), temp.c_str(), 100, 0, 1);
          HistRegPID.emplace_back(h_pv_particle[i]);

          temp = "Kalman";
          temp += name[i];
          HistRegisteredByDir.insert(std::make_pair(temp, HistRegPID));
        }

      std::vector<std::string> name_res = {"_p","_piN","_piP","_KN","_KP"};
      int id = 0;
      std::vector<TH1*> HistRegRes;
      for(auto temp_name : name_res)
	{
	  std::string temp = "mom_res";
	  temp += temp_name;
	  h_mom_res[id] = new TH2F(temp.c_str(), temp.c_str(), 25, 0., 10., 1000, 0., 0.1);
	  HistRegRes.emplace_back(h_mom_res[id]);
	  ++id;
	}
      HistRegisteredByDir.insert(std::make_pair("MomRes", HistRegRes));
      
    }

  if(EnableState[FINDING])
    {
      std::vector<TH1*> HistReg;

      h_xy = new TH2F("h_xy", "h_xy", 500, -50, 50, 500, -50, 50);
      HistReg.emplace_back(h_xy);

      h_xy_extrap = new TH2F("h_xy_extrap", "h_xy_extrap", 500, -50, 50, 500, -50, 50);
      HistReg.emplace_back(h_xy_extrap);

      h_PxPy = new TH2F("h_PxPy", "h_PxPy", 200, -1, 1, 200, -1, 1);
      HistReg.emplace_back(h_PxPy);

      h_PxPy_extrap = new TH2F("h_PxPy_extrap", "h_PxPy_extrap", 200, -1, 1, 200, -1, 1);
      HistReg.emplace_back(h_PxPy_extrap);

      h_TrackFindingStat = new TH2F("h_TrackFindingStat", "h_TrackFindingStat", 20, 0, 20, 22, -2, 20);
      HistReg.emplace_back(h_TrackFindingStat);

      HistRegisteredByDir.insert(std::make_pair("Finder", HistReg));
    }

  std::cout << " : done !" << std::endl;
}

int Ana_Hist::Write(TFile* out_file)
{

  HaveBeenWritten = true;
  out_file->cd();

  std::cout << "making directory " << std::endl;
  for(auto it : HistRegisteredByDir)
    {
      TDirectory* temp_dir = GetDir(out_file, it.first);
      temp_dir->cd();
      for(auto it_hist : it.second)
        it_hist->Write();
      out_file->cd();
    }

  return 0;
}

int Ana_Hist::WriteTemp(char* tempfile)
{
  HaveBeenWritten = true;

  TFile* ff = new TFile(tempfile, "RECREATE");
  cout << "File = " << tempfile << endl;

  ff->cd();

  for(auto it : HistRegisteredByDir)
    {
      TDirectory* temp_dir = GetDir(ff, it.first);
      temp_dir->cd();
      for(auto it_hist : it.second)
        it_hist->Write();
      ff->cd();
    }

  // TDirectory *bmtp = ff->mkdir("BetaMom_Tof+");
  // TDirectory *bmbg = ff->mkdir("BetaMom_BigTof+");
  // TDirectory *Mother[3];//
  // Mother[0] =(TDirectory *) ff->mkdir("Lambda");
  // Mother[1] =(TDirectory *) ff->mkdir("H3L");
  // Mother[2] =(TDirectory *) ff->mkdir("H4L");
  // TDirectory *invpv[3][3];
  // for(int mm=0;mm<3;mm++)
  //   {
  //     invpv[mm][0] = (TDirectory *)Mother[mm]->mkdir("pv_tight");
  //     invpv[mm][1] = (TDirectory *)Mother[mm]->mkdir("pv_mid");
  //     invpv[mm][2] = (TDirectory *)Mother[mm]->mkdir("pv_loose");
  //   }

  // //#ifdef ADD_HISTO

  // h_MassZ1BarBest->Write("",TObject::kOverwrite);

  // h_MassZ1BarBestTofs->Write("",TObject::kOverwrite);

  // hd_theWorld->Write("",TObject::kOverwrite);
  // hd_chi[0] ->Write("",TObject::kOverwrite);
  // hd_chi[1] ->Write("",TObject::kOverwrite);
  // hd_pv[0] ->Write("",TObject::kOverwrite);
  // hd_pv[1] ->Write("",TObject::kOverwrite);
  // hd_TrackFrom[0]->Write("",TObject::kOverwrite);
  // hd_TrackFrom[1]->Write("",TObject::kOverwrite);
  // hd_MassPv[0]->Write("",TObject::kOverwrite);
  // hd_MassPv[1]->Write("",TObject::kOverwrite);
  // for(auto & elem : h_MomCor)
  //   {
  //     elem->Write("",TObject::kOverwrite);
  //   }
  // for(auto & elem : hd_TimeCalTp)
  //   {
  //     bmtp->cd();

  //     elem->Write("",TObject::kOverwrite);

  //   }
  // ff->cd();
  // hd_BetaMom_Tp_all->Write("",TObject::kOverwrite);
  // hd_BetaMom_Tp_allTofs->Write("",TObject::kOverwrite);
  // hd_BetaMom_Tp_allTofs_Z2->Write("",TObject::kOverwrite);

  // hd_BetaMom_Tp_ini_all->Write("",TObject::kOverwrite);
  // hd_BetaMom_Tp_ini_allTofs->Write("",TObject::kOverwrite);
  // hd_BetaMom_Tp_kal_allTofs->Write("",TObject::kOverwrite);
  // hd_BetaMom_Tp_ini_allTofs_Z2->Write("",TObject::kOverwrite);

  // for(int dd=0;dd<18;dd++)
  //   {
  //     bmbg->cd();

  //   }
  // ff->cd();
  // hd_BetaMom_Bg_all->Write("",TObject::kOverwrite);

  // hd_BetaMom_Bg_allTofs->Write("",TObject::kOverwrite);
  // hd_BetaMom_Bg_ini_all->Write("",TObject::kOverwrite);
  // hd_BetaMom_Bg_ini_allTofs->Write("",TObject::kOverwrite);

  // ff->cd();
  // hd_IniMom_NewMom_all_tofsZ1->Write("",TObject::kOverwrite);
  // hd_IniMom_NewMom_all_tofsZ1_ini->Write("",TObject::kOverwrite);
  // //hd_Bar_DiffMom_allTofs->Write("",TObject::kOverwrite);
  // hd_IniNewMomdiff_vs_IniNewMomplus_allTofs->Write("",TObject::kOverwrite);
  // hd_charge_ini_vs_IniNewDiv->Write("",TObject::kOverwrite);
  // hd_charge_vs_IniNewDiv->Write("",TObject::kOverwrite);
  // hd_IniNewMomdiff_vs_IniNewMomplus_allTofs_Z1->Write("",TObject::kOverwrite);
  // hd_IniNewMomdiff_vs_IniNewMomplus_allTofs_ini_Z1->Write("",TObject::kOverwrite);
  // hd_best_pos_hitvstrack_tofpZ1->Write("",TObject::kOverwrite);
  // hd_best_pos_hitdifftrack_tofpZ1->Write("",TObject::kOverwrite);

  ff->Close();
  delete ff;

  return 0;
}
