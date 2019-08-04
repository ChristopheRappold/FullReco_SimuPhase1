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

  h_stats = new TH1I("stats", "stats", 20, 0, 20);
  h_statsLess3Mes = new TH2I("statsLess3", "statsLess3", 40, 0, 40,40,0,40);
  h_statsInvalid = new TH2I("statsInvalid", "statsInvalid", 40, 0, 40,40,0,40);
  h_task_exit = new TH1I("h_task_exit", "h_task_exit", 20, 0, 20);

  HistRegisteredByDir.insert(std::make_pair("stat", std::make_tuple(std::vector<TH1*>({h_stats, h_task_exit, h_statsLess3Mes, h_statsInvalid}),0)));

  std::vector<std::string> name_field = {"Bx","By","Bz"};
  std::vector<TH1*> h_fields;
  for(size_t i=0;i<name_field.size();++i)
    {
      std::string name_h_field = "FieldXY_";
      name_h_field+=name_field[i];
      FieldXY[i] = new TH2F(name_h_field.c_str(),name_h_field.c_str(),600,-150,150,600,-150,150);
      h_fields.emplace_back(FieldXY[i]);

      std::string name_h_field1 = "FieldXZ_";
      name_h_field1+=name_field[i];
      FieldXZ[i] = new TH2F(name_h_field1.c_str(),name_h_field1.c_str(),600,-150,150,1200,-300,300);
      h_fields.emplace_back(FieldXZ[i]);

      std::string name_h_field2 = "FieldYZ_";
      name_h_field2+=name_field[i];
      FieldYZ[i] = new TH2F(name_h_field2.c_str(),name_h_field2.c_str(),600,-150,150,1200,-300,300);
      h_fields.emplace_back(FieldYZ[i]);

      std::string name_h_field_n = "FieldXYn_";
      name_h_field_n+=name_field[i];
      FieldXY_n[i] = new TH2F(name_h_field_n.c_str(),name_h_field_n.c_str(),600,-150,150,600,-150,150);
      h_fields.emplace_back(FieldXY_n[i]);

      std::string name_h_field_n1 = "FieldXZn_";
      name_h_field_n1+=name_field[i];
      FieldXZ_n[i] = new TH2F(name_h_field_n1.c_str(),name_h_field_n1.c_str(),600,-150,150,1200,-300,300);
      h_fields.emplace_back(FieldXZ_n[i]);

      std::string name_h_field_n2 = "FieldYZn_";
      name_h_field_n2+=name_field[i];
      FieldYZ_n[i] = new TH2F(name_h_field_n2.c_str(),name_h_field_n2.c_str(),600,-150,150,1200,-300,300);
      h_fields.emplace_back(FieldYZ_n[i]);

    }
  
  HistRegisteredByDir.insert(std::make_pair("Field", std::make_tuple(h_fields,0)));
  
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

      h_Path = new TH1F("path_length", "path_length", 200, 0., 400.);
      HistReg.emplace_back(h_Path);
      h_Path_Back = new TH1F("path_length_back", "path_length_back", 200, 0., 400.);
      HistReg.emplace_back(h_Path_Back);
      h_MeanPath = new TH1F("Mean_path_length", "Mean_path_length", 200, 0., 400.);
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
      
      h_total_dE = new TH2F("total_dE","total_dE",30,0,30,500,0,100);
      HistReg.emplace_back(h_total_dE);
      
      HistRegisteredByDir.insert(std::make_pair("Kalman", std::make_tuple(HistReg,0)));

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
          HistRegisteredByDir.insert(std::make_pair(temp, std::make_tuple(HistRegPID,0)));
        }

      std::vector<std::string> name_res = {"_p","_piN","_piP","_KN","_KP"};
      
      std::vector<std::tuple<int,double,double> > list_bins= {std::make_tuple(500,-1.,1.), std::make_tuple(500,-0.1,0.1), std::make_tuple(500,-0.1,0.1),
						     std::make_tuple(500,-1.,1.), std::make_tuple(500,-1.,1.),
						     std::make_tuple(200,-6.,6.), std::make_tuple(200,-6.,6.), std::make_tuple(200,-6.,6.),
						     std::make_tuple(200,-6.,6.), std::make_tuple(200,-6.,6.),
      };
      std::vector<std::string> name_histres = {"h_momRes","h_upRes","h_vpRes","h_uRes","h_vRes","h_qopPull","h_upPull","h_vpPull","h_uPull","h_vPull"};

      std::vector<TH1*> HistRegRes1;
      std::vector<TH1*> HistRegRes2;
      for(size_t id = 0; id<name_res.size();++id)
	{
	  std::string tempResHist = "mom_res";
	  auto temp_nameRes = name_res[id];
	  tempResHist += temp_nameRes;
	  h_mom_res[id] = new TH2F(tempResHist.c_str(), tempResHist.c_str(), 25, 0., 5., 50, -1, 1);
	  HistRegRes1.emplace_back(h_mom_res[id]);

	  for(size_t id_h = 0; id_h < name_histres.size();++id_h)
	    {
	      std::string temp1_ResHist = name_histres[id_h];
	      temp1_ResHist+= temp_nameRes;
	      int nBin;
	      double minB, maxB;
	      std::tie(nBin,minB,maxB) = list_bins[id_h];
	      h_ResPull[id][id_h] = new TH1D(temp1_ResHist.c_str(),temp1_ResHist.c_str(),nBin, minB, maxB);
	      HistRegRes2.emplace_back(h_ResPull[id][id_h]);
	      temp1_ResHist+="2D";
	      h_ResPull_normal[id][id_h] = new TH2F(temp1_ResHist.c_str(),temp1_ResHist.c_str(),100,0,10,nBin, minB, maxB);
	      HistRegRes2.emplace_back(h_ResPull_normal[id][id_h]);
	    }
	}
      
      HistRegisteredByDir.insert(std::make_pair("MomRes", std::make_tuple(HistRegRes1,1)));
      HistRegisteredByDir.insert(std::make_pair("ResPull", std::make_tuple(HistRegRes2,0)));
      
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

      HistRegisteredByDir.insert(std::make_pair("Finder", std::make_tuple(HistReg,0)));
    }

  std::cout << " : done !" << std::endl;
}

int Ana_Hist::Write(TFile* out_file)
{

  HaveBeenWritten = true;
  out_file->cd();

  auto f_DiffResolution = [](TH2F* h)
			  {
			    TGraphErrors* g_resMean = nullptr;
			    TGraphErrors* g_resStd = nullptr;

			    if(h==nullptr)
			      return std::make_tuple(g_resMean, g_resStd);
			    
			    TString nameHist1(h->GetName());
			    nameHist1 += "Mean";
			    TString nameHist2(h->GetName());
			    nameHist2 += "Std";

			    g_resMean = new TGraphErrors;
			    g_resStd = new TGraphErrors;

			    g_resMean->SetNameTitle(nameHist1,nameHist1);
			    g_resStd->SetNameTitle(nameHist2,nameHist2);

			    for(int i_x = 1; i_x <= h->GetNbinsX(); ++i_x)
			      {
				TString nameH("Proj_");
				nameH+=i_x;
			        TH1D* h_temp = h->ProjectionY(nameH,i_x,i_x);

				double mean = h_temp->GetMean();
				double mean_err = h_temp->GetMeanError();
				double rms = h_temp->GetRMS();
				double rms_err = h_temp->GetRMSError();

				double x_c = h->GetXaxis()->GetBinCenter(i_x);
				double x_w = h->GetXaxis()->GetBinWidth(i_x);
				
				h_temp->GetXaxis()->SetRangeUser(mean-rms*2.,mean+rms*2.);
				
				mean = h_temp->GetMean();
				mean_err = h_temp->GetMeanError();
				rms = h_temp->GetRMS();
				rms_err = h_temp->GetRMSError();
				
				g_resMean->SetPoint(i_x-1, x_c, mean);
				g_resMean->SetPointError(i_x-1,x_w*0.5,mean_err);
				g_resStd->SetPoint(i_x-1, x_c, rms);
				g_resStd->SetPointError(i_x-1,x_w*0.5,rms_err);
			      }
			    return std::make_tuple(g_resMean, g_resStd);			    
			  };

  
  std::cout << "making directory " << std::endl;
  for(auto it : HistRegisteredByDir)
    {
      TDirectory* temp_dir = GetDir(out_file, it.first);
      temp_dir->cd();
      
      for(auto it_hist : std::get<0>(it.second))
        {
	  it_hist->Write();
	  if(std::get<1>(it.second)==1)
	    {
	      TGraphErrors* g1;
	      TGraphErrors* g2;
	      
	      std::tie(g1,g2) = f_DiffResolution(dynamic_cast<TH2F*>(it_hist));
	      if(g1!=nullptr)
		g1->Write();
	      if(g2!=nullptr)
		g2->Write();
	    }
	}
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
  Write(ff);
  // for(auto it : HistRegisteredByDir)
  //   {
  //     TDirectory* temp_dir = GetDir(ff, it.first);
  //     temp_dir->cd();
  //     for(auto it_hist : it.second)
  //       it_hist->Write();
  //     ff->cd();
  //   }

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
