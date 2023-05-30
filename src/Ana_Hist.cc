#include "Ana_Hist.hh"
//#include "Riostream.h"
#include "TFile.h"
#include "TGeoManager.h"

#include <sstream>
#include "spdlog/spdlog.h"

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
Ana_Hist::Ana_Hist(bool Daf, bool Vertex, bool DCproject, bool Finding, bool Riemann, bool Hough, bool Simu,
                    bool Builder, bool PrimVtx, bool PrimVtx_Si, bool DecayVtx, bool DecayVtx_piplus,
                    bool FragmentFinder, bool WASAFinder)
{
  EnableState.resize(SIZEOF_STATEHIST);
  EnableState[DAF] = Daf;
  EnableState[VERTEX] = Vertex;
  EnableState[DCPROJ] = DCproject;
  EnableState[FINDING] = Finding;
  EnableState[RIEMANN] = Riemann;
  EnableState[HOUGH] = Hough;
  EnableState[SIMU] = Simu;
  EnableState[BUILDER] = Builder;
  EnableState[PRIMVTX] = PrimVtx;
  EnableState[PRIMVTX_SI] = PrimVtx_Si;
  EnableState[DECAYVTX] = DecayVtx;
  EnableState[DECAYVTX_PIPLUS] = DecayVtx_piplus;
  EnableState[FRAGMENT] = FragmentFinder;
  EnableState[WASA] = WASAFinder;



  _logger = spdlog::get("Console");

  stringstream ss;
  ss<< "Creating Histograms ";
  for(int i = 0; i < SIZEOF_STATEHIST; ++i)
    if(EnableState[i] == true)
      ss << "Hist [" << i << "] ";
  _logger->info(ss.str());
  
  HaveBeenWritten = false;

  std::vector< std::vector<TH1*>* > Stats;
  h_stats.emplace_back(new TH1I("stats", "stats", 20, 0, 20));
  Stats.emplace_back(&h_stats.store);
  h_statsLess3Mes.emplace_back(new TH2I("statsLess3", "statsLess3", 40, 0, 40,40,0,40));
  Stats.emplace_back(&h_statsLess3Mes.store);
  h_statsInvalid.emplace_back(new TH2I("statsInvalid", "statsInvalid", 40, 0, 40,40,0,40));
  Stats.emplace_back(&h_statsInvalid.store);
  h_task_exit.emplace_back(new TH1I("h_task_exit", "h_task_exit", 20, 0, 20));
  Stats.emplace_back(&h_task_exit.store);
  
  HistRegisteredByDir.insert(std::make_pair("stat", std::make_tuple(Stats,0)));

  std::vector<std::string> name_field = {"Bx","By","Bz"};
  std::vector<std::vector<TH1*>*> h_fields;
  for(size_t i=0;i<name_field.size();++i)
    {
      std::string name_h_field = "FieldXY_";
      name_h_field+=name_field[i];
      FieldXY[i].emplace_back(new TH2F(name_h_field.c_str(),name_h_field.c_str(),300,-150,150,300,-150,150));
      h_fields.emplace_back(&FieldXY[i].store);

      std::string name_h_field1 = "FieldXZ_";
      name_h_field1+=name_field[i];
      FieldXZ[i].emplace_back(new TH2F(name_h_field1.c_str(),name_h_field1.c_str(),300,-150,150,600,-300,300));
      h_fields.emplace_back(&FieldXZ[i].store);

      std::string name_h_field2 = "FieldYZ_";
      name_h_field2+=name_field[i];
      FieldYZ[i].emplace_back(new TH2F(name_h_field2.c_str(),name_h_field2.c_str(),300,-150,150,600,-300,300));
      h_fields.emplace_back(&FieldYZ[i].store);

      std::string name_h_fieldm = "FieldXY_max";
      name_h_fieldm+=name_field[i];
      FieldXYmax[i].emplace_back(new TH2F(name_h_fieldm.c_str(),name_h_fieldm.c_str(),300,-150,150,300,-150,150));
      h_fields.emplace_back(&FieldXYmax[i].store);

      std::string name_h_field1m = "FieldXZ_max";
      name_h_field1m+=name_field[i];
      FieldXZmax[i].emplace_back(new TH2F(name_h_field1m.c_str(),name_h_field1m.c_str(),300,-150,150,600,-300,300));
      h_fields.emplace_back(&FieldXZmax[i].store);

      std::string name_h_field2m = "FieldYZ_max";
      name_h_field2m+=name_field[i];
      FieldYZmax[i].emplace_back(new TH2F(name_h_field2m.c_str(),name_h_field2m.c_str(),300,-150,150,600,-300,300));
      h_fields.emplace_back(&FieldYZmax[i].store);

      std::string name_h_field_n = "FieldXYn_";
      name_h_field_n+=name_field[i];
      FieldXY_n[i].emplace_back(new TH2F(name_h_field_n.c_str(),name_h_field_n.c_str(),300,-150,150,300,-150,150));
      h_fields.emplace_back(&FieldXY_n[i].store);

      std::string name_h_field_n1 = "FieldXZn_";
      name_h_field_n1+=name_field[i];
      FieldXZ_n[i].emplace_back(new TH2F(name_h_field_n1.c_str(),name_h_field_n1.c_str(),300,-150,150,600,-300,300));
      h_fields.emplace_back(&FieldXZ_n[i].store);

      std::string name_h_field_n2 = "FieldYZn_";
      name_h_field_n2+=name_field[i];
      FieldYZ_n[i].emplace_back(new TH2F(name_h_field_n2.c_str(),name_h_field_n2.c_str(),300,-150,150,600,-300,300));
      h_fields.emplace_back(&FieldYZ_n[i].store);

    }
  
  HistRegisteredByDir.insert(std::make_pair("Field", std::make_tuple(h_fields,0)));
  
  if(EnableState[DAF])
    {
      std::vector<std::vector<TH1*>*> HistReg;
      h_pv.emplace_back(new TH1F("total_pvalue", "total_pvalue", 200, 0, 1));
      HistReg.emplace_back(&h_pv.store);
      h_chi2.emplace_back(new TH1F("total chi2", "total chi2", 400, 0, 40));
      HistReg.emplace_back(&h_chi2.store);

      hd_pv[0].emplace_back(new TH1F("AAD_Pval_TOFp", "AAD_Pval_TOFp", 500, 0, 1));
      HistReg.emplace_back(&hd_pv[0].store);
      hd_pv[1].emplace_back(new TH1F("AAD_Pval_BigTof", "AAD_Pval_BigTof", 500, 0, 1));
      HistReg.emplace_back(&hd_pv[1].store);
      hd_chi[0].emplace_back(new TH1F("AAD_Chi2_TOFp", "AAD_Chi2_TOFp", 500, 0, 100));
      HistReg.emplace_back(&hd_chi[0].store);
      hd_chi[1].emplace_back(new TH1F("AAD_Chi2_BigTof", "AAD_Chi2_BigTof", 500, 0, 100));
      HistReg.emplace_back(&hd_chi[1].store);

      h_Path.emplace_back(new TH1F("path_length", "path_length", 200, 0., 400.));
      HistReg.emplace_back(&h_Path.store);
      h_Path_Back.emplace_back(new TH1F("path_length_back", "path_length_back", 200, 0., 400.));
      HistReg.emplace_back(&h_Path_Back.store);
      h_MeanPath.emplace_back(new TH1F("Mean_path_length", "Mean_path_length", 200, 0., 400.));
      HistReg.emplace_back(&h_MeanPath.store);
      h_dpath.emplace_back(new TH1F("Dpath", "Dpath", 200, -10., 10.));
      HistReg.emplace_back(&h_dpath.store);

      h_beta.emplace_back(new TH1F("beta", "beta", 300, 0., 3.));
      HistReg.emplace_back(&h_beta.store);
      h_beta2.emplace_back(new TH1F("beta_pvcut_0.75", "beta_pvcut_0.75", 300, 0., 3.));
      HistReg.emplace_back(&h_beta2.store);
      h_beta3.emplace_back(new TH1F("beta_pvcut_0.5", "beta_pvcut_0.5", 300, 0., 3.));
      HistReg.emplace_back(&h_beta3.store);

      h_Mass_All.emplace_back(new TH1F("mass_all", "mass_all", 500, 0, 10.));
      HistReg.emplace_back(&h_Mass_All.store);
      h_Mass_All2.emplace_back(new TH1F("mass_all_pvcut_0.75", "mass_all_pvcut_0.75", 500, 0, 10.));
      HistReg.emplace_back(&h_Mass_All2.store);
      h_Mass_All3.emplace_back(new TH1F("mass_all_pvcut_0.5", "mass_all_pvcut_0.5", 500, 0, 10.));
      HistReg.emplace_back(&h_Mass_All3.store);

      h_Mass_charge_All.emplace_back(new TH2F("mass_charge_all", "mass_charge_all", 200, 0., 10., 10, -5, 5));
      HistReg.emplace_back(&h_Mass_charge_All.store);
      h_Mass_charge_All2.emplace_back(new TH2F("mass_charge_all_pvcut_0.75", "mass_charge_all_pvcut_0.75", 200, 0., 10., 10, -5, 5));
      HistReg.emplace_back(&h_Mass_charge_All2.store);
      h_Mass_charge_All3.emplace_back(new TH2F("mass_charge_all_pvcut_0.5", "mass_charge_all_pvcut_0.5", 200, 0., 10., 10, -5, 5));
      HistReg.emplace_back(&h_Mass_charge_All.store);

      h_beta_mom.emplace_back(new TH2F("beta_mom", "beta_mom", 500, 0., 15., 200, 0., 2.5));
      HistReg.emplace_back(&h_beta_mom.store);
      h_beta_mom2.emplace_back(new TH2F("beta_mom_pvcut_0.75", "beta_mom_pvcut_0.75", 500, 0., 15., 100, 0., 2.5));
      HistReg.emplace_back(&h_beta_mom2.store);
      h_beta_mom3.emplace_back(new TH2F("beta_mom_pvcut_0.5", "beta_mom_pvcut_0.5", 500, 0., 15., 100, 0., 2.5));
      HistReg.emplace_back(&h_beta_mom3.store);
      h_beta_momcharge.emplace_back(new TH2F("beta_momcharge", "beta_momcharge", 1000, -15., 15., 200, 0., 2.5));
      HistReg.emplace_back(&h_beta_momcharge.store);
      h_beta_momcharge2.emplace_back(new TH2F("beta_momcharge_pvcut_0.75", "beta_momcharge_pvcut_0.75", 1000, -15., 15., 100, 0., 2.5));
      HistReg.emplace_back(&h_beta_momcharge2.store);
      h_beta_momcharge3.emplace_back(new TH2F("beta_momcharge_pvcut_0.5", "beta_momcharge_pvcut_0.5", 1000, -15., 15., 100, 0., 2.5));
      HistReg.emplace_back(&h_beta_momcharge3.store);

      h_pv_mom.emplace_back(new TH2F("mom_pv", "mom_pv", 200, 0., 20., 250, 0, 1));
      HistReg.emplace_back(&h_pv_mom.store);
      h_pv_beta.emplace_back(new TH2F("beta_pv", "beta_pv", 250, 0., 2.5, 250, 0, 1));
      HistReg.emplace_back(&h_pv_beta.store);
      h_pv_mass.emplace_back(new TH2F("mass_pv", "mass_pv", 200, 0., 10., 250, 0, 1));
      HistReg.emplace_back(&h_pv_mass.store);

      h_path_tof.emplace_back(new TH2F("Mean_path_TOF", "Mean_path_TOF", 500, 70., 150., 200, 0, 10.));
      HistReg.emplace_back(&h_path_tof.store);

      h_mom_tof_cut.emplace_back(new TH2F("mom_tof_pvcut_0.5", "mom_tof_pvcut_0.5", 200, 0., 20., 500, 0., 50.));
      HistReg.emplace_back(&h_mom_tof_cut.store);
      h_path_tof_cut.emplace_back(new TH2F("path_tof_pvcut_0.5", "path_tof_pvcut_0.5", 500, 15., 20., 500, 0., 50.));
      HistReg.emplace_back(&h_path_tof_cut.store);
      h_path_mom_cut.emplace_back(new TH2F("path_mom_pvcut", "path_mom_pvcut", 500, 15., 20., 200, 0, 20.));
      HistReg.emplace_back(&h_path_mom_cut.store);
      
      h_total_dE.emplace_back(new TH2F("total_dE","total_dE",30,0,30,500,0,100));
      HistReg.emplace_back(&h_total_dE.store);
      
      HistRegisteredByDir.insert(std::make_pair("Kalman", std::make_tuple(HistReg,0)));

      std::string name[4] = {"_proton", "_alpha", "_pi-", "_he3"};
      for(int i = 0; i < 4; i++)
        {
	  std::vector<std::vector<TH1*>*> HistRegPID;
          std::string temp;
          temp = "mass";
          temp += name[i];
          h_Mass[i].emplace_back(new TH1F(temp.c_str(), temp.c_str(), 10000, 0, 10.));
          HistRegPID.emplace_back(&h_Mass[i].store);

          temp = "chi2_";
          temp += name[i];
          h_chi2_particle[i].emplace_back(new TH1F(temp.c_str(), temp.c_str(), 1000, 0, 10));
          HistRegPID.emplace_back(&h_chi2_particle[i].store);

          temp = "pv_";
          temp += name[i];
          h_pv_particle[i].emplace_back(new TH1F(temp.c_str(), temp.c_str(), 100, 0, 1));
          HistRegPID.emplace_back(&h_pv_particle[i].store);

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

      std::vector<std::vector<TH1*>*> HistRegRes1;
      std::vector<std::vector<TH1*>*> HistRegRes2;
      for(size_t id = 0; id<name_res.size();++id)
	{
	  std::string tempResHist = "mom_res";
	  auto temp_nameRes = name_res[id];
	  tempResHist += temp_nameRes;
	  h_mom_res[id].emplace_back(new TH2F(tempResHist.c_str(), tempResHist.c_str(), 25, 0., 5., 50, -1, 1));
	  HistRegRes1.emplace_back(&h_mom_res[id].store);

	  for(size_t id_h = 0; id_h < name_histres.size();++id_h)
	    {
	      std::string temp1_ResHist = name_histres[id_h];
	      temp1_ResHist+= temp_nameRes;
	      int nBin;
	      double minB, maxB;
	      std::tie(nBin,minB,maxB) = list_bins[id_h];
	      h_ResPull[id][id_h].emplace_back(new TH1D(temp1_ResHist.c_str(),temp1_ResHist.c_str(),nBin, minB, maxB));
	      HistRegRes2.emplace_back(&h_ResPull[id][id_h].store);
	      temp1_ResHist+="2D";
	      h_ResPull_normal[id][id_h].emplace_back(new TH2F(temp1_ResHist.c_str(),temp1_ResHist.c_str(),100,0,10,nBin, minB, maxB));
	      HistRegRes2.emplace_back(&h_ResPull_normal[id][id_h].store);
	    }
	}
      
      HistRegisteredByDir.insert(std::make_pair("MomRes", std::make_tuple(HistRegRes1,1)));
      HistRegisteredByDir.insert(std::make_pair("ResPull", std::make_tuple(HistRegRes2,0)));
      
      std::vector<std::vector<TH1*>*> h_residual;
      for(int i=0;i<17;++i)
        for(int j=0; j<3; ++j)
	  {
	    h_ResMDC[i][j].emplace_back( new TH1F(Form("ResMDC%02d_%d",i,j),Form("ResMDC%02d_%d",i,j),300,-1,1) );
	    h_residual.emplace_back(&h_ResMDC[i][j].store);
	  }
      
      for(int i=0;i<9;++i)
	{
	  h_ResFiber[i].emplace_back( new TH1F(Form("ResFiber%d",i),Form("ResFiber%d",i),300,-1,1));
	  h_residual.emplace_back(&h_ResFiber[i].store);
	}
      for(int i=0;i<6;++i)
	{
	  h_ResMiniFiber[i].emplace_back( new TH1F(Form("ResMiniFiber%d",i),Form("ResMiniFiber%d",i),300,-1,1));
	  h_residual.emplace_back(&h_ResMiniFiber[i].store);
	}
      for(int i=0;i<2;++i)
	{
	  h_ResPSCE[i].emplace_back( new TH1F(Form("ResPSCE%d",i),Form("ResPSCE%d",i),300,-15,15));
	  h_residual.emplace_back(&h_ResPSCE[i].store);
	}

      HistRegisteredByDir.insert(std::make_pair("Residual", std::make_tuple(h_residual,0)));
    }

  if(EnableState[FINDING])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      h_xy.emplace_back(new TH2F("h_xy", "h_xy", 500, -50, 50, 500, -50, 50));
      HistReg.emplace_back(&h_xy.store);

      h_xy_extrap.emplace_back(new TH2F("h_xy_extrap", "h_xy_extrap", 500, -50, 50, 500, -50, 50));
      HistReg.emplace_back(&h_xy_extrap.store);

      h_PxPy.emplace_back(new TH2F("h_PxPy", "h_PxPy", 200, -1, 1, 200, -1, 1));
      HistReg.emplace_back(&h_PxPy.store);

      h_PxPy_extrap.emplace_back(new TH2F("h_PxPy_extrap", "h_PxPy_extrap", 200, -1, 1, 200, -1, 1));
      HistReg.emplace_back(&h_PxPy_extrap.store);

      h_TrackFindingStat.emplace_back(new TH2F("h_TrackFindingStat", "h_TrackFindingStat", 20, 0, 20, 22, -2, 20));
      HistReg.emplace_back(&h_TrackFindingStat.store);

      h_MDC_Dphi.emplace_back(new TH2F("h_MDC_Dphi","h_MDC_Dphi",25,0,25,640,-3.2,3.2));
      HistReg.emplace_back(&h_MDC_Dphi.store);
      h_MDC_NextLayer.emplace_back(new TH2F("h_MDC_NextLayer","h_MDC_NextLayer",25,0,25,40,-20,20));
      HistReg.emplace_back(&h_MDC_NextLayer.store);
      h_MDC_DiffLayer.emplace_back(new TH2F("h_MDC_DiffLayer","h_MDC_DiffLayer",25,0,25,500,-1,1));
      HistReg.emplace_back(&h_MDC_DiffLayer.store);
      h_MDC_DiffLayer2.emplace_back(new TH2F("h_MDC_DiffLayer2","h_MDC_DiffLayer2",25,0,25,500,-1,1));
      HistReg.emplace_back(&h_MDC_DiffLayer2.store);
      h_MDC_InLayer.emplace_back(new TH2F("h_MDC_InLayer","h_MDC_InLayer",25,0,25,10,0,10));
      HistReg.emplace_back(&h_MDC_InLayer.store);

      
      h_SolenoidGeo[0].emplace_back(new TH2F("h_SolenoidGeoFront","h_SolenoidGeoFront",500,-50,50,500,-50,50));
      HistReg.emplace_back(&h_SolenoidGeo[0].store);
      h_SolenoidGeo[1].emplace_back(new TH2F("h_SolenoidGeoMid"  ,"h_SolenoidGeoMid"  ,500,-50,50,500,-50,50));
      HistReg.emplace_back(&h_SolenoidGeo[1].store);
      h_SolenoidGeo[2].emplace_back(new TH2F("h_SolenoidGeoBack" ,"h_SolenoidGeoBack" ,500,-50,50,500,-50,50));
      HistReg.emplace_back(&h_SolenoidGeo[2].store);

      h_RZStats.emplace_back(new TH1F("h_RZStats","h_RZStats",10,0,10));
      HistReg.emplace_back(&h_RZStats.store);
      h_RZ.emplace_back(new TH2F("h_RZ","h_RZ",500.,0.,2.,200,0,0.5));
      HistReg.emplace_back(&h_RZ.store);
      h_RZfit_mom.emplace_back(new TH1F("h_RZfit_mom","h_RZfit_mom",500,-1.,1.));
      HistReg.emplace_back(&h_RZfit_mom.store);
      h_RZfit_Chi2.emplace_back(new TH2F("h_RZfit_Chi2","h_RZfit_Chi2",10,0,10,500,0,50.));
      HistReg.emplace_back(&h_RZfit_Chi2.store);
      h_XYfit_miniF.emplace_back(new TH2F("h_XYfit_miniF","h_XYfit_miniF",15,0,15,1000,-5.,5.));
      HistReg.emplace_back(&h_XYfit_miniF.store);
      h_MDC_Z_residu.emplace_back(new TH2F("h_MDC_Z_residu","h_MDC_Z_residu",1000.,-10.,10.,40,0,40));
      HistReg.emplace_back(&h_MDC_Z_residu.store);
      h_MDC_R_residu.emplace_back(new TH2F("h_MDC_R_residu","h_MDC_R_residu",1000.,-10.,10.,40,0,40));
      HistReg.emplace_back(&h_MDC_R_residu.store);
      h_MDC_Z_pull.emplace_back(new TH2F("h_MDC_Z_pull","h_MDC_Z_pull",1000.,-5.,5.,40,0,40));
      HistReg.emplace_back(&h_MDC_Z_pull.store);
      h_MDC_R_pull.emplace_back(new TH2F("h_MDC_R_pull","h_MDC_R_pull",1000.,-5.,5.,40,0,40));
      HistReg.emplace_back(&h_MDC_R_pull.store);

      HistRegisteredByDir.insert(std::make_pair("Finder", std::make_tuple(HistReg,0)));

      geoSolenoid.resize(17, nullptr);
    }

  if(EnableState[RIEMANN])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      h_RiemannChi2.emplace_back(new TH2F("h_RiemannChi2","h_RiemannChi2",5,0,5,120,-10,20));
      HistReg.emplace_back(&h_RiemannChi2.store);

      h_RiemannResidus.emplace_back(new TH2F("h_RiemannResidus","h_RiemannResidus",100,0,100,100,-10,10));
      HistReg.emplace_back(&h_RiemannResidus.store);

      h_PerfFinder.emplace_back(new TH2F("h_PerfFinder","h_PerfFinder",20,0,20,20,0,20));
      HistReg.emplace_back(&h_PerfFinder.store);

      h_PerfFinderLevenshtein.emplace_back(new TH2F("h_PerfFinderLevenshtein","h_PerfFinderLevenshtein",20,0,20,20,0,20));
      HistReg.emplace_back(&h_PerfFinderLevenshtein.store);

      HistRegisteredByDir.insert(std::make_pair("RiemannFinder", std::make_tuple(HistReg,0)));
    }

  if(EnableState[SIMU])
    {
      std::vector<std::vector<TH1*>*> HistReg;
      std::vector<std::vector<TH1*>*> HistRegEff;


      h_EfficiencyFiberHit.emplace_back(new TH2F("h_EfficiencyFiberHit","h_EfficiencyFiberHit",7,0,7,2,0,2));
      HistRegEff.emplace_back(&h_EfficiencyFiberHit.store);
      h_EfficiencySingleFiberHit.emplace_back(new TH2F("h_EfficiencySingleFiberHit","h_EfficiencySingleFiberHit",7,0,7,2,0,2));
      HistRegEff.emplace_back(&h_EfficiencySingleFiberHit.store);

      TString nameTempFib[] = {"UFT1", "UFT2", "UFT3", "MFT1", "MFT2", "DFT1", "DFT2"};
      std::vector<int> res_binfactor       = {1, 1, 1, 1, 1, 4, 4};
      std::vector<int> num_binfactor       = {1, 1, 3, 1, 1, 1, 1};
      std::vector<int> dvalue_binfactor    = {4, 4, 5, 1, 1, 4, 4};
      std::vector<int> dvalue_edgefactor   = {1, 1, 5, 1, 1, 1, 1};
      std::vector<int> mult_binfactor      = {1, 1, 2, 1, 1, 1, 1};
      std::vector<int> efftheta_edgefactor = {1, 1,16, 6, 6, 1, 1};
      std::vector<int> efftheta_binfactor  = {1, 1, 4, 4, 2, 1, 1};

      for(size_t i = 2; i < 7; ++i)
        {
          h_ResidualFiberHitX[i].emplace_back(new TH1F("h_ResidualFiberHitX_"+nameTempFib[i],"h_ResidualFiberHitX_"+nameTempFib[i], 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitX[i].store);
          h_ResidualFiberHitY[i].emplace_back(new TH1F("h_ResidualFiberHitY_"+nameTempFib[i],"h_ResidualFiberHitY_"+nameTempFib[i], 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitY[i].store);
          h_ResidualFiberHitR[i].emplace_back(new TH1F("h_ResidualFiberHitR_"+nameTempFib[i],"h_ResidualFiberHitR_"+nameTempFib[i], 100*res_binfactor[i], 0.,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitR[i].store);
          h_ResidualFiberHitXY[i].emplace_back(new TH2F("h_ResidualFiberHitXY_"+nameTempFib[i],"h_ResidualFiberHitXY_"+nameTempFib[i], 200*res_binfactor[i],-0.3,0.3, 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitXY[i].store);
          h_ResidualFiberHitX_Theta[i].emplace_back(new TH2F("h_ResidualFiberHitX_Theta_"+nameTempFib[i],"h_ResidualFiberHitX_Theta_"+nameTempFib[i],80,0.,40., 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitX_Theta[i].store);
          h_ResidualFiberHitY_Theta[i].emplace_back(new TH2F("h_ResidualFiberHitY_Theta_"+nameTempFib[i],"h_ResidualFiberHitY_Theta_"+nameTempFib[i],80,0.,40., 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitY_Theta[i].store);
          h_ResidualFiberHitX_HitX[i].emplace_back(new TH2F("h_ResidualFiberHitX_HitX_"+nameTempFib[i],"h_ResidualFiberHitX_HitX_"+nameTempFib[i],40,-10,10, 200,-0.3,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitX_HitX[i].store);
          h_ResidualFiberHitY_HitY[i].emplace_back(new TH2F("h_ResidualFiberHitY_HitY_"+nameTempFib[i],"h_ResidualFiberHitY_HitY_"+nameTempFib[i],40,-10,10, 200,-0.3,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitY_HitY[i].store);
          h_ResidualFiberHitR_Theta[i].emplace_back(new TH2F("h_ResidualFiberHitR_Theta_"+nameTempFib[i],"h_ResidualFiberHitR_Theta_"+nameTempFib[i],80,0.,40., 100*res_binfactor[i], 0.,0.3));
          HistReg.emplace_back(&h_ResidualFiberHitR_Theta[i].store);
          h_ResidualSingleFiberHitX[i].emplace_back(new TH1F("h_ResidualSingleFiberHitX_"+nameTempFib[i],"h_ResidualSingleFiberHitX_"+nameTempFib[i], 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualSingleFiberHitX[i].store);
          h_ResidualSingleFiberHitY[i].emplace_back(new TH1F("h_ResidualSingleFiberHitY_"+nameTempFib[i],"h_ResidualSingleFiberHitY_"+nameTempFib[i], 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualSingleFiberHitY[i].store);
          h_ResidualSingleFiberHitR[i].emplace_back(new TH1F("h_ResidualSingleFiberHitR_"+nameTempFib[i],"h_ResidualSingleFiberHitR_"+nameTempFib[i], 100*res_binfactor[i], 0.,0.3));
          HistReg.emplace_back(&h_ResidualSingleFiberHitR[i].store);
          h_ResidualSingleFiberHitX_Theta[i].emplace_back(new TH2F("h_ResidualSingleFiberHitX_Theta_"+nameTempFib[i],"h_ResidualSingleFiberHitX_Theta_"+nameTempFib[i],80,0.,40., 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualSingleFiberHitX_Theta[i].store);
          h_ResidualSingleFiberHitY_Theta[i].emplace_back(new TH2F("h_ResidualSingleFiberHitY_Theta_"+nameTempFib[i],"h_ResidualSingleFiberHitY_Theta_"+nameTempFib[i],80,0.,40., 200*res_binfactor[i],-0.3,0.3));
          HistReg.emplace_back(&h_ResidualSingleFiberHitY_Theta[i].store);
          h_ResidualSingleFiberHitR_Theta[i].emplace_back(new TH2F("h_ResidualSingleFiberHitR_Theta_"+nameTempFib[i],"h_ResidualSingleFiberHitR_Theta_"+nameTempFib[i],80,0.,40., 100*res_binfactor[i], 0.,0.3));
          HistReg.emplace_back(&h_ResidualSingleFiberHitR_Theta[i].store);
          h_EfficiencyFiberHit_Theta[i].emplace_back(new TH2F("h_EfficiencyFiberHit_Theta_"+nameTempFib[i],"h_EfficiencyFiberHit_Theta_"+nameTempFib[i],10*efftheta_binfactor[i],0,5*efftheta_edgefactor[i],2,0,2));
          HistRegEff.emplace_back(&h_EfficiencyFiberHit_Theta[i].store);
          h_EfficiencyFiberHit_dvalue[i].emplace_back(new TH2F("h_EfficiencyFiberHit_dvalue_"+nameTempFib[i],"h_EfficiencyFiberHit_dvalue_"+nameTempFib[i],100,-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i],2,0,2));
          HistRegEff.emplace_back(&h_EfficiencyFiberHit_dvalue[i].store);
          h_EfficiencyFiberHit_mult[i].emplace_back(new TH2F("h_EfficiencyFiberHit_mult_"+nameTempFib[i],"h_EfficiencyFiberHit_mult_"+nameTempFib[i],8*mult_binfactor[i],0,8*mult_binfactor[i],2,0,2));
          HistRegEff.emplace_back(&h_EfficiencyFiberHit_mult[i].store);
          h_EfficiencySingleFiberHit_Theta[i].emplace_back(new TH2F("h_EfficiencySingleFiberHit_Theta_"+nameTempFib[i],"h_EfficiencySingleFiberHit_Theta_"+nameTempFib[i],10*efftheta_binfactor[i],0,5*efftheta_edgefactor[i],2,0,2));
          HistRegEff.emplace_back(&h_EfficiencySingleFiberHit_Theta[i].store);
          h_EfficiencySingleFiberHit_dvalue[i].emplace_back(new TH2F("h_EfficiencySingleFiberHit_dvalue_"+nameTempFib[i],"h_EfficiencySingleFiberHit_dvalue_"+nameTempFib[i],100,-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i],2,0,2));
          HistRegEff.emplace_back(&h_EfficiencySingleFiberHit_dvalue[i].store);
          h_NumFiberHit_GoodReco[i].emplace_back(new TH2F("h_NumFiberHit_GoodReco_"+nameTempFib[i],"h_NumFiberHit_GoodReco_"+nameTempFib[i],7*num_binfactor[i],0,7*num_binfactor[i],7*num_binfactor[i],0,7*num_binfactor[i]));
          HistReg.emplace_back(&h_NumFiberHit_GoodReco[i].store);
          h_NumFiberHit_Ghost[i].emplace_back(new TH2F("h_NumFiberHit_Ghost_"+nameTempFib[i],"h_NumFiberHit_Ghost_"+nameTempFib[i],7*num_binfactor[i],0,7*num_binfactor[i],14*num_binfactor[i],0,14*num_binfactor[i]));
          HistReg.emplace_back(&h_NumFiberHit_Ghost[i].store);
          h_FiberHit_dvalue[i].emplace_back(new TH1F("h_FiberHit_dvalue_"+nameTempFib[i],"h_FiberHit_dvalue_"+nameTempFib[i], 100*dvalue_binfactor[i],-4.,4.));
          HistReg.emplace_back(&h_FiberHit_dvalue[i].store);
          h_FiberHitSingle_dvalue[i].emplace_back(new TH1F("h_FiberHitSingle_dvalue_"+nameTempFib[i],"h_FiberHitSingle_dvalue_"+nameTempFib[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitSingle_dvalue[i].store);
          h_FiberHitReal_dvalue[i].emplace_back(new TH1F("h_FiberHitReal_dvalue_"+nameTempFib[i],"h_FiberHitReal_dvalue_"+nameTempFib[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue[i].store);
          h_FiberHitReal_dvalue_Theta[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_Theta_"+nameTempFib[i],"h_FiberHitReal_dvalue_Theta_"+nameTempFib[i], 80,0.,40., 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_Theta[i].store);
          h_FiberHitReal_dvalue_Phi[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_Phi_"+nameTempFib[i],"h_FiberHitReal_dvalue_Phi_"+nameTempFib[i], 180,-180.,180., 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_Phi[i].store);
          h_FiberHitReal_dvalue_Theta03_Phi[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_Theta03_Phi_"+nameTempFib[i],"h_FiberHitReal_dvalue_Theta03_Phi_"+nameTempFib[i], 180,-180.,180., 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_Theta03_Phi[i].store);
          h_FiberHitReal_dvalue_Theta310_Phi[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_Theta310_Phi_"+nameTempFib[i],"h_FiberHitReal_dvalue_Theta310_Phi_"+nameTempFib[i], 180,-180.,180., 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_Theta310_Phi[i].store);
          h_FiberHitReal_dvalue_Theta1020_Phi[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_Theta1020_Phi_"+nameTempFib[i],"h_FiberHitReal_dvalue_Theta1020_Phi_"+nameTempFib[i], 180,-180.,180., 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_Theta1020_Phi[i].store);
          h_FiberHitReal_dvalue_HitX[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_HitX_"+nameTempFib[i],"h_FiberHitReal_dvalue_HitX_"+nameTempFib[i], 200,-10.,10., 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_HitX[i].store);
          h_FiberHitReal_dvalue_HitY[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_HitY_"+nameTempFib[i],"h_FiberHitReal_dvalue_HitY_"+nameTempFib[i], 200,-10.,10., 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_HitY[i].store);
          h_FiberHitReal_dvalue_PosX[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_PosX_"+nameTempFib[i],"h_FiberHitReal_dvalue_PosX_"+nameTempFib[i], 200,-10.,10., 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_PosX[i].store);
          h_FiberHitReal_dvalue_tanThetacosPhi[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_tanThetacosPhi_"+nameTempFib[i],"h_FiberHitReal_dvalue_tanThetacosPhi_"+nameTempFib[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_tanThetacosPhi[i].store);
          h_FiberHitReal_dvalue_dfunction[i].emplace_back(new TH2F("h_FiberHitReal_dvalue_dfunction_"+nameTempFib[i],"h_FiberHitReal_dvalue_dfunction_"+nameTempFib[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHitReal_dvalue_dfunction[i].store);
          h_FiberHit_Residualdvalue[i].emplace_back(new TH1F("h_FiberHit_Residualdvalue_"+nameTempFib[i],"h_FiberHit_Residualdvalue_"+nameTempFib[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHit_Residualdvalue[i].store);
          h_FiberHit_Residualdvalue_Realdvalue[i].emplace_back(new TH2F("h_FiberHit_Residualdvalue_Realdvalue_"+nameTempFib[i],"h_FiberHit_Residualdvalue_Realdvalue_"+nameTempFib[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i], 100*dvalue_binfactor[i],-1.*dvalue_edgefactor[i],1.*dvalue_edgefactor[i]));
          HistReg.emplace_back(&h_FiberHit_Residualdvalue_Realdvalue[i].store);
        }

      h_ResidualFiberX.emplace_back(new TH2F("h_ResidualFiberX","h_ResidualFiberX",20,0,20,1000,-1,1));
      HistReg.emplace_back(&h_ResidualFiberX.store);

      h_ResidualFiberY.emplace_back(new TH2F("h_ResidualFiberY","h_ResidualFiberY",20,0,20,1000,-1,1));
      HistReg.emplace_back(&h_ResidualFiberY.store);

      TString nameTemp[] = {"FiberD3_u","FiberD3_v","MiniFiberD1_x","MiniFiberD1_u","MiniFiberD1_v",
			    "MiniFiberD2_x","MiniFiberD2_u","MiniFiberD2_v","PSFE","PSCE","PSBE"};
      TString nameTemp1[] = {"_qN","_qP"};

      for(int i=0;i<11;++i)
	for(int j=0;j<2;++j)
	  {
	    h_ResidualFiberX_Angle[i][j].emplace_back(new TH2F("h_ResidualFiberX_Angle_"+nameTemp[i]+nameTemp1[j],"h_ResidualFiberX_Angle_"+nameTemp[i]+nameTemp1[j],180,0.,90, 200,-1,1));
	    HistReg.emplace_back(&h_ResidualFiberX_Angle[i][j].store);

	    h_ResidualFiberY_Angle[i][j].emplace_back(new TH2F("h_ResidualFiberY_Angle_"+nameTemp[i]+nameTemp1[j],"h_ResidualFiberY_Angle_"+nameTemp[i]+nameTemp1[j],180,0.,90,200,-1,1));
	    HistReg.emplace_back(&h_ResidualFiberY_Angle[i][j].store);
	  }

      h_ResidualFiberDzDphi.emplace_back(new TH2F("h_ResidualDzDphi","h_ResidualDzDphi",300,0,300,600,-30,30));
      HistReg.emplace_back(&h_ResidualFiberDzDphi.store);

      h_ResidualFiberDzDtheta.emplace_back(new TH2F("h_ResidualDzDtheta","h_ResidualDzDtheta",300,0,300,500,-10,10));
      HistReg.emplace_back(&h_ResidualFiberDzDtheta.store);

      HistRegisteredByDir.insert(std::make_pair("Simu", std::make_tuple(HistReg,0)));
      HistRegisteredByDir.insert(std::make_pair("Simu_Eff", std::make_tuple(HistRegEff,2)));
    }


  if(EnableState[BUILDER])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      h_Builderstats.emplace_back(new TH1I("Builderstats", "Builderstats", 10, 0, 10));
      HistReg.emplace_back(&h_Builderstats.store);

      int range10[7] = {256, 256, 384, 256, 256, 256, 256};
      int num_f13[7] = {256, 256, 384, 256, 256, 256, 256};
      for(int i=0; i<7; ++i)
        {
          for(int j=0; j<3; ++j)
            {
              h10[i][j].emplace_back(new TH2D(Form("h10[%d][%d]",i,j),Form("h10[%d][%d]",i,j),range10[i],-0.5,range10[i]-0.5, 300, -1700, -1400));
              HistReg.emplace_back(&h10[i][j].store);
              h11[i][j].emplace_back(new TH1D(Form("h11[%d][%d]",i,j),Form("h11[%d][%d]",i,j),200,-100,100));
              HistReg.emplace_back(&h11[i][j].store);
              h12[i][j].emplace_back(new TH1D(Form("h12[%d][%d]",i,j),Form("h12[%d][%d]",i,j),200,-100,100));
              HistReg.emplace_back(&h12[i][j].store);
              h13[i][j].emplace_back(new TH1D(Form("h13[%d][%d]",i,j),Form("h13[%d][%d]",i,j),10,-0.5,9.5));
              HistReg.emplace_back(&h13[i][j].store);
              h14[i][j].emplace_back(new TH1D(Form("h14[%d][%d]",i,j),Form("h14[%d][%d]",i,j),10,-0.5,9.5));
              HistReg.emplace_back(&h14[i][j].store);
              h15[i][j].emplace_back(new TH1D(Form("h15[%d][%d]",i,j),Form("h15[%d][%d]",i,j),8,-0.5,7.5));
              HistReg.emplace_back(&h15[i][j].store);
              h75[i][j].emplace_back(new TH2D(Form("h75[%d][%d]",i,j),Form("h75[%d][%d]",i,j),400,0,400, 300, -150, 150));
              HistReg.emplace_back(&h75[i][j].store);
              hfiber_13_0[i][j].emplace_back(new TH1D(Form("hfiber_13_0[%d][%d]",i,j),Form("Fiber %d-%d Leading Time" ,i,j), 300, -1700-0.5, -1400-0.5));
              HistReg.emplace_back(&hfiber_13_0[i][j].store);
              hfiber_13_1[i][j].emplace_back(new TH1D(Form("hfiber_13_1[%d][%d]",i,j),Form("Fiber %d-%d Time"         ,i,j), 200, -20, 20));
              HistReg.emplace_back(&hfiber_13_1[i][j].store);
              hfiber_13_2[i][j].emplace_back(new TH1D(Form("hfiber_13_2[%d][%d]",i,j),Form("Fiber %d-%d Time (w/o T0)",i,j), 200, -20, 20));
              HistReg.emplace_back(&hfiber_13_2[i][j].store);
              hfiber_13_3[i][j].emplace_back(new TH2D(Form("hfiber_13_3[%d][%d]",i,j),Form("Fiber %d-%d LeadingTime : Fiber" ,i,j), num_f13[i], -0.5, num_f13[i]-0.5, 300, -1700, -1400));
              HistReg.emplace_back(&hfiber_13_3[i][j].store);
              hfiber_13_4[i][j].emplace_back(new TH2D(Form("hfiber_13_4[%d][%d]",i,j),Form("Fiber %d-%d Time : Fiber"        ,i,j), num_f13[i], -0.5, num_f13[i]-0.5, 200, -30  , 30 ));
              HistReg.emplace_back(&hfiber_13_4[i][j].store);
            }
          h16[i].emplace_back(new TH2D(Form("h16[%d]",i),Form("h16[%d]",i),200,-100,100, 200, -100, 100));
          HistReg.emplace_back(&h16[i].store);
          h17[i].emplace_back(new TH1D(Form("h17[%d]",i),Form("h17[%d]",i),10 ,-0.5 ,9.5));
          HistReg.emplace_back(&h17[i].store);
          h17_2[i].emplace_back(new TH1D(Form("h17_2[%d]",i),Form("h17_2[%d]",i),200,-30  ,30 ));
          HistReg.emplace_back(&h17_2[i].store);
        }
      h18_3_1.emplace_back(new TH1D("h18_3_1","DFT12 ntrack",10, -0.5, 9.5));
      HistReg.emplace_back(&h18_3_1.store);
      h18_3_2.emplace_back(new TH1D("h18_3_2","DFT12 chi2"  ,200,   0, 100));
      HistReg.emplace_back(&h18_3_2.store);
      h18_3_3.emplace_back(new TH2D("h18_3_3","DFT12 XA"    ,100, -30, 30, 100, -30, 30));
      HistReg.emplace_back(&h18_3_3.store);
      h18_3_4.emplace_back(new TH2D("h18_3_4","DFT12 YB"    ,100, -30, 30, 100, -30, 30));
      HistReg.emplace_back(&h18_3_4.store);
      h18_3_5.emplace_back(new TH1D("h18_3_5","DFT12 X"     ,200, -30, 30));
      HistReg.emplace_back(&h18_3_5.store);
      h18_3_6.emplace_back(new TH1D("h18_3_6","DFT12 Y"     ,200, -30, 30));
      HistReg.emplace_back(&h18_3_6.store);
      h18_3_7.emplace_back(new TH1D("h18_3_7","DFT12 A"     ,200, -30, 30));
      HistReg.emplace_back(&h18_3_7.store);
      h18_3_8.emplace_back(new TH1D("h18_3_8","DFT12 B"     ,200, -30, 30));
      HistReg.emplace_back(&h18_3_8.store);
      int num51[3] = {384, 256, 256};
      for(int i=0; i<3 ; ++i)
        {
          for(int j=0; j<3; ++j)
            {
              for(int k=0; k<2; ++k)
                {
                  h51[i][j][k].emplace_back(new TH2D(Form("h51[%d][%d][%d]",i,j,k),Form("Fiber %d-%d %d",i,j,k) ,num51[i], -0.5, num51[i]-0.5, 100, 0, 100 ));
                  HistReg.emplace_back(&h51[i][j][k].store);
                }
            }
        }

      //MFT12
      hfiber_1_1.emplace_back(new TH1D("hfiber_1_1","MFT12 combi"        ,101, -0.5, 100.5));
      HistReg.emplace_back(&hfiber_1_1.store);
      hfiber_1_2.emplace_back(new TH1D("hfiber_1_2","MFT12 combi [k]"    ,500,  0  , 100  ));
      HistReg.emplace_back(&hfiber_1_2.store);
      hfiber_1_3.emplace_back(new TH1D("hfiber_1_3","MFT12 combi [k]"    ,500,  0  , 1000 ));
      HistReg.emplace_back(&hfiber_1_3.store);
      hfiber_1_4.emplace_back(new TH1D("hfiber_1_4","MFT12 combi [M]"    ,500,  0  , 10   ));
      HistReg.emplace_back(&hfiber_1_4.store);
      hfiber_1_5.emplace_back(new TH1D("hfiber_1_5","MFT12 combi all [k]",500,  0  , 100  ));
      HistReg.emplace_back(&hfiber_1_5.store);
      hfiber_1_6.emplace_back(new TH1D("hfiber_1_6","MFT12 combi all [k]",500,  0  , 1000 ));
      HistReg.emplace_back(&hfiber_1_6.store);
      hfiber_1_7.emplace_back(new TH1D("hfiber_1_7","MFT12 combi all [M]",500,  0  , 10   ));
      HistReg.emplace_back(&hfiber_1_7.store);
      hfiber_1_9.emplace_back(new TH1D("hfiber_1_9","MFT12 pair combi [k]",500, 0  , 100));
      HistReg.emplace_back(&hfiber_1_9.store);
      hfiber_2_1_1.emplace_back(new TH1D("hfiber_2_1_1","MFT12 combi pre"    ,101, -0.5, 100.5));
      HistReg.emplace_back(&hfiber_2_1_1.store);
      hfiber_2_1_2.emplace_back(new TH1D("hfiber_2_1_2","MFT12 combi pre [k]",500,  0  , 100  ));
      HistReg.emplace_back(&hfiber_2_1_2.store);
      hfiber_2_2_1.emplace_back(new TH1D("hfiber_2_2_1","MFT12 combi"        ,101, -0.5, 100.5));
      HistReg.emplace_back(&hfiber_2_2_1.store);
      hfiber_2_2_2.emplace_back(new TH1D("hfiber_2_2_2","MFT12 combi [k]"    ,500,  0  , 100  ));
      HistReg.emplace_back(&hfiber_2_2_2.store);
      hfiber_2_3.emplace_back(new TH1D("hfiber_2_3"  ,"MFT12 combi buf1"   ,101, -0.5, 100.5));
      HistReg.emplace_back(&hfiber_2_3.store);
      hfiber_3_0.emplace_back(new TH1D("hfiber_3_0"  ,"MFT12 ntrack"           ,11, -0.5, 10.5));
      HistReg.emplace_back(&hfiber_3_0.store);
      hfiber_3_0_2.emplace_back(new TH1D("hfiber_3_0_2","MFT12 ntrack (XUV)"     ,11, -0.5, 10.5));
      HistReg.emplace_back(&hfiber_3_0_2.store);
      hfiber_6_1.emplace_back(new TH1D("hfiber_6_1","MFT PSB diff phi (combi)" , 200, -0.7, 0.7));
      HistReg.emplace_back(&hfiber_6_1.store);
      hfiber_6_2.emplace_back(new TH1D("hfiber_6_2","MFT PSB diff Z (combi)  " , 200, -300, 300));
      HistReg.emplace_back(&hfiber_6_2.store);
      hfiber_6_3.emplace_back(new TH2D("hfiber_6_3","MFT PSB phi (combi)"      , 100, -180, 180, 100, -180, 180));
      HistReg.emplace_back(&hfiber_6_3.store);
      hfiber_6_4.emplace_back(new TH2D("hfiber_6_4","MFT PSB Z (combi)  "      , 100, -500, 500, 100, -500, 500));
      HistReg.emplace_back(&hfiber_6_4.store);
      hfiber_12_1_1.emplace_back(new TH2D("hfiber_12_1_1","MFT X1 X2 (all)" , 100, -100, 100, 100, -100 ,100));
      HistReg.emplace_back(&hfiber_12_1_1.store);
      hfiber_12_2_1.emplace_back(new TH2D("hfiber_12_2_1","MFT U1 U2 (all)" , 100, -100, 100, 100, -100 ,100));
      HistReg.emplace_back(&hfiber_12_2_1.store);
      hfiber_12_3_1.emplace_back(new TH2D("hfiber_12_3_1","MFT V1 V2 (all)" , 100, -100, 100, 100, -100 ,100));
      HistReg.emplace_back(&hfiber_12_3_1.store);
      hfiber_12_1_2.emplace_back(new TH2D("hfiber_12_1_2","MFT X1 X2 (pair)", 100, -100, 100, 100, -100 ,100));
      HistReg.emplace_back(&hfiber_12_1_2.store);
      hfiber_12_2_2.emplace_back(new TH2D("hfiber_12_2_2","MFT U1 U2 (pair)", 100, -100, 100, 100, -100 ,100));
      HistReg.emplace_back(&hfiber_12_2_2.store);
      hfiber_12_3_2.emplace_back(new TH2D("hfiber_12_3_2","MFT V1 V2 (pair)", 100, -100, 100, 100, -100 ,100));
      HistReg.emplace_back(&hfiber_12_3_2.store);

      //DFT12
      hfiber_4_1.emplace_back(new TH1D("hfiber_4_1","DFT12 ntrack", 10, -0.5, 9.5));
      HistReg.emplace_back(&hfiber_4_1.store);
      hfiber_4_2_1.emplace_back(new TH2D("hfiber_4_2_1", "DFT12 XY single" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_2_1.store);
      hfiber_4_3_1.emplace_back(new TH2D("hfiber_4_3_1", "DFT12 XA single" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_3_1.store);
      hfiber_4_4_1.emplace_back(new TH2D("hfiber_4_4_1", "DFT12 YB single" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_4_1.store);
      hfiber_4_5_1.emplace_back(new TH1D("hfiber_4_5_1", "DFT12 TOT single", 150, 0, 150));
      HistReg.emplace_back(&hfiber_4_5_1.store);
      hfiber_4_2_2.emplace_back(new TH2D("hfiber_4_2_2", "DFT12 XY multi" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_2_2.store);
      hfiber_4_3_2.emplace_back(new TH2D("hfiber_4_3_2", "DFT12 XA multi" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_3_2.store);
      hfiber_4_4_2.emplace_back(new TH2D("hfiber_4_4_2", "DFT12 YB multi" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_4_2.store);
      hfiber_4_5_2.emplace_back(new TH1D("hfiber_4_5_2", "DFT12 TOT multi", 150, 0, 150));
      HistReg.emplace_back(&hfiber_4_5_2.store);
      hfiber_4_1_3.emplace_back(new TH1D("hfiber_4_1_3","DFT12 ntrack", 10, -0.5, 9.5));
      HistReg.emplace_back(&hfiber_4_1_3.store);
      hfiber_4_2_3.emplace_back(new TH2D("hfiber_4_2_3", "DFT12 XY all" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_2_3.store);
      hfiber_4_3_3.emplace_back(new TH2D("hfiber_4_3_3", "DFT12 XA all" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_3_3.store);
      hfiber_4_4_3.emplace_back(new TH2D("hfiber_4_4_3", "DFT12 YB all" , 100, -100, 100, 100, -100, 100));
      HistReg.emplace_back(&hfiber_4_4_3.store);
      hfiber_4_5_3.emplace_back(new TH1D("hfiber_4_5_3", "DFT12 TOT all", 150, 0, 150));
      HistReg.emplace_back(&hfiber_4_5_3.store);
      hfiber_5_1.emplace_back(new TH1D("hfiber_5_1","DFT12 combi"        ,101, -0.5, 100.5));
      HistReg.emplace_back(&hfiber_5_1.store);
      hfiber_5_2.emplace_back(new TH1D("hfiber_5_2","DFT12 combi [k]"    ,500,  0  , 100  ));
      HistReg.emplace_back(&hfiber_5_2.store);
      hfiber_5_3.emplace_back(new TH1D("hfiber_5_3","DFT12 combi [k]"    ,500,  0  , 1000 ));
      HistReg.emplace_back(&hfiber_5_3.store);
      hfiber_5_4.emplace_back(new TH1D("hfiber_5_4","DFT12 combi [M]"    ,500,  0  , 10   ));
      HistReg.emplace_back(&hfiber_5_4.store);
      hfiber_5_5.emplace_back(new TH1D("hfiber_5_5","DFT12 combi all [k]",500,  0  , 100  ));
      HistReg.emplace_back(&hfiber_5_5.store);
      hfiber_5_6.emplace_back(new TH1D("hfiber_5_6","DFT12 combi all [k]",500,  0  , 1000 ));
      HistReg.emplace_back(&hfiber_5_6.store);
      hfiber_5_7.emplace_back(new TH1D("hfiber_5_7","DFT12 combi all [M]",500,  0  , 10   ));
      HistReg.emplace_back(&hfiber_5_7.store);

      hpsb_0_1.emplace_back(new TH1D("hpsb_0_1","tu cut" ,200, -21000, -19000));
      HistReg.emplace_back(&hpsb_0_1.store);
      hpsb_0_2.emplace_back(new TH1D("hpsb_0_2","nhit"   ,30, -0.5, 29.5));
      HistReg.emplace_back(&hpsb_0_2.store);
      hpsb_0_3.emplace_back(new TH1D("hpsb_0_3","phi"    ,200, -3.2*2, 3.2*2));
      HistReg.emplace_back(&hpsb_0_3.store);
      hpsb_0_4.emplace_back(new TH2D("hpsb_0_4","phi vs seg" ,46, -0.5,  45.5,  200, -3.2, 3.2));
      HistReg.emplace_back(&hpsb_0_4.store);
      hpsb_1_1.emplace_back(new TH1D("hpsb_1_1","pos z" ,200, -500, 500));
      HistReg.emplace_back(&hpsb_1_1.store);
      for(int i=0; i<46; ++i)
      {
        hpsb_2[i].emplace_back(new TH1D(Form("hpsb_2[%d]",i),Form("PSB Seg %d",i),200, -15, 15));
        HistReg.emplace_back(&hpsb_2[i].store);
        hpsb_3[i].emplace_back(new TH2D(Form("hpsb_3[%d]",i),Form("PSB Seg %d",i),200, -8, 8, 200, -400, 400));
        HistReg.emplace_back(&hpsb_3[i].store);
        hpsb_4[i].emplace_back(new TH1D(Form("hpsb_4[%d]",i),Form("PSB Seg %d",i),200, -500, 500));
        HistReg.emplace_back(&hpsb_4[i].store);
      }
      h76.emplace_back(new TH2D("h76","HitID Phi" ,50, 0, 50, 300, -3.5, 3.5));
      HistReg.emplace_back(&h76.store);

      hpsfe_0_1.emplace_back(new TH1D("hpsfe_0_1","t cut" ,200, -26500, -23500));
      HistReg.emplace_back(&hpsfe_0_1.store);
      hpsfe_0_2.emplace_back(new TH1D("hpsfe_0_2","nhit"   ,30, -0.5, 29.5));
      HistReg.emplace_back(&hpsfe_0_2.store);
      hpsfe_0_3.emplace_back(new TH1D("hpsfe_0_3","phi"    ,200, -3.2*2, 3.2*2));
      HistReg.emplace_back(&hpsfe_0_3.store);
      hpsfe_0_4.emplace_back(new TH2D("hpsfe_0_4","phi vs seg" ,46, -0.5,  45.5,  200, -3.2, 3.2));
      HistReg.emplace_back(&hpsfe_0_4.store);

      hpsbe_0_1.emplace_back(new TH1D("hpsbe_0_1","t cut" ,200, -26500, -23500));
      HistReg.emplace_back(&hpsbe_0_1.store);
      hpsbe_0_2.emplace_back(new TH1D("hpsbe_0_2","nhit"   ,30, -0.5, 29.5));
      HistReg.emplace_back(&hpsbe_0_2.store);
      hpsbe_0_3.emplace_back(new TH1D("hpsbe_0_3","phi"    ,200, -3.2*2, 3.2*2));
      HistReg.emplace_back(&hpsbe_0_3.store);
      hpsbe_0_4.emplace_back(new TH2D("hpsbe_0_4","phi vs seg" ,46, -0.5,  45.5,  200, -3.2, 3.2));
      HistReg.emplace_back(&hpsbe_0_4.store);
      hpsbe_1_0.emplace_back(new TH2D("hpsbe_1_0","t all" ,38, -0.5, 37.5, 200, -26500, -23500));
      HistReg.emplace_back(&hpsbe_1_0.store);

      ht0_0_1.emplace_back(new TH1D("ht0_0_1","tu cut" ,200, -26500, -25500));
      HistReg.emplace_back(&ht0_0_1.store);
      ht0_0_2.emplace_back(new TH1D("ht0_0_2","td cut" ,200, -26500, -25500));
      HistReg.emplace_back(&ht0_0_2.store);
      ht0_0_3.emplace_back(new TH1D("ht0_0_3","nhit"   ,30, -0.5, 29.5));
      HistReg.emplace_back(&ht0_0_3.store);
      ht0_0_4.emplace_back(new TH1D("ht0_0_4","seg"    ,28, -0.5, 27.5));
      HistReg.emplace_back(&ht0_0_4.store);
      for(int i=0; i<28; ++i)
      {
        ht0_1[i].emplace_back(new TH1D(Form("ht0_1[%d]",i),Form("T0 Seg %d",i),200, -15, 15));
        HistReg.emplace_back(&ht0_1[i].store);
      }

      hmdc_0_1.emplace_back(new TH1D("hmdc_0_1","mdc leading"     ,200, -2000, 1000));
      HistReg.emplace_back(&hmdc_0_1.store);
      hmdc_0_2.emplace_back(new TH2D("hmdc_0_2","mdc leading tot" ,200, -2000, 1000, 200, 0, 1000));
      HistReg.emplace_back(&hmdc_0_2.store);
      hmdc_0_3.emplace_back(new TH1D("hmdc_0_3","mdc r"           ,200, 0, 500));
      HistReg.emplace_back(&hmdc_0_3.store);
      hmdc_0_4.emplace_back(new TH1D("hmdc_0_4","mdc phi"         ,200, -10, 10));
      HistReg.emplace_back(&hmdc_0_4.store);
      hmdc_0_5.emplace_back(new TH1D("hmdc_0_5","mdc leading layer 16"     ,200, -2000, 1000));
      HistReg.emplace_back(&hmdc_0_5.store);
      hmdc_0_6.emplace_back(new TH2D("hmdc_0_6","mdc leading tot layer 16" ,200, -2000, 1000, 200, 0, 1000));
      HistReg.emplace_back(&hmdc_0_6.store);
      hmdc_0_9.emplace_back(new TH2D("hmdc_0_9","mdc wire vs phi layer 16" ,300, 0, 300, 360, -180, 180));
      HistReg.emplace_back(&hmdc_0_9.store);
      for(int i=0; i<17; ++i)
      {
        hmdc_1[i].emplace_back(new TH1D(Form("hmdc_1[%d]",i),Form("TL layer %d",i), 300,-2000,-500));
        HistReg.emplace_back(&hmdc_1[i].store);
        hmdc_2[i].emplace_back(new TH1D(Form("hmdc_2[%d]",i),Form("DT layer %d",i), 125,-50,200));
        HistReg.emplace_back(&hmdc_2[i].store);
        hmdc_3[i].emplace_back(new TH1D(Form("hmdc_3[%d]",i),Form("DL layer %d",i), 100,-1,5));
        HistReg.emplace_back(&hmdc_3[i].store);
        hmdc_3_2[i].emplace_back(new TH1D(Form("hmdc_3_2[%d]",i),Form("DL layer %d cut",i), 100,-1,5));
        HistReg.emplace_back(&hmdc_3_2[i].store);
        hmdc_3_3[i].emplace_back(new TH1D(Form("hmdc_3_3[%d]",i),Form("DL layer %d cut",i), 100,-1,5));
        HistReg.emplace_back(&hmdc_3_3[i].store);
      }

      hmwdc_1_1.emplace_back(new TH1D("hmwdc_1_1","MWDC X"   ,500, -200 ,200));
      HistReg.emplace_back(&hmwdc_1_1.store);
      hmwdc_1_2.emplace_back(new TH1D("hmwdc_1_2","MWDC Y"   ,500, -50  ,50 ));
      HistReg.emplace_back(&hmwdc_1_2.store);
      hmwdc_1_3.emplace_back(new TH1D("hmwdc_1_3","MWDC A"   ,500, -30  ,30));
      HistReg.emplace_back(&hmwdc_1_3.store);
      hmwdc_1_4.emplace_back(new TH1D("hmwdc_1_4","MWDC B"   ,500, -30  ,30 ));
      HistReg.emplace_back(&hmwdc_1_4.store);
      hmwdc_1_5.emplace_back(new TH1D("hmwdc_1_5","MWDC Chi2",500, 0    ,10));
      HistReg.emplace_back(&hmwdc_1_5.store);
      hmwdc_1_6.emplace_back(new TH1D("hmwdc_1_6","MWDC nt"  ,10 , 0    ,10));
      HistReg.emplace_back(&hmwdc_1_6.store);

      hs4sci_1_1.emplace_back(new TH1D("hs4sci_1_1","S4 SC31 dE"     ,500, 0 ,3000));
      HistReg.emplace_back(&hs4sci_1_1.store);
      hs4sci_1_2.emplace_back(new TH1D("hs4sci_1_2","S4 SC41 dE"     ,500, 0 ,3000));
      HistReg.emplace_back(&hs4sci_1_2.store);
      hs4sci_1_3.emplace_back(new TH1D("hs4sci_1_3","S4 SC42 High dE",500, 0 ,3000));
      HistReg.emplace_back(&hs4sci_1_3.store);
      hs4sci_1_4.emplace_back(new TH1D("hs4sci_1_4","S4 SC42 Low dE",500, 0 ,3000));
      HistReg.emplace_back(&hs4sci_1_4.store);

      hs4sci_2_1.emplace_back(new TH2D("hs4sci_2_1","S4 TOF31-41 : SC31 dE"     ,200, 30 ,100, 200, 0, 3000));
      HistReg.emplace_back(&hs4sci_2_1.store);
      hs4sci_2_2.emplace_back(new TH2D("hs4sci_2_2","S4 TOF31-41 : SC41 dE"     ,200, 30 ,100, 200, 0, 3000));
      HistReg.emplace_back(&hs4sci_2_2.store);
      hs4sci_2_3.emplace_back(new TH2D("hs4sci_2_3","S4 TOF31-41 : SC42 High dE",200, 30 ,100, 200, 0, 3000));
      HistReg.emplace_back(&hs4sci_2_3.store);
      hs4sci_2_4.emplace_back(new TH2D("hs4sci_2_4","S4 TOF31-41 : SC42 Low  dE",200, 30 ,100, 200, 0, 3000));
      HistReg.emplace_back(&hs4sci_2_4.store);

      htrig_0.emplace_back(new TH1D("htrig_0","Trig (pre)"     ,16, -0.5, 15.5));
      HistReg.emplace_back(&htrig_0.store);
      htrig_1.emplace_back(new TH2D("htrig_1","V775 (pre)"     ,32, -0.5, 31.5, 200, 0, 4000));
      HistReg.emplace_back(&htrig_1.store);
      htrig_2.emplace_back(new TH2D("htrig_2","V775 cut (pre)" ,32, -0.5, 31.5, 200, 0, 4000));
      HistReg.emplace_back(&htrig_2.store);
      htrig_3.emplace_back(new TH1D("htrig_3","Trig"     ,16, -0.5, 15.5));
      HistReg.emplace_back(&htrig_3.store);
      htrig_4.emplace_back(new TH2D("htrig_4","V775 cut" ,32, -0.5, 31.5, 200, 0, 4000));
      HistReg.emplace_back(&htrig_4.store);

      HistRegisteredByDir.insert(std::make_pair("Builder", std::make_tuple(HistReg,0)));
    }


  if(EnableState[FRAGMENT])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      hopt_1_1.emplace_back(new TH1D("hopt_1_1","Optics A2 reco" ,500, -100 ,100));
      HistReg.emplace_back(&hopt_1_1.store);
      hopt_1_2.emplace_back(new TH1D("hopt_1_2","Optics B2 reco" ,500, -100 ,100));
      HistReg.emplace_back(&hopt_1_2.store);
      hopt_1_3.emplace_back(new TH1D("hopt_1_3","Optics Mom He3" ,500, 7.0  ,8.0));
      HistReg.emplace_back(&hopt_1_3.store);
      hopt_1_4.emplace_back(new TH1D("hopt_1_4","Optics Num"     ,  5, -0.5 ,4.5));
      HistReg.emplace_back(&hopt_1_4.store);

      hopt_2_1.emplace_back(new TH2D("hopt_2_1","Optics A2 reco : A2" ,500, -100 ,100, 500, -100, 100));
      HistReg.emplace_back(&hopt_2_1.store);
      hopt_2_2.emplace_back(new TH2D("hopt_2_2","Optics B2 reco : B2" ,500, -100 ,100, 500, -100, 100));
      HistReg.emplace_back(&hopt_2_2.store);
      hopt_2_3.emplace_back(new TH1D("hopt_2_3","Optics A2 reco - A2" ,500, -30 ,30 ));
      HistReg.emplace_back(&hopt_2_3.store);
      hopt_2_4.emplace_back(new TH1D("hopt_2_4","Optics B2 reco - B2" ,500, -30 ,30 ));
      HistReg.emplace_back(&hopt_2_4.store);

      HistRegisteredByDir.insert(std::make_pair("FragmentFinder", std::make_tuple(HistReg,0)));
    }

  if(EnableState[WASA])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      h23_1.emplace_back(new TH2D("h23_1","Fiber PSB Phi"     ,200, -3.2, 3.2, 50 , -3.21, 3.34));
      HistReg.emplace_back(&h23_1.store);
      h23_2.emplace_back(new TH2D("h23_2","Fiber PSB Z"       ,200, -500, 500, 200, -500 , 500 ));
      HistReg.emplace_back(&h23_2.store);
      h24_1.emplace_back(new TH2D("h24_1","Fiber MDC Phi"     ,200, -3.2, 3.2, 200 , -3.2, 3.2));
      HistReg.emplace_back(&h24_1.store);
      h24_2.emplace_back(new TH2D("h24_2","Fiber MDC Phi cut" ,200, -3.2, 3.2, 200 , -3.2, 3.2));
      HistReg.emplace_back(&h24_2.store);
      h24_2_1.emplace_back(new TH1D("h24_2_1","Fiber MDC Phi"     ,200, -0.5, 0.5));
      HistReg.emplace_back(&h24_2_1.store);
      h24_2_2.emplace_back(new TH1D("h24_2_2","Fiber MDC Phi cut" ,200, -0.5, 0.5));
      HistReg.emplace_back(&h24_2_2.store);
      for(int i = 0; i < 17; ++i)
        {
          h24_3[i].emplace_back(new TH2D(Form("h24_3[%d]",i),Form("Fiber MDC Phi layer %d",i)     ,200, -3.2, 3.2, 200 , -3.2, 3.2));
          HistReg.emplace_back(&h24_3[i].store);
          h24_4[i].emplace_back(new TH2D(Form("h24_4[%d]",i),Form("Fiber MDC Phi layer %d cut",i) ,200, -3.2, 3.2, 200 , -3.2, 3.2));
          HistReg.emplace_back(&h24_4[i].store);
          h24_2_3[i].emplace_back(new TH1D(Form("h24_2_3[%d]",i),Form("Fiber MDC Phi layer %d",i)     ,200, -0.5, 0.5));
          HistReg.emplace_back(&h24_2_3[i].store);
          h24_2_4[i].emplace_back(new TH1D(Form("h24_2_4[%d]",i),Form("Fiber MDC Phi layer %d cut",i) ,200, -0.5, 0.5));
          HistReg.emplace_back(&h24_2_4[i].store);
          hmdc_2_2[i].emplace_back(new TH1D(Form("hmdc_2_2[%d]",i),Form("DT layer %d cut",i), 125,-50,200));
          HistReg.emplace_back(&hmdc_2_2[i].store);
          hmdc_3_2[i].emplace_back(new TH1D(Form("hmdc_3_2[%d]",i),Form("DL layer %d cut",i), 100,-1,5));
          HistReg.emplace_back(&hmdc_3_2[i].store);
        }

      HistRegisteredByDir.insert(std::make_pair("WasaFinder", std::make_tuple(HistReg,0)));
    }

if(EnableState[PRIMVTX])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      TString nameTempFib[] = {"UFT1", "UFT2", "UFT3", "MFT1", "MFT2", "DFT1", "DFT2"};

      h_NFiberMult.emplace_back(new TH2F("h_NFiberMult", "h_NFiberMult", 20, 0, 20, 15, 0, 15));
      HistReg.emplace_back(&h_NFiberMult.store);
      h_NCombsXUV_UFT12.emplace_back(new TH2F("h_NCombsXUV_UFT12", "h_NCombsXUV_UFT12", 10, 0, 10, 10, 0, 10));
      HistReg.emplace_back(&h_NCombsXUV_UFT12.store);

      for(size_t i = 0; i < 5; ++i)
        {
          h_NCombsXUV[i].emplace_back(new TH2F("h_NCombsXUV_"+nameTempFib[i], "h_NCombsXUV_"+nameTempFib[i], 50, 0, 50, 2, 0, 2));
          HistReg.emplace_back(&h_NCombsXUV[i].store);
        }
      for(size_t i = 2; i < 5; ++i)
        {
          h_CombsXUV_dvalue_theta[i].emplace_back(new TH2F("h_CombsXUV_dvalue_theta_"+nameTempFib[i], "h_CombsXUV_dvalue_theta_"+nameTempFib[i], 90, 0, 90, 500, -5, 5));
          HistReg.emplace_back(&h_CombsXUV_dvalue_theta[i].store);
          h_CombsXUV_dvalue_phi[i].emplace_back(new TH2F("h_CombsXUV_dvalue_phi_"+nameTempFib[i], "h_CombsXUV_dvalue_phi_"+nameTempFib[i], 180, -180, 180, 500, -5, 5));
          HistReg.emplace_back(&h_CombsXUV_dvalue_phi[i].store);
          h_CombsXUV_dvalue_phi_theta5[i].emplace_back(new TH2F("h_CombsXUV_dvalue_phi_theta5_"+nameTempFib[i], "h_CombsXUV_dvalue_phi_theta5_"+nameTempFib[i], 180, -180, 180, 500, -5, 5));
          HistReg.emplace_back(&h_CombsXUV_dvalue_phi_theta5[i].store);
          h_CombsXUV_dvalue_phi_theta10[i].emplace_back(new TH2F("h_CombsXUV_dvalue_phi_theta10_"+nameTempFib[i], "h_CombsXUV_dvalue_phi_theta10_"+nameTempFib[i], 180, -180, 180, 500, -5, 5));
          HistReg.emplace_back(&h_CombsXUV_dvalue_phi_theta10[i].store);
        }

      h_NHits_PrimaryTracks.emplace_back(new TH2F("h_NHits_PrimaryTracks", "h_NHits_PrimaryTracks", 10, 0, 10, 5, 0, 5));
      HistReg.emplace_back(&h_NHits_PrimaryTracks.store);

      h_nTrackCandidates.emplace_back(new TH2F("h_nTrackCandidates","h_nTrackCandidates", 100000, 0, 100000, 8, 0, 8));
      HistReg.emplace_back(&h_nTrackCandidates.store);
      h_DistanceBeamTracks.emplace_back(new TH2F("h_DistanceBeamTracks","h_DistanceBeamTracks", 400, 0, 20, 3, 0, 3));
      HistReg.emplace_back(&h_DistanceBeamTracks.store);
      h_PosZBeamTracks.emplace_back(new TH2F("h_PosZBeamTracks","h_PosZBeamTracks", 800, 100, 300, 3, 0, 3));
      HistReg.emplace_back(&h_PosZBeamTracks.store);
      h_thetaTracks.emplace_back(new TH2F("h_thetaTracks","h_thetaTracks", 200, 0, 50, 5, 0, 5));
      HistReg.emplace_back(&h_thetaTracks.store);
      h_chi2ndfTracks.emplace_back(new TH2F("h_chi2ndfTracks","h_chi2ndfTracks", 1000, 0, 1000, 5, 0, 5));
      HistReg.emplace_back(&h_chi2ndfTracks.store);

      h_fvalues.emplace_back(new TH1F("h_fvalues", "h_fvalues", 200, -0.5, 1.5));
      HistReg.emplace_back(&h_fvalues.store);

      h_InteractionPointPosX.emplace_back(new TH1F("h_InteractionPointPosX","h_InteractionPointPosX", 100, -2, 2));
      HistReg.emplace_back(&h_InteractionPointPosX.store);
      h_InteractionPointPosY.emplace_back(new TH1F("h_InteractionPointPosY","h_InteractionPointPosY", 100, -2, 2));
      HistReg.emplace_back(&h_InteractionPointPosY.store);
      h_InteractionPointPosZ.emplace_back(new TH1F("h_InteractionPointPosZ","h_InteractionPointPosZ", 100, 192, 200));
      HistReg.emplace_back(&h_InteractionPointPosZ.store);

      h_InteractionPointDistance_V_value.emplace_back(new TH2F("h_InteractionPointDistance_V_value","h_InteractionPointDistance_V_value", 100, 0, 5, 200, 0, 10));
      HistReg.emplace_back(&h_InteractionPointDistance_V_value.store);

      h_InteractionPointDistanceX.emplace_back(new TH1F("h_InteractionPointDistanceX","h_InteractionPointDistanceX", 600, -3, 3));
      HistReg.emplace_back(&h_InteractionPointDistanceX.store);
      h_InteractionPointDistanceY.emplace_back(new TH1F("h_InteractionPointDistanceY","h_InteractionPointDistanceY", 600, -3, 3));
      HistReg.emplace_back(&h_InteractionPointDistanceY.store);
      h_InteractionPointDistanceZ.emplace_back(new TH1F("h_InteractionPointDistanceZ","h_InteractionPointDistanceZ", 500, -5, 5));
      HistReg.emplace_back(&h_InteractionPointDistanceZ.store);

      h_InteractionPointDistanceX_pull.emplace_back(new TH1F("h_InteractionPointDistanceX_pull","h_InteractionPointDistanceX_pull", 2000, -5, 5));
      HistReg.emplace_back(&h_InteractionPointDistanceX_pull.store);
      h_InteractionPointDistanceY_pull.emplace_back(new TH1F("h_InteractionPointDistanceY_pull","h_InteractionPointDistanceY_pull", 2000, -5, 5));
      HistReg.emplace_back(&h_InteractionPointDistanceY_pull.store);
      h_InteractionPointDistanceZ_pull.emplace_back(new TH1F("h_InteractionPointDistanceZ_pull","h_InteractionPointDistanceZ_pull", 2000, -5, 5));
      HistReg.emplace_back(&h_InteractionPointDistanceZ_pull.store);

      h_CovarianceSigmaX.emplace_back(new TH1F("h_SqrtCovariance_SigmaX", "h_SqrtCovariance_SigmaX", 1000, 0, 0.1));
      HistReg.emplace_back(&h_CovarianceSigmaX.store);
      h_CovarianceSigmaY.emplace_back(new TH1F("h_SqrtCovariance_SigmaY", "h_SqrtCovariance_SigmaY", 1000, 0, 0.1));
      HistReg.emplace_back(&h_CovarianceSigmaY.store);
      h_CovarianceSigmaZ.emplace_back(new TH1F("h_SqrtCovariance_SigmaZ", "h_SqrtCovariance_SigmaZ", 2000, 0, 0.2));
      HistReg.emplace_back(&h_CovarianceSigmaZ.store);

      h_PrimStatus.emplace_back(new TH2F("PrimVtxStatus","PrimVtxStatus",20,0,20,20,0,20));
      HistReg.emplace_back(&h_PrimStatus.store);
      h_PrimVtxstats.emplace_back(new TH1F("PrimVtxstats", "PrimVtxstats", 10, 0, 10));
      HistReg.emplace_back(&h_PrimVtxstats.store);
       
      HistRegisteredByDir.insert(std::make_pair("PrimaryVtx", std::make_tuple(HistReg,0)));
    }

  if(EnableState[PRIMVTX_SI])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      h_HitMultiplicity_Si1.emplace_back(new TH1F("h_HitMultiplicity_Si1","h_HitMultiplicity_Si1", 20, 0, 20));
      HistReg.emplace_back(&h_HitMultiplicity_Si1.store);
      h_HitMultiplicityRecons_Si1.emplace_back(new TH1F("h_HitMultiplicityRecons_Si1","h_HitMultiplicityRecons_Si1", 40, 0, 40));
      HistReg.emplace_back(&h_HitMultiplicityRecons_Si1.store);
      h_HitMultiplicityDiff_Si1.emplace_back(new TH1F("h_HitMultiplicityDiff_Si1","h_HitMultiplicityDiff_Si1", 50, -5, 45));
      HistReg.emplace_back(&h_HitMultiplicityDiff_Si1.store);
      h_HitMultiplicityDiffNHits_Si1.emplace_back(new TH2F("h_HitMultiplicityDiffNHits_Si1","h_HitMultiplicityDiffNHits_Si1", 50, -5, 45, 20, 0, 20));
      HistReg.emplace_back(&h_HitMultiplicityDiffNHits_Si1.store);

      h_EnergyDiffStrips_Si1.emplace_back(new TH1F("h_EnergyDiffStrips_Si1", "h_EnergyDiffStrips_Si1", 1000, -1, 1));
      HistReg.emplace_back(&h_EnergyDiffStrips_Si1.store);
      h_nEventsGoodrecons_Si1.emplace_back(new TH1F("h_nEventsGoodrecons_Si1", "h_nEventsGoodrecons_Si1", 20, 0, 20));
      HistReg.emplace_back(&h_nEventsGoodrecons_Si1.store);
      h_nEventsGhost_Si1.emplace_back(new TH1F("h_nEventsGhost_Si1", "h_nEventsGhost_Si1", 40, 0, 40));
      HistReg.emplace_back(&h_nEventsGhost_Si1.store);
      h_nEventsGoodreconsGhost_Si1.emplace_back(new TH2F("h_nEventsGoodreconsGhost_Si1", "h_nEventsGoodreconsGhost_Si1", 20, 0, 20, 40, 0, 40));
      HistReg.emplace_back(&h_nEventsGoodreconsGhost_Si1.store);
      h_nEventsRealGoodrecons_Si1.emplace_back(new TH2F("h_nEventsRealGoodrecons_Si1", "h_nEventsRealGoodrecons_Si1", 20, 0, 20, 20, 0, 20));
      HistReg.emplace_back(&h_nEventsRealGoodrecons_Si1.store);
      h_nEventsRealRejectPad_Si1.emplace_back(new TH2F("h_nEventsRealRejectPad_Si1", "h_nEventsRealRejectPad_Si1", 20, 0, 20, 300, 0, 300));
      HistReg.emplace_back(&h_nEventsRealRejectPad_Si1.store);


      h_HitMultiplicity_Si2.emplace_back(new TH1F("h_HitMultiplicity_Si2","h_HitMultiplicity_Si2", 20, 0, 20));
      HistReg.emplace_back(&h_HitMultiplicity_Si2.store);
      h_HitMultiplicityRecons_Si2.emplace_back(new TH1F("h_HitMultiplicityRecons_Si2","h_HitMultiplicityRecons_Si2", 40, 0, 40));
      HistReg.emplace_back(&h_HitMultiplicityRecons_Si2.store);
      h_HitMultiplicityDiff_Si2.emplace_back(new TH1F("h_HitMultiplicityDiff_Si2","h_HitMultiplicityDiff_Si2", 50, -5, 45));
      HistReg.emplace_back(&h_HitMultiplicityDiff_Si2.store);
      h_HitMultiplicityDiffNHits_Si2.emplace_back(new TH2F("h_HitMultiplicityDiffNHits_Si2","h_HitMultiplicityDiffNHits_Si2", 50, -5, 45, 20, 0, 20));
      HistReg.emplace_back(&h_HitMultiplicityDiffNHits_Si2.store);
      h_EnergyDiffStrips_Si2.emplace_back(new TH1F("h_EnergyDiffStrips_Si2", "h_EnergyDiffStrips_Si2", 1000, -1, 1));
      HistReg.emplace_back(&h_EnergyDiffStrips_Si2.store);

      h_nEventsGoodrecons_Si2.emplace_back(new TH1F("h_nEventsGoodrecons_Si2","h_nEventsGoodrecons_Si2", 20, 0, 20));
      HistReg.emplace_back(&h_nEventsGoodrecons_Si2.store);
      h_nEventsGhost_Si2.emplace_back(new TH1F("h_nEventsGhost_Si2","h_nEventsGhost_Si2", 40, 0, 40));
      HistReg.emplace_back(&h_nEventsGhost_Si2.store);
      h_nEventsGoodreconsGhost_Si2.emplace_back(new TH2F("h_nEventsGoodreconsGhost_Si2","h_nEventsGoodreconsGhost_Si2", 20, 0, 20, 40, 0, 40));
      HistReg.emplace_back(&h_nEventsGoodreconsGhost_Si2.store);
      h_nEventsRealGoodrecons_Si2.emplace_back(new TH2F("h_nEventsRealGoodrecons_Si2","h_nEventsRealGoodrecons_Si2", 20, 0, 20, 20, 0, 20));
      HistReg.emplace_back(&h_nEventsRealGoodrecons_Si2.store);
      h_nEventsRealRejectPad_Si2.emplace_back(new TH2F("h_nEventsRealRejectPad_Si2", "h_nEventsRealRejectPad_Si2", 20, 0, 20, 300, 0, 300));
      HistReg.emplace_back(&h_nEventsRealRejectPad_Si2.store);

      h_MFCheck_Theta_MomSi1MomSi2.emplace_back(new TH1F("h_MFCheck_Theta_MomSi1MomSi2", "h_MFCheck_Theta_MomSi1MomSi2", 2000, 0, 5));
      HistReg.emplace_back(&h_MFCheck_Theta_MomSi1MomSi2.store);
      h_MFCheck_Dist_MomSi1HitSi2.emplace_back(new TH1F("h_MFCheck_Dist_MomSi1HitSi2", "h_MFCheck_Dist_MomSi1HitSi2", 2000, 0, 1));
      HistReg.emplace_back(&h_MFCheck_Dist_MomSi1HitSi2.store);
      h_MFCheck_Dist_MomSi2HitSi1.emplace_back(new TH1F("h_MFCheck_Dist_MomSi2HitSi1", "h_MFCheck_Dist_MomSi2HitSi1", 2000, 0, 1));
      HistReg.emplace_back(&h_MFCheck_Dist_MomSi2HitSi1.store);
      h_MFCheck_Dist_MomSi1HitIP.emplace_back(new TH1F("h_MFCheck_Dist_MomSi1HitIP", "h_MFCheck_Dist_MomSi1HitIP", 2000, 0, 1));
      HistReg.emplace_back(&h_MFCheck_Dist_MomSi1HitIP.store);
      h_MFCheck_Dist_MomSi2HitIP.emplace_back(new TH1F("h_MFCheck_Dist_MomSi2HitIP", "h_MFCheck_Dist_MomSi2HitIP", 2000, 0, 1));
      HistReg.emplace_back(&h_MFCheck_Dist_MomSi2HitIP.store);

      h_EnergyStripEnergyTotalReal.emplace_back(new TH2F("h_EnergyStripEnergyTotalReal", "h_EnergyStripEnergyTotalReal", 800, 0, 1.6, 2000, 0, 4));
      HistReg.emplace_back(&h_EnergyStripEnergyTotalReal.store);
      h_EnergyStripEnergyTotal.emplace_back(new TH2F("h_EnergyStripEnergyTotal", "h_EnergyStripEnergyTotal", 800, 0, 1.6, 2000, 0, 4));
      HistReg.emplace_back(&h_EnergyStripEnergyTotal.store);
      h_EnergyDiffSilicons.emplace_back(new TH2F("h_EnergyDiffSilicons","h_EnergyDiffSilicons", 10000, -1, 1, 3, 0, 3));
      HistReg.emplace_back(&h_EnergyDiffSilicons.store);

      h_EnergyDepositionMother.emplace_back(new TH2F("h_EnergyDepositionMother","h_EnergyDepositionMother", 1000, 0, 2, 5, 0, 5));
      HistReg.emplace_back(&h_EnergyDepositionMother.store);
      h_EnergyDepositionDaughters.emplace_back(new TH2F("h_EnergyDepositionDaughters","h_EnergyDepositionDaughters", 1500, 0, 3, 5, 0, 5));
      HistReg.emplace_back(&h_EnergyDepositionDaughters.store);
        
      h_nTrackCandidates.emplace_back(new TH2F("h_nTrackCandidates","h_nTrackCandidates", 400, 0, 400, 5, 0, 5));
      HistReg.emplace_back(&h_nTrackCandidates.store);
      h_DistanceBeamTracks.emplace_back(new TH2F("h_DistanceBeamTracks","h_DistanceBeamTracks", 10000, 0, 0.5, 3, 0, 3));
      HistReg.emplace_back(&h_DistanceBeamTracks.store);
      h_PosZBeamTracks.emplace_back(new TH2F("h_PosZBeamTracks","h_PosZBeamTracks", 1000, 22, 29, 3, 0, 3));
      HistReg.emplace_back(&h_PosZBeamTracks.store);
      h_thetaTracks.emplace_back(new TH2F("h_thetaTracks","h_thetaTracks", 1000, 0, 90, 3, 0, 3));
      HistReg.emplace_back(&h_thetaTracks.store);
      h_thetaResol.emplace_back(new TH1F("h_thetaResol","h_thetaResol", 10000, -5, 5));
      HistReg.emplace_back(&h_thetaResol.store);

      h_Acc_ThetaCandidates.emplace_back(new TH1F("h_Acc_ThetaCandidates", "h_Acc_ThetaCandidates", 1800, 0, 90));
      HistReg.emplace_back(&h_Acc_ThetaCandidates.store);
      h_Acc_ThetaAllReal.emplace_back(new TH1F("h_Acc_ThetaAllReal", "h_Acc_ThetaAllReal", 1800, 0, 90));
      HistReg.emplace_back(&h_Acc_ThetaAllReal.store);

      h_nCandidatesRealTracks.emplace_back(new TH2F("h_nCandidatesRealTracks","h_nCandidatesRealTracks", 100, 0, 100, 100, 0, 100));
      HistReg.emplace_back(&h_nCandidatesRealTracks.store);
      h_nCandidatesRealTracks_IfRecons.emplace_back(new TH2F("h_nCandidatesRealTracks_IfRecons","h_nCandidatesRealTracks_IfRecons", 100, 0, 100, 100, 0, 100));
      HistReg.emplace_back(&h_nCandidatesRealTracks_IfRecons.store);

      h_nHypernucleiTrack.emplace_back(new TH1F("h_nHypernucleiTrack", "h_nHypernucleiTrack", 9, 0, 9));
      HistReg.emplace_back(&h_nHypernucleiTrack.store);
      h_fvalues.emplace_back(new TH1F("h_fvalues", "h_fvalues", 2000, -0.5, 1.5));
      HistReg.emplace_back(&h_fvalues.store);

      h_InteractionPointDistance.emplace_back(new TH1F("h_InteractionPointDistance","h_InteractionPointDistance", 1000, 0, 1));
      HistReg.emplace_back(&h_InteractionPointDistance.store);
      h_InteractionPointDistanceX.emplace_back(new TH1F("h_InteractionPointDistanceX","h_InteractionPointDistanceX", 2000, -1, 1));
      HistReg.emplace_back(&h_InteractionPointDistanceX.store);
      h_InteractionPointDistanceY.emplace_back(new TH1F("h_InteractionPointDistanceY","h_InteractionPointDistanceY", 2000, -1, 1));
      HistReg.emplace_back(&h_InteractionPointDistanceY.store);
      h_InteractionPointDistanceZ.emplace_back(new TH1F("h_InteractionPointDistanceZ","h_InteractionPointDistanceZ", 2000, -1, 1));
      HistReg.emplace_back(&h_InteractionPointDistanceZ.store);

      h_InteractionPointDistanceX_pull.emplace_back(new TH1F("h_InteractionPointDistanceX_pull","h_InteractionPointDistanceX_pull", 2000, -5, 5));
      HistReg.emplace_back(&h_InteractionPointDistanceX_pull.store);
      h_InteractionPointDistanceY_pull.emplace_back(new TH1F("h_InteractionPointDistanceY_pull","h_InteractionPointDistanceY_pull", 2000, -5, 5));
      HistReg.emplace_back(&h_InteractionPointDistanceY_pull.store);
      h_InteractionPointDistanceZ_pull.emplace_back(new TH1F("h_InteractionPointDistanceZ_pull","h_InteractionPointDistanceZ_pull", 2000, -5, 5));
      HistReg.emplace_back(&h_InteractionPointDistanceZ_pull.store);

      h_CovarianceSigmaX.emplace_back(new TH1F("h_SqrtCovariance_SigmaX", "h_SqrtCovariance_SigmaX", 1000, 0, 0.1));
      HistReg.emplace_back(&h_CovarianceSigmaX.store);
      h_CovarianceSigmaY.emplace_back(new TH1F("h_SqrtCovariance_SigmaY", "h_SqrtCovariance_SigmaY", 1000, 0, 0.1));
      HistReg.emplace_back(&h_CovarianceSigmaY.store);
      h_CovarianceSigmaZ.emplace_back(new TH1F("h_SqrtCovariance_SigmaZ", "h_SqrtCovariance_SigmaZ", 2000, 0, 0.2));
      HistReg.emplace_back(&h_CovarianceSigmaZ.store);

      h_IP_DecayDistance.emplace_back(new TH1F("h_IP_DecayDistance","h_IP_DecayDistance", 6000, 0, 24));
      HistReg.emplace_back(&h_IP_DecayDistance.store);
      h_IP_DecayDistanceX.emplace_back(new TH1F("h_IP_DecayDistanceX","h_IP_DecayDistanceX", 1000, -2, 2));
      HistReg.emplace_back(&h_IP_DecayDistanceX.store);
      h_IP_DecayDistanceY.emplace_back(new TH1F("h_IP_DecayDistanceY","h_IP_DecayDistanceY", 1000, -2, 2));
      HistReg.emplace_back(&h_IP_DecayDistanceY.store);
      h_IP_DecayDistanceZ.emplace_back(new TH1F("h_IP_DecayDistanceZ","h_IP_DecayDistanceZ", 6000, -2, 22));
      HistReg.emplace_back(&h_IP_DecayDistanceZ.store);

      h_DecayPositionDistance.emplace_back(new TH1F("h_DecayPositionDistance","h_DecayPositionDistance", 6000, 0, 24));
      HistReg.emplace_back(&h_DecayPositionDistance.store);
      h_DecayPositionDistanceX.emplace_back(new TH1F("h_DecayPositionDistanceX","h_DecayPositionDistanceX", 1000, -2, 2));
      HistReg.emplace_back(&h_DecayPositionDistanceX.store);
      h_DecayPositionDistanceY.emplace_back(new TH1F("h_DecayPositionDistanceY","h_DecayPositionDistanceY", 1000, -2, 2));
      HistReg.emplace_back(&h_DecayPositionDistanceY.store);
      h_DecayPositionDistanceZ.emplace_back(new TH1F("h_DecayPositionDistanceZ","h_DecayPositionDistanceZ", 6000, -2, 22));
      HistReg.emplace_back(&h_DecayPositionDistanceZ.store);

      h_PrimStatus.emplace_back(new TH2F("PrimVtxStatus","PrimVtxStatus",20,0,20,20,0,20));
      HistReg.emplace_back(&h_PrimStatus.store);
      h_PrimVtxstats.emplace_back(new TH1F("PrimVtxstats", "PrimVtxstats", 10, 0, 10));
      HistReg.emplace_back(&h_PrimVtxstats.store);
       
      HistRegisteredByDir.insert(std::make_pair("PrimaryVtx_Si", std::make_tuple(HistReg,0)));
    }

  if(EnableState[DECAYVTX])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      TString namePionType[] = {"_pi-", "_pi+"};

      for(size_t i = 0; i < 1 + EnableState[DECAYVTX_PIPLUS]; ++i)
        {
          h_P_fragments[i].emplace_back(new TH1F("h_P_fragments"+namePionType[i], "h_P_fragments"+namePionType[i], 500, 5, 10));
          HistReg.emplace_back(&h_P_fragments[i].store);
          h_Pt_fragments[i].emplace_back(new TH1F("h_Pt_fragments"+namePionType[i], "h_Pt_fragments"+namePionType[i], 100, 0, 1));
          HistReg.emplace_back(&h_Pt_fragments[i].store);
          h_Pz_fragments[i].emplace_back(new TH1F("h_Pz_fragments"+namePionType[i], "h_Pz_fragments"+namePionType[i], 500, 5, 10));
          HistReg.emplace_back(&h_Pz_fragments[i].store);
          h_Dist_FragmentTrackPrimVtx[i].emplace_back(new TH1F("h_Dist_FragmentTrackPrimVtx"+namePionType[i], "h_Dist_FragmentTrackPrimVtx"+namePionType[i], 100, 0, 10));
          HistReg.emplace_back(&h_Dist_FragmentTrackPrimVtx[i].store);

          h_P_pions[i].emplace_back(new TH1F("h_P_pions"+namePionType[i], "h_P_pions"+namePionType[i], 100, 0, 10));
          HistReg.emplace_back(&h_P_pions[i].store);
          h_Pt_pions[i].emplace_back(new TH1F("h_Pt_pions"+namePionType[i], "h_Pt_pions"+namePionType[i], 200, 0, 2));
          HistReg.emplace_back(&h_Pt_pions[i].store);
          h_Pz_pions[i].emplace_back(new TH1F("h_Pz_pions"+namePionType[i], "h_Pz_pions"+namePionType[i], 150, -5, 10));
          HistReg.emplace_back(&h_Pz_pions[i].store);
          h_Chi2ndf_pions[i].emplace_back(new TH1F("h_Chi2ndf_pions"+namePionType[i], "h_Chi2ndf_pions"+namePionType[i], 500, 0, 50));
          HistReg.emplace_back(&h_Chi2ndf_pions[i].store);

          h_Pt_cutpions[i].emplace_back(new TH1F("h_Pt_cutpions"+namePionType[i], "h_Pt_cutpions"+namePionType[i], 200, 0, 2));
          HistReg.emplace_back(&h_Pt_cutpions[i].store);
          h_Pz_cutpions[i].emplace_back(new TH1F("h_Pz_cutpions"+namePionType[i], "h_Pz_cutpions"+namePionType[i], 150, -5, 10));
          HistReg.emplace_back(&h_Pz_cutpions[i].store);

          h_Nrealpions[i].emplace_back(new TH1F("h_Nrealpions"+namePionType[i], "h_Nrealpions"+namePionType[i], 10, 0, 10));
          HistReg.emplace_back(&h_Nrealpions[i].store);
          h_Ncutpions[i].emplace_back(new TH1F("h_Ncutpions"+namePionType[i], "h_Ncutpions"+namePionType[i], 10, 0, 10));
          HistReg.emplace_back(&h_Ncutpions[i].store);
          h_Npions[i].emplace_back(new TH1F("h_Npions"+namePionType[i], "h_Npions"+namePionType[i], 10, 0, 10));
          HistReg.emplace_back(&h_Npions[i].store);


          h_Closedist_Distance[i].emplace_back(new TH1F("h_Closedist_Distance"+namePionType[i],"h_Closedist_Distance"+namePionType[i], 200, 0, 20));
          HistReg.emplace_back(&h_Closedist_Distance[i].store);
          h_Closedist_PosZ[i].emplace_back(new TH1F("h_Closedist_PosZ"+namePionType[i],"h_Closedist_PosZ"+namePionType[i], 500, -10, 40));
          HistReg.emplace_back(&h_Closedist_PosZ[i].store);
          h_Dist_DecayTrackPrimVtx[i].emplace_back(new TH2F("h_Dist_DecayTrackPrimVtx"+namePionType[i], "h_Dist_DecayTrackPrimVtx"+namePionType[i], 400, 0, 20, 3, 0, 3));
          HistReg.emplace_back(&h_Dist_DecayTrackPrimVtx[i].store);

          h_Closedist_cutDistance[i].emplace_back(new TH1F("h_Closedist_cutDistance"+namePionType[i],"h_Closedist_cutDistance"+namePionType[i], 200, 0, 20));
          HistReg.emplace_back(&h_Closedist_cutDistance[i].store);
          h_Closedist_cutPosZ[i].emplace_back(new TH1F("h_Closedist_cutPosZ"+namePionType[i],"h_Closedist_cutPosZ"+namePionType[i], 500, -10, 40));
          HistReg.emplace_back(&h_Closedist_cutPosZ[i].store);
          h_Dist_cutDecayTrackPrimVtx[i].emplace_back(new TH2F("h_Dist_cutDecayTrackPrimVtx"+namePionType[i], "h_Dist_cutDecayTrackPrimVtx"+namePionType[i], 400, 0, 40, 3, 0, 3));
          HistReg.emplace_back(&h_Dist_cutDecayTrackPrimVtx[i].store);


          h_DecayVertexDistance[i].emplace_back(new TH1F("h_DecayVertexDistance"+namePionType[i], "h_DecayVertexDistance"+namePionType[i], 400, 0, 40));
          HistReg.emplace_back(&h_DecayVertexDistance[i].store);
          h_DecayVertexDistanceX[i].emplace_back(new TH1F("h_DecayVertexDistanceX"+namePionType[i], "h_DecayVertexDistanceX"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceX[i].store);
          h_DecayVertexDistanceY[i].emplace_back(new TH1F("h_DecayVertexDistanceY"+namePionType[i], "h_DecayVertexDistanceY"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceY[i].store);
          h_DecayVertexDistanceZ[i].emplace_back(new TH1F("h_DecayVertexDistanceZ"+namePionType[i], "h_DecayVertexDistanceZ"+namePionType[i], 200, -10, 10));
          HistReg.emplace_back(&h_DecayVertexDistanceZ[i].store);

          h_DecayVertexDistance_centroid[i].emplace_back(new TH1F("h_DecayVertexDistance_centroid"+namePionType[i], "h_DecayVertexDistance_centroid"+namePionType[i], 400, 0, 40));
          HistReg.emplace_back(&h_DecayVertexDistance_centroid[i].store);
          h_DecayVertexDistanceX_centroid[i].emplace_back(new TH1F("h_DecayVertexDistanceX_centroid"+namePionType[i], "h_DecayVertexDistanceX_centroid"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceX_centroid[i].store);
          h_DecayVertexDistanceY_centroid[i].emplace_back(new TH1F("h_DecayVertexDistanceY_centroid"+namePionType[i], "h_DecayVertexDistanceY_centroid"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceY_centroid[i].store);
          h_DecayVertexDistanceZ_centroid[i].emplace_back(new TH1F("h_DecayVertexDistanceZ_centroid"+namePionType[i], "h_DecayVertexDistanceZ_centroid"+namePionType[i], 200, -10, 10));
          HistReg.emplace_back(&h_DecayVertexDistanceZ_centroid[i].store);

          h_DecayVertexDistance_KFPart[i].emplace_back(new TH1F("h_DecayVertexDistance_KFPart"+namePionType[i], "h_DecayVertexDistance_KFPart"+namePionType[i], 400, 0, 40));
          HistReg.emplace_back(&h_DecayVertexDistance_KFPart[i].store);
          h_DecayVertexDistanceX_KFPart[i].emplace_back(new TH1F("h_DecayVertexDistanceX_KFPart"+namePionType[i], "h_DecayVertexDistanceX_KFPart"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceX_KFPart[i].store);
          h_DecayVertexDistanceY_KFPart[i].emplace_back(new TH1F("h_DecayVertexDistanceY_KFPart"+namePionType[i], "h_DecayVertexDistanceY_KFPart"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceY_KFPart[i].store);
          h_DecayVertexDistanceZ_KFPart[i].emplace_back(new TH1F("h_DecayVertexDistanceZ_KFPart"+namePionType[i], "h_DecayVertexDistanceZ_KFPart"+namePionType[i], 200, -10, 10));
          HistReg.emplace_back(&h_DecayVertexDistanceZ_KFPart[i].store);

          h_DecayVertexDistance_KFPart_PrimVtx[i].emplace_back(new TH1F("h_DecayVertexDistance_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexDistance_KFPart_PrimVtx"+namePionType[i], 400, 0, 40));
          HistReg.emplace_back(&h_DecayVertexDistance_KFPart_PrimVtx[i].store);
          h_DecayVertexDistanceX_KFPart_PrimVtx[i].emplace_back(new TH1F("h_DecayVertexDistanceX_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexDistanceX_KFPart_PrimVtx"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceX_KFPart_PrimVtx[i].store);
          h_DecayVertexDistanceY_KFPart_PrimVtx[i].emplace_back(new TH1F("h_DecayVertexDistanceY_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexDistanceY_KFPart_PrimVtx"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceY_KFPart_PrimVtx[i].store);
          h_DecayVertexDistanceZ_KFPart_PrimVtx[i].emplace_back(new TH1F("h_DecayVertexDistanceZ_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexDistanceZ_KFPart_PrimVtx"+namePionType[i], 200, -10, 10));
          HistReg.emplace_back(&h_DecayVertexDistanceZ_KFPart_PrimVtx[i].store);

          h_DecayVertexDistance_KFPart_PrimVtx_Mass[i].emplace_back(new TH1F("h_DecayVertexDistance_KFPart_PrimVtx_Mass"+namePionType[i], "h_DecayVertexDistance_KFPart_PrimVtx_Mass"+namePionType[i], 400, 0, 40));
          HistReg.emplace_back(&h_DecayVertexDistance_KFPart_PrimVtx_Mass[i].store);
          h_DecayVertexDistanceX_KFPart_PrimVtx_Mass[i].emplace_back(new TH1F("h_DecayVertexDistanceX_KFPart_PrimVtx_Mass"+namePionType[i], "h_DecayVertexDistanceX_KFPart_PrimVtx_Mass"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceX_KFPart_PrimVtx_Mass[i].store);
          h_DecayVertexDistanceY_KFPart_PrimVtx_Mass[i].emplace_back(new TH1F("h_DecayVertexDistanceY_KFPart_PrimVtx_Mass"+namePionType[i], "h_DecayVertexDistanceY_KFPart_PrimVtx_Mass"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceY_KFPart_PrimVtx_Mass[i].store);
          h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass[i].emplace_back(new TH1F("h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass"+namePionType[i], "h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass"+namePionType[i], 200, -10, 10));
          HistReg.emplace_back(&h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass[i].store);

          h_DecayVertexcutDistance[i].emplace_back(new TH1F("h_DecayVertexcutDistance"+namePionType[i], "h_DecayVertexcutDistance"+namePionType[i], 400, 0, 40));
          HistReg.emplace_back(&h_DecayVertexcutDistance[i].store);
          h_DecayVertexcutDistanceX[i].emplace_back(new TH1F("h_DecayVertexcutDistanceX"+namePionType[i], "h_DecayVertexcutDistanceX"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexcutDistanceX[i].store);
          h_DecayVertexcutDistanceY[i].emplace_back(new TH1F("h_DecayVertexcutDistanceY"+namePionType[i], "h_DecayVertexcutDistanceY"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexcutDistanceY[i].store);
          h_DecayVertexcutDistanceZ[i].emplace_back(new TH1F("h_DecayVertexcutDistanceZ"+namePionType[i], "h_DecayVertexcutDistanceZ"+namePionType[i], 200, -10, 10));
          HistReg.emplace_back(&h_DecayVertexcutDistanceZ[i].store);

    /*
          h_DecayVertexcutDistance_KFPart[i].emplace_back(new TH1F("h_DecayVertexcutDistance_KFPart"+namePionType[i], "h_DecayVertexcutDistance_KFPart"+namePionType[i], 4000, 0, 40));
          HistReg.emplace_back(&h_DecayVertexcutDistance_KFPart[i].store);
          h_DecayVertexcutDistanceX_KFPart[i].emplace_back(new TH1F("h_DecayVertexcutDistanceX_KFPart"+namePionType[i], "h_DecayVertexcutDistanceX_KFPart"+namePionType[i], 1000, -5, 5));
          HistReg.emplace_back(&h_DecayVertexcutDistanceX_KFPart[i].store);
          h_DecayVertexcutDistanceY_KFPart[i].emplace_back(new TH1F("h_DecayVertexcutDistanceY_KFPart"+namePionType[i], "h_DecayVertexcutDistanceY_KFPart"+namePionType[i], 1000, -5, 5));
          HistReg.emplace_back(&h_DecayVertexcutDistanceY_KFPart[i].store);
          h_DecayVertexcutDistanceZ_KFPart[i].emplace_back(new TH1F("h_DecayVertexcutDistanceZ_KFPart"+namePionType[i], "h_DecayVertexcutDistanceZ_KFPart"+namePionType[i], 2000, -10, 10));
          HistReg.emplace_back(&h_DecayVertexcutDistanceZ_KFPart[i].store);
    */

          h_DecayVertexcutDistance_KFPart_PrimVtx[i].emplace_back(new TH1F("h_DecayVertexcutDistance_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexcutDistance_KFPart_PrimVtx"+namePionType[i], 400, 0, 40));
          HistReg.emplace_back(&h_DecayVertexcutDistance_KFPart_PrimVtx[i].store);
          h_DecayVertexcutDistanceX_KFPart_PrimVtx[i].emplace_back(new TH1F("h_DecayVertexcutDistanceX_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexcutDistanceX_KFPart_PrimVtx"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexcutDistanceX_KFPart_PrimVtx[i].store);
          h_DecayVertexcutDistanceY_KFPart_PrimVtx[i].emplace_back(new TH1F("h_DecayVertexcutDistanceY_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexcutDistanceY_KFPart_PrimVtx"+namePionType[i], 100, -5, 5));
          HistReg.emplace_back(&h_DecayVertexcutDistanceY_KFPart_PrimVtx[i].store);
          h_DecayVertexcutDistanceZ_KFPart_PrimVtx[i].emplace_back(new TH1F("h_DecayVertexcutDistanceZ_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexcutDistanceZ_KFPart_PrimVtx"+namePionType[i], 200, -10, 10));
          HistReg.emplace_back(&h_DecayVertexcutDistanceZ_KFPart_PrimVtx[i].store);

          h_DecayVertexPosZ_real[i].emplace_back(new TH1F("h_DecayVertexPosZ_real"+namePionType[i], "h_DecayVertexPosZ_real"+namePionType[i], 300, -20, 40));
          HistReg.emplace_back(&h_DecayVertexPosZ_real[i].store);
          h_DecayVertexPosZ_vfunction[i].emplace_back(new TH1F("h_DecayVertexPosZ_vfunction"+namePionType[i], "h_DecayVertexPosZ_vfunction"+namePionType[i], 300, -20, 40));
          HistReg.emplace_back(&h_DecayVertexPosZ_vfunction[i].store);
          h_DecayVertexPosZ_centroid[i].emplace_back(new TH1F("h_DecayVertexPosZ_centroid"+namePionType[i], "h_DecayVertexPosZ_centroid"+namePionType[i], 300, -20, 40));
          HistReg.emplace_back(&h_DecayVertexPosZ_centroid[i].store);
          h_DecayVertexPosZ_KFPart[i].emplace_back(new TH1F("h_DecayVertexPosZ_KFPart"+namePionType[i], "h_DecayVertexPosZ_KFPart"+namePionType[i], 300, -20, 40));
          HistReg.emplace_back(&h_DecayVertexPosZ_KFPart[i].store);
          h_DecayVertexPosZ_AllVfunc[i].emplace_back(new TH1F("h_DecayVertexPosZ_AllVfunc"+namePionType[i], "h_DecayVertexPosZ_AllVfunc"+namePionType[i], 300, -20, 40));
          HistReg.emplace_back(&h_DecayVertexPosZ_AllVfunc[i].store);
          h_DecayVertexPosZ_AllCentroid[i].emplace_back(new TH1F("h_DecayVertexPosZ_AllCentroid"+namePionType[i], "h_DecayVertexPosZ_AllCentroid"+namePionType[i], 300, -20, 40));
          HistReg.emplace_back(&h_DecayVertexPosZ_AllCentroid[i].store);
          h_DecayVertexPosZ_AllKFPart[i].emplace_back(new TH1F("h_DecayVertexPosZ_AllKFPart"+namePionType[i], "h_DecayVertexPosZ_AllKFPart"+namePionType[i], 300, -20, 40));
          HistReg.emplace_back(&h_DecayVertexPosZ_AllKFPart[i].store);

          h_N_MotherTracks[i].emplace_back(new TH2F("h_N_MotherTracks"+namePionType[i], "h_N_MotherTracks"+namePionType[i], 10, 0, 10, 500, 0, 5));
          HistReg.emplace_back(&h_N_MotherTracks[i].store);
          h_Dist_DaughterTracks[i].emplace_back(new TH2F("h_Dist_DaughterTracks"+namePionType[i], "h_Dist_DaughterTracks"+namePionType[i], 400, 0, 4, 500, 0, 5));
          HistReg.emplace_back(&h_Dist_DaughterTracks[i].store);
          h_Angle_MotherFragment[i].emplace_back(new TH2F("h_Angle_MotherFragment"+namePionType[i], "h_Angle_MotherFragment"+namePionType[i], 100, 0, 10, 500, 0, 5));
          HistReg.emplace_back(&h_Angle_MotherFragment[i].store);
          h_Angle_MotherPion[i].emplace_back(new TH2F("h_Angle_MotherPion"+namePionType[i], "h_Angle_MotherPion"+namePionType[i], 400, 0, 40, 500, 0, 5));
          HistReg.emplace_back(&h_Angle_MotherPion[i].store);
          h_Chi2ndf_MotherTracks[i].emplace_back(new TH2F("h_Chi2ndf_MotherTracks"+namePionType[i], "h_Chi2ndf_MotherTracks"+namePionType[i], 200, 0, 1000, 500, 0, 5));
          HistReg.emplace_back(&h_Chi2ndf_MotherTracks[i].store);
          h_Dist_MotherTrackPrimVtx[i].emplace_back(new TH2F("h_Dist_MotherTrackPrimVtx"+namePionType[i], "h_Dist_MotherTrackPrimVtx"+namePionType[i], 100, 0, 0.2, 500, 0, 5));
          HistReg.emplace_back(&h_Dist_MotherTrackPrimVtx[i].store);
          h_Theta_MotherTrackPrimVtx[i].emplace_back(new TH2F("h_Theta_MotherTrackPrimVtx"+namePionType[i], "h_Theta_MotherTrackPrimVtx"+namePionType[i], 100, 0, 1, 500, 0, 5));
          HistReg.emplace_back(&h_Theta_MotherTrackPrimVtx[i].store);
          h_DecayVertexPosZ_KFPart_PrimVtx[i].emplace_back(new TH2F("h_DecayVertexPosZ_KFPart_PrimVtx"+namePionType[i], "h_DecayVertexPosZ_KFPart_PrimVtx"+namePionType[i], 500, -10, 40, 500, 0, 5));
          HistReg.emplace_back(&h_DecayVertexPosZ_KFPart_PrimVtx[i].store);
          h_DecayFragmentMomZ_KFPart_PrimVtx[i].emplace_back(new TH2F("h_DecayFragmentMomZ_KFPart_PrimVtx"+namePionType[i], "h_DecayFragmentMomZ_KFPart_PrimVtx"+namePionType[i], 500, 5, 10, 500, 0, 5));
          HistReg.emplace_back(&h_DecayFragmentMomZ_KFPart_PrimVtx[i].store);
          h_DecayPionMomZ_KFPart_PrimVtx[i].emplace_back(new TH2F("h_DecayPionMomZ_KFPart_PrimVtx"+namePionType[i], "h_DecayPionMomZ_KFPart_PrimVtx"+namePionType[i], 150, -5, 10, 500, 0, 5));
          HistReg.emplace_back(&h_DecayPionMomZ_KFPart_PrimVtx[i].store);
          h_Hyp_ArmenterosPodolanski[i].emplace_back(new TH2F("h_Hyp_ArmenterosPodolanski"+namePionType[i], "h_Hyp_ArmenterosPodolanski"+namePionType[i], 100, 0, 1, 100, 0, 1));
          HistReg.emplace_back(&h_Hyp_ArmenterosPodolanski[i].store);
          h_Hyp_CutArmenterosPodolanski[i].emplace_back(new TH2F("h_Hyp_CutArmenterosPodolanski"+namePionType[i], "h_Hyp_CutArmenterosPodolanski"+namePionType[i], 100, 0, 1, 100, 0, 1));
          HistReg.emplace_back(&h_Hyp_CutArmenterosPodolanski[i].store);

          h_HypInvariantMass[i].emplace_back(new TH1F("h_HypInvariantMass"+namePionType[i], "h_HypInvariantMass"+namePionType[i], 500, 0, 5));
          HistReg.emplace_back(&h_HypInvariantMass[i].store);
          h_HypInvariantMass_Z05[i].emplace_back(new TH1F("h_HypInvariantMass_Z05"+namePionType[i], "h_HypInvariantMass_Z05"+namePionType[i], 500, 0, 5));
          HistReg.emplace_back(&h_HypInvariantMass_Z05[i].store);
          h_HypInvariantMass_Z10[i].emplace_back(new TH1F("h_HypInvariantMass_Z10"+namePionType[i], "h_HypInvariantMass_Z10"+namePionType[i], 500, 0, 5));
          HistReg.emplace_back(&h_HypInvariantMass_Z10[i].store);
          h_HypInvariantMass_Z15[i].emplace_back(new TH1F("h_HypInvariantMass_Z15"+namePionType[i], "h_HypInvariantMass_Z15"+namePionType[i], 500, 0, 5));
          HistReg.emplace_back(&h_HypInvariantMass_Z15[i].store);
          h_HypInvariantMass_Z20[i].emplace_back(new TH1F("h_HypInvariantMass_Z20"+namePionType[i], "h_HypInvariantMass_Z20"+namePionType[i], 500, 0, 5));
          HistReg.emplace_back(&h_HypInvariantMass_Z20[i].store);
          h_HypErrorInvariantMass[i].emplace_back(new TH1F("h_HypErrorInvariantMass"+namePionType[i], "h_HypErrorInvariantMass"+namePionType[i], 1000, 0, 0.5));
          HistReg.emplace_back(&h_HypErrorInvariantMass[i].store);
          

          h_Hyp_RealLifeTime[i].emplace_back(new TH1F("h_Hyp_RealLifeTime"+namePionType[i], "h_Hyp_RealLifeTime"+namePionType[i], 100, 0, 1000));
          HistReg.emplace_back(&h_Hyp_RealLifeTime[i].store);
          h_HypLifeTime_PrimVtx[i].emplace_back(new TH1F("h_HypLifeTime_PrimVtx"+namePionType[i], "h_HypLifeTime_PrimVtx"+namePionType[i], 100, 0, 1000));
          HistReg.emplace_back(&h_HypLifeTime_PrimVtx[i].store);
          h_HypErrorLifeTime_PrimVtx[i].emplace_back(new TH1F("h_HypErrorLifeTime_PrimVtx"+namePionType[i], "h_HypErrorLifeTime_PrimVtx"+namePionType[i], 200, 0, 200));
          HistReg.emplace_back(&h_HypErrorLifeTime_PrimVtx[i].store);
          h_HypcutLifeTime_PrimVtx[i].emplace_back(new TH1F("h_HypcutLifeTime_PrimVtx"+namePionType[i], "h_HypcutLifeTime_PrimVtx"+namePionType[i], 100, 0, 1000));
          HistReg.emplace_back(&h_HypcutLifeTime_PrimVtx[i].store);


          h_HypInvariantMassCheck[i].emplace_back(new TH2F("h_HypInvariantMassCheck"+namePionType[i], "h_HypInvariantMassCheck"+namePionType[i], 1000, -10, 10, 2, 0, 2));
          HistReg.emplace_back(&h_HypInvariantMassCheck[i].store);
          h_HypInvariantErrorMassCheck[i].emplace_back(new TH2F("h_HypInvariantErrorMassCheck"+namePionType[i], "h_HypInvariantErrorMassCheck"+namePionType[i], 150, -300, 1200, 2, 0, 2));
          HistReg.emplace_back(&h_HypInvariantErrorMassCheck[i].store);


          h_HypInvariantMass_LorentzVect[i].emplace_back(new TH1F("h_HypInvariantMass_LorentzVect"+namePionType[i], "h_HypInvariantMass_LorentzVect"+namePionType[i], 500, 0, 5));
          HistReg.emplace_back(&h_HypInvariantMass_LorentzVect[i].store);
          h_HypInvariantMass_CutLorentzVect[i].emplace_back(new TH1F("h_HypInvariantMass_CutLorentzVect"+namePionType[i], "h_HypInvariantMass_CutLorentzVect"+namePionType[i], 500, 0, 5));
          HistReg.emplace_back(&h_HypInvariantMass_CutLorentzVect[i].store);


          h_EffPosZ_real[i].emplace_back(new TH1F("h_EffPosZ_real"+namePionType[i], "h_EffPosZ_real"+namePionType[i], 600, -10, 50));
          HistReg.emplace_back(&h_EffPosZ_real[i].store);
          h_EffPosZ_preKF[i].emplace_back(new TH1F("h_EffPosZ_preKF"+namePionType[i], "h_EffPosZ_preKF"+namePionType[i], 600, -10, 50));
          HistReg.emplace_back(&h_EffPosZ_preKF[i].store);
          h_EffPosZ_postKF[i].emplace_back(new TH1F("h_EffPosZ_postKF"+namePionType[i], "h_EffPosZ_postKF"+namePionType[i], 600, -10, 50));
          HistReg.emplace_back(&h_EffPosZ_postKF[i].store);
          h_EffPosZ_preKFPart[i].emplace_back(new TH1F("h_EffPosZ_preKFPart"+namePionType[i], "h_EffPosZ_preKFPart"+namePionType[i], 600, -10, 50));
          HistReg.emplace_back(&h_EffPosZ_preKFPart[i].store);
          h_EffPosZ_postKFPart[i].emplace_back(new TH1F("h_EffPosZ_postKFPart"+namePionType[i], "h_EffPosZ_postKFPart"+namePionType[i], 600, -10, 50));
          HistReg.emplace_back(&h_EffPosZ_postKFPart[i].store);

          h_EffPosZPosR_real[i].emplace_back(new TH2F("h_EffPosZPosR_real"+namePionType[i], "h_EffPosZPosR_real"+namePionType[i], 600, -10, 50, 50, 0, 5));
          HistReg.emplace_back(&h_EffPosZPosR_real[i].store);
          h_EffPosZPosR_postKFPart[i].emplace_back(new TH2F("h_EffPosZPosR_postKFPart"+namePionType[i], "h_EffPosZPosR_postKFPart"+namePionType[i], 600, -10, 50, 50, 0, 5));
          HistReg.emplace_back(&h_EffPosZPosR_postKFPart[i].store);

    /*
          h_N_SiHits_ReconsTracks[i].emplace_back(new TH2F("h_N_SiHits_ReconsTracks"+namePionType[i], "h_N_SiHits_ReconsTracks"+namePionType[i], 20, 0, 20, 6, 0, 6));
          HistReg.emplace_back(&h_N_SiHits_ReconsTracks[i].store);

          h_N_Si_MotherTracks[i].emplace_back(new TH1F("h_N_Si_MotherTracks"+namePionType[i], "h_N_Si_MotherTracks"+namePionType[i], 10, 0, 10));
          HistReg.emplace_back(&h_N_Si_MotherTracks[i].store);


          h_DecayVertexDistance_AllVfunc[i].emplace_back(new TH1F("h_DecayVertexDistance_AllVfunc"+namePionType[i], "h_DecayVertexDistance_AllVfunc"+namePionType[i], 4000, 0, 40));
          HistReg.emplace_back(&h_DecayVertexDistance_AllVfunc[i].store);
          h_DecayVertexDistanceX_AllVfunc[i].emplace_back(new TH1F("h_DecayVertexDistanceX_AllVfunc"+namePionType[i], "h_DecayVertexDistanceX_AllVfunc"+namePionType[i], 1000, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceX_AllVfunc[i].store);
          h_DecayVertexDistanceY_AllVfunc[i].emplace_back(new TH1F("h_DecayVertexDistanceY_AllVfunc"+namePionType[i], "h_DecayVertexDistanceY_AllVfunc"+namePionType[i], 1000, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceY_AllVfunc[i].store);
          h_DecayVertexDistanceZ_AllVfunc[i].emplace_back(new TH1F("h_DecayVertexDistanceZ_AllVfunc"+namePionType[i], "h_DecayVertexDistanceZ_AllVfunc"+namePionType[i], 2000, -10, 10));
          HistReg.emplace_back(&h_DecayVertexDistanceZ_AllVfunc[i].store);

          h_DecayVertexDistance_AllCentroid[i].emplace_back(new TH1F("h_DecayVertexDistance_AllCentroid"+namePionType[i], "h_DecayVertexDistance_AllCentroid"+namePionType[i], 4000, 0, 40));
          HistReg.emplace_back(&h_DecayVertexDistance_AllCentroid[i].store);
          h_DecayVertexDistanceX_AllCentroid[i].emplace_back(new TH1F("h_DecayVertexDistanceX_AllCentroid"+namePionType[i], "h_DecayVertexDistanceX_AllCentroid"+namePionType[i], 1000, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceX_AllCentroid[i].store);
          h_DecayVertexDistanceY_AllCentroid[i].emplace_back(new TH1F("h_DecayVertexDistanceY_AllCentroid"+namePionType[i], "h_DecayVertexDistanceY_AllCentroid"+namePionType[i], 1000, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceY_AllCentroid[i].store);
          h_DecayVertexDistanceZ_AllCentroid[i].emplace_back(new TH1F("h_DecayVertexDistanceZ_AllCentroid"+namePionType[i], "h_DecayVertexDistanceZ_AllCentroid"+namePionType[i], 2000, -10, 10));
          HistReg.emplace_back(&h_DecayVertexDistanceZ_AllCentroid[i].store);

          h_DecayVertexDistance_AllKFPart[i].emplace_back(new TH1F("h_DecayVertexDistance_AllKFPart"+namePionType[i], "h_DecayVertexDistance_AllKFPart"+namePionType[i], 4000, 0, 40));
          HistReg.emplace_back(&h_DecayVertexDistance_AllKFPart[i].store);
          h_DecayVertexDistanceX_AllKFPart[i].emplace_back(new TH1F("h_DecayVertexDistanceX_AllKFPart"+namePionType[i], "h_DecayVertexDistanceX_AllKFPart"+namePionType[i], 1000, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceX_AllKFPart[i].store);
          h_DecayVertexDistanceY_AllKFPart[i].emplace_back(new TH1F("h_DecayVertexDistanceY_AllKFPart"+namePionType[i], "h_DecayVertexDistanceY_AllKFPart"+namePionType[i], 1000, -5, 5));
          HistReg.emplace_back(&h_DecayVertexDistanceY_AllKFPart[i].store);
          h_DecayVertexDistanceZ_AllKFPart[i].emplace_back(new TH1F("h_DecayVertexDistanceZ_AllKFPart"+namePionType[i], "h_DecayVertexDistanceZ_AllKFPart"+namePionType[i], 2000, -10, 10));
          HistReg.emplace_back(&h_DecayVertexDistanceZ_AllKFPart[i].store);
    */
          h_DecayVtxstats[i].emplace_back(new TH1F("h_DecayVtxstats"+namePionType[i], "h_DecayVtxstats"+namePionType[i], 10, 0, 10));
          HistReg.emplace_back(&h_DecayVtxstats[i].store);
        }

      HistRegisteredByDir.insert(std::make_pair("DecayVtx", std::make_tuple(HistReg,0)));               
    }


  _logger->info( " : done !");
}

void Ana_Hist::DebugHists()
{
  _logger->debug("Histo debuging ");
  _logger->debug("h_stats: {} / {}",fmt::ptr(h_stats.h), h_stats.h->GetEntries());
  for(auto& hh : h_stats.store)
    _logger->debug("h_stats: storing : {} / {}",fmt::ptr(hh), hh->GetEntries());

  return;  
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

  auto Merging = [&](std::vector<TH1*>& vecH) {
    //_logger->debug("merging function: {}",vecH.size());
    if(vecH.size()==1)
      return;

    //_logger->debug("histo vec {}, name {}, n {}",vecH.size(),vecH[0]->GetName(), vecH[0]->GetEntries());
    TList list;
    for(size_t i = 1; i<vecH.size();++i)
      {
	//_logger->debug("#{}, {}",i,vecH[i]->GetEntries());
	list.Add(vecH[i]);
      }
    
    vecH[0]->Reset();
    vecH[0]->Merge(&list);
    //_logger->debug("merged : {}",vecH[0]->GetEntries());
    return;
  };
  
    auto f_dEffdTheta = [](TH2F* h)
  {
    TString nameAll("Proj_All");
    TH1D* h_All = h->ProjectionX(nameAll,1,1);
    h_All->Sumw2();

    TString nameAcc("Proj_Acc");
    TH1D* h_Acc = h->ProjectionX(nameAcc,2,2);
    h_Acc->Sumw2();

    TString nameH = h->GetName();
    Int_t nBins = h->GetXaxis()->GetNbins();
    Double_t minedge = h->GetXaxis()->GetXmin();
    Double_t maxedge = h->GetXaxis()->GetXmax();

    TH1F* h_dEffdTheta = new TH1F(nameH, nameH, nBins, minedge, maxedge);
    h_dEffdTheta->Divide(h_Acc, h_All);

    return h_dEffdTheta;
  };

  _logger->info( "making directory ");
  for(auto it : HistRegisteredByDir)
    {
      TDirectory* temp_dir = GetDir(out_file, it.first);
      temp_dir->cd();
      //_logger->debug("dir: {}",it.first);
      for(auto& it_hist : std::get<0>(it.second))
        {
	  //_logger->debug("vecHist[0]: {} / {} {}",it_hist->size(),fmt::ptr(it_hist->at(0)),it_hist->at(0)->GetName());
	  Merging(*it_hist);

    if(std::get<1>(it.second)==2)
      {
        TH1F* h_dEff;
        h_dEff = f_dEffdTheta(dynamic_cast<TH2F*>(it_hist->at(0)));
        h_dEff->Write();
        continue;
      }

	  it_hist->at(0)->Write();
	  if(std::get<1>(it.second)==1)
	    {
	      TGraphErrors* g1;
	      TGraphErrors* g2;
	      
	      std::tie(g1,g2) = f_DiffResolution(dynamic_cast<TH2F*>(it_hist->at(0)));
	      if(g1!=nullptr)
		g1->Write();
	      if(g2!=nullptr)
		g2->Write();
	    }
	}
      out_file->cd();
    }

  TDirectory* temp_dir = GetDir(out_file, "Finder");
  temp_dir->cd();
  for(auto el : geoSolenoid)
    if(el!=nullptr)
      el->Write();
  
  return 0;
}

int Ana_Hist::WriteTemp(char* tempfile)
{
  HaveBeenWritten = true;

  TFile* ff = new TFile(tempfile, "RECREATE");
  _logger->info( "File = {}", tempfile);

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
