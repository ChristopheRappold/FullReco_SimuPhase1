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
Ana_Hist::Ana_Hist(bool Daf, bool Vertex, bool DCproject, bool Finding, bool Hough, bool Simu, bool PrimVtx, bool DecayVtx)
{
  EnableState.resize(SIZEOF_STATEHIST);
  EnableState[DAF] = Daf;
  EnableState[VERTEX] = Vertex;
  EnableState[DCPROJ] = DCproject;
  EnableState[FINDING] = Finding;
  EnableState[HOUGH] = Hough;
  EnableState[SIMU] = Simu;
  EnableState[PRIMVTX] = PrimVtx;
  EnableState[DECAYVTX] = DecayVtx;



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
      FieldXY[i].emplace_back(new TH2F(name_h_field.c_str(),name_h_field.c_str(),600,-150,150,600,-150,150));
      h_fields.emplace_back(&FieldXY[i].store);

      std::string name_h_field1 = "FieldXZ_";
      name_h_field1+=name_field[i];
      FieldXZ[i].emplace_back(new TH2F(name_h_field1.c_str(),name_h_field1.c_str(),600,-150,150,1200,-300,300));
      h_fields.emplace_back(&FieldXZ[i].store);

      std::string name_h_field2 = "FieldYZ_";
      name_h_field2+=name_field[i];
      FieldYZ[i].emplace_back(new TH2F(name_h_field2.c_str(),name_h_field2.c_str(),600,-150,150,1200,-300,300));
      h_fields.emplace_back(&FieldYZ[i].store);

      std::string name_h_field_n = "FieldXYn_";
      name_h_field_n+=name_field[i];
      FieldXY_n[i].emplace_back(new TH2F(name_h_field_n.c_str(),name_h_field_n.c_str(),600,-150,150,600,-150,150));
      h_fields.emplace_back(&FieldXY_n[i].store);

      std::string name_h_field_n1 = "FieldXZn_";
      name_h_field_n1+=name_field[i];
      FieldXZ_n[i].emplace_back(new TH2F(name_h_field_n1.c_str(),name_h_field_n1.c_str(),600,-150,150,1200,-300,300));
      h_fields.emplace_back(&FieldXZ_n[i].store);

      std::string name_h_field_n2 = "FieldYZn_";
      name_h_field_n2+=name_field[i];
      FieldYZ_n[i].emplace_back(new TH2F(name_h_field_n2.c_str(),name_h_field_n2.c_str(),600,-150,150,1200,-300,300));
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

      h_pv_mom.emplace_back(new TH2F("mom_pv", "mom_pv", 200, 0., 20., 250, 0, 1));
      HistReg.emplace_back(&h_pv_mom.store);
      h_pv_beta.emplace_back(new TH2F("beta_pv", "beta_pv", 250, 0., 2.5, 250, 0, 1));
      HistReg.emplace_back(&h_pv_beta.store);
      h_pv_mass.emplace_back(new TH2F("mass_pv", "mass_pv", 200, 0., 10., 250, 0, 1));
      HistReg.emplace_back(&h_pv_mass.store);

      h_path_tof.emplace_back(new TH2F("Mean_path_TOF", "Mean_path_TOF", 500, 15., 20., 90, 10, 40));
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
	  h_ResFiber[i].emplace_back( new TH1F(Form("ResFiber%d",i),Form("ResFiber%d",i),300,-0.2,0.2));
	  h_residual.emplace_back(&h_ResFiber[i].store);
	}
      for(int i=0;i<6;++i)
	{
	  h_ResMiniFiber[i].emplace_back( new TH1F(Form("ResMiniFiber%d",i),Form("ResMiniFiber%d",i),300,-0.2,0.2));
	  h_residual.emplace_back(&h_ResMiniFiber[i].store);
	}
      for(int i=0;i<2;++i)
	{
	  h_ResPSCE[i].emplace_back( new TH1F(Form("ResPSCE%d",i),Form("ResPSCE%d",i),300,-5,5));
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

  if(EnableState[PRIMVTX])
    {
      std::vector<std::vector<TH1*>*> HistReg;

      h_HitMultiplicity_Si1.emplace_back(new TH1F("h_HitMultiplicity_Si1","h_HitMultiplicity_Si1", 20, 0, 20));
      HistReg.emplace_back(&h_HitMultiplicity_Si1.store);
      h_HitMultiplicityRecons_Si1.emplace_back(new TH1F("h_HitMultiplicityRecons_Si1","h_HitMultiplicityRecons_Si1", 30, 0, 30));
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



       h_HitMultiplicity_Si2.emplace_back(new TH1F("h_HitMultiplicity_Si2","h_HitMultiplicity_Si2", 20, 0, 20));
       HistReg.emplace_back(&h_HitMultiplicity_Si2.store);
       h_HitMultiplicityRecons_Si2.emplace_back(new TH1F("h_HitMultiplicityRecons_Si2","h_HitMultiplicityRecons_Si2", 30, 0, 30));
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


       h_EnergyStripEnergyTotalReal.emplace_back(new TH2F("h_EnergyStripEnergyTotalReal", "h_EnergyStripEnergyTotalReal", 400, 0, 0.8, 2000, 0, 4));
       HistReg.emplace_back(&h_EnergyStripEnergyTotalReal.store);
       h_EnergyStripEnergyTotal.emplace_back(new TH2F("h_EnergyStripEnergyTotal", "h_EnergyStripEnergyTotal", 400, 0, 0.8, 2000, 0, 4));
       HistReg.emplace_back(&h_EnergyStripEnergyTotal.store);
       h_EnergyDiffSilicons.emplace_back(new TH2F("h_EnergyDiffSilicons","h_EnergyDiffSilicons", 10000, -1, 1, 3, 0, 3));
       HistReg.emplace_back(&h_EnergyDiffSilicons.store);

       h_EnergyDepositionMother.emplace_back(new TH2F("h_EnergyDepositionMother","h_EnergyDepositionMother", 1000, 0, 1, 5, 0, 5));
       HistReg.emplace_back(&h_EnergyDepositionMother.store);
       h_EnergyDepositionDaughters.emplace_back(new TH2F("h_EnergyDepositionDaughters","h_EnergyDepositionDaughters", 1500, 0, 1.5, 5, 0, 5));
       HistReg.emplace_back(&h_EnergyDepositionDaughters.store);

        
       h_nTrackCandidates.emplace_back(new TH2F("h_nTrackCandidates","h_nTrackCandidates", 400, 0, 400, 5, 0, 5));
       HistReg.emplace_back(&h_nTrackCandidates.store);
       h_DistanceBeamTracks.emplace_back(new TH2F("h_DistanceBeamTracks","h_DistanceBeamTracks", 10000, 0, 0.5, 3, 0, 3));
       HistReg.emplace_back(&h_DistanceBeamTracks.store);
       h_PosZBeamTracks.emplace_back(new TH2F("h_PosZBeamTracks","h_PosZBeamTracks", 1000, 22, 29, 3, 0, 3));
       HistReg.emplace_back(&h_PosZBeamTracks.store);
       h_thetaTracks.emplace_back(new TH2F("h_thetaTracks","h_thetaTracks", 1000, 0, 90, 3, 0, 3));
       HistReg.emplace_back(&h_thetaTracks.store);

       //h_nHypTrackReal.emplace_back(new TH1F("h_nHypTrackReal", "h_nHypTrackReal", 3, 0, 3));
       //HistReg.emplace_back(&h_nHypTrackReal.store);
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
       
       HistRegisteredByDir.insert(std::make_pair("PrimaryVtx", std::make_tuple(HistReg,0)));
    }

  if(EnableState[DECAYVTX])
    {
      std::vector<std::vector<TH1*>*> HistReg;
      
      h_Pt_fragments.emplace_back(new TH1F("h_Pt_fragments", "h_Pt_fragments", 1000, 0, 1));
      HistReg.emplace_back(&h_Pt_fragments.store);
      h_Pz_fragments.emplace_back(new TH1F("h_Pz_fragments", "h_Pz_fragments", 1500, 0, 15));
      HistReg.emplace_back(&h_Pz_fragments.store);
      h_Dist_FragmentTrackPrimVtx.emplace_back(new TH1F("h_Dist_FragmentTrackPrimVtx", "h_Dist_FragmentTrackPrimVtx", 10000, 0, 10));
      HistReg.emplace_back(&h_Dist_FragmentTrackPrimVtx.store);

      h_Pt_pions.emplace_back(new TH1F("h_Pt_pions", "h_Pt_pions", 1000, 0, 1));
      HistReg.emplace_back(&h_Pt_pions.store);
      h_Pz_pions.emplace_back(new TH1F("h_Pz_pions", "h_Pz_pions", 1000, -3, 5));
      HistReg.emplace_back(&h_Pz_pions.store);
      h_Chi2ndf_pions.emplace_back(new TH1F("h_Chi2ndf_pions", "h_Chi2ndf_pions", 1000, 0, 10));
      HistReg.emplace_back(&h_Chi2ndf_pions.store);

      h_Pt_cutpions.emplace_back(new TH1F("h_Pt_cutpions", "h_Pt_cutpions", 1000, 0, 1));
      HistReg.emplace_back(&h_Pt_cutpions.store);
      h_Pz_cutpions.emplace_back(new TH1F("h_Pz_cutpions", "h_Pz_cutpions", 1000, -3, 5));
      HistReg.emplace_back(&h_Pz_cutpions.store);

      h_Nrealpions.emplace_back(new TH1F("h_Nrealpions", "h_Nrealpions", 10, 0, 10));
      HistReg.emplace_back(&h_Nrealpions.store);
      h_Ncutpions.emplace_back(new TH1F("h_Ncutpions", "h_Ncutpions", 10, 0, 10));
      HistReg.emplace_back(&h_Ncutpions.store);
      h_Npions.emplace_back(new TH1F("h_Npions", "h_Npions", 10, 0, 10));
      HistReg.emplace_back(&h_Npions.store);


      h_Closedist_Distance.emplace_back(new TH1F("h_Closedist_Distance","h_Closedist_Distance", 1000, 0, 5));
      HistReg.emplace_back(&h_Closedist_Distance.store);
      h_Closedist_PosZ.emplace_back(new TH1F("h_Closedist_PosZ","h_Closedist_PosZ", 3000, 20, 50));
      HistReg.emplace_back(&h_Closedist_PosZ.store);
      h_Dist_DecayTrackPrimVtx.emplace_back(new TH2F("h_Dist_DecayTrackPrimVtx", "h_Dist_DecayTrackPrimVtx", 2000, 0, 20, 3, 0, 3));
      HistReg.emplace_back(&h_Dist_DecayTrackPrimVtx.store);

      h_Closedist_cutDistance.emplace_back(new TH1F("h_Closedist_cutDistance","h_Closedist_cutDistance", 1000, 0, 5));
      HistReg.emplace_back(&h_Closedist_cutDistance.store);
      h_Closedist_cutPosZ.emplace_back(new TH1F("h_Closedist_cutPosZ","h_Closedist_cutPosZ", 3000, 20, 50));
      HistReg.emplace_back(&h_Closedist_cutPosZ.store);
      h_Dist_cutDecayTrackPrimVtx.emplace_back(new TH2F("h_Dist_cutDecayTrackPrimVtx", "h_Dist_cutDecayTrackPrimVtx", 2000, 0, 20, 3, 0, 3));
      HistReg.emplace_back(&h_Dist_cutDecayTrackPrimVtx.store);


      h_DecayVertexDistance.emplace_back(new TH1F("h_DecayVertexDistance", "h_DecayVertexDistance", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexDistance.store);
      h_DecayVertexDistanceX.emplace_back(new TH1F("h_DecayVertexDistanceX", "h_DecayVertexDistanceX", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceX.store);
      h_DecayVertexDistanceY.emplace_back(new TH1F("h_DecayVertexDistanceY", "h_DecayVertexDistanceY", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceY.store);
      h_DecayVertexDistanceZ.emplace_back(new TH1F("h_DecayVertexDistanceZ", "h_DecayVertexDistanceZ", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexDistanceZ.store);

      h_DecayVertexDistance_centroid.emplace_back(new TH1F("h_DecayVertexDistance_centroid", "h_DecayVertexDistance_centroid", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexDistance_centroid.store);
      h_DecayVertexDistanceX_centroid.emplace_back(new TH1F("h_DecayVertexDistanceX_centroid", "h_DecayVertexDistanceX_centroid", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceX_centroid.store);
      h_DecayVertexDistanceY_centroid.emplace_back(new TH1F("h_DecayVertexDistanceY_centroid", "h_DecayVertexDistanceY_centroid", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceY_centroid.store);
      h_DecayVertexDistanceZ_centroid.emplace_back(new TH1F("h_DecayVertexDistanceZ_centroid", "h_DecayVertexDistanceZ_centroid", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexDistanceZ_centroid.store);

      h_DecayVertexDistance_KFPart.emplace_back(new TH1F("h_DecayVertexDistance_KFPart", "h_DecayVertexDistance_KFPart", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexDistance_KFPart.store);
      h_DecayVertexDistanceX_KFPart.emplace_back(new TH1F("h_DecayVertexDistanceX_KFPart", "h_DecayVertexDistanceX_KFPart", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceX_KFPart.store);
      h_DecayVertexDistanceY_KFPart.emplace_back(new TH1F("h_DecayVertexDistanceY_KFPart", "h_DecayVertexDistanceY_KFPart", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceY_KFPart.store);
      h_DecayVertexDistanceZ_KFPart.emplace_back(new TH1F("h_DecayVertexDistanceZ_KFPart", "h_DecayVertexDistanceZ_KFPart", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexDistanceZ_KFPart.store);

      h_DecayVertexDistance_KFPart_PrimVtx.emplace_back(new TH1F("h_DecayVertexDistance_KFPart_PrimVtx", "h_DecayVertexDistance_KFPart_PrimVtx", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexDistance_KFPart_PrimVtx.store);
      h_DecayVertexDistanceX_KFPart_PrimVtx.emplace_back(new TH1F("h_DecayVertexDistanceX_KFPart_PrimVtx", "h_DecayVertexDistanceX_KFPart_PrimVtx", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceX_KFPart_PrimVtx.store);
      h_DecayVertexDistanceY_KFPart_PrimVtx.emplace_back(new TH1F("h_DecayVertexDistanceY_KFPart_PrimVtx", "h_DecayVertexDistanceY_KFPart_PrimVtx", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceY_KFPart_PrimVtx.store);
      h_DecayVertexDistanceZ_KFPart_PrimVtx.emplace_back(new TH1F("h_DecayVertexDistanceZ_KFPart_PrimVtx", "h_DecayVertexDistanceZ_KFPart_PrimVtx", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexDistanceZ_KFPart_PrimVtx.store);

      h_DecayVertexDistance_KFPart_PrimVtx_Mass.emplace_back(new TH1F("h_DecayVertexDistance_KFPart_PrimVtx_Mass", "h_DecayVertexDistance_KFPart_PrimVtx_Mass", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexDistance_KFPart_PrimVtx_Mass.store);
      h_DecayVertexDistanceX_KFPart_PrimVtx_Mass.emplace_back(new TH1F("h_DecayVertexDistanceX_KFPart_PrimVtx_Mass", "h_DecayVertexDistanceX_KFPart_PrimVtx_Mass", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceX_KFPart_PrimVtx_Mass.store);
      h_DecayVertexDistanceY_KFPart_PrimVtx_Mass.emplace_back(new TH1F("h_DecayVertexDistanceY_KFPart_PrimVtx_Mass", "h_DecayVertexDistanceY_KFPart_PrimVtx_Mass", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceY_KFPart_PrimVtx_Mass.store);
      h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass.emplace_back(new TH1F("h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass", "h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexDistanceZ_KFPart_PrimVtx_Mass.store);

      h_DecayVertexcutDistance.emplace_back(new TH1F("h_DecayVertexcutDistance", "h_DecayVertexcutDistance", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexcutDistance.store);
      h_DecayVertexcutDistanceX.emplace_back(new TH1F("h_DecayVertexcutDistanceX", "h_DecayVertexcutDistanceX", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexcutDistanceX.store);
      h_DecayVertexcutDistanceY.emplace_back(new TH1F("h_DecayVertexcutDistanceY", "h_DecayVertexcutDistanceY", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexcutDistanceY.store);
      h_DecayVertexcutDistanceZ.emplace_back(new TH1F("h_DecayVertexcutDistanceZ", "h_DecayVertexcutDistanceZ", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexcutDistanceZ.store);

/*
      h_DecayVertexcutDistance_KFPart.emplace_back(new TH1F("h_DecayVertexcutDistance_KFPart", "h_DecayVertexcutDistance_KFPart", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexcutDistance_KFPart.store);
      h_DecayVertexcutDistanceX_KFPart.emplace_back(new TH1F("h_DecayVertexcutDistanceX_KFPart", "h_DecayVertexcutDistanceX_KFPart", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexcutDistanceX_KFPart.store);
      h_DecayVertexcutDistanceY_KFPart.emplace_back(new TH1F("h_DecayVertexcutDistanceY_KFPart", "h_DecayVertexcutDistanceY_KFPart", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexcutDistanceY_KFPart.store);
      h_DecayVertexcutDistanceZ_KFPart.emplace_back(new TH1F("h_DecayVertexcutDistanceZ_KFPart", "h_DecayVertexcutDistanceZ_KFPart", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexcutDistanceZ_KFPart.store);
*/

      h_DecayVertexcutDistance_KFPart_PrimVtx.emplace_back(new TH1F("h_DecayVertexcutDistance_KFPart_PrimVtx", "h_DecayVertexcutDistance_KFPart_PrimVtx", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexcutDistance_KFPart_PrimVtx.store);
      h_DecayVertexcutDistanceX_KFPart_PrimVtx.emplace_back(new TH1F("h_DecayVertexcutDistanceX_KFPart_PrimVtx", "h_DecayVertexcutDistanceX_KFPart_PrimVtx", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexcutDistanceX_KFPart_PrimVtx.store);
      h_DecayVertexcutDistanceY_KFPart_PrimVtx.emplace_back(new TH1F("h_DecayVertexcutDistanceY_KFPart_PrimVtx", "h_DecayVertexcutDistanceY_KFPart_PrimVtx", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexcutDistanceY_KFPart_PrimVtx.store);
      h_DecayVertexcutDistanceZ_KFPart_PrimVtx.emplace_back(new TH1F("h_DecayVertexcutDistanceZ_KFPart_PrimVtx", "h_DecayVertexcutDistanceZ_KFPart_PrimVtx", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexcutDistanceZ_KFPart_PrimVtx.store);

      h_DecayVertexPosZ_real.emplace_back(new TH1F("h_DecayVertexPosZ_real", "h_DecayVertexPosZ_real", 300, 20, 50));
      HistReg.emplace_back(&h_DecayVertexPosZ_real.store);
      h_DecayVertexPosZ_vfunction.emplace_back(new TH1F("h_DecayVertexPosZ_vfunction", "h_DecayVertexPosZ_vfunction", 300, 20, 50));
      HistReg.emplace_back(&h_DecayVertexPosZ_vfunction.store);
      h_DecayVertexPosZ_centroid.emplace_back(new TH1F("h_DecayVertexPosZ_centroid", "h_DecayVertexPosZ_centroid", 300, 20, 50));
      HistReg.emplace_back(&h_DecayVertexPosZ_centroid.store);
      h_DecayVertexPosZ_KFPart.emplace_back(new TH1F("h_DecayVertexPosZ_KFPart", "h_DecayVertexPosZ_KFPart", 300, 20, 50));
      HistReg.emplace_back(&h_DecayVertexPosZ_KFPart.store);
      h_DecayVertexPosZ_AllVfunc.emplace_back(new TH1F("h_DecayVertexPosZ_AllVfunc", "h_DecayVertexPosZ_AllVfunc", 300, 20, 50));
      HistReg.emplace_back(&h_DecayVertexPosZ_AllVfunc.store);
      h_DecayVertexPosZ_AllCentroid.emplace_back(new TH1F("h_DecayVertexPosZ_AllCentroid", "h_DecayVertexPosZ_AllCentroid", 300, 20, 50));
      HistReg.emplace_back(&h_DecayVertexPosZ_AllCentroid.store);
      h_DecayVertexPosZ_AllKFPart.emplace_back(new TH1F("h_DecayVertexPosZ_AllKFPart", "h_DecayVertexPosZ_AllKFPart", 300, 20, 50));
      HistReg.emplace_back(&h_DecayVertexPosZ_AllKFPart.store);

      h_N_MotherTracks.emplace_back(new TH2F("h_N_MotherTracks", "h_N_MotherTracks", 10, 0, 10, 3000, 2, 5));
      HistReg.emplace_back(&h_N_MotherTracks.store);
      h_Dist_DaughterTracks.emplace_back(new TH2F("h_Dist_DaughterTracks", "h_Dist_DaughterTracks", 1000, 0, 2, 3000, 2, 5));
      HistReg.emplace_back(&h_Dist_DaughterTracks.store);
      h_Angle_MotherFragment.emplace_back(new TH2F("h_Angle_MotherFragment", "h_Angle_MotherFragment", 1000, 0, 5, 3000, 2, 5));
      HistReg.emplace_back(&h_Angle_MotherFragment.store);
      h_Angle_MotherPion.emplace_back(new TH2F("h_Angle_MotherPion", "h_Angle_MotherPion", 1000, 0, 30, 3000, 2, 5));
      HistReg.emplace_back(&h_Angle_MotherPion.store);
      h_Chi2ndf_MotherTracks.emplace_back(new TH2F("h_Chi2ndf_MotherTracks", "h_Chi2ndf_MotherTracks", 2000, 0, 40, 3000, 2, 5));
      HistReg.emplace_back(&h_Chi2ndf_MotherTracks.store);
      h_Dist_MotherTrackPrimVtx.emplace_back(new TH2F("h_Dist_MotherTrackPrimVtx", "h_Dist_MotherTrackPrimVtx", 1000, 0, 0.2, 3000, 2, 5));
      HistReg.emplace_back(&h_Dist_MotherTrackPrimVtx.store);
      h_Theta_MotherTrackPrimVtx.emplace_back(new TH2F("h_Theta_MotherTrackPrimVtx", "h_Theta_MotherTrackPrimVtx", 1000, 0, 2, 3000, 2, 5));
      HistReg.emplace_back(&h_Theta_MotherTrackPrimVtx.store);
      h_DecayVertexPosZ_KFPart_PrimVtx.emplace_back(new TH2F("h_DecayVertexPosZ_KFPart_PrimVtx", "h_DecayVertexPosZ_KFPart_PrimVtx", 1000, 20, 60, 3000, 2, 5));
      HistReg.emplace_back(&h_DecayVertexPosZ_KFPart_PrimVtx.store);
      h_Hyp_ArmenterosPodolanski.emplace_back(new TH2F("h_Hyp_ArmenterosPodolanski", "h_Hyp_ArmenterosPodolanski", 1000, 0, 1, 1000, 0, 1));
      HistReg.emplace_back(&h_Hyp_ArmenterosPodolanski.store);
      h_Hyp_CutArmenterosPodolanski.emplace_back(new TH2F("h_Hyp_CutArmenterosPodolanski", "h_Hyp_CutArmenterosPodolanski", 1000, 0, 1, 1000, 0, 1));
      HistReg.emplace_back(&h_Hyp_CutArmenterosPodolanski.store);

      h_HypInvariantMass.emplace_back(new TH1F("h_HypInvariantMass", "h_HypInvariantMass", 3000, 2, 5));
      HistReg.emplace_back(&h_HypInvariantMass.store);
      h_HypErrorInvariantMass.emplace_back(new TH1F("h_HypErrorInvariantMass", "h_HypErrorInvariantMass", 1000, 0, 0.5));
      HistReg.emplace_back(&h_HypErrorInvariantMass.store);
      

      h_Hyp_RealLifeTime.emplace_back(new TH1F("h_Hyp_RealLifeTime", "h_Hyp_RealLifeTime", 1000, 0, 1000));
      HistReg.emplace_back(&h_Hyp_RealLifeTime.store);
      h_HypLifeTime_PrimVtx.emplace_back(new TH1F("h_HypLifeTime_PrimVtx", "h_HypLifeTime_PrimVtx", 1000, 0, 1000));
      HistReg.emplace_back(&h_HypLifeTime_PrimVtx.store);
      h_HypErrorLifeTime_PrimVtx.emplace_back(new TH1F("h_HypErrorLifeTime_PrimVtx", "h_HypErrorLifeTime_PrimVtx", 2000, 0, 200));
      HistReg.emplace_back(&h_HypErrorLifeTime_PrimVtx.store);
      h_HypcutLifeTime_PrimVtx.emplace_back(new TH1F("h_HypcutLifeTime_PrimVtx", "h_HypcutLifeTime_PrimVtx", 1000, 0, 1000));
      HistReg.emplace_back(&h_HypcutLifeTime_PrimVtx.store);


      h_HypInvariantMassCheck.emplace_back(new TH2F("h_HypInvariantMassCheck", "h_HypInvariantMassCheck", 1000, -10, 10, 2, 0, 2));
      HistReg.emplace_back(&h_HypInvariantMassCheck.store);
      h_HypInvariantErrorMassCheck.emplace_back(new TH2F("h_HypInvariantErrorMassCheck", "h_HypInvariantErrorMassCheck", 1500, -300, 1200, 2, 0, 2));
      HistReg.emplace_back(&h_HypInvariantErrorMassCheck.store);


      h_HypInvariantMass_LorentzVect.emplace_back(new TH1F("h_HypInvariantMass_LorentzVect", "h_HypInvariantMass_LorentzVect", 3000, 2, 5));
      HistReg.emplace_back(&h_HypInvariantMass_LorentzVect.store);
      h_HypInvariantMass_CutLorentzVect.emplace_back(new TH1F("h_HypInvariantMass_CutLorentzVect", "h_HypInvariantMass_CutLorentzVect", 3000, 2, 5));
      HistReg.emplace_back(&h_HypInvariantMass_CutLorentzVect.store);


      h_EffPosZ_real.emplace_back(new TH1F("h_EffPosZ_real", "h_EffPosZ_real", 30, 20, 60));
      HistReg.emplace_back(&h_EffPosZ_real.store);
      h_EffPosZ_preKF.emplace_back(new TH1F("h_EffPosZ_preKF", "h_EffPosZ_preKF", 30, 20, 60));
      HistReg.emplace_back(&h_EffPosZ_preKF.store);
      h_EffPosZ_postKF.emplace_back(new TH1F("h_EffPosZ_postKF", "h_EffPosZ_postKF", 30, 20, 60));
      HistReg.emplace_back(&h_EffPosZ_postKF.store);
      h_EffPosZ_preKFPart.emplace_back(new TH1F("h_EffPosZ_preKFPart", "h_EffPosZ_preKFPart", 30, 20, 60));
      HistReg.emplace_back(&h_EffPosZ_preKFPart.store);
      h_EffPosZ_postKFPart.emplace_back(new TH1F("h_EffPosZ_postKFPart", "h_EffPosZ_postKFPart", 30, 20, 60));
      HistReg.emplace_back(&h_EffPosZ_postKFPart.store);
/*
      h_N_Si_MotherTracks.emplace_back(new TH1F("h_N_Si_MotherTracks", "h_N_Si_MotherTracks", 10, 0, 10));
      HistReg.emplace_back(&h_N_Si_MotherTracks.store);


      h_DecayVertexDistance_AllVfunc.emplace_back(new TH1F("h_DecayVertexDistance_AllVfunc", "h_DecayVertexDistance_AllVfunc", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexDistance_AllVfunc.store);
      h_DecayVertexDistanceX_AllVfunc.emplace_back(new TH1F("h_DecayVertexDistanceX_AllVfunc", "h_DecayVertexDistanceX_AllVfunc", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceX_AllVfunc.store);
      h_DecayVertexDistanceY_AllVfunc.emplace_back(new TH1F("h_DecayVertexDistanceY_AllVfunc", "h_DecayVertexDistanceY_AllVfunc", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceY_AllVfunc.store);
      h_DecayVertexDistanceZ_AllVfunc.emplace_back(new TH1F("h_DecayVertexDistanceZ_AllVfunc", "h_DecayVertexDistanceZ_AllVfunc", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexDistanceZ_AllVfunc.store);

      h_DecayVertexDistance_AllCentroid.emplace_back(new TH1F("h_DecayVertexDistance_AllCentroid", "h_DecayVertexDistance_AllCentroid", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexDistance_AllCentroid.store);
      h_DecayVertexDistanceX_AllCentroid.emplace_back(new TH1F("h_DecayVertexDistanceX_AllCentroid", "h_DecayVertexDistanceX_AllCentroid", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceX_AllCentroid.store);
      h_DecayVertexDistanceY_AllCentroid.emplace_back(new TH1F("h_DecayVertexDistanceY_AllCentroid", "h_DecayVertexDistanceY_AllCentroid", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceY_AllCentroid.store);
      h_DecayVertexDistanceZ_AllCentroid.emplace_back(new TH1F("h_DecayVertexDistanceZ_AllCentroid", "h_DecayVertexDistanceZ_AllCentroid", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexDistanceZ_AllCentroid.store);

      h_DecayVertexDistance_AllKFPart.emplace_back(new TH1F("h_DecayVertexDistance_AllKFPart", "h_DecayVertexDistance_AllKFPart", 4000, 0, 40));
      HistReg.emplace_back(&h_DecayVertexDistance_AllKFPart.store);
      h_DecayVertexDistanceX_AllKFPart.emplace_back(new TH1F("h_DecayVertexDistanceX_AllKFPart", "h_DecayVertexDistanceX_AllKFPart", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceX_AllKFPart.store);
      h_DecayVertexDistanceY_AllKFPart.emplace_back(new TH1F("h_DecayVertexDistanceY_AllKFPart", "h_DecayVertexDistanceY_AllKFPart", 1000, -5, 5));
      HistReg.emplace_back(&h_DecayVertexDistanceY_AllKFPart.store);
      h_DecayVertexDistanceZ_AllKFPart.emplace_back(new TH1F("h_DecayVertexDistanceZ_AllKFPart", "h_DecayVertexDistanceZ_AllKFPart", 2000, -10, 10));
      HistReg.emplace_back(&h_DecayVertexDistanceZ_AllKFPart.store);
*/
      h_DecayVtxstats.emplace_back(new TH1F("h_DecayVtxstats", "h_DecayVtxstats", 10, 0, 10));
      HistReg.emplace_back(&h_DecayVtxstats.store);
       
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
