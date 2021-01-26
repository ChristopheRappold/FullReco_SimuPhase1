#ifndef HYPHIRHISTO_H
#define HYPHIRHISTO_H

//#include "ShadowAna_Hist.hh"

//#include <iostream>

#include <string>
#include <unordered_map>
#include <vector>

#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TGraphErrors.h"
#include "TEllipse.h"

//#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMath.h"

#include "spdlog/logger.h"

enum StateHist : int
{
  DAF = 0,
  VERTEX,
  DCPROJ,
  FINDING,
  HOUGH,
  SIMU,
  SIZEOF_STATEHIST
};

template <typename T>
struct Hist {
  T* h;
  std::vector<TH1*> store;
  void emplace_back(T* _h) {
    h = _h;
    h->SetDirectory(0);
    store.emplace_back(_h);
  }
};

class Ana_Hist
{
  std::shared_ptr<spdlog::logger> _logger;
  
  public:
  bool HaveBeenWritten;
  std::vector<bool> EnableState;

  Hist<TH1I> h_stats;
  Hist<TH2I> h_statsLess3Mes;
  Hist<TH2I> h_statsInvalid;
  Hist<TH1I> h_task_exit;
  //Field
  Hist<TH2F> FieldXY[3];
  Hist<TH2F> FieldXZ[3];
  Hist<TH2F> FieldYZ[3];
  Hist<TH2F> FieldXY_n[3];
  Hist<TH2F> FieldXZ_n[3];
  Hist<TH2F> FieldYZ_n[3];
  // Finder
  Hist<TH2F> h_xy;
  Hist<TH2F> h_PxPy;
  Hist<TH2F> h_xy_extrap;
  Hist<TH2F> h_PxPy_extrap;
  Hist<TH2F> h_TrackFindingStat;
  Hist<TH2F> h_MDC_Dphi;
  Hist<TH2F> h_MDC_NextLayer;
  Hist<TH2F> h_MDC_DiffLayer;
  Hist<TH2F> h_MDC_DiffLayer2;
  Hist<TH2F> h_MDC_InLayer;
  
  Hist<TH2F> h_SolenoidGeo[3];
  std::vector<TEllipse*> geoSolenoid;

  Hist<TH1F> h_RZStats;
  Hist<TH2F> h_RZ;
  Hist<TH1F> h_RZfit_mom;
  Hist<TH2F> h_RZfit_Chi2;
  Hist<TH2F> h_XYfit_miniF;
  Hist<TH2F> h_MDC_Z_residu;
  Hist<TH2F> h_MDC_R_residu;
  Hist<TH2F> h_MDC_Z_pull;
  Hist<TH2F> h_MDC_R_pull;
  // K>lman:
  Hist<TH1F> h_pv;
  Hist<TH1F> h_chi2;
  Hist<TH1F> hd_pv[2];
  Hist<TH1F> hd_chi[2];
  
  Hist<TH1F> h_Path;
  Hist<TH1F> h_Path_Back;
  Hist<TH1F> h_MeanPath;
  Hist<TH1F> h_dpath;
  
  Hist<TH1F> h_beta;
  Hist<TH1F> h_beta2;
  Hist<TH1F> h_beta3;
  
  Hist<TH1F> h_Mass_All;
  Hist<TH1F> h_Mass_All2;
  Hist<TH1F> h_Mass_All3;
  
  Hist<TH2F> h_Mass_charge_All;
  Hist<TH2F> h_Mass_charge_All2;
  Hist<TH2F> h_Mass_charge_All3;
  
  Hist<TH2F> h_beta_mom;
  Hist<TH2F> h_beta_mom2;
  Hist<TH2F> h_beta_mom3;
  
  Hist<TH2F> h_pv_mom;
  Hist<TH2F> h_pv_beta;
  Hist<TH2F> h_pv_mass;
  
  Hist<TH2F> h_path_tof;
  
  Hist<TH2F> h_mom_tof_cut;
  Hist<TH2F> h_path_mom_cut;
  Hist<TH2F> h_path_tof_cut;
  
  Hist<TH1F> h_Mass[4];
  Hist<TH1F> h_chi2_particle[4];
  Hist<TH1F> h_pv_particle[4];
  
  Hist<TH2F> h_mom_res[5];
  Hist<TH1D> h_ResPull[5][10];
  Hist<TH2F> h_ResPull_normal[5][10];

  Hist<TH2F> h_total_dE;

  //residual
  Hist<TH1F> h_ResFiber[9];
  Hist<TH1F> h_ResMiniFiber[6];
  Hist<TH1F> h_ResMDC[17][3];
  Hist<TH1F> h_ResPSCE[2];
  
  std::unordered_map<std::string, std::tuple<std::vector<std::vector<TH1*>*>, int> > HistRegisteredByDir;

  Ana_Hist(bool Daf = true, bool Vertex = true, bool DCproject = true, bool Finding = true, bool Hough = true, bool Simu = false);
  ~Ana_Hist();

  int Write(TFile*);
  int WriteTemp(char* tempfile);

  void DebugHists();
  
  template <typename T>
  T* CloneAndRegister(Hist<T>& hist) {
    T* newHist = dynamic_cast<T*>(hist.h->Clone());
    hist.store.emplace_back(newHist);
    newHist->SetDirectory(0);
    return newHist;
  }

  TDirectory* GetDir(TFile* out_file, const std::string& name_dir)
  {
    if(!out_file->GetDirectory(name_dir.c_str()))
      return out_file->mkdir(name_dir.c_str());
    else
      return out_file->GetDirectory(name_dir.c_str());
  }
  
};

#endif
