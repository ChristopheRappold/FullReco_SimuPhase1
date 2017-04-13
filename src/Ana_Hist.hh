#ifndef HYPHIRHISTO_H
#define HYPHIRHISTO_H

//#include "ShadowAna_Hist.hh"

#include <iostream>

#include <string>
#include <unordered_map>
#include <vector>

#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2I.h"

//#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMath.h"

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

class Ana_Hist
{

  public:
  bool HaveBeenWritten;
  std::vector<bool> EnableState;

  TH1I* h_stats;
  TH1I* h_task_exit;
  // Finder
  TH2F* h_xy;
  TH2F* h_PxPy;
  TH2F* h_xy_extrap;
  TH2F* h_PxPy_extrap;
  TH2F* h_TrackFindingStat;
  // Kalman:

  TH1F* h_pv;
  TH1F* h_chi2;
  TH1F* hd_pv[2];
  TH1F* hd_chi[2];

  TH1F* h_Path;
  TH1F* h_Path_Back;
  TH1F* h_MeanPath;
  TH1F* h_dpath;

  TH1F* h_beta;
  TH1F* h_beta2;
  TH1F* h_beta3;

  TH1F* h_Mass_All;
  TH1F* h_Mass_All2;
  TH1F* h_Mass_All3;

  TH2F* h_Mass_charge_All;
  TH2F* h_Mass_charge_All2;
  TH2F* h_Mass_charge_All3;

  TH2F* h_beta_mom;
  TH2F* h_beta_mom2;
  TH2F* h_beta_mom3;

  TH2F* h_pv_mom;
  TH2F* h_pv_beta;
  TH2F* h_pv_mass;

  TH2F* h_path_tof;

  TH2F* h_mom_tof_cut;
  TH2F* h_path_mom_cut;
  TH2F* h_path_tof_cut;

  TH1F* h_Mass[4];
  TH1F* h_chi2_particle[4];
  TH1F* h_pv_particle[4];

  TH2F* h_mom_res[5];

  
  std::unordered_map<std::string, std::vector<TH1*> > HistRegisteredByDir;

  Ana_Hist(bool Daf = true, bool Vertex = true, bool DCproject = true, bool Finding = true, bool Hough = true, bool Simu = false);
  ~Ana_Hist();

  int Write(TFile*);
  int WriteTemp(char* tempfile);
  TDirectory* GetDir(TFile* out_file, const std::string& name_dir)
  {
    if(!out_file->GetDirectory(name_dir.c_str()))
      return out_file->mkdir(name_dir.c_str());
    else
      return out_file->GetDirectory(name_dir.c_str());
  }
};

#endif
