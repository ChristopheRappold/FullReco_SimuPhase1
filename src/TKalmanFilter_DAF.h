#ifndef TKALMANFILTERDAF
#define TKALMANFILTERDAF

#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "TDataProcess.h"
//#include "Ana_EventNew.hh"
//#include "Ana_Event/Ana_EventNew_v10.hh"
//#include "Ana_Event/Ana_EventNew_v16.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "THyphiAttributes.h"
//#include "Genfit/HypLSLTrackFieldMapRep_v2.h"
//#include "Genfit_v3/core/GFTools.h"
//#include "Genfit_v3/core/GFDaf.h"
//#include "Genfit_v3/core/GFTrack.h"
//#include "Genfit_v3/RKTrackRep/RKTrackRep.h"

//#include "GFTools.h"
#include "AbsKalmanFitter.h"
#include "DAF.h"
#include "KalmanFitter.h"
#include "KalmanFitterRefTrack.h"

#include "EventDisplay.h"
#include "RKTrackRep.h"
#include "Track.h"

//#include "GFRaveVertexFactory.h"

//#include "Math/ProbFunc.h"

using TDataProcessInterface = TDataProcess<FullRecoEvent, MCAnaEventG4Sol>;
// namespace MathKalman
// {
//   struct Prob {
//     double operator() (double chi2,double ndf)
//     {
//       if (ndf <= 0)
//  return 0; // Set CL to zero in case ndf<=0

//       if (chi2 <= 0)
//  {
//    if (chi2 < 0)
//      return 0;
//    else
//      return 1;
//  }

//       return ::ROOT::Math::chisquared_cdf_c(chi2,ndf);
//     }
//   };

//}

class TKalmanFilter_DAF : public TDataProcessInterface
{
  public:
    const THyphiAttributes& att;

    explicit TKalmanFilter_DAF(const THyphiAttributes& attr);
    ~TKalmanFilter_DAF() final;

  // int Init(Ana_Hist* h);
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) final;
  void InitMT() final;
  private:
  int Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) final;
  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

  int Kalman_Filter_FromTrack(FullRecoEvent& RecoEvent);

    // genfit::DAF* Fitter;
    genfit::AbsKalmanFitter* Fitter;
    genfit::AbsKalmanFitter* Fitter_rescue;
    genfit::Track* Vtracks;
    genfit::RKTrackRep* rep;
    genfit::EventDisplay* display;
    //genfit::GFRaveVertexFactory* vertexFactory;

    std::vector<genfit::DetPlane*> list_Plane;

    TVector3 Plane_time;

  int Nb_event = 0;
  struct LocalHists {
    TH1I* h_stats;
    TH2I* h_statsLess3Mes;
    TH2I* h_statsInvalid;

    TH1F* h_pv;
    TH1F* h_chi2;
    TH1F* hd_chi[2];
    TH1F* hd_pv[2];

    TH1F* h_Path;
    TH1F* h_MeanPath;

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
    TH1D* h_ResPull[5][10];
    TH2F* h_ResPull_normal[5][10];

    TH2F* h_total_dE;

  };
  LocalHists LocalHisto;
};

#endif
