#ifndef TRPHIZTRACKMDC
#define TRPHIZTRACKMDC

#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "THyphiAttributes.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "TFile.h"
//#include "TKalmanFilter_DAF_ZMQ.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH2I.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

#include "TMultiGraph.h"

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent,Out>;

namespace MDC_PRZ {

  int f_checkDiscontinuity(TGraphErrors* g, double limit);
  std::tuple<int,double,double> f_calInterset(std::array<double,4>& x, std::array<double,4> y);
  void f_MeanCov(const std::vector<std::tuple<double,double>>& ResZPR, TVectorD& mean, TMatrixD& cov);
  std::tuple<int,double,double,double,double,double,double,double,double> f_interset(const TList* listPZ,const TList* listRZ, int id_1, int id_2);
  std::tuple<int,double,double,double,double,double,double> f_intersetTrack(const TList* listRZ, int id_1, double p0, double p1, const TMatrixD& covP01);
};

template<class Out>
class TRPhiZTrackMDC final :  public TDataProcessInterface<Out>
{
  public :
  const THyphiAttributes& att;

  TRPhiZTrackMDC(const THyphiAttributes& attr);
  ~TRPhiZTrackMDC();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,Out* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;
  int CheckTrackFinding(FullRecoEvent& RecoEvent);

  bool OutputEvents = false;

  TString namefilePhiZ;
  TFile* f_phiZ = nullptr;
  TTree* t_phiZ = nullptr;

  TMultiGraph* mg_trackPhiZ;
  TMultiGraph* mg_trackRZ;

  bool RZfit = false;
  int MDCWireType = 0;

  int tempEvent;
  struct LocalHists
  {
    TH2F* h_ResidualMDC_dZ1;
    TH2F* h_ResidualMDC_dZ2;
    TH2F* h_RPhiZMDC_Sigma;
    TH2F* h_ResidualMDC_dZ_PSB;
    TH2F* h_ResidualMDC_dZ_PSBE;
    TH2F* h_ResidualMDC_dZ_PSFE;
    TH2F* h_ResidualMDC_dZ_More6;

  };
  LocalHists LocalHisto;
};


#endif
