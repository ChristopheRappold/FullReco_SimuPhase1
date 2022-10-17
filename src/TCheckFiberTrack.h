#ifndef TCHECKFIBERTRACK
#define TCHECKFIBERTRACK

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

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent,Out>;

template<class Out>
class TCheckFiberTrack final :  public TDataProcessInterface<Out>
{
  public :
  const THyphiAttributes& att;

  TCheckFiberTrack(const THyphiAttributes& attr);
  ~TCheckFiberTrack();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,Out* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;
  int CheckTrackFinding(const FullRecoEvent& RecoEvent);

  struct LocalHists
  {
    TH2F* h_ResidualFiberX;
    TH2F* h_ResidualFiberY;
  };
  LocalHists LocalHisto;
};


#endif
