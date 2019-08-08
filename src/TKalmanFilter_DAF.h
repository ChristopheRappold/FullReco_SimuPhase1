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

//#include "Math/ProbFunc.h"

using TDataProcessInterface = TDataProcess<FullRecoEvent, MCAnaEventG4Sol>;
// namespace MathKalman
// {
//   struct Prob {
//     double operator() (double chi2,double ndf)
//     {
//       if (ndf <= 0)
// 	return 0; // Set CL to zero in case ndf<=0

//       if (chi2 <= 0)
// 	{
// 	  if (chi2 < 0)
// 	    return 0;
// 	  else
// 	    return 1;
// 	}

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
  int operator()(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) final;

  private:
  int Exec(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) final;
  int SoftExit(int) final;
  int Kalman_Filter_FromTrack(FullRecoEvent& RecoEvent);

  // genfit::DAF* Fitter;
  genfit::AbsKalmanFitter* Fitter;
  genfit::AbsKalmanFitter* Fitter_rescue;
  genfit::Track* Vtracks;
  genfit::RKTrackRep* rep;
  genfit::EventDisplay* display;

  std::vector<genfit::DetPlane*> list_Plane;

  TVector3 Plane_time;

  int Nb_event = 0;
};

#endif
