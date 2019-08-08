#ifndef TKALMANFILTERFRS
#define TKALMANFILTERFRS
#include "TDataProcess.h"
#include "FullRecoEvent.hh"
//#include "Ana_EventNew.hh"
//#include "Ana_Event/Ana_EventNew_v10.hh"
#include "Ana_Event/MCAnaEvent.hh"
#include "THyphiAttributes.h"
//#include "Genfit/HypLSLTrackFieldMapRep_v2.h"
//#include "Genfit_v3/core/GFTools.h"
//#include "Genfit_v3/core/GFDaf.h"
//#include "Genfit_v3/core/GFTrack.h"
//#include "Genfit_v3/RKTrackRep/RKTrackRep.h"

#include "GFTools.h"
#include "GFDaf.h"
#include "GFTrack.h"
#include "RKTrackRep.h"

#include "Math/ProbFunc.h"

using TDataProcessInterfaceMC = TDataProcess<FullRecoEvent,MCAnaEvent> ;
namespace MathKalmanFRS
{
  struct Prob {
    double operator() (double chi2,double ndf)
    {
      if (ndf <= 0) 
	return 0; // Set CL to zero in case ndf<=0
      
      if (chi2 <= 0) 
	{
	  if (chi2 < 0) 
	    return 0;
	  else
	    return 1;
	}
      
      return ::ROOT::Math::chisquared_cdf_c(chi2,ndf); 
    }
  };

}
class TKalmanFilter_FRS :  public TDataProcessInterfaceMC
{
  public :
  const THyphiAttributes& att;


  TKalmanFilter_FRS(const THyphiAttributes& attr);
  ~TKalmanFilter_FRS();

  //int Init(Ana_Hist* h);
  int operator() (FullRecoEvent& RecoEvent,MCAnaEvent* OutTree);
 private:
  int Exec(FullRecoEvent& RecoEvent,MCAnaEvent* OutTree);
  int SoftExit(int);
  int Kalman_Filter_FromTrack(FullRecoEvent& RecoEvent);

  GFDaf Daf;
  GFKalman Kalman;
  GFTrack* Vtracks;
  RKTrackRep* rep;
  RKTrackRep* rep_length;

  std::vector<GFDetPlane*> list_Plane;

  TVector3 Plane_time;

};


#endif
