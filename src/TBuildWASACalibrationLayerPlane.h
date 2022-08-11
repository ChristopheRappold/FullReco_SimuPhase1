#ifndef TDATABUILDWASACALIBPL
#define TDATABUILDWASACALIBPL

#include "TDataBuilder.h"

#include "THyphiAttributes.h"
#include "TRandom3.h"

#include "FullRecoEvent.hh"

#include "Ana_Event/Ana_WasaEvent.hh"
#include "EventWASAUnpack/WASAUnpackBranch.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"


#include "Debug.hh"

#include "TGeoManager.h"
#include "TGeoElement.h"

#include <sstream>

#include "spdlog/spdlog.h"

//#include "TClustering.h"
#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif




class TBuildWASACalibrationLayerPlane final : public TDataBuilder
{

public:

  const THyphiAttributes& att;

  explicit TBuildWASACalibrationLayerPlane(const THyphiAttributes& att);
  ~TBuildWASACalibrationLayerPlane() final;

#ifdef ROOT6
  ReturnRes::InfoM operator() (const EventWASAUnpack& event, FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree);
#else
  ReturnRes::InfoM operator() (const EventWASAUnpack& event, FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree);
#endif
  ReturnRes::InfoM operator() (const EventWASAUnpack& event, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree) {}

private :

#ifdef ROOT6
  int Exec(const EventWASAUnpack& event, FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree);
#else
  int Exec(const EventWASAUnpack& event, FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree);
#endif
  
  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

private:
  


  struct LocalHists
  {
    TH1I* h_stats;
  };
  LocalHists LocalHisto;
};

#endif
