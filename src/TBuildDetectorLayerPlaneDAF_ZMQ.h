#ifndef TDATABUILDLPDAF_ZMQ
#define TDATABUILDLPDAF_ZMQ

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Debug.hh"
#include "EventG4Sol/TG4Sol_Event.hh"
#include "EventG4Sol/TG4Sol_Hit.hh"
#include "FullRecoEventZMQ.hh"
#include "TDataBuilder.h"
#include "THyphiAttributes.h"

//#include "MathematicalTools.hh"
#include "TFile.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "spdlog/spdlog.h"

#include <sstream>

//#include "TClustering.h"
#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif


class TBuildDetectorLayerPlaneDAF_ZMQ : public TDataBuilder
{

public:
  const THyphiAttributes& att;

  explicit TBuildDetectorLayerPlaneDAF_ZMQ(const THyphiAttributes& att);
  ~TBuildDetectorLayerPlaneDAF_ZMQ() final;

#ifdef ROOT6
  ReturnRes::InfoM operator()(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
                              ZMQ::DataBuilderOut& RecoEvent);
#else
  ReturnRes::InfoM operator()(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits,
                              ZMQ::DataBuilderOut& RecoEvent);
#endif

private:
  // int Exec(THyphiEvent_Phys_new *event,Ana_Event* OutTree);
#ifdef ROOT6
  int Exec(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits,
           ZMQ::DataBuilderOut& RecoEvent);
#else
  int Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, ZMQ::DataBuilderOut& RecoEvent);
#endif

  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

private:
  // std::vector<std::vector<GFAbsRecoHit*> > ListAllHits;
  std::unordered_map<int, int> orderDetectors;
  ZMQ::PDG_fromName pid_fromName;
  struct LocalHists
  {
    TH1I* h_stats;
  };
  LocalHists LocalHisto;
};

#endif
