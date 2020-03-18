#ifndef TMERGEROUTPUT_ZMQ
#define TMERGEROUTPUT_ZMQ

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Debug.hh"
#include "FullRecoEventZMQ.hh"
#include "TDataMerger.h"
#include "THyphiAttributes.h"

#include "TRandom3.h"

#include "spdlog/spdlog.h"

#include <sstream>

class TMergerOutput_ZMQ : public TDataMerger
{
public:
  const THyphiAttributes& att;

  explicit TMergerOutput_ZMQ(const THyphiAttributes& att);
  ~TMergerOutput_ZMQ() final;

  ReturnRes::InfoM operator()(ZMQ::DataFitterOut& In, MCAnaEventG4Sol* OutTree);

private:
  int Exec(ZMQ::DataFitterOut& In, MCAnaEventG4Sol* OutTree);

  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

  void PushToTreeParticles(ZMQ::DataFitterOut& RecoEvent, MCAnaEventG4Sol* OutTree);
  void PushToTreeHits(ZMQ::DataFitterOut& RecoEvent, MCAnaEventG4Sol* OutTree);
  void PushToTreeTracks(ZMQ::DataFitterOut& RecoEvent, MCAnaEventG4Sol* OutTree);

  struct LocalHists {
    TH1I* h_stats;
  };
  LocalHists LocalHisto;
};

#endif
