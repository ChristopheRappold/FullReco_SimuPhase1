#ifndef TMERGEROUTPUT_MT
#define TMERGEROUTPUT_MT

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "TDataMerger.h"
#include "THyphiAttributes.h"

#include "TRandom3.h"

#include "spdlog/spdlog.h"

#include <sstream>

class TMergerOutput_MT final : public TDataMerger
{
public:
  const THyphiAttributes& att;

  explicit TMergerOutput_MT(const THyphiAttributes& att);
  ~TMergerOutput_MT() final;

  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, ReturnRes::InfoM Status, MCAnaEventG4Sol* OutTree);

private:
  int Exec(FullRecoEvent& RecoEvent, ReturnRes::InfoM Status, MCAnaEventG4Sol* OutTree);

  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

  void PushToTreeParticles(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
  void PushToTreeHits(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
  void PushToTreeTracks(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);

  struct LocalHists {
    TH1I* h_stats;
  };
  LocalHists LocalHisto;
};

#endif
