#ifndef TGNNFINDER
#define TGNNFINDER

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"

#include "HitAna/ConstantParameter.hh"
#include "HitAna/FiberHitAna.hh"
#include "HitAna/PSBHitAna.hh"
#include "HitAna/PSFEHitAna.hh"
#include "HitAna/FiberTrackAna.hh"
#include "HitAna/FiberAnalyzer.hh"
#include "HitAna/TrackHit.hh"

#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "THyphiAttributes.h"

#include "TVector3.h"
#include "TRandom3.h"

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent, Out>;


template<class Out>
class TGNNFinder final : public TDataProcessInterface<Out> //TDataProcess<FullRecoEvent, Out> //TDataProcessInterface
{
public:
  const THyphiAttributes& att;

  TGNNFinder(const THyphiAttributes& attr);
  ~TGNNFinder();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, Out* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderGNN(FullRecoEvent& RecoEvent);


};

#endif
