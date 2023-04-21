#ifndef TWASAFINDER
#define TWASAFINDER

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

//typedef TDataProcess<FullRecoEvent, MCAnaEventG4Sol> TDataProcessInterface;
// template<class Out>
// struct TDataProcessInterfaceImp { using type = TDataProcess<FullRecoEvent, Out>; } ;

// template<class Out>
// using TDataProcessInterface = typename TDataProcessInterfaceImp<Out>::type;

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent, Out>;


template<class Out>
class TWASAFinder final : public TDataProcessInterface<Out> //TDataProcess<FullRecoEvent, Out> //TDataProcessInterface
{
public:
  const THyphiAttributes& att;

  TWASAFinder(const THyphiAttributes& attr);
  ~TWASAFinder();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, Out* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderWASA(FullRecoEvent& RecoEvent);

  double GetPSB_R(int _seg);
  double GetPSB_Phi(int _seg);
  //double GetFirstMFT_Z(TrackHit* Track);
  void CloseDist(FragmentTrack FragTrack, TrackHit* WASATrack, double& distance, TVector3& centroid);


  struct LocalHists
  {
    TH2D* h23_1;
    TH2D* h23_2;
    TH2D* h24_1;
    TH2D* h24_2;
    TH2D* h24_3[17];
    TH2D* h24_4[17];
    TH1D* h24_2_1;
    TH1D* h24_2_2;
    TH1D* h24_2_3[17];
    TH1D* h24_2_4[17];
    TH1D* hmdc_2_2[17];
    TH1D* hmdc_3_2[17];
  };
  
  LocalHists LocalHisto;
};

#endif
