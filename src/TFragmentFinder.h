#ifndef TFRAGMENTFINDER
#define TFRAGMENTFINDER

#include "Ana_Event/MCAnaEventG4Sol.hh"
#include "Ana_Event/Ana_WasaEvent.hh"

#include "HitAna/OpticsMom.hh"

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
class TFragmentFinder final : public TDataProcessInterface<Out> //TDataProcess<FullRecoEvent, Out> //TDataProcessInterface
{
public:
  const THyphiAttributes& att;

  TFragmentFinder(const THyphiAttributes& attr);
  ~TFragmentFinder();

  // int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator()(FullRecoEvent& RecoEvent, Out* OutTree) override;

private:
  int Exec(FullRecoEvent& RecoEvent, Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int FinderFragment(FullRecoEvent& RecoEvent);







  struct LocalHists
  {
    TH1D* hopt_1_1;
    TH1D* hopt_1_2;
    TH1D* hopt_1_3;
    TH1D* hopt_1_4;
    TH2D* hopt_2_1;
    TH2D* hopt_2_2;
    TH1D* hopt_2_3;
    TH1D* hopt_2_4;
  };
  LocalHists LocalHisto;
};

#endif
