#ifndef TDATABUILDWASACALIBPL
#define TDATABUILDWASACALIBPL

#include "TDataBuilder.h"

#include "THyphiAttributes.h"
#include "TRandom3.h"

#include "FullRecoEvent.hh"

#include "HitAna/ConstantParameter.hh"
#include "HitAna/FiberAnalyzer.hh"
#include "HitAna/FiberHitAna.hh"
#include "HitAna/FiberHitXUV.hh"
#include "HitAna/FiberTrackAna.hh"
#include "HitAna/MDCHitAna.hh"
#include "HitAna/ParaManager.hh"
#include "HitAna/PSBEHitAna.hh"
#include "HitAna/PSBHitAna.hh"
#include "HitAna/PSFEHitAna.hh"
#include "HitAna/T0HitAna.hh"

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
  ReturnRes::InfoM operator() (const EventWASAUnpack& event, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);

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
    TH1I* h_Builderstats;

    TH2D* h10[7][3];
    TH1D* h11[7][3];
    TH1D* h12[7][3];
    TH1D* h13[7][3];
    TH1D* h14[7][3];
    TH1D* h15[7][3];
    TH2D* h16[7];
    TH1D* h17[7];
    TH1D* h17_2[7];

    TH1D* hpsb_0_1;
    TH1D* hpsb_0_2;
    TH1D* hpsb_0_3;
    TH2D* hpsb_0_4;
    TH1D* hpsb_1_1;
    TH1D* hpsb_2[46];
    TH2D* hpsb_3[46];
    TH1D* hpsb_4[46];

    TH1D* hpsfe_0_1;
    TH1D* hpsfe_0_2;
    TH1D* hpsfe_0_3;
    TH2D* hpsfe_0_4;

    TH1D* hpsbe_0_1;
    TH1D* hpsbe_0_2;
    TH1D* hpsbe_0_3;
    TH2D* hpsbe_0_4;
    TH2D* hpsbe_1_0;

    TH1D* ht0_0_1;
    TH1D* ht0_0_2;
    TH1D* ht0_0_3;
    TH1D* ht0_0_4;
    TH1D* ht0_1[28];

    TH1D* hmdc_0_1;
    TH2D* hmdc_0_2;
    TH1D* hmdc_0_3;
    TH1D* hmdc_0_4;
    TH1D* hmdc_0_5;
    TH2D* hmdc_0_6;
    TH1D* hmdc_1[17];
    TH1D* hmdc_2[17];
    TH1D* hmdc_2_2[17];
    TH1D* hmdc_2_3[17];
    TH1D* hmdc_3[17];
    TH1D* hmdc_3_2[17];
    TH1D* hmdc_3_3[17];
  };

  LocalHists LocalHisto;

public:

  std::unique_ptr<ParaManager> par;

};

#endif
