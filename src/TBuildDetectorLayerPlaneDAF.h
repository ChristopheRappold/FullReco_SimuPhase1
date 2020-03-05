#ifndef TDATABUILDLPDAF
#define TDATABUILDLPDAF

#include "TDataBuilder.h"

#include "THyphiAttributes.h"
#include "TRandom3.h"

#include "FullRecoEvent.hh"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "EventG4Sol/TG4Sol_Event.hh"
#include "EventG4Sol/TG4Sol_Hit.hh"

//#include "MathematicalTools.hh"
#include "Debug.hh"

#include "TProfile.h"
#include "TFile.h"

#include "TGeoManager.h"
#include "TGeoElement.h"

#include <sstream>

#include "spdlog/spdlog.h"

//#include "TClustering.h"
#ifdef ROOT6
#include "TTreeReaderArray.h"
#endif


constexpr bool IsPlanar(G4Sol::SolDet idDet) {
  switch (idDet) {
  case G4Sol::InSi0 : ;
  case G4Sol::InSi1 : ;
  case G4Sol::InSi2 : ;
  case G4Sol::InSi3 : ;
  case G4Sol::TR1 : ;
  case G4Sol::TR2 : ;
  case G4Sol::PSFE : ;
  case G4Sol::PSBE : ;
  case G4Sol::TrFwd0 : ;
  case G4Sol::TrFwd1 : ;
  case G4Sol::TrFwd2 : ;
  case G4Sol::RPC_l : ;
  case G4Sol::RPC_h : ;
  case G4Sol::FMF2Stop0 : ;
  case G4Sol::FMF2Stop1 : ;
  case G4Sol::FMF2Stop2 :
    return true;
  default:
    return false;
  };  
}; 



class TBuildDetectorLayerPlaneDAF : public TDataBuilder
{

 public:
  
  const THyphiAttributes& att;

  explicit TBuildDetectorLayerPlaneDAF(const THyphiAttributes& att);
  ~TBuildDetectorLayerPlaneDAF() final;

#ifdef ROOT6
  ReturnRes::InfoM operator() (const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#else
  ReturnRes::InfoM operator() (const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#endif
  
  private :
  //int Exec(THyphiEvent_Phys_new *event,Ana_Event* OutTree);
#ifdef ROOT6
  int Exec(const TG4Sol_Event& event, const std::vector<TTreeReaderArray<TG4Sol_Hit>*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#else  
  int Exec(const TG4Sol_Event& event, const std::vector<TClonesArray*>& hits, FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree);
#endif
  
  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

private:
  
  //std::vector<std::vector<GFAbsRecoHit*> > ListAllHits;
  std::unordered_map<int,int> orderDetectors;
  PDG_fromName pid_fromName;
  struct LocalHists
  {
    TH1I* h_stats;
  };
  LocalHists LocalHisto;
};

#endif
