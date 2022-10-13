#ifndef TRIEMANNFINDER
#define TRIEMANNFINDER

#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "THyphiAttributes.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TH2I.h"

#include "tricktrack/FKDPoint.h"
#include "tricktrack/FKDTree.h"
#include "tricktrack/RiemannFit.h"


//typedef TDataProcess<FullRecoEvent,MCAnaEventG4Sol> TDataProcessInterface;
template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent,Out>;

using RPhiHit = tricktrack::FKDPoint<double, 3>;

struct HTrack {

  std::vector<double> Par = std::vector<double>(5);
  TMatrixD Cov = TMatrixD(5,5);
  int q = 0;
  double chi2_circle = -1.;
  double chi2_line = -1.;

};

struct RTrack {
  std::set<int> hits;
  std::vector<int> sortedHits;
  double x0;
  double y0;
  double r0;
  double chi2;
  double minPhi;
  double maxPhi;
  int q;
  bool toRefit = false;
  HTrack helix;
  RTrack(const std::set<int>& h, double x,double y,double r,double c2, double minP, double maxP, int Q_):hits(h),x0(x),y0(y),r0(r),chi2(c2),minPhi(minP),maxPhi(maxP),q(Q_) {};
};


template<class Out>
class TRiemannFinder final :  public TDataProcessInterface<Out>
{
  public :
  const THyphiAttributes& att;

  TRiemannFinder(const THyphiAttributes& attr);
  ~TRiemannFinder();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,Out* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int BuildKDTree(const FullRecoEvent& RecoEvent, tricktrack::FKDTree<RPhiHit, double, 3>& KDtree, std::vector<RPhiHit>& TempHits);
  void BuildTrackCand(const FullRecoEvent& RecoEvent, tricktrack::FKDTree<RPhiHit, double, 3>& KDtree, const std::vector<RPhiHit>& TempHits, std::vector<RTrack>& newTracksCand);
  void AddStereoWire(const FullRecoEvent& RecoEvent, std::vector<RTrack>& newTracksCand);
  void AddEndCap(const FullRecoEvent& RecoEvent, std::vector<RTrack>& newTracksCand);
  void AddOtherMDC(const FullRecoEvent& RecoEvent, std::vector<RTrack>& newTracksCand);
  void AddOtherPSCE(const FullRecoEvent& RecoEvent, std::vector<RTrack>& newTracksCand);
  void SortAndPosZ(std::vector<RTrack>& newTracksCand);
  void FitterRiemann(std::vector<RTrack>& newTracksCand);

  int FinderTrack(FullRecoEvent& RecoEvent,  std::vector<RTrack>& newTracksCand);

  std::vector<std::tuple<int, int, int> > TempHitToAllHits;
  std::vector<TVector3> TempHitXYZ;
  std::vector<std::array<double,9>> TempCovXYZ;

  std::set<int> redo_layer;
  std::unordered_map<int, std::unordered_map<int,int>> redo_layerToTrackID;

  std::vector<std::tuple<int, int, int> > TempHitPSBE;

  tricktrack::helix_fit res_h;

  bool OutputEvents = false;

  TString namefileFinder;
  TFile* f_finder = nullptr;
  TTree* t_finder = nullptr;

  TH2I* h_XY;
  TH2I* h_XYSim;
  TH2I* h_XYReco;
  TH2I* h_RPhi;
  TH2I* h_RPhiSim;
  TH2I* h_ZPhiSim;
  TH2I* h_XYConformal;
  TH2I* h_XYRiemann;
  TH2I* h_XYTracks;
  TH2I* h_XYRiemannTracks;

  constexpr double dl_max(int idDet) {
    switch(idDet){
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
      return 0.2;
      break;
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
      return 0.3;
      break;
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
      return 0.4;
      break;
    default:
      return 0.;
      break;
    }
  };

  constexpr int ch_max(int idDet) {
    int ch[17] = {52,64,76,88,100,76,86,96,106,116,126,104,112,120,130,138,148};
    return ch[idDet];
  };

  G4Sol::SolDet LastFrontWall;

  struct LocalHists
  {
    TH2F* h_RiemannChi2;
    TH2F* h_RiemannResidus;
    //TH2F* h_Stats;
    //TH2F* h_mom;

  };
  LocalHists LocalHisto;
};


#endif
