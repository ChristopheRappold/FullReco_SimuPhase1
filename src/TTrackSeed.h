#ifndef TTRACKSEED
#define TTRACKSEED

#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "THyphiAttributes.h"

#include "tricktrack/FKDPoint.h"
#include "tricktrack/RiemannFit.h"


//typedef TDataProcess<FullRecoEvent,MCAnaEventG4Sol> TDataProcessInterface;
template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent,Out>;

namespace Seed {
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
}

template<class Out>
class TTrackSeed final :  public TDataProcessInterface<Out>
{
  public :
  const THyphiAttributes& att;

  TTrackSeed(const THyphiAttributes& attr);
  ~TTrackSeed();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,Out* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  void SortAndPosZ(std::vector<Seed::RTrack>& newTracksCand);
  void FitterRiemann(std::vector<Seed::RTrack>& newTracksCand);

  int TrackFit(FullRecoEvent& RecoEvent,  std::vector<Seed::RTrack>& newTracksCand);

  std::tuple<int, double, double, double> extract_WireXY(const std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, int id_det, int id_hit, bool DoCov = false, TMatrixD* cov = nullptr, double* dl = nullptr);
  std::tuple<double,double,double> extract_PSBXY(const  std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, int id_det, int id_hit, bool DoCov = false, TMatrixD* cov = nullptr);

  constexpr int ch_max(int idDet) {
    int ch[17] = {52,64,76,88,100,76,86,96,106,116,126,104,112,120,130,138,148};
    return ch[idDet];
  };

  std::vector<std::tuple<int, int, int> > TempHitToAllHits;
  std::vector<TVector3> TempHitXYZ;
  std::vector<std::array<double,9>> TempCovXYZ;

  tricktrack::helix_fit res_h;

  int LastFrontWall;

  struct LocalHists
  {
    TH2F* h_SeedRiemannChi2;
    TH2F* h_SeedRiemannResidus;
    //TH2F* h_Stats;
    //TH2F* h_mom;

  };
  LocalHists LocalHisto;
};


#endif
