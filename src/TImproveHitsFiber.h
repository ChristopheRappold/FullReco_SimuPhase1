#ifndef TIMPROVEHITSFIBER
#define TIMPROVEHITSFIBER

#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "THyphiAttributes.h"

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent,Out>;

namespace Improve {
  struct miniF
  {
    int det = -1;
    int hit = -1;
    std::array<TVector3,2> xyz;
    std::array<TVector3,2> xyz_low;
    double res = -1;
  };
}


template<class Out>
class TImproveHitsFiber final :  public TDataProcessInterface<Out>
{
  public :
  const THyphiAttributes& att;

  TImproveHitsFiber(const THyphiAttributes& attr);
  ~TImproveHitsFiber();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,Out* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;

  int ImproveFiber(FullRecoEvent& RecoEvent);

  //TVector3 ImprovePSFBE(const std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, int id_det, int id_hit, const ParamFitsRZ& fitRZ);
  //Improve::miniF ImproveFiber(const std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, int id_det, int id_hit, const ParamFitsRZ& fitRZ);
  Improve::miniF ImproveFiber2(const std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > >& ListHits, int id_det, int id_hit);
  std::tuple<double,double> Intersection(const Improve::miniF& fiber1, const Improve::miniF& fiber2);
  std::tuple<double,double> Inter(const TVector3& p1, const TVector3& d1, const TVector3& p2, const TVector3& d2);
  int LastFrontWall;

  struct LocalHists
  {
    TH2F* h_ImproveHitsResidus;
    TH2F* h_ImproveHitsResidus2;
  };
  LocalHists LocalHisto;
};


#endif
