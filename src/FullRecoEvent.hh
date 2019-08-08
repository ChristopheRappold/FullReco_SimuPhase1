#ifndef FullRecoEvent_h
#define FullRecoEvent_h

#include <iterator>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"

//#include "Genfit/Track.h"
//#include "Genfit/HypRecoHit.h"
#include "AbsMeasurement.h"
#include "PlanarMeasurement.h"
#include "SpacepointMeasurement.h"
#include "ProlateSpacepointMeasurement.h"

#include <cassert>

class MomRef;
class CompMomRef;

using map_mom3 = std::map<MomRef, std::vector<double>, CompMomRef> ;

namespace G4Sol
{
  enum SolDet : int
    {
      InSi0 = 0,
      InSi1,
      InSi2,
      InSi3 /*3*/,
      TR1  /*4*/,
      TR2  /*5*/,
      PSFE /*6*/,
      MG01 /*7*/,
      MG02,
      MG03,
      MG04,
      MG05,
      MG06,
      MG07,
      MG08,
      MG09,
      MG10,
      MG11,
      MG12,
      MG13,
      MG14,
      MG15,
      MG16,
      MG17 /*23*/,
      PSCE /*24*/,
      PSBE /*25*/,
      CDC_layer0 /*4 - 26*/,
      CDC_layer1,
      CDC_layer2,
      CDC_layer3,
      CDC_layer4,
      CDC_layer5,
      CDC_layer6,
      CDC_layer7,
      CDC_layer8,
      CDC_layer9,
      CDC_layer10,
      CDC_layer11,
      CDC_layer12,
      CDC_layer13,
      CDC_layer14 /*18 - 40*/,
      CDHBar /*19 - 41*/,
      TrFwd0 /*20 - 42*/,
      TrFwd1,
      TrFwd2 /*22 - 44*/,
      RPC_l /*23 - 44*/,
      RPC_h /*24 - 45*/,
      FMF2Stop0 /*25 - 46*/,
      FMF2Stop1,
      FMF2Stop2 /*27 - 48*/,
      SIZEOF_G4SOLDETTYPE
    };
  
  constexpr auto nameLiteralDet = {  "InSi0", "InSi1", "InSi2", "InSi3", "TR1", "TR2", "PSFE", "MG01", "MG02", "MG03", "MG04",
				     "MG05", "MG06", "MG07", "MG08", "MG09", "MG10", "MG11", "MG12", "MG13", "MG14", "MG15",
				     "MG16", "MG17", "PSCE", "PSBE", "CDC0", "CDC1", "CDC2", "CDC3", "CDC4",
				     "CDC5", "CDC6", "CDC7", "CDC8", "CDC9", "CDC10", "CDC11", "CDC12",
				     "CDC13", "CDC14", "CDHBar", "TrFwd0", "TrFwd1", "TrFwd2", "RPCl", "RPCh",
				     "FMF2S0", "FMF2S1", "FMF2S2"};
  

  
  template <typename T, T... args>
  struct EnumIter : public std::iterator<std::input_iterator_tag, T, std::ptrdiff_t, const T*, const T&>
  {
    static constexpr T values[] = {args...};
    static constexpr std::size_t size = sizeof...(args);

    int pos;
    EnumIter() // No value is end
      : pos(size)
    {    }
    EnumIter(T val) : pos(std::distance(&values[0], std::find(&values[0], &values[size], val))) {}

    const T& operator*() const { return values[pos]; }
    EnumIter& operator++()
    {
      ++pos;
      return *this;
    }
    EnumIter operator++(int)
    {
      EnumIter r(*this);
      this->operator++();
      return r;
    }
    bool operator==(EnumIter const& rhs) { return pos == rhs.pos; }
    bool operator!=(EnumIter const& rhs) { return pos != rhs.pos; }
  };
  template <typename T, T... args>
  constexpr T EnumIter<T, args...>::values[];
  using SolDetIter = struct EnumIter<SolDet, InSi0, InSi1, InSi2, InSi3 /*3*/, TR1 /*4*/, TR2 /*5*/, PSFE /*6*/, MG01 /*7*/, MG02, MG03, MG04, MG05, MG06, MG07, MG08,
				     MG09, MG10, MG11, MG12, MG13, MG14, MG15, MG16, MG17 /*23*/, PSCE /*24*/, PSBE /*25*/, CDC_layer0 /*4 - 26*/, CDC_layer1,
				     CDC_layer2, CDC_layer3, CDC_layer4, CDC_layer5, CDC_layer6, CDC_layer7, CDC_layer8, CDC_layer9, CDC_layer10,
				     CDC_layer11, CDC_layer12, CDC_layer13, CDC_layer14 /*18 - 40*/, CDHBar /*19 - 41*/, TrFwd0 /*20 - 42*/, TrFwd1,
				     TrFwd2 /*22 - 43*/, RPC_l /*23 - 44*/, RPC_h /*24 - 45*/, FMF2Stop0 /*25 - 46*/, FMF2Stop1, FMF2Stop2> ;
  
}

class MomRef
{
  public:
  int id_track;

  double p_value;
  double mass;
  double x_TOF;
  double y_TOF;

  MomRef() : id_track(-1), p_value(-1.), mass(-1.), x_TOF(-1.), y_TOF(-1.){};
  MomRef(int id_, double pv_, double mass_, double x_, double y_) : id_track(id_), p_value(pv_), mass(mass_), x_TOF(x_), y_TOF(y_){};
  ~MomRef() = default;
};

class CompMomRef
{
  public:
  bool operator()(const MomRef& a, const MomRef& b) const
  {
    if(TMath::Abs(a.p_value - b.p_value) < 1e-7)
      return a.id_track < b.id_track;
    else
      return a.p_value < b.p_value;
  }
};

class CompMomRef_fromXtof
{
  public:
  bool operator()(const MomRef& a, const MomRef& b) const
  {
    if(TMath::Abs(a.x_TOF - b.x_TOF) < 1e-4)
      {
        if(TMath::Abs(a.p_value - b.p_value) < 1e-7)
          return a.id_track < b.id_track;
        else
          return a.p_value < b.p_value;
      }
    else
      return a.x_TOF < b.x_TOF;
  }
};

class TVertex
{
  public:
  std::string type;
  TVector3 vertex;
  double Chi2;
  int Ndf;
  TMatrixT<double> Cov;

  TVertex(std::string t_, const TVector3& v_, double chi2_, int ndf_, const TMatrixT<double>& c_)
    : type(std::move(t_)), vertex(v_), Chi2(chi2_), Ndf(ndf_), Cov(c_)
  {
  }
  explicit TVertex(std::string t_) : type(std::move(t_)), vertex(TVector3(0., 0., 0.)), Chi2(0.), Ndf(0), Cov(TMatrixT<double>(3, 3)) {}

  ~TVertex() = default;
};

struct SimHit
{
  int layerID = -1;
  double hitX = -999.;
  double hitY = -999.;
  double hitZ = -999.;

  double momX = -999.;
  double momY = -999.;
  double momZ = -999.;

  int pdg = 0;
  double mass = 0.;

  double Eloss = 0.;
  double time = 0.;
  double length = 0.;
};

struct InfoPar
{
  int pdg = 0;
  double momX = -999.;
  double momY = -999.;
  double momZ = -999.;
  double mass = -999.;
  double Eloss = 0.;
  double time = 0.;
  double length = 0.;
};

struct ResSolDAF
{
  int charge = 0;
  int pdg_ini = 0;
  int pdg_guess = 0;

  double pvalue = -1.;
  double chi2 = -1.;
  double ndf = 0.;
  int fitter = -1;
  
  double momX_init = -999.;
  double momY_init = -999.;
  double momZ_init = -999.;

  double momX = -999.;
  double momY = -999.;
  double momZ = -999.;

  double posX = -999.;
  double posY = -999.;
  double posZ = -999.;

  double mass = -999.;
  double mass2 = -999.;
  double beta = -999.;
  double beta2 = -999.;
  double tof = -999.;
  double tof2 = -999.;

  double path_length = -999.;
  double path_length2 = -999.;
  double path_time = -999.;

  double dE = -999.;

  int firstHit = -1;
  int lastHit = -1;
  int Ncentral = 0;
  
  double cov_matrix[6][6] = {{-999., -999., -999., -999., -999., -999.}, {-999., -999., -999., -999., -999., -999.},
                             {-999., -999., -999., -999., -999., -999.}, {-999., -999., -999., -999., -999., -999.},
                             {-999., -999., -999., -999., -999., -999.}, {-999., -999., -999., -999., -999., -999.}};
  struct HitInfo
  {
    double hitX = -999.;
    double hitY = -999.;
    double hitZ = -999.;
  };

  HitInfo Allhit[G4Sol::SIZEOF_G4SOLDETTYPE];
};

class FullRecoEvent
{
  public:
  // map_mom3 Mom_Particle;
  std::unordered_map<int, ResSolDAF> DAF_results;

  std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > > ListHits;

  // std::vector< std::vector<std::vector<genfit::AbsMeasurement*> > > ListHitsDAF;
  // std::vector< std::vector<std::vector<std::vector<int> > > > ListIdHitsDAFCluster;
  std::unordered_map<int, std::vector<int> > TrackDAF;
  std::unordered_map<int, std::vector<SimHit> > TrackDAFSim;
  std::unordered_map<int, std::vector<InfoPar> > TrackInfo;
  std::unordered_map<int, std::tuple<int,double,double,double,double> > TrackMother;

  FullRecoEvent();
  ~FullRecoEvent();

  void Clear(int toclean = 0);
};

#endif
