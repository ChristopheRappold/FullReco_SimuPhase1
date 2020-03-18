#ifndef FullRecoEventZMQ_h
#define FullRecoEventZMQ_h

#include "TDatabasePDG.h"
#include "TGeoElement.h"
#include "TGeoManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"

#include <array>
#include <iterator>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

//#include "Genfit/Track.h"
//#include "Genfit/HypRecoHit.h"
#include "AbsMeasurement.h"
#include "PlanarMeasurement.h"
#include "ProlateSpacepointMeasurement.h"
#include "ReturnRes.hh"
#include "SpacepointMeasurement.h"
#include "msgpack.hpp"
#include "spdlog/spdlog.h"

#include <cassert>

namespace ZMQ
{
namespace G4Sol
{
enum SolDet : int
{
  InSi0 = 0,
  InSi1,
  InSi2,
  InSi3 /*3*/,
  TR1 /*4*/,
  TR2 /*5*/,
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


constexpr bool IsPlanar(ZMQ::G4Sol::SolDet idDet)
{
  switch(idDet)
    {
    case ZMQ::G4Sol::InSi0:;
    case ZMQ::G4Sol::InSi1:;
    case ZMQ::G4Sol::InSi2:;
    case ZMQ::G4Sol::InSi3:;
    case ZMQ::G4Sol::TR1:;
    case ZMQ::G4Sol::TR2:;
    case ZMQ::G4Sol::PSFE:;
    case ZMQ::G4Sol::PSBE:;
    case ZMQ::G4Sol::TrFwd0:;
    case ZMQ::G4Sol::TrFwd1:;
    case ZMQ::G4Sol::TrFwd2:;
    case ZMQ::G4Sol::RPC_l:;
    case ZMQ::G4Sol::RPC_h:;
    case ZMQ::G4Sol::FMF2Stop0:;
    case ZMQ::G4Sol::FMF2Stop1:;
    case ZMQ::G4Sol::FMF2Stop2:
      return true;
    default:
      return false;
    };
};

  
constexpr auto nameLiteralDet = {"InSi0", "InSi1", "InSi2",  "InSi3",  "TR1",   "TR2",    "PSFE",   "MG01",   "MG02",
                                 "MG03",  "MG04",  "MG05",   "MG06",   "MG07",  "MG08",   "MG09",   "MG10",   "MG11",
                                 "MG12",  "MG13",  "MG14",   "MG15",   "MG16",  "MG17",   "PSCE",   "PSBE",   "CDC0",
                                 "CDC1",  "CDC2",  "CDC3",   "CDC4",   "CDC5",  "CDC6",   "CDC7",   "CDC8",   "CDC9",
                                 "CDC10", "CDC11", "CDC12",  "CDC13",  "CDC14", "CDHBar", "TrFwd0", "TrFwd1", "TrFwd2",
                                 "RPCl",  "RPCh",  "FMF2S0", "FMF2S1", "FMF2S2"};

template <typename T, T... args>
struct EnumIter : public std::iterator<std::input_iterator_tag, T, std::ptrdiff_t, const T*, const T&>
{
  static constexpr T values[]       = {args...};
  static constexpr std::size_t size = sizeof...(args);

  int pos;
  EnumIter() // No value is end
      : pos(size)
  {
  }
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
using SolDetIter =
    struct EnumIter<SolDet, InSi0, InSi1, InSi2, InSi3 /*3*/, TR1 /*4*/, TR2 /*5*/, PSFE /*6*/, MG01 /*7*/, MG02, MG03,
                    MG04, MG05, MG06, MG07, MG08, MG09, MG10, MG11, MG12, MG13, MG14, MG15, MG16, MG17 /*23*/,
                    PSCE /*24*/, PSBE /*25*/, CDC_layer0 /*4 - 26*/, CDC_layer1, CDC_layer2, CDC_layer3, CDC_layer4,
                    CDC_layer5, CDC_layer6, CDC_layer7, CDC_layer8, CDC_layer9, CDC_layer10, CDC_layer11, CDC_layer12,
                    CDC_layer13, CDC_layer14 /*18 - 40*/, CDHBar /*19 - 41*/, TrFwd0 /*20 - 42*/, TrFwd1,
                    TrFwd2 /*22 - 43*/, RPC_l /*23 - 44*/, RPC_h /*24 - 45*/, FMF2Stop0 /*25 - 46*/, FMF2Stop1,
                    FMF2Stop2>;

} // namespace G4Sol

class PDG_fromName
{

  static const std::string ElName2[];

  TGeoElementTable* tableRN;
  std::unordered_map<std::string, double> cache_NucleiPID;

  const double u        = 0.931494061;
  const double m_eminus = 5.10998909999999971e-04; // GeV

  int lastPDGIon = 10004;

public:
  PDG_fromName() : tableRN(nullptr) {}
  ~PDG_fromName()
  {
    std::stringstream ss;
    ss << "PDG_fromName cache \n";
    for(const auto& it : cache_NucleiPID)
      ss << "[" << it.first << ", " << it.second << "] ";
    spdlog::get("Console")->info(ss.str());
  }

  int operator()(const std::string& name)
  {
    // std::cout<<"!> PID from name:"<<name<<"\n";

    if(name == "proton")
      return 2212;
    if(name == "He3")
      return 10003;
    if(name == "pi-")
      return -211;
    if(name == "pi+")
      return 211;
    if(name == "neutron")
      return 2112;
    if(name == "kaon0")
      return 311;
    if(name == "kaon+")
      return 321;
    if(name == "kaon-")
      return -321;
    if(name == "H3L")
      return 20001;
    if(name == "pi0")
      return 111;
    if(name == "triton")
      return 10001;
    if(name == "deuteron")
      return 10000;
    if(name == "alpha")
      return 10002;
    if(name == "He3")
      return 10003;
    if(name == "Li6")
      return 10004;

    auto it_name = cache_NucleiPID.find(name);
    if(it_name != cache_NucleiPID.end())
      return it_name->second;

    // std::cout<<"ADD NEW ELEM TO PDG DATABASE: "<<name<<"\n";

    std::size_t posCharge = name.find_first_of("0123456789");

    if(posCharge == std::string::npos)
      return 0;

    // std::cout<<"found Charge at"<<posCharge;
    int AtomMass            = std::stoi(name.substr(posCharge, std::string::npos));
    std::string nameElement = name.substr(0, posCharge);

    // std::cout<<" Unpack into:"<<nameElement<<" "<<AtomMass<<"\n";

    int id_Elem = -1;
    for(size_t i = 0; i < 111; ++i)
      if(nameElement == ElName2[i])
        {
          id_Elem = i;
          break;
        }
    if(id_Elem <= 0)
      return 0;

    // std::cout<<"Element found :"<<id_Elem<<" "<<ElName2[id_Elem]<<"\n";

    if(tableRN == nullptr)
      tableRN = gGeoManager->GetElementTable();

    auto* TempElement = tableRN->GetElementRN(AtomMass, id_Elem);
    double Dmass      = 0.;
    if(TempElement == nullptr)
      {
        spdlog::get("Console")->info("E> no element ! {} {} {}", AtomMass, id_Elem, fmt::ptr(TempElement));
        if(AtomMass == 23 && id_Elem == 14)
          Dmass = 23.073 * 1e-3; // MeV -> GeV
        else
          return 0;
      }
    else
      Dmass = TempElement->MassEx() * 1e-3; // MeV -> GeV

    double Mass = Dmass + AtomMass * u - id_Elem * m_eminus;

    TDatabasePDG::Instance()->AddParticle(name.c_str(), name.c_str(), Mass, kFALSE, 0., id_Elem * 3., "Ions",
                                          ++lastPDGIon);

    cache_NucleiPID.insert(std::make_pair(name, lastPDGIon));

    return lastPDGIon;

    // return 0;
  }
};

class MomRef
{
public:
  int id_track;

  double p_value;
  double mass;
  double x_TOF;
  double y_TOF;

  MomRef() : id_track(-1), p_value(-1.), mass(-1.), x_TOF(-1.), y_TOF(-1.){};
  MomRef(int id_, double pv_, double mass_, double x_, double y_)
      : id_track(id_), p_value(pv_), mass(mass_), x_TOF(x_), y_TOF(y_){};
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
  explicit TVertex(std::string t_)
      : type(std::move(t_)), vertex(TVector3(0., 0., 0.)), Chi2(0.), Ndf(0), Cov(TMatrixT<double>(3, 3))
  {
  }

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

  int pdg     = 0;
  double mass = 0.;

  double Eloss  = 0.;
  double time   = 0.;
  double length = 0.;
};

struct InfoPar
{
  int pdg       = 0;
  double momX   = -999.;
  double momY   = -999.;
  double momZ   = -999.;
  double mass   = -999.;
  double Eloss  = 0.;
  double time   = 0.;
  double length = 0.;
};

struct ResSolDAF
{
  int charge    = 0;
  int pdg_ini   = 0;
  int pdg_guess = 0;

  double pvalue = -1.;
  double chi2   = -1.;
  double ndf    = 0.;
  int fitter    = -1;

  double momX_init = -999.;
  double momY_init = -999.;
  double momZ_init = -999.;

  double momX = -999.;
  double momY = -999.;
  double momZ = -999.;

  double posX = -999.;
  double posY = -999.;
  double posZ = -999.;

  double mass  = -999.;
  double mass2 = -999.;
  double beta  = -999.;
  double beta2 = -999.;
  double tof   = -999.;
  double tof2  = -999.;

  double path_length  = -999.;
  double path_length2 = -999.;
  double path_time    = -999.;

  double dE = -999.;

  int firstHit = -1;
  int lastHit  = -1;
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

struct OutParticle
{
  std::string type;
  int Mc_id;
  int Mother_id;
  int Pdg;
  int Charge;
  std::array<double, 4> MomMass;
  std::array<double, 4> Vtx;
  double Weigth;
  bool GeoAcc;

  MSGPACK_DEFINE(type, Mc_id, Mother_id, Pdg, Charge, MomMass, Vtx, Weigth, GeoAcc);
};

struct OutHit
{
  std::string name;
  int LayerID;
  int HitID;
  std::array<double, 3> MCHit;
  std::array<double, 3> Hit;
  double Energy;
  double MCTime;
  double Time;
  double Length;
  int MC_id;
  int Charge;
  int Pdg;
  double Brho;
  double MagnetInteraction;
  std::array<double, 4> MCparticle;

  MSGPACK_DEFINE(name, LayerID, HitID, MCHit, Hit, Energy, MCTime, Time, Length, MC_id, Charge, Pdg, Brho,
                 MagnetInteraction, MCparticle);
};

struct OutTrack
{
  std::string type;
  int MC_status;
  double Chi2;
  double Chi2_X;
  double Chi2_Y;
  double Mass;
  int pdgcode;
  std::array<double, 4> MomMass;
  std::array<double, 3> Mom;

  int Charge;
  int BarId;
  double dE;
  double Beta;
  std::array<double,3> RefPoint;
  double Pval2;
  int TofsBar;
  double PathLength;
  double TOF;

  std::array<double, 3> MomIni;
  float RChiIni;
  double PathLengthIni;
  double TOFIni;
  double BetaIni;
  double MassIni;

  std::array<double, 4> Sim2Vtx;

  double State[6];
  double Cov[6][6];

  MSGPACK_DEFINE(type, MC_status, Chi2, Chi2_X, Chi2_Y, Mass, pdgcode, MomMass, Mom, Charge, BarId, dE, Beta, RefPoint, Pval2,
                 TofsBar, PathLength, TOF, MomIni, RChiIni, PathLengthIni, TOFIni, BetaIni, MassIni, Sim2Vtx, State,
                 Cov);
};

struct EventStatus
{
  ReturnRes::InfoM status;
  Long64_t BeginEvent;
  Long64_t EndEvent;

  MSGPACK_DEFINE(status, BeginEvent, EndEvent);
};

struct DataBuilderOut
{
  ReturnRes::InfoM status;
  Long64_t idEvent;
  std::vector<OutParticle> DumpParticles;
  std::vector<std::vector<OutHit> > DumpHits;

  MSGPACK_DEFINE(status, idEvent, DumpParticles, DumpHits);
};

struct DataFitterOut
{

  ReturnRes::InfoM status;
  Long64_t idEvent;
  DataBuilderOut previousStep;
  std::vector<OutTrack> DumpTracks;

  MSGPACK_DEFINE(status, idEvent, previousStep, DumpTracks);
};
}
#endif
