#ifndef FullRecoEvent_h
#define FullRecoEvent_h

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
#include "WireMeasurement.h"
#include "WirePointMeasurement.h"
//#include "WireMeasurementNew.h"
#include "ProlateSpacepointMeasurement.h"
#include "SpacepointMeasurement.h"
#include "spdlog/spdlog.h"

#include "HitAna/FiberHitXUV.hh"
#include "HitAna/MWDCTracking.hh"
#include "HitAna/FiberTrackAna.hh"
#include "HitAna/FiberHitAna.hh"

#include <cassert>

class MomRef;
class CompMomRef;

using map_mom3 = std::map<MomRef, std::vector<double>, CompMomRef>;

namespace G4Sol
{
enum SolDet : int
{
  InSi0 = 0,
  InSi1,
  InSi2,
  InSi3, /*3*/
  TR1, /*4*/
  TR2, /*5*/
  Si1x, /*6*/
  Si1y,
  Si2x,
  Si2y, /*9*/
  Si1x_SD, /*10*/
  Si1y_SD,
  Si2x_SD,
  Si2y_SD, /*13*/
  Si1x_SD_pad, /*14*/
  Si1y_SD_pad,
  Si2x_SD_pad,
  Si2y_SD_pad, /*17*/
  TO_Counter, /*18*/
  FiberD1_x, /*19*/
  FiberD1_u,
  FiberD1_v,
  FiberD2_x, /*22*/
  FiberD2_u,
  FiberD2_v,
  FiberD3_x, /*25*/
  FiberD3_u,
  FiberD3_v,
  MiniFiberD1_x, /*28*/
  MiniFiberD1_u,
  MiniFiberD1_v,
  MiniFiberD2_x,
  MiniFiberD2_v,
  MiniFiberD2_u,
  FiberD4_v, /*34*/
  FiberD4_u,
  FiberD4_x,
  FiberD5_x, /*37*/
  FiberD5_u,
  FiberD5_v,
  PSFE, /*40*/ //PSB_F
  MG01, /*41*/
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
  MG17, /*57*/
  PSCE, /*58*/   //PSB
  PSBE, /*59*/   //PSB_B
  MWDC01, /*60*/
  MWDC02,
  MWDC03,
  MWDC04,
  MWDC05,
  MWDC06,
  MWDC07,
  MWDC08,
  MWDC09,
  MWDC10,
  MWDC11,
  MWDC12,
  MWDC13,
  MWDC14,
  MWDC15,
  MWDC16, /*75*/
  FiberD1_xy,
  FiberD2_xy,
  FiberD3_xy,
  FiberD4_xy,
  FiberD5_xy,
  MiniFiberD1_xy,
  MiniFiberD2_xy,
  CDC_layer0, /*60*/
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
  CDC_layer14, /*74*/
  CDHBar, /*75*/
  TrFwd0, /*76*/
  TrFwd1,
  TrFwd2, /*77*/
  RPC_l, /*78*/
  RPC_h, /*79*/
  FMF2Stop0, /*80*/
  FMF2Stop1,
  FMF2Stop2, /*82*/
  SIZEOF_G4SOLDETTYPE
};

constexpr auto nameLiteralDet = {
    "InSi0"         ,          "InSi1",          "InSi2",          "InSi3",            "TR1",            "TR2",
    "Si1x"          ,           "Si1y",           "Si2x",           "Si2y",
    "Si1x_SD"          ,           "Si1y_SD",           "Si2x_SD",           "Si2y_SD",
    "Si1x_SD_pad"          ,       "Si1y_SD_pad",       "Si2x_SD_pad",       "Si2y_SD_pad",       "T0_Counter",
    "FiberD1_x",           "FiberD1_u",      "FiberD1_v",      "FiberD2_x",      "FiberD2_u",      "FiberD2_v",
    "FiberD3_x",           "FiberD3_u",      "FiberD3_v",      
    "MiniFiberD1_x", "MiniFiberD1_u", "MiniFiberD1_v", "MiniFiberD2_x", "MiniFiberD2_v", "MiniFiberD2_u",
    "FiberD4_v",      "FiberD4_u",      "FiberD4_x",           "FiberD5_x",      "FiberD5_u",      "FiberD5_v",
    "PSFE",     "MG01",           "MG02",
    "MG03",                     "MG04",           "MG05",           "MG06",           "MG07",           "MG08",
    "MG09",                     "MG10",           "MG11",           "MG12",           "MG13",           "MG14",
    "MG15",                     "MG16",           "MG17",           "PSCE",           "PSBE",
    "MWDC01", "MWDC02", "MWDC03", "MWDC04", "MWDC05", "MWDC06", "MWDC07", "MWDC08", "MWDC09", "MWDC10", "MWDC11",
    "MWDC12", "MWDC13", "MWDC14", "MWDC15", "MWDC16",
    "FiberD1_xy", "FiberD2_xy", "FiberD3_xy", "FiberD4_xy",  "FiberD5_xy",  "MiniFiberD1_xy", "MiniFiberD2_xy",
    "CDC0",
    "CDC1",                     "CDC2",           "CDC3",           "CDC4",           "CDC5",           "CDC6",
    "CDC7",                     "CDC8",           "CDC9",          "CDC10",          "CDC11",          "CDC12",
    "CDC13",                   "CDC14",         "CDHBar",         "TrFwd0",         "TrFwd1",         "TrFwd2",
    "RPCl",                     "RPCh",         "FMF2S0",         "FMF2S1",         "FMF2S2"};

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
  using SolDetIter = struct EnumIter<
    SolDet, InSi0, InSi1, InSi2, InSi3 /*3*/, TR1 /*4*/, TR2 /*5*/,
    Si1x, Si1y, Si2x, Si2y,
    Si1x_SD, Si1y_SD, Si2x_SD, Si2y_SD, TO_Counter,
    FiberD1_x, FiberD1_u, FiberD1_v, FiberD2_x, FiberD2_u, FiberD2_v, FiberD3_x, FiberD3_u, FiberD3_v,
    MiniFiberD1_x, MiniFiberD1_u, MiniFiberD1_v, MiniFiberD2_x, MiniFiberD2_v, MiniFiberD2_u,
    FiberD4_v, FiberD4_u, FiberD4_x, FiberD5_x, FiberD5_u, FiberD5_v,
    PSFE /*6*/, MG01 /*7*/, MG02, MG03, MG04, MG05, MG06, MG07,
    MG08, MG09, MG10, MG11, MG12, MG13, MG14, MG15, MG16, MG17 /*23*/, PSCE /*24*/, PSBE /*25*/,
    FiberD1_xy,FiberD2_xy,FiberD3_xy,FiberD4_xy,FiberD5_xy,MiniFiberD1_xy,MiniFiberD2_xy,
    CDC_layer0 /*4 - 26*/,
    CDC_layer1, CDC_layer2, CDC_layer3, CDC_layer4, CDC_layer5, CDC_layer6, CDC_layer7, CDC_layer8, CDC_layer9,
    CDC_layer10, CDC_layer11, CDC_layer12, CDC_layer13, CDC_layer14 /*18 - 40*/, CDHBar /*19 - 41*/, TrFwd0 /*20 - 42*/,
    TrFwd1, TrFwd2 /*22 - 43*/, RPC_l /*23 - 44*/, RPC_h /*24 - 45*/, FMF2Stop0 /*25 - 46*/, FMF2Stop1, FMF2Stop2>;

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

struct InfoInit
{
  int charge  = -999 ;
  double time = -999.;
  double posX = -999.;
  double posY = -999.;
  double posZ = -999.;
  double momX = -999.;
  double momY = -999.;
  double momZ = -999.;
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
  double TOT    = 0.;
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
  int Nmfiber = 0;
  int iterNum  = -1;

  double cov_matrix[6][6] = {{-999., -999., -999., -999., -999., -999.}, {-999., -999., -999., -999., -999., -999.},
                             {-999., -999., -999., -999., -999., -999.}, {-999., -999., -999., -999., -999., -999.},
                             {-999., -999., -999., -999., -999., -999.}, {-999., -999., -999., -999., -999., -999.}};
  struct HitInfo
  {
    double hitX = -999.;
    double hitY = -999.;
    double hitZ = -999.;
  };

  double ResMDC[17][3]    = {{-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.},
                          {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.},
                          {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.},
                          {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.},
                          {-999., -999., -999.}};
  double WeightMDC[17][3] = {{-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.},
                             {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.},
                             {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.},
                             {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.}, {-999., -999., -999.},
                             {-999., -999., -999.}};
  double ResFiber[9]      = {-999., -999., -999., -999., -999., -999., -999., -999., -999.};
  double ResMiniFiber[6]   = {-999., -999., -999., -999., -999., -999.};
  double WeightFiber[9]   = {-999., -999., -999., -999., -999., -999., -999., -999., -999.};
  double WeightMiniFiber[6]   = {-999., -999., -999., -999., -999., -999.};
  double ResPSCE[2]       = {-999., -999.};
  double WeightPSCE[2]    = {-999., -999.};

  HitInfo Allhit[G4Sol::SIZEOF_G4SOLDETTYPE];
};


struct DecayTrackInfo
{
  int Pdg       = 0;
  int Id        = 0;
  double Chi2   = -1.;
  TVector3 Hit_Pos = TVector3(0.,0.,0.);
  TLorentzVector Hit_MomEnergy = TLorentzVector(0.,0.,0.,0.);
};

struct IdHit
{
  int id_det;
  int id_hit;
  int id_track;
  IdHit(int iD,int iH, int iT):id_det(iD),id_hit(iH),id_track(iT) { }
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
};

struct OutHit
{
  std::string name;
  int LayerID;
  int HitID;
  std::array<double, 3> MCHit;
  std::array<double, 3> Hit;

  int MC_id;
  int Charge;
  int Pdg;
  double Brho;
  double MagnetInteraction;
  std::array<double, 4> MCparticle;
};


class PrimaryVtxTrack 
{
  public:
    PrimaryVtxTrack() { };
    ~PrimaryVtxTrack() = default;

    void SetX(double _x) { x = _x; };
    void SetY(double _y) { y = _y; };
    void SetA(double _a) { a = _a; };
    void SetB(double _b) { b = _b; };
    void SetChi2NDF(double _chi2ndf) { chi2ndf = _chi2ndf; };
    void SetFTHit(size_t i_ft, double hit_pos) { ft_hits.emplace_back(std::make_tuple(i_ft, hit_pos)); };

    double GetX() { return x; };
    double GetY() { return y; };
    double GetA() { return a; };
    double GetB() { return b; };
    double GetChi2NDF() { return chi2ndf; };
    std::vector<std::tuple<size_t,double>> GetFTHits() { return ft_hits; };

    //bool IsFTHit(size_t i_ft);
    double GetTheta();
    double GetPhi();
    size_t GetNHits() { return ft_hits.size(); };

  private:
    double x = -999.;
    double y = -999.; //z = mid of target;
    double a = -999.;
    double b = -999.;
    double chi2ndf = -999.;
    std::vector<std::tuple<size_t,double>> ft_hits = {};
};


class FragmentTrack 
{
  public:
    FragmentTrack() { };
    ~FragmentTrack() = default;

    void SetPos(TVector3 _pos) { pos = _pos; };
    void SetMom(TVector3 _mom) { mom = _mom; };
    void SetCovMatrix(std::vector<double> _cov_matrix) { cov_matrix = _cov_matrix; };
    void SetTOT(double _tot) { tot = _tot; };
    void SetTime(double _time) { time = _time; };
    void SetChi2NDF(double _chi2ndf) { chi2ndf = _chi2ndf; };
    void SetIsBest(bool _isbest) { isbest = _isbest; };
    void SetPID(int _pid) { pid = _pid; };


    TVector3 GetPos() { return pos; };
    TVector3 GetMom() { return mom; };
    std::vector<double> GetCovMatrix() { return cov_matrix; };
    double GetTOT() { return tot; };
    double GetTime() { return time; };
    double GetChi2NDF() { return chi2ndf; };
    bool GetIsBest() { return isbest; };
    int GetPID() { return pid; };

  private:
    TVector3 pos = TVector3(-999., -999., -999.);
    TVector3 mom = TVector3(-999., -999., -999.);
    std::vector<double> cov_matrix = {};
    double tot = -999.;
    double time = -999.;
    double chi2ndf = -999.;
    bool isbest = false;
    int pid = -999;
};


struct Hyp
{
  Int_t Pattern; /// 1 = Simulation / 2 = KFParticle_real / 3 = KFParticle_cut / 4 = KFParticle / 5 = KFParticle_Mass / 6 = LorentzVector
  
  //Mother:
  Int_t PDG;
  Int_t N_Mother;
  Double32_t Chi2ndf;
  Int_t NDF;
  TLorentzVector MomE;
  TVector3 PrimVtx;
  TVector3 DecayVtx;
  TVector3 Dist_RealReconsVtx;
  Double32_t Dist_MotherPrimVtx;
  Double32_t Angle_MotherPrimVtx;
  Double32_t InvMass;
  Double32_t ErrInvMass;
  Int_t ErrGetMass;
  Double32_t LifeTime;
  Double32_t ErrLifeTime;
  Int_t ErrGetLifeTime;
  Int_t Mother_IsFromHyp;

  //Daughters:
  Int_t PDG_Fragment;
  Int_t Id_Fragment;
  TLorentzVector MomE_Fragment;
  Double32_t Chi2ndf_Fragment;
  Int_t NDF_Fragment;
  Double32_t Pvalue_Fragment;
  Double32_t Angle_MotherFragment;
  Int_t Fragment_IsFromHyp; // 0-> No; 1-> Yes

  Int_t PDG_Pion;
  Int_t Id_Pion;
  TLorentzVector MomE_Pion;
  Double32_t Chi2ndf_Pion;
  Int_t NDF_Pion;
  Double32_t Pvalue_Pion;
  Double32_t Angle_MotherPion;
  Int_t NHitsMDC_Pion;
  Int_t NHitsMinifiber_Pion;
  Int_t N_Pion;
  Int_t Pion_IsFromHyp; // 0-> No; 1-> Yes

  Double32_t Dist_Daughters;
  Double32_t ArmPod_Qt;
  Double32_t ArmPod_Alfa;
};

struct KFFitInfo
{
  Double32_t Pvalue = -1.;
  Int_t NHitsMDC = -1;
  Int_t NHitsMinifiber = -1;
};

struct CandTrack
{
  std::vector<int> orderedHitIds;
  bool Rfitted = false;
  int q = 0;
  std::vector<double> par = std::vector<double>(5);
  TMatrixD Cov = TMatrixD(5,5);
  double chi2_circle = -1.;
  double chi2_line = -1.;

  CandTrack() = default;
  CandTrack(const std::vector<int>& Hits, bool R, int Q, const std::vector<double>& P, const TMatrixD& C, double chi2_1, double chi2_2):orderedHitIds(Hits),Rfitted(R),q(Q),par(P),Cov(C),chi2_circle(chi2_1),chi2_line(chi2_2) {};
};

// struct OutTrack
// {
//   std::string type;
//   int MC_status;
//   double Chi2;
//   double Chi2_X;
//   double Chi2_Y;
//   double Mass;
//   int pdgcode;
//   std::array<double,4> MomMass;
//   std::array<double,3> Mom;

//   int Charge;
//   int BarId;
//   double dE;
//   double Beta;
//   double Pval2;
//   int TofsBar;
//   double PathLength;
//   double TOF;

//   std::array<double,3> MomIni;
//   float RChiIni;
//   double PathLengthIni;
//   double TOFIni;
//   double BetaIni;
//   double MassIni;

//   std::array<double,4> Sim2Vtx;

//   double State[6];
//   double Cov[6][6];
// };


struct MeasurementInfo
{
  double TOT = -9999.;
  double time = -9999.;
  double dE = -9999.;
  //size_t layer = 100000;

  MeasurementInfo() = default;
  MeasurementInfo(double TOT_, double time_, double dE_):TOT(TOT_),time(time_),dE(dE_) {};
  void SetInfo(double TOT_, double time_, double dE_) { TOT = TOT_; time = time_; dE = dE_;};
};


class FullRecoEvent
{
public:
  unsigned int idThread;
  Long64_t idEvent;

  std::vector<OutParticle> ToDumpParticles;
  std::vector<std::vector<OutHit> > ToDumpHits;
  // std::vector<OutTrack> ToDumpTracks;

  // map_mom3 Mom_Particle;
  std::unordered_map<int, ResSolDAF> DAF_results;

  std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > > ListHits;
  std::vector<std::vector<std::unique_ptr<genfit::AbsMeasurement> > > OldListHits;
  std::vector<std::vector<MeasurementInfo> > ListHitsInfo;

  std::vector<std::vector<int> > ListHitsToTracks;
  // std::vector< std::vector<std::vector<genfit::AbsMeasurement*> > > ListHitsDAF;
  // std::vector< std::vector<std::vector<std::vector<int> > > > ListIdHitsDAFCluster;
  std::unordered_map<int, std::vector<int> > TrackDAF;
  std::unordered_map<int, std::vector<std::vector<SimHit> > > TrackDAFSim;
  std::unordered_map<int, InfoInit> TrackDAFInit;
  std::unordered_map<int, InfoInit> TrackDAFInitSim;
  std::unordered_map<int, std::vector<InfoPar> > TrackInfo;

  std::vector<CandTrack> TracksFound;
  std::vector<IdHit> IdHitsToMeasurement;

  std::unordered_map<int, std::tuple<int, double, double, double, double> > TrackMother;
  std::unordered_map<int, InfoInit> DaughtersTrackDAFInit;

  std::vector<MWDCTracking> MWDCTracks;

  std::vector<PrimaryVtxTrack> BeamTracks;
  std::vector<PrimaryVtxTrack> PrimaryTracks;

  std::vector<FragmentTrack> FragmentTracks;
  int FragmentPID;

  std::map< std::string, std::vector<FiberTrackAna*> > FiberTrackCont;
  std::vector<std::vector<std::vector<FiberHitAna*> > > FiberHitClCont;

/*
  std::vector<std::unordered_map<size_t, double > > Si_HitsEnergyLayer;
  
  std::vector<std::vector<double> > Hits_Si1{};
  std::vector<std::vector<double> > Hits_Si2{};
  std::vector<std::vector<double> > Hits_Si3{};
  std::vector<std::vector<double> > Hits_Si4{};
  
  std::vector<std::tuple<double, size_t> > HitsX_Si1{};
  std::vector<std::tuple<double, size_t> > HitsY_Si1{};
  std::vector<std::tuple<double, size_t> > HitsX_Si2{};
  std::vector<std::tuple<double, size_t> > HitsY_Si2{};
*/

  TLorentzVector Mother_MomE;
  std::array<double,3> InteractionPoint;
  std::array<double,3> DecayVertex;
  double Hyp_LifeTime;

  TVector3 PrimVtxRecons;
  std::array<double,6> CovMatrix_IP; //lower triangle

  TVector3 DecayVtxRecons;
  std::array<double,6> CovMatrix_SV; //lower triangle

  std::vector<Hyp> Hyp_Vect;

  
  explicit FullRecoEvent(unsigned int idTh = 1);
  ~FullRecoEvent();

    void Clear(int toclean = 0);
};

#endif
