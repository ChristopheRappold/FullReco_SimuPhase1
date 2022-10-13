#include "TDataHit.hh"
//#include "Riostream.h"
//#define DEBUG_MEM

using namespace std;

// ClassImp(TDataHit)

// TDataHit::TDataHit():name("default_hit"),LayerID(-1),HitID(-1),TrackID(-1),Hit(-999.,-999.,-999.),Charge(0),Pdg(0),Brho(0.),MagnetInteraction(0.),Time(0.),Energy(0.),TrackLength(0.)
// {

// }

// TDataHit::~TDataHit()
// {
// #ifdef DEBUG_MEM
//   cout<<"HitDet Destructor :"<<name<<endl;
// #endif
// }

// TDataHit::TDataHit(const TDataHit&
// M):LayerID(M.LayerID),HitID(M.HitID),TrackID(M.TrackID),Hit(M.Hit),Charge(M.Charge),Pdg(M.Pdg),Brho(M.Brho),MagnetInteraction(M.MagnetInteraction),Time(M.Time),Energy(M.Energy),TrackLength(M.TrackLength)
// {
//   name = "const_copy";
//   name+=M.name;
// #ifdef DEBUG_MEM
//   cout<<"HitDet const Copy Constructor :"<<name<<endl;
// #endif

//   //HitPos = M.HitPos;
//   //Id = M.Id;
//   //Qdc = M.Qdc;
//   //TimeCh=M.TimeCh;
//   //TimeNs = M.TimeNs;
// }

// TDataHit& TDataHit::operator=(const TDataHit& M)
// {
//   name = M.name;
//   name+="op=";
//   LayerID = M.LayerID;
//   HitID = M.HitID;
//   TrackID = M.TrackID;
//   Hit.SetXYZ(M.Hit.X(),M.Hit.Y(),M.Hit.Z());
//   Charge = M.Charge;
//   Pdg = M.Pdg;
//   Brho = M.Brho;
//   MagnetInteraction = M.MagnetInteraction;
//   Time = M.Time;
//   Energy = M.Energy;
//   TrackLength = M.TrackLength;

//   return *this;
// }

// void TDataHit::Clear(Option_t * /*option*/)
// {
// #ifdef DEBUG_MEM
//   cout<<"HitDet Clearing:"<<name<<endl;
// #endif
//   name="";
//   LayerID = -1;
//   HitID = -1;
//   TrackID = -1;
//   Hit.SetXYZ(-999.,-999.,-999.);
//   Charge = 0;
//   Pdg = 0;
//   Brho = 0;
//   MagnetInteraction = 0.;
//   Time = 0.;
//   Energy = 0.;
//   TrackLength = 0.;
// }

void TFrsTpcHit::Clear(Option_t* /*option*/)
{

  LayerID  = -1;
  tpc_dt_s = {-9999, -9999, -9999, -9999};
  tpc_lt_s = {-9999, -9999};
  tpc_rt_s = {-9999, -9999};
  tpc_xraw = {-9999, -9999};
  tpc_yraw = {-9999, -9999, -9999, -9999};
  tpc_csum = {-9999, -9999, -9999, -9999};

  tpc_x    = -9999.;
  tpc_y    = -9999.;
  tpc_de   = -9999.;
  tpc_dx12 = -9999.;

  tpc_timeref_s = -9999;
};

void TS4MwdcHit::Clear(Option_t* /*option*/)
{

  LayerID = -1;
  HitID   = -1;

  t_leading  = -9999.;
  t_trailing = -9999.;
};

void TWfdHit::Clear(Option_t* /*option*/)
{

  TypeID = NWFD;
  HitID  = -1;

  Wfd_qdc       = -9999.;
  Wfd_max       = -9999.;
  Wfd_min       = -9999.;
  Wfd_max_index = -9999.;
  Wfd_min_index = -9999.;
};

void TWasaPSHit::Clear(Option_t* /*option*/)
{

  TypeID = NPTYPE; // Back or Front
  HitID  = -1;
  Time   = -9999.;
  PosZ   = -9999.;
  dE     = -9999.;
};

void TMDCHit::Clear(Option_t* /*option*/)
{
  LayerID = -1;
  HitID   = -1;

  t_leading  = -9999;
  t_trailing = -9999;
};

void TFiberHit::Clear(Option_t* /*option*/)
{
  TypeID = NFTYPE;
  HitID  = -1;

  t_leading  = -9999;
  t_trailing = -9999;
};

void TCsIHit::Clear(Option_t* /*option*/)
{

  LR             = -1;
  longitudinalId = -1;
  phiId          = -1;
  adc.clear();
  ampl  = -999.;
  trapz = -999.;
};

void TSciHit::Clear(Option_t* /*option*/)
{
  TypeID = NSC;
  dT     = -9999.;
  dE     = -9999.;
};
