#ifndef THYPHIWASADATAHIT_H
#define THYPHIWASADATAHIT_H

#include "Riostream.h"
#include "TLorentzVector.h"
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

#include <array>

// class TDataHit : public TObject
// {
// public :

//   TString name; //!

//   Int_t LayerID;
//   Int_t HitID;
//   Int_t TrackID;

//   TVector3 Hit;

//   Int_t Charge;
//   Int_t Pdg;
//   Double_t Brho;
//   Double_t MagnetInteraction;

//   Double_t Time;
//   Double_t Energy;
//   Double_t TrackLength;

//   TDataHit();

//   TDataHit(const TDataHit& M);
//   TDataHit& operator=(const TDataHit&);

//   virtual  void Clear(Option_t* ="");
//   ~TDataHit();

//   ClassDef(TDataHit,3)

// };

class TFrsTpcHit : public TObject
{
public:
  Int_t LayerID                 = -1;
  std::array<Int_t, 4> tpc_dt_s = {-9999, -9999, -9999, -9999};
  std::array<Int_t, 2> tpc_lt_s = {-9999, -9999};
  std::array<Int_t, 2> tpc_rt_s = {-9999, -9999};
  std::array<Int_t, 2> tpc_xraw = {-9999, -9999};
  std::array<Int_t, 4> tpc_yraw = {-9999, -9999, -9999, -9999};
  std::array<Int_t, 4> tpc_csum = {-9999, -9999, -9999, -9999};

  Float_t tpc_x    = -9999.;
  Float_t tpc_y    = -9999.;
  Float_t tpc_de   = -9999.;
  Float_t tpc_dx12 = -9999.;

  Int_t tpc_timeref_s = -9999;

  TFrsTpcHit(){};

  TFrsTpcHit(const TFrsTpcHit& M) = default;
  TFrsTpcHit& operator=(const TFrsTpcHit&) = default;

  virtual void Clear(Option_t* = "");
  ~TFrsTpcHit(){};

  ClassDef(TFrsTpcHit, 1)
};

class TS4MwdcHit : public TObject
{
public:
  Int_t LayerID = -1;
  Int_t HitID   = -1;

  Int_t t_leading  = -9999.;
  Int_t t_trailing = -9999.;

  TS4MwdcHit(){};

  TS4MwdcHit(const TS4MwdcHit& M) = default;
  TS4MwdcHit& operator=(const TS4MwdcHit&) = default;

  virtual void Clear(Option_t* = "");
  ~TS4MwdcHit(){};

  ClassDef(TS4MwdcHit, 1)
};

class TWfdHit : public TObject
{
public:
  enum WType
  {
    SC31_L = 0,
    SC31_R,
    SC41_L,
    SC41_R,
    SC42_L,
    SC42_R,
    SC43_L,
    SC43_R,
    PSB_L,
    PSB_R,
    PSBE,
    PSFE,
    NWFD
  };

  WType TypeID = NWFD;
  Int_t HitID  = -1;

  Int_t Wfd_qdc       = -9999.;
  Int_t Wfd_max       = -9999.;
  Int_t Wfd_min       = -9999.;
  Int_t Wfd_max_index = -9999.;
  Int_t Wfd_min_index = -9999.;

  TWfdHit(){};

  TWfdHit(const TWfdHit& M) = default;
  TWfdHit& operator=(const TWfdHit&) = default;

  virtual void Clear(Option_t* = "");
  ~TWfdHit(){};

  ClassDef(TWfdHit, 1)
};

class TWasaPSHit : public TObject
{
public:
  enum PType
  {
    PSB,
    PSFE,
    PSBE,
    T0,
    NPTYPE
  };

  PType TypeID = NPTYPE; // Back or Front
  Int_t HitID  = -1;

  // Short_t Qdc;
  // Int_t Nhit;
  // std::vector<Int_t> Tdc;
  Float_t Time = -9999.;
  Float_t PosZ = -9999.;
  Float_t dE   = -9999.;

  TWasaPSHit(){};

  TWasaPSHit(const TWasaPSHit& M) = default;
  TWasaPSHit& operator=(const TWasaPSHit&) = default;

  virtual void Clear(Option_t* = "");
  ~TWasaPSHit(){};

  ClassDef(TWasaPSHit, 1)
};

class TMDCHit : public TObject
{
public:
  Int_t LayerID = -1;
  Int_t HitID   = -1;

  Int_t t_leading  = -9999;
  Int_t t_trailing = -9999;

  TMDCHit(){};

  TMDCHit(const TMDCHit& M) = default;
  TMDCHit& operator=(const TMDCHit&) = default;

  virtual void Clear(Option_t* = "");
  ~TMDCHit(){};

  ClassDef(TMDCHit, 1)
};

class TFiberHit : public TObject
{
public:
  enum FType
  {
    UFT1_0 = 0,
    UFT1_1,
    UFT1_2,
    UFT2_0,
    UFT2_1,
    UFT2_2,
    UFT3_0,
    UFT3_1,
    UFT3_2,
    MFT_0,
    MFT_1,
    MFT_2,
    MFT_3,
    MFT_4,
    MFT_5,
    DFT1_0,
    DFT1_1,
    DFT1_2,
    DFT2_0,
    DFT2_1,
    DFT2_2,
    NFTYPE
  };

  FType TypeID = NFTYPE;
  Int_t HitID  = -1;

  Int_t t_leading  = -9999;
  Int_t t_trailing = -9999;

  TFiberHit(){};

  TFiberHit(const TFiberHit& M) = default;
  TFiberHit& operator=(const TFiberHit&) = default;

  virtual void Clear(Option_t* = "");
  ~TFiberHit(){};

  ClassDef(TFiberHit, 1)
};

class TCsIHit : public TObject
{
public:
  Int_t LR             = -1;
  Int_t longitudinalId = -1;
  Int_t phiId          = -1;
  std::vector<Short_t> adc;
  Float_t ampl  = -999.;
  Float_t trapz = -999.;

  TCsIHit(){};

  TCsIHit(const TCsIHit& M) = default;
  TCsIHit& operator=(const TCsIHit&) = default;

  virtual void Clear(Option_t* = "");
  ~TCsIHit(){};

  ClassDef(TCsIHit, 1)
};

class TSciHit : public TObject
{
public:
  enum SType
  {
    SC31,
    SC41_0,
    SC41_1,
    SC41_2,
    SC42,
    SC42LOWGAIN,
    SC43_0,
    SC43_1,
    SC43_2,
    SC43_3,
    NSC
  };
  SType TypeID = NSC;

  Float_t dT = -9999.;
  Float_t dE = -9999.;

  TSciHit();

  TSciHit(const TSciHit& M) = default;
  TSciHit& operator=(const TSciHit&) = default;

  virtual void Clear(Option_t* = "");
  ~TSciHit(){};

  ClassDef(TSciHit, 1)
};

#endif
