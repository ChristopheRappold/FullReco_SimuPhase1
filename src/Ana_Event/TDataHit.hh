#ifndef THYPHIDATAHIT_H
#define THYPHIDATAHIT_H

#include "TObject.h"
#include "TString.h"

#include "Riostream.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class TDataHit : public TObject
{
public :

  TString name; //!

  Int_t LayerID;
  Int_t HitID;
  Int_t TrackID;

  TVector3 Hit;

  Int_t Charge;
  Int_t Pdg;
  Double_t Brho;
  Double_t MagnetInteraction;

  Double_t Time;
  Double_t Energy;
  Double_t TrackLength;

  TDataHit();

  TDataHit(const TDataHit& M); 
  TDataHit& operator=(const TDataHit&);

  virtual  void Clear(Option_t* ="");
  ~TDataHit();

 
  ClassDef(TDataHit,3)

};

class TFrsTpcHit : public TObject
{
public :

  Int_t LayerID;
  std::array<Int_t, 4> tpc_dt_s;
  std::array<Int_t, 2> tpc_lt_s;
  std::array<Int_t, 2> tpc_rt_s;
  std::array<Int_t, 2> tpc_xraw;
  std::array<Int_t, 4> tpc_yraw;
  std::array<Int_t, 4> tpc_csum;

  Float_t tpc_x;
  Float_t tpc_y;
  Float_t tpc_de;
  Float_t tpc_dx12;

  Int_t tpc_timeref_s;

  TFrsTpcHit();

  TFrsTpcHit(const TFrsTpcHit& M);
  TFrsTpcHit& operator=(const TFrsTpcHit&);

  virtual  void Clear(Option_t* ="");
  ~TFrsTpcHit();


  ClassDef(TFrsTpcHit,1)

};


class S4MwdcHit : public TObject
{
public :

  Int_t LayerID;
  Int_t HitID;

  Int_t t_leading;
  Int_t t_trailing;

  S4MwdcHit();

  S4MwdcHit(const S4MwdcHit& M);
  S4MwdcHit& operator=(const S4MwdcHit&);

  virtual  void Clear(Option_t* ="");
  ~S4MwdcHit();


  ClassDef(S4MwdcHit,1)

};

class TWfdHit : public TObject
{
public :
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

  WType TypeID; //!
  Int_t HitID;

  Int_t Wfd_qdc;
  Int_t Wfd_max;
  Int_t Wfd_min;
  Int_t Wfd_max_index;
  Int_t Wfd_min_index;

  TWfdHit();

  TWfdHit(const TWfdHit& M);
  TWfdHit& operator=(const TWfdHit&);

  virtual  void Clear(Option_t* ="");
  ~TWfdHit();

  ClassDef(TWfdHit,1)

};

class TWasaPSHit : public TObject
{
public :
  enum PType
    {
      PSB,
      PSFE,
      PSBE,
      T0,
      NPTYPE
    };


  PType TypeID; // Back or Front
  Int_t HitID;

  // Short_t Qdc;
  // Int_t Nhit;
  // std::vector<Int_t> Tdc;
  Float_t Time;
  Float_t PosZ;
  Float_t dE;


  TPSBHit();

  TPSBHit(const TPSBHit& M);
  TPSBHit& operator=(const TPSBHit&);

  virtual  void Clear(Option_t* ="");
  ~TPSBHit();


  ClassDef(TPSBHit,1)

};

class TMDCHit : public TObject
{
public :

  TString name; //!

  Int_t LayerID;
  Int_t HitID;

  Int_t t_leading;
  Int_t t_trailing;

  TMDCHit();

  TMDCHit(const TMDCHit& M);
  TMDCHit& operator=(const TMDCHit&);

  virtual  void Clear(Option_t* ="");
  ~TMDCHit();


  ClassDef(TMDCHit,1)

};

class TFiberHit : public TObject
{
public :

  enum Ftype
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

  FType TypeID;
  Int_t HitID;

  Int_t t_leading  = -9999;
  Int_t t_trailing = -9999;


  TFiberHit();

  TFiberHit(const TFiberHit& M);
  TFiberHit& operator=(const TFiberHit&);

  virtual  void Clear(Option_t* ="");
  ~TFiberHit();

  ClassDef(TFiberHit,1)

};

class TCsIHit : public TObject
{
public :

  Int_t LR             = -1;
  Int_t longitudinalId = -1;
  Int_t phiId          = -1;
  std::vector<Short_t> adc;
  Float_t ampl         = -999.;
  Float_t trapz        = -999.;


  TCsIHit();

  TCsIHit(const TCsIHit& M);
  TCsIHit& operator=(const TCsIHit&);

  virtual  void Clear(Option_t* ="");
  ~TCsIHit();


  ClassDef(TCsIHit,1)

};


class TSciHit : public TObject
{
public :

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
  SType TypeID;

  Float_t dT;
  Float_t dE;


  TSciHit();

  TSciHit(const TSciHit& M);
  TSciHit& operator=(const TSciHit&);

  virtual  void Clear(Option_t* ="");
  ~TSciHit();


  ClassDef(TSciHit,1)

};



#endif
