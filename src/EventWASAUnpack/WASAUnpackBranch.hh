#pragma once
#include <TObject.h>
#include <vector>
#include <array>
#include <TString.h>

class LMDHeader : public TObject
{
  public:
  LMDHeader();
  ~LMDHeader();

  TString lmd_filename;
  TString lmd_open_time;

  ClassDef(LMDHeader, 1);
};

class S4TQ : public TObject
{
  template <class T, std::size_t N, std::size_t M>
  using array2 = std::array<std::array<T, M>, N>; // T val[N][M]
  public:
  S4TQ();
  ~S4TQ();
  void Init();
  ULong64_t ts = 0;
  Int_t eventnr = -999;
  Int_t mtdc_ec = -9999;
  Int_t v792_ec = -9999;
  Int_t v775_ec = -9999;
  Int_t v830_ec = -9999;
  Int_t mtdc_ec_diff = -9999;
  Int_t v792_ec_diff = -9999;
  Int_t v775_ec_diff = -9999;
  Int_t v830_ec_diff = -9999;
  std::array<Int_t, 32> v775;
  std::array<Bool_t, 16> trig;
  std::array<UInt_t, 32> v830;
  std::array<Int_t, 32> v830diff;
  std::array<Int_t, 34> mtdc_nhit; //32->34
  std::vector<std::vector<Int_t>> mtdc; // mtdc[34][]
  std::array<Short_t, 32> v792;
  Bool_t v792_ok = false;
  Bool_t v775_ok = false;
  Bool_t mtdc_ok = false;
  Bool_t v830_ok = false;
  Bool_t all_ok  = false;
  array2<Int_t, 3, 2> tdc_sc41;
  std::array<Int_t, 2> tdc_sc42;
  array2<Int_t, 4, 2> tdc_sc43;
  std::array<Int_t, 2> tdc_sc31;
  array2<Int_t, 3, 2> qdc_sc41;
  std::array<Int_t, 2> qdc_sc42;
  std::array<Int_t, 2> qdc_sc42_lowgain;
  array2<Int_t, 4, 2> qdc_sc43;
  std::array<Int_t, 2> qdc_sc31;
  Float_t t_sc41_online  = -99999;
  Float_t t_sc42_online  = -99999;
  Float_t t_sc43_online  = -99999;
  Float_t t_sc31_online  = -99999;
  Float_t tof4143_online = -99999;
  Float_t tof3141_online = -99999;
  Float_t de_sc41_online = -9999;
  Float_t de_sc42_online = -9999;
  Float_t de_sc43_online = -9999;
  Float_t de_sc31_online = -9999;
  std::array<Int_t, 3> tdc_trig;
  std::array<Int_t, 3> tdc_s2tref;
  ClassDef(S4TQ, 1);
};

class MWDCHit : public TObject
{
  public:
  MWDCHit();
  ~MWDCHit();
  Int_t i_ctdc;
  Int_t ch_ctdc;
  Int_t i_plane;
  Int_t i_wire;
  Int_t t_leading;
  Int_t t_trailing;
  ClassDef(MWDCHit, 1);
};

class FRSTPC : public TObject
{
  template <class T, std::size_t N, std::size_t M>
  using array2 = std::array<std::array<T, M>, N>; // T val[N][M]
  public:
  FRSTPC();
  ~FRSTPC();
  void Init();
  ULong64_t ts  = 0;            // from unpack
  Int_t eventnr = -999;         // from unpack
  Bool_t all_ok = false;        // from unpack
  array2<Int_t, 7, 4> tpc_dt_s; // the rest is copy from TFRSCalibrEvent
  array2<Int_t, 7, 2> tpc_lt_s;
  array2<Int_t, 7, 2> tpc_rt_s;
  array2<Int_t, 7, 2> tpc_xraw;
  array2<Int_t, 7, 4> tpc_yraw;
  array2<Int_t, 7, 4> tpc_csum;
  std::array<Float_t, 7> tpc_x;
  std::array<Float_t, 7> tpc_y;
  array2<Bool_t, 7, 4> b_tpc_csum;
  std::array<Bool_t, 7> b_tpc_xy;
  std::array<Float_t, 7> tpc_de;
  std::array<Float_t, 7> tpc_dx12;
  std::array<Bool_t, 7> b_tpc_de;
  std::array<Bool_t, 8> b_tpc_timeref;
  std::array<Int_t, 8> tpc_timeref_s;
  ClassDef(FRSTPC, 1);
};

class S4MWDC : public TObject
{
  public:
  S4MWDC();
  ~S4MWDC();
  void Init();
  ULong64_t ts;
  Int_t eventnr;
  Bool_t all_ok;
  std::vector<MWDCHit> mwdchit;
  std::vector<Int_t> tref;
  Float_t x_online = -9999;
  Float_t a_online = -9999;
  Float_t y_online = -9999;
  Float_t b_online = -9999;
  Float_t chi2_online = -9999;
  std::array<Float_t, 16> res_online;
  std::array<Float_t, 16> dl_online;
  std::array<Int_t, 16> t_online;
  std::array<Int_t, 16> wire_online;
  ClassDef(S4MWDC, 1);
};

class S4WFD : public TObject
{
  template <class T, std::size_t N, std::size_t M>
  using array2 = std::array<std::array<T, M>, N>; // T val[N][M]
  public:
  S4WFD();
  ~S4WFD();
  void Init();
  ULong64_t ts = 0;
  Int_t eventnr = -999;
  std::vector<std::vector<Short_t>> v1742; // v1742[36][]
  std::array<Int_t, 36> v1742_min;
  std::array<Int_t, 36> v1742_max;
  std::array<Int_t, 36> v1742_min_index;
  std::array<Int_t, 36> v1742_max_index;
  std::array<Int_t, 36> v1742_qdc;
  std::array<UInt_t, 4> v1742_ts     = {0, 0, 0, 0};
  std::array<Int_t, 4> v1742_ts_diff = {0, 0, 0, 0};
  Int_t v1742_ec      = -9999;
  Int_t v1742_ec_diff = -9999;
  Int_t v1742_gr_mask = -9999;
  std::array<Int_t, 4> v1742_size;
  std::array<Int_t, 4> v1742_freq;
  std::array<Int_t, 4> v1742_trn;
  std::array<Int_t, 4> v1742_start_index;
  Bool_t v1742_ok = false;
  Int_t v1290n_ec = -9999;
  Int_t v1290n_ec_diff = -9999;
  std::vector<std::vector<Int_t>> v1290n;
  Bool_t v1290n_ok = false;
  Bool_t all_ok    = false;
  std::array<Bool_t, 2> has_wfd_sc31;
  std::array<Int_t, 2> sc31_wfd_qdc;
  std::array<Int_t, 2> sc31_wfd_max;
  std::array<Int_t, 2> sc31_wfd_min;
  std::array<Int_t, 2> sc31_wfd_max_index;
  std::array<Int_t, 2> sc31_wfd_min_index;
  array2<Bool_t, 3, 2> has_wfd_sc41;
  array2<Int_t, 3, 2> sc41_wfd_qdc;
  array2<Int_t, 3, 2> sc41_wfd_max;
  array2<Int_t, 3, 2> sc41_wfd_min;
  array2<Int_t, 3, 2> sc41_wfd_max_index;
  array2<Int_t, 3, 2> sc41_wfd_min_index;
  std::array<Bool_t, 2> has_wfd_sc42;
  std::array<Int_t, 2> sc42_wfd_qdc;
  std::array<Int_t, 2> sc42_wfd_max;
  std::array<Int_t, 2> sc42_wfd_min;
  std::array<Int_t, 2> sc42_wfd_max_index;
  std::array<Int_t, 2> sc42_wfd_min_index;
  array2<Bool_t, 4, 2> has_wfd_sc43;
  array2<Int_t, 4, 2> sc43_wfd_qdc;
  array2<Int_t, 4, 2> sc43_wfd_max;
  array2<Int_t, 4, 2> sc43_wfd_min;
  array2<Int_t, 4, 2> sc43_wfd_max_index;
  array2<Int_t, 4, 2> sc43_wfd_min_index;

  void Draw(int ich);
  void ShowWFDExistCh();
  ClassDef(S4WFD, 1);
};

class S2TQ1 : public TObject
{
  template <class T, std::size_t N, std::size_t M>
  using array2 = std::array<std::array<T, M>, N>; // T val[N][M]
  public:
  S2TQ1();
  ~S2TQ1();
  void Init();
  ULong64_t ts = 0;
  Int_t eventnr;
  Short_t v792[3][32];
  std::vector<std::vector<std::vector<Int_t>>> v1290; // [3][32][]
  Int_t v1290_nhit[3][32];
  UInt_t v830[32];
  Int_t v830diff[32];
  Int_t v792_ec[3];
  Int_t v1290_ec[3];
  Int_t v792_ec_diff[3];
  Int_t v1290_ec_diff[3];
  std::vector<std::vector<Int_t>> v1290_tref; // [3][] : [mod][nhit]
  Int_t v830_ec   = -9999;
  Int_t v830_ec_diff   = -9999;
  std::array<Bool_t, 3> v792_ok  = {false, false, false};
  std::array<Bool_t, 3> v1290_ok = {false, false, false};
  Bool_t v830_ok  = false;
  Bool_t all_ok   = false;
  std::vector<std::vector<std::vector<Int_t>>> tdc_psb; // [46][2][]
  array2<Int_t, 46, 2> nhit_psb;
  std::array<Int_t, 4> tdc_dummy;
  Int_t nhit_dummy = -1;
  array2<Short_t, 46, 2> qdc_psb;
  Short_t qdc_dummy = -9999;
  std::array<Bool_t, 46> has_hit_psb_online;
  array2<Bool_t, 46, 2> has_hit_psb_online_ud;
  std::array<Float_t, 46> t_psb_online;
  std::array<Float_t, 46> timedif_psb_online;
  std::array<Float_t, 46> z_psb_online;
  std::array<Float_t, 46> de_psb_online;
  ClassDef(S2TQ1, 1);
};

class MDCHit : public TObject
{
  public:
  MDCHit();
  ~MDCHit();
  void Init();
  Int_t i_ctdc = -1;
  Int_t ch_ctdc = -1;
  Int_t i_layer = -1;
  Int_t i_wire = -1;
  Int_t t_leading = -9999;
  Int_t t_trailing = -9999;
  ClassDef(MDCHit, 1);
};

class S2MDC : public TObject
{
  public:
  S2MDC();
  ~S2MDC();
  void Init();
  ULong64_t ts = 0;
  Int_t eventnr = -999;
  std::vector<MDCHit> mdchit;
  Bool_t all_ok = false;
  std::vector<Int_t> tref;
  ClassDef(S2MDC, 1);
};

class S2WFD123 : public TObject
{
  template <class T, std::size_t N, std::size_t M>
  using array2 = std::array<std::array<T, M>, N>; // T val[N][M]
  public:
  S2WFD123();
  ~S2WFD123();
  void Init();
  std::array<ULong64_t, 3> ts;
  std::array<Int_t, 3> eventnr = {-999, -999, -999};
  std::vector<std::vector<std::vector<Short_t>>> v1742; // v1742[3][36][]
  array2<Int_t, 3, 36> v1742_min;
  array2<Int_t, 3, 36> v1742_max;
  array2<Int_t, 3, 36> v1742_min_index;
  array2<Int_t, 3, 36> v1742_max_index;
  array2<Int_t, 3, 36> v1742_qdc;
  array2<UInt_t, 3, 4> v1742_ts;
  array2<Int_t, 3, 4> v1742_ts_diff;
  std::array<Int_t, 3> v1742_ec;
  std::array<Int_t, 3> v1742_ec_diff;
  std::array<Int_t, 3> v1742_gr_mask;
  array2<Int_t, 3, 4> v1742_size;
  array2<Int_t, 3, 4> v1742_freq;
  array2<Int_t, 3, 4> v1742_trn;
  array2<Int_t, 3, 4> v1742_start_index;
  std::array<Bool_t, 3> v1742_ok;
  Int_t v1290n_ec = -9999;
  Int_t v1290n_ec_diff = -9999;
  std::vector<std::vector<Int_t>> v1290n; // v1290n[16][]
  Bool_t v1290n_ok = false;
  Bool_t all_ok = false;
  array2<Bool_t, 46, 2> has_wfd_psb;
  array2<Int_t, 46, 2> psb_wfd_qdc;
  array2<Int_t, 46, 2> psb_wfd_max;
  array2<Int_t, 46, 2> psb_wfd_min;
  array2<Int_t, 46, 2> psb_wfd_max_index;
  array2<Int_t, 46, 2> psb_wfd_min_index;
  void Draw(int imod, int ich);
  void ShowWFDExistCh();
  ClassDef(S2WFD123, 1);
};

class S2TQ2 : public TObject
{
  template <class T, std::size_t N, std::size_t M>
  using array2 = std::array<std::array<T, M>, N>; // T val[N][M]
  template <class T, std::size_t N1, std::size_t N2, std::size_t N3>
  using array3 = std::array<array2<T, N2, N3>, N1>;

  public:
  S2TQ2();
  ~S2TQ2();
  void Init();
  ULong64_t ts = 0;
  Int_t eventnr = -999;
  array2<Short_t, 2, 32> v792;
  // Int_t v1290[5][32][4];
  std::vector<std::vector<std::vector<Int_t>>> v1290;
  array2<Int_t, 5, 32> v1290_nhit; // orginally v1290_nhit[5][32][4]
  std::vector<std::vector<Int_t>> v1290_tref; // [5][] : [mod][nhit]
  std::array<Int_t, 2> v792_ec;
  std::array<Int_t, 2> v792_ec_diff;
  std::array<Int_t, 5> v1290_ec;
  std::array<Int_t, 5> v1290_ec_diff;
  std::array<Bool_t, 2> v792_ok;
  std::array<Bool_t, 5> v1290_ok;
  Bool_t all_ok = false;
  array2<Int_t, 44, 4> tdc_psfe;
  array2<Int_t, 38, 4> tdc_psbe;
  array3<Int_t, 28, 2, 4> tdc_t0;
  std::array<Int_t, 44> nhit_psfe;
  std::array<Int_t, 38> nhit_psbe;
  array2<Int_t, 28, 2> nhit_t0;
  std::array<Short_t, 44> qdc_psfe;
  std::array<Short_t, 38> qdc_psbe;
  array2<Short_t, 28, 2> qdc_t0;
  std::array<Bool_t, 44> has_hit_psfe_online;
  std::array<Float_t, 44> t_psfe_online;
  std::array<Bool_t, 38> has_hit_psbe_online;
  std::array<Float_t, 38> t_psbe_online;
  std::array<Bool_t, 28> has_hit_t0_online;
  std::array<Float_t, 28> t_t0_online;
  std::array<Float_t, 28> de_t0_online;
  ClassDef(S2TQ2, 1);
};

class S2WFD45 : public TObject
{
  template <class T, std::size_t N, std::size_t M>
  using array2 = std::array<std::array<T, M>, N>; // T val[N][M]

  public:
  S2WFD45();
  ~S2WFD45();
  void Init();
  ULong64_t ts;
  Int_t eventnr;
  std::vector<std::vector<std::vector<Short_t>>> v1742; // v[2][36][]
  array2<Int_t, 2, 36> v1742_min;
  array2<Int_t, 2, 36> v1742_max;
  array2<Int_t, 2, 36> v1742_min_index;
  array2<Int_t, 2, 36> v1742_max_index;
  array2<Int_t, 2, 36> v1742_qdc;
  array2<UInt_t, 2, 4> v1742_ts;
  std::array<Int_t, 2> v1742_ec;
  array2<Int_t, 2, 4> v1742_ts_diff;
  std::array<Int_t, 2> v1742_ec_diff;
  std::array<Int_t, 2> v1742_gr_mask;   // 3->2??
  array2<Int_t, 2, 4> v1742_size;
  array2<Int_t, 2, 4> v1742_freq;
  array2<Int_t, 2, 4> v1742_trn;
  array2<Int_t, 2, 4> v1742_start_index;
  std::array<Bool_t, 2> v1742_ok;
  Bool_t all_ok = false;
  std::array<Bool_t, 44> has_wfd_psfe;
  std::array<Int_t, 44> psfe_wfd_qdc;
  std::array<Int_t, 44> psfe_wfd_max;
  std::array<Int_t, 44> psfe_wfd_min;
  std::array<Int_t, 44> psfe_wfd_max_index;
  std::array<Int_t, 44> psfe_wfd_min_index;
  std::array<Bool_t, 38> has_wfd_psbe;
  std::array<Int_t, 38> psbe_wfd_qdc;
  std::array<Int_t, 38> psbe_wfd_max;
  std::array<Int_t, 38> psbe_wfd_min;
  std::array<Int_t, 38> psbe_wfd_max_index;
  std::array<Int_t, 38> psbe_wfd_min_index;
  void Draw(int imod, int ich);
  void ShowWFDExistCh();
  ClassDef(S2WFD45, 1);
};

class FiberHit : public TObject
{
  public:
  FiberHit();
  ~FiberHit();
  void Init();
  Int_t i_sfp      = -1;
  Int_t i_ctdc     = -1;
  Int_t i_ch       = -1;
  Int_t i_detector = -1;
  Int_t i_layer    = -1;
  Int_t i_fiber    = -1;
  Int_t t_leading  = -9999;
  Int_t t_trailing = -9999;
  ClassDef(FiberHit, 1);
};

class FiberTrack : public TObject
{
  public:
  FiberTrack();
  ~FiberTrack();
  void Init();
  Float_t x    = -9999;
  Float_t y    = -9999;
  Float_t a    = -9999;
  Float_t b    = -9999;
  Float_t chi2 = -9999;
  ClassDef(FiberTrack, 1);
};

class S2Fiber : public TObject
{
  public:
  S2Fiber();
  ~S2Fiber();
  void Init();
  ULong64_t ts = 0;
  Int_t eventnr = -999;
  std::vector<FiberHit> fiberhit;
  Bool_t all_ok = false;
  std::array<Int_t, 8> tref;
  std::vector<FiberTrack> tr_uft12_online;
  std::vector<FiberTrack> tr_mft12_online;
  std::vector<FiberTrack> tr_dft12_online;
  std::vector<FiberTrack> tr_uft3mft12_online;
  std::vector<FiberTrack> tr_uft123dft12_online;
  ClassDef(S2Fiber, 1);
};

class CsIHit : public TObject
{
  public:
  CsIHit();
  ~CsIHit();
  void Init();
  Int_t i_pc    = -1;
  Int_t i_sfp   = -1;
  Int_t i_febex = -1;
  Int_t i_ch    = -1;
  Int_t i_lr    = -1;
  Int_t i_theta = -1;
  Int_t i_phi   = -1;
  std::vector<Short_t> adc;
  Float_t trapez_online = -9999;
  ClassDef(CsIHit, 1);
};

class S2CsI : public TObject
{
  template <class T, std::size_t N, std::size_t M>
  using array2 = std::array<std::array<T, M>, N>;

  public:
  S2CsI();
  ~S2CsI();
  void Init();
  std::array<ULong64_t, 2> ts;
  std::array<Int_t, 2> eventnr;
  std::vector<CsIHit> csihit;
  std::vector<std::vector<Short_t>> adc_trig; // adc_trig[2][]
  std::array<Bool_t, 2> all_ok;
  array2<Int_t, 2, 4> polarity;
  array2<Int_t, 2, 4> n_slave;
  void Draw(int ihit);
  void DrawADC(int ipc);
  ClassDef(S2CsI, 1);
};

namespace rse {
  template <class T, class U>
  void InitCArray(T &arr, U init_val) {
    for (auto &r : arr) {
      r = init_val;
    }
  }

  template <class T, class U>
  void InitCArray2d(T &arr, U init_val) {
    for (auto &r : arr) {
      for (auto &element : r) {
        element = init_val;
      }
    }
  }
}
