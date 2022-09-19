#include "WASAUnpackBranch.hh"
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <memory>
#include <iostream>

ClassImp(MWDCHit);
ClassImp(FRSTPC);
ClassImp(S4MWDC);
ClassImp(S4TQ);
ClassImp(S4WFD);
ClassImp(S2TQ1);
ClassImp(MDCHit);
ClassImp(S2MDC);
ClassImp(S2TQ2);
ClassImp(S2WFD123);
ClassImp(S2WFD45);
ClassImp(FiberHit);
ClassImp(FiberTrack);
ClassImp(S2Fiber);
ClassImp(CsIHit);
ClassImp(S2CsI);
ClassImp(LMDHeader);

LMDHeader::LMDHeader() : lmd_filename(""), lmd_open_time("") {}

LMDHeader::~LMDHeader() {}

////// ---- MWDCHit : Sub Class ---- ////
MWDCHit::MWDCHit() : i_ctdc(-1), ch_ctdc(-1), i_plane(-1), i_wire(-1), t_leading(-9999), t_trailing(-9999) {}

MWDCHit::~MWDCHit() {}

////// ---- S4MWDC ---- ////
S4MWDC::S4MWDC() { Init(); }

void S4MWDC::Init() {
  ts     = 0;
  eventnr = -999;
  all_ok = false;
  mwdchit.clear();
  tref.clear();
  x_online    = -9999;
  a_online    = -9999;
  y_online    = -9999;
  b_online    = -9999;
  chi2_online = -9999;
  rse::InitCArray(res_online, -9999);
  rse::InitCArray(dl_online, -9999);
  rse::InitCArray(t_online, -9999);
  rse::InitCArray(wire_online, -9999);
}

S4MWDC::~S4MWDC() {}

//// ---- S4TQ1 ---- ////
S4TQ::S4TQ() : mtdc(34) { Init(); }

void S4TQ::Init() {
  ts             = 0;
  eventnr        = -999;
  mtdc_ec        = -9999;
  v792_ec        = -9999;
  v775_ec        = -9999;
  v830_ec        = -9999;
  mtdc_ec_diff   = -9999;
  v792_ec_diff   = -9999;
  v775_ec_diff   = -9999;
  v830_ec_diff   = -9999;
  v792_ok        = false;
  v775_ok        = false;
  mtdc_ok        = false;
  v830_ok        = false;
  all_ok         = false;
  t_sc41_online  = -99999;
  t_sc42_online  = -99999;
  t_sc43_online  = -99999;
  t_sc31_online  = -99999;
  tof4143_online = -99999;
  tof3141_online = -99999;
  de_sc41_online = -9999;
  de_sc42_online = -9999;
  de_sc43_online = -9999;
  de_sc31_online = -9999;
  rse::InitCArray(v775, 0);
  rse::InitCArray(trig, false);
  rse::InitCArray(v830, 0);
  rse::InitCArray(v830diff, -9999);
  rse::InitCArray(v792, -9999);
  rse::InitCArray(tdc_sc42, -9999);
  rse::InitCArray(tdc_sc31, -9999);
  rse::InitCArray(qdc_sc42, -9999);
  rse::InitCArray(qdc_sc42_lowgain, -9999);
  rse::InitCArray(qdc_sc31, -9999);
  rse::InitCArray(tdc_trig, -9999);
  rse::InitCArray(tdc_s2tref, -9999);
  rse::InitCArray(mtdc_nhit, 0);
  rse::InitCArray2d(tdc_sc41, -9999);
  rse::InitCArray2d(tdc_sc43, -9999);
  rse::InitCArray2d(qdc_sc41, -9999);
  rse::InitCArray2d(qdc_sc43, -9999);
  for (auto &r : mtdc) {
    r.clear();
  }
}

S4TQ::~S4TQ() {}

//// ---- S4WFD ---- ////
S4WFD::S4WFD() : v1742(36), v1290n(16) {
  Init();
  v1742.shrink_to_fit();
  v1290n.shrink_to_fit();
}

void S4WFD::ShowWFDExistCh() {
  std::cout << "---- S2WFD45: Waveform existing channel ---- \n";
  int count    = 0;
  int ch_count = 0;
  for (int ich = 0; ich < (int)v1742.size(); ++ich) {
      if (v1742[ich].size() != 0) {
        std::cout << "ch = " << ich << "\n";
        ++count;
      }
      ++ch_count;
  }
  std::cout << "-- Num. of channels: " << count << "/" << ch_count << std::endl;
}

void S4WFD::Init() {
  ts            = 0;
  eventnr       = -999;
  rse::InitCArray(v1742_ts, 0);
  rse::InitCArray(v1742_ts_diff, -9999);
  v1742_ec      = -9999;
  v1742_ec_diff = -9999;
  v1742_gr_mask = -9999;
  v1742_ok      = false;
  v1290n_ec     = -9999;
  v1290n_ec_diff= -9999;
  v1290n_ok     = false;
  all_ok        = false;
  for (auto &r : v1742) {
    r.clear();
  }
  for (auto &r : v1290n) {
    r.clear();
  }
  rse::InitCArray(v1742_size, -9999);
  rse::InitCArray(v1742_freq, -9999);
  rse::InitCArray(v1742_trn, -9999);
  rse::InitCArray(v1742_start_index, -9999);
  rse::InitCArray(v1742_min, -9999);
  rse::InitCArray(v1742_max, 9999);
  rse::InitCArray(v1742_min_index, -1);
  rse::InitCArray(v1742_max_index, -1);
  rse::InitCArray(v1742_qdc, -1);
  rse::InitCArray(has_wfd_sc42, false);
  rse::InitCArray(has_wfd_sc31, false);
  rse::InitCArray2d(has_wfd_sc41, false);
  rse::InitCArray2d(has_wfd_sc43, false);

  rse::InitCArray(sc31_wfd_qdc, -1);
  rse::InitCArray(sc31_wfd_max, 9999);
  rse::InitCArray(sc31_wfd_min, -9999);
  rse::InitCArray(sc31_wfd_max_index, -1);
  rse::InitCArray(sc31_wfd_min_index, -1);
  rse::InitCArray(sc42_wfd_qdc, -1);
  rse::InitCArray(sc42_wfd_max, 9999);
  rse::InitCArray(sc42_wfd_min, -9999);
  rse::InitCArray(sc42_wfd_max_index, -1);
  rse::InitCArray(sc42_wfd_min_index, -1);
  rse::InitCArray2d(sc41_wfd_qdc, -1);
  rse::InitCArray2d(sc41_wfd_max, 9999);
  rse::InitCArray2d(sc41_wfd_min, -9999);
  rse::InitCArray2d(sc41_wfd_max_index, -1);
  rse::InitCArray2d(sc41_wfd_min_index, -1);
  rse::InitCArray2d(sc43_wfd_qdc, -1);
  rse::InitCArray2d(sc43_wfd_max, 9999);
  rse::InitCArray2d(sc43_wfd_min, -9999);
  rse::InitCArray2d(sc43_wfd_max_index, -1);
  rse::InitCArray2d(sc43_wfd_min_index, -1);

}

void S4WFD::DrawHist(int ich) {
  static std::shared_ptr<TCanvas> c;
  c.reset();
  c = std::make_shared<TCanvas>("c_s4wfd", Form("s4wfd123[%d]", ich), 700, 500);
  static std::shared_ptr<TGraph> graph;
  graph.reset();

  graph = std::make_shared<TGraph>();
  graph->SetName(Form("s4wfd_ch%d", ich));
  for (int ipoint = 0; ipoint < (int)v1742[ich].size(); ++ipoint) {
    graph->SetPoint(ipoint, ipoint, v1742[ich][ipoint]);
  }
  graph->GetXaxis()->SetTitle("sample");
  graph->SetMarkerSize(0.5);
  graph->SetMarkerStyle(20);
  graph->Draw("AP");
}

S4WFD::~S4WFD() {}

//// ---- S2TQ1 ---- ////
S2TQ1::S2TQ1() : v1290(3, std::vector<std::vector<int>>(32)), v1290_tref(3), tdc_psb(46, std::vector<std::vector<int>>(2)) {
  Init();
}

void S2TQ1::Init() {
  ts         = 0;
  eventnr    = -999;
  v830_ec    = -9999;
  v830_ec_diff    = -9999;
  v792_ok[0] = false;
  v792_ok[1] = false;
  v792_ok[2] = false;
  v1290_ok[0]= false;
  v1290_ok[1]= false;
  v1290_ok[2]= false;
  v830_ok    = false;
  all_ok     = false;
  nhit_dummy = -1;
  qdc_dummy  = -9999;

  rse::InitCArray(v830, 0);
  rse::InitCArray(v830diff, -1);
  rse::InitCArray(v792_ec, -9999);
  rse::InitCArray(v1290_ec, -9999);
  rse::InitCArray(v792_ec_diff, -9999);
  rse::InitCArray(v1290_ec_diff, -9999);
  rse::InitCArray(tdc_dummy, -9999);
  rse::InitCArray(has_hit_psb_online, false);
  rse::InitCArray2d(has_hit_psb_online_ud, false);
  rse::InitCArray(t_psb_online, false);
  rse::InitCArray(timedif_psb_online, -9999);
  rse::InitCArray(z_psb_online, -9999);
  rse::InitCArray(de_psb_online, -9999);
  rse::InitCArray2d(v792, -9999);
  rse::InitCArray2d(v1290_nhit, 0);
  rse::InitCArray2d(nhit_psb, 0);
  rse::InitCArray2d(qdc_psb, -9999);
  for (auto &mod : v1290) {
    for (auto &hit_ch : mod) {
      hit_ch.clear();
    }
  }
  for (auto &mod : tdc_psb) {
    for (auto &multitdc : mod) {
      multitdc.clear();
    }
  }
  for(auto &tref : v1290_tref) {
    tref.clear();
  }
}

S2TQ1::~S2TQ1() {}

//// ---- MDCHit : Sub Class ---- ////
MDCHit::MDCHit() : i_ctdc(-1), ch_ctdc(-1), i_layer(-1), i_wire(-1), t_leading(-9999), t_trailing(-9999) {
}

void MDCHit::Init() {
  i_ctdc     = -1;
  ch_ctdc    = -1;
  i_layer    = -1;
  i_wire     = -1;
  t_leading  = -9999;
  t_trailing = -9999;
}

MDCHit::~MDCHit() {}

//// ---- S2MDC ---- ////
S2MDC::S2MDC() {
  Init();
}

void S2MDC::Init() {
  ts = 0;
  eventnr = -999;
  mdchit.clear();
  all_ok = false;
  tref.clear();
}

S2MDC::~S2MDC(){};


//// ---- FRS TPC---- ////
FRSTPC::FRSTPC() {
  Init();
}

void FRSTPC::Init() {
  ts = 0;
  eventnr = -999;
  all_ok = false;

  for(int i=0; i<8; i++){
    tpc_timeref_s[i] = -999;
    b_tpc_timeref[i] = kFALSE;
  }
  for (int i = 0; i<7;i++){
    tpc_x[i] = 9999.9;
    tpc_y[i] = 9999.9;
    tpc_de[i] = 0.0;
    tpc_dx12[i]=-9999.;
    b_tpc_xy[i] = kFALSE;
    for (int j=0;j<4;j++){
      tpc_csum[i][j] =0;
      b_tpc_csum[i][j] = kFALSE;
      b_tpc_de[i]= kFALSE;
      tpc_dt_s[i][j] = -999;
      tpc_yraw[i][j] = -9999999;
    }
    for (int j=0;j<2;j++){
      tpc_lt_s[i][j] = -999;
      tpc_rt_s[i][j] = -999;
      tpc_xraw[i][j] = -9999999;
    }
  }
}

FRSTPC::~FRSTPC(){};


//// ---- S2WFD123 ---- ;;
S2WFD123::S2WFD123() : v1742(3, std::vector<std::vector<int16_t>>(36)), v1290n(16) {
  v1742.shrink_to_fit();
  for (auto &r : v1742) {
    r.shrink_to_fit();
  }
  v1290n.shrink_to_fit();
  for (auto &r : v1290n) {
    r.shrink_to_fit();
  }
  Init();
}

void S2WFD123::DrawHist(int imod, int ich) {
  static std::shared_ptr<TCanvas> c;
  c.reset();
  c = std::make_shared<TCanvas>("c_wfd123", Form("s2wfd123[%d][%d]", imod, ich), 700, 500);
  static std::shared_ptr<TGraph> graph;
  graph.reset();

  graph = std::make_shared<TGraph>();
  graph->SetName(Form("s2wfd123_%d_%d", imod, ich));
  for (int ipoint = 0; ipoint < (int)v1742[imod][ich].size(); ++ipoint) {
    graph->SetPoint(ipoint, ipoint, v1742[imod][ich][ipoint]);
  }
  graph->GetXaxis()->SetTitle("sample");
  graph->SetMarkerSize(0.5);
  graph->SetMarkerStyle(20);
  graph->Draw("AP");
}

void S2WFD123::ShowWFDExistCh() {
  std::cout << "---- S2WFD123: Waveform existing channel ---- \n";
  int count    = 0;
  int ch_count = 0;
  for (int imod = 0; imod < (int)v1742.size(); ++imod) {
    for (int ich = 0; ich < (int)v1742[imod].size(); ++ich) {
      if (v1742[imod][ich].size() != 0) {
        std::cout << "(imod, ich) = (" << imod << ", " << ich << ")\n";
        ++count;
      }
      ++ch_count;
    }
  }
  std::cout << "-- Num. of channels: " << count << "/" << ch_count << std::endl;
}

void S2WFD123::Init() {
  v1290n_ec = -9999;
  v1290n_ec_diff = -9999;
  v1290n_ok = false;
  all_ok    = false;
  rse::InitCArray(ts, 0);
  rse::InitCArray(eventnr, -999);
  rse::InitCArray2d(v1742_ts, 0);
  rse::InitCArray(v1742_ec, -9999);
  rse::InitCArray2d(v1742_ts_diff, -9999);
  rse::InitCArray(v1742_ec_diff, -9999);
  rse::InitCArray(v1742_gr_mask, -1);
  rse::InitCArray(v1742_ok, false);
  rse::InitCArray2d(v1742_min, -9999);
  rse::InitCArray2d(v1742_max, 9999);
  rse::InitCArray2d(v1742_min_index, -1);
  rse::InitCArray2d(v1742_max_index, -1);
  rse::InitCArray2d(v1742_qdc, -1);
  rse::InitCArray2d(v1742_size, 0);
  rse::InitCArray2d(v1742_freq, 0);
  rse::InitCArray2d(v1742_trn, -9999);
  rse::InitCArray2d(v1742_start_index, -1);
  rse::InitCArray2d(has_wfd_psb, false);
  rse::InitCArray2d(psb_wfd_qdc, -1);
  rse::InitCArray2d(psb_wfd_max, 9999);
  rse::InitCArray2d(psb_wfd_min, -9999);
  rse::InitCArray2d(psb_wfd_max_index, -1);
  rse::InitCArray2d(psb_wfd_min_index, -1);

  for (auto &r : v1290n) {
    r.clear();
  }
  for (auto &r : v1742) {
    for (auto &s : r) {
      s.clear();
    }
  }
}

S2WFD123::~S2WFD123() {}

//// ---- S2TQ2 ---- ////
S2TQ2::S2TQ2() : v1290(5, std::vector<std::vector<int>>(32)), v1290_tref(5) { Init(); }

void S2TQ2::Init() {
  ts     = 0;
  eventnr = -999;
  all_ok = false;
  rse::InitCArray(v792_ec, -9999);
  rse::InitCArray(v1290_ec, -9999);
  rse::InitCArray(v792_ec_diff, -9999);
  rse::InitCArray(v1290_ec_diff, -9999);
  rse::InitCArray(v792_ok, false);
  rse::InitCArray(v1290_ok, false);
  rse::InitCArray(nhit_psfe, 0);
  rse::InitCArray(nhit_psbe, 0);
  rse::InitCArray(qdc_psfe, -9999);
  rse::InitCArray(qdc_psbe, -9999);
  rse::InitCArray(has_hit_psfe_online, false);
  rse::InitCArray(t_psfe_online, -9999);
  rse::InitCArray(has_hit_psbe_online, false);
  rse::InitCArray(t_psbe_online, -9999);
  rse::InitCArray(has_hit_t0_online, false);
  rse::InitCArray(t_t0_online, -9999);
  rse::InitCArray(de_t0_online, -9999);
  rse::InitCArray2d(v792, -9999);
  rse::InitCArray2d(v1290_nhit, -9999); //corrected definition of s2tq2 v1290 nhit
  rse::InitCArray2d(tdc_psfe, -9999);
  rse::InitCArray2d(tdc_psbe, -9999);
  rse::InitCArray2d(nhit_t0, 0);
  rse::InitCArray2d(qdc_t0, -9999);
  for (auto &r : v1290) {
    for (auto &s : r) {
      s.clear();
    }
  }
  for (auto &r : tdc_t0) {
    rse::InitCArray2d(r, -9999);
  }
  for(auto &r : v1290_tref) {
    r.clear();
  }
}

S2TQ2::~S2TQ2() {}

//// ---- S2WFD45 ---- ////
S2WFD45::S2WFD45() : v1742(2, std::vector<std::vector<int16_t>>(36)) {
  v1742.shrink_to_fit();
  for (auto &r : v1742)
    r.shrink_to_fit();
  Init();
}

void S2WFD45::ShowWFDExistCh() {
  std::cout << "---- S2WFD45: Waveform existing channel ---- \n";
  int count    = 0;
  int ch_count = 0;
  for (int imod = 0; imod < (int)v1742.size(); ++imod) {
    for (int ich = 0; ich < (int)v1742[imod].size(); ++ich) {
      if (v1742[imod][ich].size() != 0) {
        std::cout << "(imod, ich) = (" << imod << ", " << ich << ")\n";
        ++count;
      }
      ++ch_count;
    }
  }
  std::cout << "-- Num. of channels: " << count << "/" << ch_count << std::endl;
}

void S2WFD45::Init() {
  ts = 0;
  eventnr = -999;
  rse::InitCArray2d(v1742_ts, 0);
  rse::InitCArray(v1742_ec, -9999);
  rse::InitCArray2d(v1742_ts_diff, -9999);
  rse::InitCArray(v1742_ec_diff, -9999);
  rse::InitCArray(v1742_gr_mask, -1);
  rse::InitCArray(v1742_ok, false);
  rse::InitCArray(has_wfd_psfe, false);
  rse::InitCArray(has_wfd_psbe, false);
  rse::InitCArray2d(v1742_min, -9999);
  rse::InitCArray2d(v1742_max, 9999);
  rse::InitCArray2d(v1742_min_index, -1);
  rse::InitCArray2d(v1742_max_index, -1);
  rse::InitCArray2d(v1742_qdc, -1);
  rse::InitCArray2d(v1742_size, -1);
  rse::InitCArray2d(v1742_freq, -9999);
  rse::InitCArray2d(v1742_trn, -1);
  rse::InitCArray2d(v1742_start_index, -1);
  for (auto &r : v1742) {
    for (auto &s : r) {
      s.clear();
    }
  }
  rse::InitCArray(psbe_wfd_qdc, -1);
  rse::InitCArray(psbe_wfd_max, 9999);
  rse::InitCArray(psbe_wfd_min, -9999);
  rse::InitCArray(psbe_wfd_max_index, -1);
  rse::InitCArray(psbe_wfd_min_index, -1);
  rse::InitCArray(psfe_wfd_qdc, -1);
  rse::InitCArray(psfe_wfd_max, 9999);
  rse::InitCArray(psfe_wfd_min, -9999);
  rse::InitCArray(psfe_wfd_max_index, -1);
  rse::InitCArray(psfe_wfd_min_index, -1);
}

void S2WFD45::DrawHist(int imod, int ich) {
  static std::shared_ptr<TCanvas> c;
  c.reset();
  c = std::make_shared<TCanvas>("c_s2wfd45", Form("s2wfd45[%d][%d]", imod, ich), 700, 500);
  static std::shared_ptr<TGraph> graph;
  graph.reset();

  graph = std::make_shared<TGraph>();
  graph->SetName(Form("s2wfd45_%d_%d", imod, ich));
  for (int ipoint = 0; ipoint < (int)v1742[imod][ich].size(); ++ipoint) {
    graph->SetPoint(ipoint, ipoint, v1742[imod][ich][ipoint]);
  }
  graph->GetXaxis()->SetTitle("sample");
  graph->SetMarkerSize(0.5);
  graph->SetMarkerStyle(20);
  graph->Draw("AP");
}

S2WFD45::~S2WFD45() {
}

//// ---- FiberHit : Sub Class ---- ////
FiberHit::FiberHit() {}
FiberHit::~FiberHit() {}

void FiberHit::Init() {
  i_sfp      = -1;
  i_ctdc     = -1;
  i_ch       = -1;
  i_detector = -1;
  i_layer    = -1;
  i_fiber    = -1;
  t_leading  = -9999;
  t_trailing = -9999;
}

//// ---- FiberTrack : Sub Class ---- ////
FiberTrack::FiberTrack() {
}

void FiberTrack::Init() {
  x    = -9999;
  y    = -9999;
  a    = -9999;
  b    = -9999;
  chi2 = -9999;
}

FiberTrack::~FiberTrack() {
}

//// ---- S2Fiber ---- ////
S2Fiber::S2Fiber() {
  Init();
}

void S2Fiber::Init() {
  ts      = 0;
  eventnr = -999;
  all_ok  = false;
  rse::InitCArray(tref, -9999);
  fiberhit.clear();
  tr_uft12_online.clear();
  tr_mft12_online.clear();
  tr_dft12_online.clear();
  tr_uft3mft12_online.clear();
  tr_uft123dft12_online.clear();
}

S2Fiber::~S2Fiber() {}

//// ---- CsIHit : Sub Class ---- ////
CsIHit::CsIHit() { Init(); }

void CsIHit::Init() {
  i_pc    = -1;
  i_sfp   = -1;
  i_febex = -1;
  i_ch    = -1;
  i_theta = -1;
  i_phi   = -1;
  adc.clear();
  trapez_online = -9999;
}

CsIHit::~CsIHit() {}

//// ---- S2CsI ---- ////
S2CsI::S2CsI() : adc_trig(2) {
  adc_trig.shrink_to_fit();
  Init();
}

void S2CsI::Init() {
  rse::InitCArray(ts, 0);
  rse::InitCArray(eventnr, -999);
  csihit.clear();
  for (auto &r : adc_trig) {
    r.clear();
  }
  for(int i=0; i<2; i++){
    all_ok[i]=kFALSE;
    for(int j=0; j<4; j++){
       polarity[i][j]=-9999;
       n_slave[i][j]=0;
    }
  }
  rse::InitCArray(all_ok, false);
}

oid S2CsI::DrawHist(int ihit){
  static std::shared_ptr<TCanvas> c;
  c.reset();
  c = std::make_shared<TCanvas>("c_csi", "", 700, 500);
  static std::shared_ptr<TGraph> graph;
  graph.reset();
  graph = std::make_shared<TGraph>();
  for(int i = 0; i < (int)csihit[ihit].adc.size(); ++i){
    graph->SetPoint(i, i, csihit[ihit].adc[i]);
  }
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.5);
  graph->Draw("AP");
  graph->SetTitle(Form("pc%d sfp%d febex%d ch%d", csihit[ihit].i_pc, csihit[ihit].i_sfp, csihit[ihit].i_febex, csihit[ihit].i_ch));
}

void S2CsI::DrawADC(int ipc){
    static std::shared_ptr<TCanvas> c;
  c.reset();
  c = std::make_shared<TCanvas>("c_adc_trig", "", 700, 500);
  static std::shared_ptr<TGraph> graph;
  graph.reset();
  graph = std::make_shared<TGraph>();
  for(int i = 0; i < (int)adc_trig[ipc].size(); ++i){
    graph->SetPoint(i, i, adc_trig[ipc][i]);
  }
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.5);
  graph->Draw("AP");
  graph->SetTitle(Form("pc%d adc trig.", ipc));
}

S2CsI::~S2CsI() {
}
