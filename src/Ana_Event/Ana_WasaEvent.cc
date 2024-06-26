#include "Ana_WasaEvent.hh"

using namespace std;

ClassImp(Ana_WasaEvent)


Ana_WasaEvent::Ana_WasaEvent()
{
  UFT = new TClonesArray("TFiberHit",20);
  MFT = new TClonesArray("TFiberHit",20);
  DFT = new TClonesArray("TFiberHit",20);

  FRSTPC = new TClonesArray("TFrsTpcHit",20);

  MDC = new TClonesArray("TMDCHit",20);
  CSI = new TClonesArray("TCsIHit",20);

  T0_COUNTER = new TClonesArray("TWasaPSHit",20);
  PSBE = new TClonesArray("TWasaPSHit",20);
  PSFE = new TClonesArray("TWasaPSHit",20);
  PSCE = new TClonesArray("TWasaPSHit",20);

  SCI = new TClonesArray("TSciHit",20);
  MWDC = new TClonesArray("TS4MwdcHit",20);

  WFD = new TClonesArray("TWfdHit",20);

  TrackCand = new TClonesArray("TTrackCand",100);
  fTrack = new TClonesArray("THyphiTrack",20);
  fHyp = new TClonesArray("THypernucleus",20);


  trigger = 0;
  Field = 0.;
  ReducFactor = 1.;

  NtrackCand = 0;
  Ntrack=0;
  Nhyp=0;

  NUft = 0;
  NMft = 0;
  NDft = 0;

  NFrstpc = 0;

  NMdc = 0;
  NCsi = 0;

  NT0_counter = 0;
  NPsfe = 0;
  NPsbe = 0;
  NPsce = 0;

  NSci = 0;
  NMwdc = 0;

  NWfd = 0;
  
  Setup();
}

Ana_WasaEvent::~Ana_WasaEvent()
{
  Clear();
  Reset();
}

void Ana_WasaEvent::Clear(Option_t *option)
{
  UFT->Clear("C");
  MFT->Clear("C");
  DFT->Clear("C");

  FRSTPC->Clear("C");

  MDC->Clear("C");
  CSI->Clear("C");

  T0_COUNTER->Clear("C");
  PSBE->Clear("C");
  PSFE->Clear("C");
  PSCE->Clear("C");

  SCI->Clear("C");
  MWDC->Clear("C");

  WFD->Clear("C");

  TrackCand->Clear("C");
  fTrack->Clear("C");
  fHyp->Clear("C");

  Setup();
}

void Ana_WasaEvent::Reset()
{
  delete UFT; UFT = 0;
  delete MFT; MFT = 0;
  delete DFT; DFT = 0;

  delete FRSTPC; FRSTPC = 0;

  delete MDC; MDC = 0;
  delete CSI; CSI = 0;

  delete T0_COUNTER; T0_COUNTER = 0;
  delete PSFE; PSFE = 0;
  delete PSBE; PSBE = 0;
  delete PSCE; PSCE = 0;

  delete SCI; SCI = 0;
  delete MWDC; MWDC = 0;

  delete WFD; WFD = 0;

  delete TrackCand; TrackCand = 0;
  delete fTrack; fTrack = 0;
  delete fHyp; fHyp = 0;
}


int Ana_WasaEvent::Setup()
{
  trigger = 0;
  Field = 0;
  ReducFactor =1.;

  NtrackCand =0;
  Ntrack=0;
  Nhyp=0;

  NUft = 0;
  NMft = 0;
  NDft = 0;

  NFrstpc = 0;

  NMdc = 0;
  NCsi = 0;

  NT0_counter = 0;
  NPsfe = 0;
  NPsbe = 0;
  NPsce = 0;

  NSci = 0;
  NMwdc = 0;

  NWfd = 0;

  return 0;
}
