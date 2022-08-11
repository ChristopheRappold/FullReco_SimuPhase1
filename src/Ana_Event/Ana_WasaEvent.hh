#ifndef Ana_WasaEvent_h
#define Ana_WasaEvent_h

#include "TVector3.h"
#include "TClonesArray.h"
#include <string>
#include <vector>
#include "TObject.h"
#include "THypernucleus.hh"
#include "THyphiHitDet.hh"
#include "THyphiTrack_v4.hh"
#include "TTrackCand.hh"
#include "TDataHit.hh"


class THypernucleus;
class THyphiHitDet;
class THyphiTrack;
class TDataHit;

class Ana_WasaEvent : public TObject{
public :

  Ana_WasaEvent();
  ~Ana_WasaEvent();
  
  void Clear(Option_t *option ="");
  int Setup();
  void Reset();

  Int_t trigger;

  ULong64_t timestamp;

  TString lmd_filename;
  TString lmd_open_time;

  Double32_t Field;
  Double32_t ReducFactor;

  Int_t NUft;
  Int_t NMft;
  Int_t NDft;

  Int_t NFrstpc;

  Int_t NMdc;
  Int_t NCsi;

  Int_t NT0_counter;
  Int_t NPsfe;
  Int_t NPsbe;
  Int_t NPsce;

  Int_t NSci;
  Int_t NMwdc;

  Int_t NWfd;


  TClonesArray* UFT; //->
  TClonesArray* MFT; //->
  TClonesArray* DFT; //->

  TClonesArray* FRSTPC; //->

  TClonesArray* MDC; //->
  TClonesArray* CSI; //->

  TClonesArray* T0_COUNTER; //->
  TClonesArray* PSFE; //->
  TClonesArray* PSBE; //->
  TClonesArray* PSCE; //->

  TClonesArray* SCI; //->
  TClonesArray* MWDC; //->

  TClonesArray* WFD; //->

  Int_t NtrackCand;
  TClonesArray* TrackCand; //->  

  Int_t Ntracks;
  TClonesArray* fTrack; //->

  Int_t Nhyp;
  TClonesArray* fHyp; //->

  ClassDef(Ana_WasaEvent,1)
};

#endif
