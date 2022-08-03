#ifndef Ana_Event_vData_h
#define Ana_Event_vData_h

#include "TVector3.h"
#include "TClonesArray.h"
#include <string>
#include <vector>
#include "TObject.h"
#include "THypernucleus.hh"
#include "THyphiHitDet.hh"
#include "TMcParticle.hh"
#include "THyphiTrack_v4.hh"
#include "TTrackCand.hh"
#include "TDataHit.hh"


class THypernucleus;
class THyphiHitDet;
class TMcParticle;
class THyphiTrack;
class TDataHit;

class Ana_Event : public TObject{
public :

  Ana_Event();
  ~Ana_Event();
  
  void Clear(Option_t *option ="");
  int Setup();
  static void Reset();
  int Add_MC(const TMcParticle& M);
  int Add_Hit(const THyphiHitDet& H,TString detector);
  int Add_Track(const THyphiTrack& T);
  int Add_Hyp(const THypernucleus& H);

  Int_t trigger;

  Int_t Nmc;
  TClonesArray* fMC_Particle; //->

  Int_t NT0_counter;
  Int_t NUft;
  Int_t NMft;
  Int_t NDft;

  Int_t NMdc;
  Int_t NCsi;
  Int_t NPsfe;
  Int_t NPsbe;
  Int_t NPsce;

  Int_t NSci;
  Int_t NMwdc;

  TClonesArray* T0_COUNTER; //->
  TClonesArray* UFT; //->
  TClonesArray* MFT; //->
  TClonesArray* DFT; //->

  TClonesArray* MDC; //->
  TClonesArray* CSI; //->
  TClonesArray* PSFE; //->
  TClonesArray* PSBE; //->
  TClonesArray* PSCE; //->

  TClonesArray* SCI; //->
  TClonesArray* MWDC; //->

  Int_t NtrackCand;
  TClonesArray* TrackCand; //->  

  Int_t Ntracks;
  TClonesArray* fTrack; //->

  Int_t Nhyp;
  TClonesArray* fHyp; //->


  static TClonesArray* gMC_Particle; //!
  static TClonesArray* gfTrack; //! 
  static TClonesArray* gfHyp; //!

  ClassDef(Ana_Event,17)
};

#endif
