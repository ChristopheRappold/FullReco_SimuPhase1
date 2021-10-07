#ifndef MCAnaEventG4Sol_h
#define MCAnaEventG4Sol_h

#include "TClonesArray.h"
#include <string>

#include "TObject.h"

#include "TMcHit.hh"
#include "TMcParticle.hh"
#include "TTrackCand.hh"
#include "THyphiTrack_v4.hh"
#include "THypernucleus.hh"

// class TMcHit;
// class TMcParticle;

class MCAnaEventG4Sol : public TObject{
public :

  MCAnaEventG4Sol();
  ~MCAnaEventG4Sol();
  
  void Clear(Option_t *option ="");
  int Setup();
  void   Reset();
  // int Add_MC(const TMcParticle& M);
  // int Add_Hit(const TMcHit& H,TString detector);
  // int Add_Track(const THyphiTrack& T);

  Int_t trigger;

  Double32_t Field;
  Double32_t ReducFactor;

  Int_t Nmc;
  TClonesArray* fMC_Particle; //->

  Int_t NInSi;
  Int_t NTr;
  Int_t NFiber;
  Int_t NCdc;
  Int_t NCdh;
  Int_t NFwdtracker;
  Int_t NRpc;
  Int_t NFmf2;
  Int_t NPsfe;
  Int_t NPsbe;
  Int_t NPsce;
  
  TClonesArray* InSi; //->
  TClonesArray* TR; //->
  TClonesArray* Fiber; //->
  TClonesArray* CDC; //->
  TClonesArray* CDH; //->
  TClonesArray* FwdTracker; //->
  TClonesArray* RPC; //->
  TClonesArray* FMF2; //->

  TClonesArray* PSFE; //->
  TClonesArray* PSBE; //->
  TClonesArray* PSCE; //->
  
  Int_t NtrackCand;
  TClonesArray* TrackCand; //->

  Int_t Ntrack;
  TClonesArray* fTrack; //->

  Int_t Nhyp;
  TClonesArray* fHyp; //->

  // static TClonesArray* gMC_Particle; //!

  // static TClonesArray* gInSi; //!
  // static TClonesArray* gTR; //!
  // static TClonesArray* gFiber; //!
  // static TClonesArray* gCDC; //!
  // static TClonesArray* gCDH; //!
  // static TClonesArray* gFwdTracker; //!
  // static TClonesArray* gRPC; //!
  // static TClonesArray* gFMF2; //!

  // static TClonesArray* gPSFE; //!
  // static TClonesArray* gPSBE; //!
  // static TClonesArray* gPSCE; //!

  // static TClonesArray* gTrackCand; //!

  // static TClonesArray* gfTrack; //!
  // static TClonesArray* gfHyp; //!
  
  ClassDef(MCAnaEventG4Sol,6)
};

#endif
