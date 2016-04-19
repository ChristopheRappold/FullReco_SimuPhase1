#ifndef MCAnaEventG4Sol_h
#define MCAnaEventG4Sol_h

#include "TClonesArray.h"
#include <string>

#include "TObject.h"

#include "TMcHit.hh"
#include "TMcParticle.hh"
#include "THyphiTrack_v4.hh"

// class TMcHit;
// class TMcParticle;

class MCAnaEventG4Sol : public TObject{
public :

  MCAnaEventG4Sol();
  ~MCAnaEventG4Sol();
  
  void Clear(Option_t *option ="");
  int Setup();
  static void   Reset();
  // int Add_MC(const TMcParticle& M);
  // int Add_Hit(const TMcHit& H,TString detector);
  // int Add_Track(const THyphiTrack& T);

  Int_t trigger;

  Double32_t Field;
  Double32_t ReducFactor;

  Int_t Nmc;
  TClonesArray* fMC_Particle; //->

  Int_t NInSi;
  Int_t NCdc;
  Int_t NCdh;
  Int_t NFwdtracker;
  Int_t NRpc;
  Int_t NFmf2;
  TClonesArray* InSi; //->
  TClonesArray* CDC; //->
  TClonesArray* CDH; //->
  TClonesArray* FwdTracker; //->
  TClonesArray* RPC; //->
  TClonesArray* FMF2; //->

  Int_t Ntrack;
  TClonesArray* fTrack; //->

  static TClonesArray* gMC_Particle; //!
  
  static TClonesArray* gInSi; //!
  static TClonesArray* gCDC; //!
  static TClonesArray* gCDH; //!
  static TClonesArray* gFwdTracker; //!
  static TClonesArray* gRPC; //!
  static TClonesArray* gFMF2; //!

  static TClonesArray* gfTrack; //! 
  
  ClassDef(MCAnaEventG4Sol,1)
};

#endif
