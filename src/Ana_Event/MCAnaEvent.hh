#ifndef MCAnaEvent_h
#define MCAnaEvent_h

#include "TClonesArray.h"
#include <string>

#include "TObject.h"
//#include "THypernucleus.hh"
//#include "TTrackSimple_v2.hh"
//#include "THyphiSimpleHit.hh"
#include "TMcHit.hh"
#include "TMcParticle.hh"
#include "THyphiTrack_v3.hh"

//class THypernucleus;
//class THyphiSimpleHit;
class TMcHit;
//class TTrackSimple;
class TMcParticle;
//class THyphiTrack;

class MCAnaEvent : public TObject{
public :

  MCAnaEvent();
  ~MCAnaEvent();
  
  void Clear(Option_t *option ="");
  int Setup();
  static void   Reset();
  int Add_MC(const TMcParticle& M);
  //  int Add_HitTr0(const THyphiSimpleHit& H);
  //  int Add_HitTr0(const TMcHit& H);
  int Add_Hit(const TMcHit& H,TString detector);
  //int Add_Track(const TTrackSimple& T);
  int Add_Track(const THyphiTrack& T);
  //int Add_Hyp(const THypernucleus& H);

  Int_t trigger;

  Int_t Nsystematic;
  Double32_t SysX_shift;
  Double32_t SysZ_shift;
  Double32_t SysY_angle;
  Double32_t Field;

  Double32_t SysX_shift2; // cm
  Double32_t SysY_shift2; // cm
  Double32_t SysZ_angle2; // degree
  Double32_t Field2; // Telsa
  
  Int_t Nmc;
  TClonesArray* fMC_Particle; //->

  Int_t Nhit_tr0;
  Int_t Nhit_tr1;
  Int_t Nhit_tr2;
  Int_t Nhit_dc1;
  Int_t Nhit_dc2;
  Int_t Nhit_dc3;
  Int_t Nhit_dc2stop;
  Int_t Nhit_tofp;
  Int_t Nhit_stop;
  Int_t Nhit_stop2;



  TClonesArray* TR0; //->
  TClonesArray* TR1; //->
  TClonesArray* TR2; //->

  TClonesArray* DC1; //->
  TClonesArray* DC2; //->
  TClonesArray* DC3; //->
  TClonesArray* DC2stop; //->
  TClonesArray* TOFp; //->
  TClonesArray* STOP; //->
  TClonesArray* STOP2; //->

  

  Int_t Ntracks;
  TClonesArray* fTrack; //->
  // Int_t Nhyp;
  // TClonesArray* fHyp; //->

  static TClonesArray* gMC_Particle; //!
  //  static TClonesArray* gTR0; //!

  static TClonesArray* gTR0; //!
  static TClonesArray* gTR1; //!
  static TClonesArray* gTR2; //!

  static TClonesArray* gDC1; //!
  static TClonesArray* gDC2; //!
  static TClonesArray* gDC3; //!
  static TClonesArray* gDC2stop; //!
  static TClonesArray* gTOFp; //!
  static TClonesArray* gSTOP; //!
  static TClonesArray* gSTOP2; //!

   static TClonesArray* gfTrack; //! 
  // static TClonesArray* gfHyp; //!

  ClassDef(MCAnaEvent,3)
};

#endif
