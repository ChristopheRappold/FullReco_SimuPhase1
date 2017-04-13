#ifndef Ana_Event_v16_h
#define Ana_Event_v16_h

#include "TVector3.h"
#include "TClonesArray.h"
#include <string>
#include <vector>
#include "TObject.h"
#include "THypernucleus.hh"
//#include "TTrackSimple_v2.hh"
//#include "THyphiSimpleHit.hh"
#include "THyphiHitDet.hh"
#include "TMcParticle.hh"
#include "THyphiTrack_v4.hh"

class THypernucleus;
//class THyphiSimpleHit;
class THyphiHitDet;
//class TTrackSimple;
class TMcParticle;
class THyphiTrack;

class Ana_Event : public TObject{
public :

  Ana_Event();
  ~Ana_Event();
  
  void Clear(Option_t *option ="");
  int Setup();
  static void   Reset();
  int Add_MC(const TMcParticle& M);
  //  int Add_HitTr0(const THyphiSimpleHit& H);
  //  int Add_HitTr0(const THyphiHitDet& H);
  int Add_Hit(const THyphiHitDet& H,TString detector);
  //  int Add_Track(const TTrackSimple& T);
  int Add_Track(const THyphiTrack& T);
  int Add_Hyp(const THypernucleus& H);

  Int_t trigger;

  Int_t Nmc;
  TClonesArray* fMC_Particle; //->

  Int_t Nhit_tr0x;
  Int_t Nhit_tr0y;
  Int_t Nhit_tr1x;
  Int_t Nhit_tr1y;
  Int_t Nhit_tr2x;
  Int_t Nhit_tr2y;
  Int_t Nhit_dc1x;
  Int_t Nhit_dc1xp;
  Int_t Nhit_dc1u;
  Int_t Nhit_dc1up;
  Int_t Nhit_dc1v;
  Int_t Nhit_dc1vp;
  Int_t Nhit_dc2x;
  Int_t Nhit_dc2xp;
  Int_t Nhit_dc2y;
  Int_t Nhit_dc2yp;
  Int_t Nhit_dc2u;
  Int_t Nhit_tofs;
  Int_t Nhit_tofp;
  Int_t Nhit_bg;



  TClonesArray* TR0x; //->
  TClonesArray* TR0y; //->
  TClonesArray* TR1x; //->
  TClonesArray* TR1y; //->
  TClonesArray* TR2x; //->
  TClonesArray* TR2y; //->

  TClonesArray* DC1x; //->
  TClonesArray* DC1xp; //->
  TClonesArray* DC1u; //->
  TClonesArray* DC1up; //->
  TClonesArray* DC1v; //->
  TClonesArray* DC1vp; //->
  TClonesArray* DC2x; //->
  TClonesArray* DC2xp; //->
  TClonesArray* DC2y; //->
  TClonesArray* DC2yp; //->
  TClonesArray* DC2u; //->
  TClonesArray* TOFs; //->
  TClonesArray* TOFp; //->
  TClonesArray* Bg; //->

  

  Int_t Ntracks;
  TClonesArray* fTrack; //->
  Int_t Nhyp;
  TClonesArray* fHyp; //->

  static TClonesArray* gMC_Particle; //!
  //  static TClonesArray* gTR0; //!

  static TClonesArray* gTR0x; //!
  static TClonesArray* gTR0y; //!
  static TClonesArray* gTR1x; //!
  static TClonesArray* gTR1y; //!
  static TClonesArray* gTR2x; //!
  static TClonesArray* gTR2y; //!

  static TClonesArray* gDC1x; //!
  static TClonesArray* gDC1xp; //!
  static TClonesArray* gDC1u; //!
  static TClonesArray* gDC1up; //!
  static TClonesArray* gDC1v; //!
  static TClonesArray* gDC1vp; //!
  static TClonesArray* gDC2x; //!
  static TClonesArray* gDC2xp; //!
  static TClonesArray* gDC2y; //!
  static TClonesArray* gDC2yp; //!
  static TClonesArray* gDC2u; //!
  static TClonesArray* gTOFs; //!
  static TClonesArray* gTOFp; //!
  static TClonesArray* gBg; //!

  static TClonesArray* gfTrack; //! 
  static TClonesArray* gfHyp; //!

  ClassDef(Ana_Event,16)
};

#endif
