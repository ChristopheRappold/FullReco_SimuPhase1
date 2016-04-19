#ifndef Ana_Event_h
#define Ana_Event_h

#include "TVector3.h"
#include "TClonesArray.h"
#include <string>
#include <vector>
#include "TObject.h"
#include "THypernucleus.hh"
#include "TTrackSimple.hh"
#include "THyphiSimpleHit.hh"

//using namespace std;

// class Hypernucleus : public TObject{

// public :

//   TString type;
//   Int_t pattern; ///  1=no cut / 2=mass_cut / 4 = best_mass_cut / 10 = vtx cut / 100 = best vtx cut
//   Int_t Ndecay;
//   //Double32_t Chi2;
//   Double32_t Pvalue;
//   Double32_t InvMass;
//   Double32_t Dist;
//   TLorentzVector MomMass;
//   TLorentzVector Vtx;
//   TLorentzVector MomMassD1; 
//   TLorentzVector MomMassD2; 
//   TLorentzVector MomMassD3;
   
//   Hypernucleus();
//   Hypernucleus(const Hypernucleus& H);
//   ~Hypernucleus();
  
//   virtual void Clear(Option_t *option ="");
  
//   ClassDef(Hypernucleus,2)
   
// };

class THypernucleus;
class THyphiSimpleHit;
class TTrackSimple;

class Ana_Event : public TObject{
public :

  Ana_Event();
  ~Ana_Event();
  
  void Clear(Option_t *option ="");
  int Setup();
  static void   Reset();
  int Add_MC(TVector3 mom,TVector3 pos);
  int Add_HitTr0(const THyphiSimpleHit& H);
  int Add_Track(const TTrackSimple& T);
  int Add_Hyp(const THypernucleus& H);

  std::vector<std::string> fMC_name;
  Int_t Nmc;
  TClonesArray* fMC_Mom; //->
  TClonesArray* fMC_Vtx; //->


  Int_t Nhit_tr0;
  TClonesArray* TR0;


  Int_t Ntracks;
  TClonesArray* fTrack; //->
  Int_t Nhyp;
  TClonesArray* fHyp; //->

  static TClonesArray* gMC_Mom; //!
  static TClonesArray* gMC_Vtx; //!
  static TClonesArray* gTR0; //!
  static TClonesArray* gfTrack; //! 
  static TClonesArray* gfHyp; //!


  ClassDef(Ana_Event,4)
};

#endif
