#ifndef Ana_Event_h
#define Ana_Event_h

#include "TVector3.h"
#include "TClonesArray.h"
#include <string>
#include <vector>
#include "TObject.h"
#include "THypernucleus.hh"

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

class Ana_Event : public TObject{
public :

  Ana_Event();
  ~Ana_Event();
  
  void Clear(Option_t *option ="");
  int Setup();
  static void   Reset();
  int Add_MC(TVector3 mom,TVector3 pos);
  //int Add_Seed(TVector3 mom);
  int Add_Track(TVector3 mom,TVector3 pos);
  //int Add_HypAndD(TVector3 mom,TVector3 vtx,TVector3 D1mom, TVector3 D2mom);
  //int Add_HypAndD(TVector3 mom,TVector3 vtx,TVector3 D1mom, TVector3 D2mom, TVector3 D3mom);
  int Add_Hyp(const THypernucleus& H);

  std::vector<std::string> fMC_name;
  //vector<TVector3> MC_Mom;
  //vector<TVector3> MC_Vtx;
  Int_t Nmc;
  TClonesArray* fMC_Mom; //->
  TClonesArray* fMC_Vtx; //->


  Int_t NId;
  std::vector<Int_t> fId;
  std::vector<std::string> fIdName;


  std::vector<Int_t> fTR0id;
  //std::vector<Double_t> fTR0_X;
  //std::vector<Double_t> fTR0_Y;
  //std::vector<Double_t> fTR0_Z;
  //std::vector<Double_t> fTR0_t;
  std::vector<Double_t> fTR0_E;

  //std::vector<Int_t> fTR1id;
  //std::vector<Double_t> fTR1_X;
  //std::vector<Double_t> fTR1_Y;
  //std::vector<Double_t> fTR1_Z;
  //vector<Double_t> fTR1_t;
  //vector<Double_t> fTR1_E;

  //std::vector<Int_t> fTR2id;
  //std::vector<Double_t> fTR2_X;
  //std::vector<Double_t> fTR2_Y;
  //std::vector<Double_t> fTR2_Z;
  //vector<Double_t> fTR2_t;
  //vector<Double_t> fTR2_E;

  //std::vector<Int_t> fTOFid;
  //std::vector<Double_t> fTOF_X;
  //std::vector<Double_t> fTOF_Y;
  //std::vector<Double_t> fTOF_Z;
  //std::vector<Double_t> fTOF_t;
  //std::vector<Double_t> fTOF_E;

  //std::vector<Int_t> fDC1id;
  //std::vector<Double_t> fDC1_X;
  //std::vector<Double_t> fDC1_Y;
  //std::vector<Double_t> fDC1_Z;
  //vector<Double_t> fDC1_t;
  //vector<Double_t> fDC1_E;

  //std::vector<Int_t> fDC2id;
  //std::vector<Double_t> fDC2_X;
  //std::vector<Double_t> fDC2_Y;
  //std::vector<Double_t> fDC2_Z;
  //vector<Double_t> fDC2_t;
  //vector<Double_t> fDC2_E;

  //Int_t Nseed;
  //std::vector<Int_t> fId_bfm;
  //std::vector<Int_t> fId_dc;
  //std::vector<Int_t> fId_tof;
  //TClonesArray* fSeedMom; //

  Int_t Ntracks;
  TClonesArray* fTrackMom; //->
  //TClonesArray* fTrackPos; //->
  //std::vector<Int_t> fTrackCharge;
  std::vector<Double_t> fTrackMass;
  std::vector<Double_t> fTrackChi2;
  //std::vector<Double_t> fTrackPV;

  Int_t Nhyp;
  TClonesArray* fHyp; //->
  //std::vector<std::string> fHypName;
  //std::vector<std::string> fHypD1Name;
  //std::vector<std::string> fHypD2Name;
  //std::vector<std::string> fHypD3Name;

  //std::vector<Int_t> fHypPV;
  //std::vector<Double_t> fHypDist;
  //std::vector<Double_t> fHypInvMass;
  //TClonesArray* fHypMom; //->
  //TClonesArray* fHypVtx; //->
  //TClonesArray* fHypD1Mom; //->
  //TClonesArray* fHypD2Mom; //->
  //TClonesArray* fHypD3Mom; //->


  static TClonesArray* gMC_Mom; //!
  static TClonesArray* gMC_Vtx; //!
  //static TClonesArray* gfSeedMom; //!
  static TClonesArray* gfTrackMom; //! 
  //static TClonesArray* gfTrackPos; //!
  static TClonesArray* gfHyp; //!
  //static TClonesArray* gfHypMom; //!
  //static TClonesArray* gfHypVtx; //!
 // static TClonesArray* gfHypD1Mom; //!
 // static TClonesArray* gfHypD2Mom; //!
 // static TClonesArray* gfHypD3Mom; //!


  ClassDef(Ana_Event,2)
};

#endif
