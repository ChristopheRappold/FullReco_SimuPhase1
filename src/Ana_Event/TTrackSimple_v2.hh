#ifndef TRACKSIMPLE_h
#define TRACKSIMPLE_h

#include "TObject.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TVector3.h"

class TTrackSimple : public TObject{

public :
  
  TString type;
  Int_t MC_status;
  Double32_t Chi2;
  Double32_t Chi2_X;
  Double32_t Chi2_Y;
  Double32_t Mass;
  Int_t pdgcode;
  TLorentzVector MomMass;
  TVector3 Mom;


  TTrackSimple();
  TTrackSimple(const TTrackSimple& H);
  ~TTrackSimple();

  virtual void Clear(Option_t *option ="");
  void SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec);
  ClassDef(TTrackSimple,4)
    
};

#endif
