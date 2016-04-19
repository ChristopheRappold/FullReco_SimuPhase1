#ifndef MCPARTICLE_h
#define MCPARTICLE_h

#include "TVector3.h"
#include <string>
#include <vector>
#include "TObject.h"
#include "TLorentzVector.h"

//using namespace std;

class TMcParticle : public TObject{

public :
  
  TString type;
  Int_t Mc_id;
  Int_t Mother_id;
  Int_t Pdg;
  Int_t Charge;
  TLorentzVector MomMass;
  TLorentzVector Vtx;
  Double_t Weigth;
  Bool_t GeoAcc;

  TMcParticle();
  TMcParticle(const TMcParticle& H);
  ~TMcParticle();

  virtual void Clear(Option_t* option ="");
  virtual void Print(Option_t* option = "") const;

  ClassDef(TMcParticle,3)
    
};

#endif
