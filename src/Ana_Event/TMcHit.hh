#ifndef THYPHIMCHIT_H
#define THYPHIMCHIT_H

#include "TObject.h"
#include "TString.h"

#include "Riostream.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class TMcHit : public TObject
{
public :

  TString name; //!

  Int_t LayerID;
  Int_t HitID;
  Int_t TrackID;

  TVector3 MCHit;
  TVector3 Hit;

  Int_t MC_id;
  Int_t Charge;
  Int_t Pdg;
  Double_t Brho;
  Double_t MagnetInteraction;
  TLorentzVector MCparticle;

  Double_t Time;
  Double_t Energy;
  Double_t TrackLength;
  Double_t HitLength;

  TMcHit();

  TMcHit(const TMcHit& M); 
  TMcHit& operator=(const TMcHit&);
  //friend TBuffer &operator<<(TBuffer &b, const TMcHit*);

  //virtual TMcHit* CloneD(const char*);
  virtual  void Clear(Option_t* ="");
  ~TMcHit();

 
  ClassDef(TMcHit,3)

};

#endif
