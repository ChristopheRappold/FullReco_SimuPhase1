#ifndef THYPHIDATAHIT_H
#define THYPHIDATAHIT_H

#include "TObject.h"
#include "TString.h"

#include "Riostream.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class TDataHit : public TObject
{
public :

  TString name; //!

  Int_t LayerID;
  Int_t HitID;
  Int_t TrackID;

  TVector3 Hit;

  Int_t Charge;
  Int_t Pdg;
  Double_t Brho;
  Double_t MagnetInteraction;

  Double_t Time;
  Double_t Energy;
  Double_t TrackLength;

  TDataHit();

  TDataHit(const TDataHit& M); 
  TDataHit& operator=(const TDataHit&);

  virtual  void Clear(Option_t* ="");
  ~TDataHit();

 
  ClassDef(TDataHit,3)

};

#endif
