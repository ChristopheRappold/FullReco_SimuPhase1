#ifndef THYPHISIMPLEHIT_H
#define THYPHISIMPLEHIT_H

#include "TObject.h"
#include "TString.h"
#include "TClass.h"
#include "Riostream.h"

class THyphiSimpleHit : public TObject
{
public :

  TString name; //!
  Double32_t x;
  Double32_t y;
  Double32_t z;
  Double32_t t;
  Double32_t E;
  Int_t PID;


  THyphiSimpleHit();
  THyphiSimpleHit(TString name,double x,double y,double z,double t,double E,int pid);
  THyphiSimpleHit(THyphiSimpleHit& M);
  THyphiSimpleHit(const THyphiSimpleHit& M); 
  THyphiSimpleHit& operator=(const THyphiSimpleHit&);
  friend TBuffer &operator<<(TBuffer &b, const THyphiSimpleHit*);

  //virtual THyphiSimpleHit* CloneD(const char*);
  virtual  void Clear(Option_t* ="");
  ~THyphiSimpleHit();

 
  ClassDef(THyphiSimpleHit,1)

};

#endif
