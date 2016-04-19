#ifndef THYPHIHITDET_H
#define THYPHIHITDET_H

#include "TObject.h"
#include "TString.h"
#include "TClass.h"
#include "Riostream.h"
#include "TVector3.h"
#include "TVector2.h"

class THyphiHitDet : public TObject
{
public :

  TString name; //!

  TVector3 HitPos;
  std::vector<int> Id;
  std::vector<double> Qdc;
  std::vector<int> TimeCh;
  double TimeNs;
  THyphiHitDet();
  // THyphiHitDet(TString name,double x,double y,double z,double t,double E,int pid);
  THyphiHitDet(THyphiHitDet& M);
  THyphiHitDet(const THyphiHitDet& M); 
  THyphiHitDet& operator=(const THyphiHitDet&);
  friend TBuffer &operator<<(TBuffer &b, const THyphiHitDet*);

  //virtual THyphiHitDet* CloneD(const char*);
  virtual  void Clear(Option_t* ="");
  ~THyphiHitDet();

 
  ClassDef(THyphiHitDet,1)

};

#endif
