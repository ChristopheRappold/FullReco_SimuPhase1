#ifndef TRACKSIMPLE_h
#define TRACKSIMPLE_h

#include "TObject.h"
#include "TLorentzVector.h"

//using namespace std;

class TTrackSimple : public TObject{

public :
  
  Double32_t Chi2;
  Double32_t Mass;
  TLorentzVector MomMass;
  
  TTrackSimple();
  TTrackSimple(const TTrackSimple& H);
  ~TTrackSimple();

  virtual void Clear(Option_t *option ="");

  ClassDef(TTrackSimple,2)
    
};

#endif
