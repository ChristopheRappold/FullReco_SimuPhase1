#ifndef HYP_h
#define HYP_h

#include "TVector3.h"
#include <string>
#include <vector>
#include "TObject.h"
#include "TLorentzVector.h"

//using namespace std;

class THypernucleus : public TObject{

public :
  
  TString type;
  Int_t pattern; ///  1=no cut / 2=mass_cut / 4 = best_mass_cut / 10 = vtx cut / 100 = best vtx cut
  Int_t Ndecay;
  //Double32_t Chi2;
  Double32_t Pvalue;
  Double32_t InvMass;
  Double32_t Dist;
  TLorentzVector MomMass;
  TLorentzVector Vtx;
  TLorentzVector MomMassD1;
  TLorentzVector MomMassD2;
  TLorentzVector MomMassD3;

  THypernucleus();
  THypernucleus(const THypernucleus& H);
  ~THypernucleus();

  virtual void Clear(Option_t *option ="");

  ClassDef(THypernucleus,2)
    
};

#endif
