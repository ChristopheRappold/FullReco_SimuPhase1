#ifndef HYPHITRACKV2_h
#define HYPHITRACKV2_h

#include "TObject.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TVector3.h"
#include "TMatrixT.h"
#include "TClass.h"

class THyphiTrack : public TObject{

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

  Int_t Charge;
  Int_t BarId;
  Double32_t dE;
  Double32_t Beta;
  TVector3 RefPoint;
  TVector3 HitTR1;
  TVector3 HitTR2;
  TVector3 HitDC1;
  TVector3 HitDC1p;
  TVector3 NumDC1;
  TVector3 HitDC2;
  TVector3 NumDC2;
  TVector3 HitTOF;
  Double32_t Pval2;
  Int_t TofsBar;

  Double32_t PathLength;
  Double32_t TOF;


  TVector3 MomIni;
  Double32_t RChiIni;
  //  Double32_t PathLengthIni;
  Double32_t BetaIni;
  

  ////////////////////////
  Double32_t State[5];
  Double32_t Cov[5][5];
  ////////////////////////
  TVector3 MomTof;  
  ///////////////////////

//   TMatrixT<double> state;
//   TMatrixT<double> cov;
  
  THyphiTrack();
  THyphiTrack(const THyphiTrack& H);
  ~THyphiTrack();

  virtual void Clear(Option_t *option ="");
  void SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec);
  void SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec,Int_t charge,Int_t barid,Double_t beta,const TVector3& hit_tr1,const TVector3& hit_tr2, const TVector3 hit_dc2, const TVector3 num_dc2,Double_t pval2,Int_t tofsbar );

//   inline void setState(const TMatrixT<double>& aState) {
//     state = aState;
//   }
//   inline void setCov(const TMatrixT<double>& aCov) {
//     cov = aCov;
//   }
  ClassDef(THyphiTrack,2)
    
};

#endif
