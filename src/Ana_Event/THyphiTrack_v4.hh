#ifndef HYPHITRACKV3_h
#define HYPHITRACKV3_h

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
  Double_t Chi2;
  Float_t Chi2_X;
  Float_t Chi2_Y;
  Double_t Mass;
  Int_t pdgcode;
  TLorentzVector MomMass;
  TVector3 Mom;

  Int_t Charge;
  Int_t BarId;
  Float_t dE;
  Float_t Beta;
  TVector3 RefPoint;
  // TVector3 HitTR1;
  // TVector3 HitTR2;
  // TVector3 HitDC1;
  // TVector3 HitDC1p;
  // TVector3 NumDC1;
  // TVector3 HitDC2;
  // TVector3 NumDC2;
  // TVector3 HitTOF;
  Double_t Pval2;
  Int_t TofsBar;

  Double_t PathLength;
  Double_t TOF;


  TVector3 MomIni;
  Float_t RChiIni;
  //  Double32_t PathLengthIni;
  Float_t BetaIni;
  

  ////////////////////////
  Double_t State[6];
  Double_t Cov[6][6];
  ////////////////////////
  // TVector3 MomTof;  
  // ///////////////////////
  // Double32_t Residu_X_TR1;
  // Double32_t Residu_X_TR2;
  // Double32_t Residu_X_DC2;
  // Double32_t Residu_X_TOF;

  // Double32_t Residu_Y_TR1;
  // Double32_t Residu_Y_TR2;
  // Double32_t Residu_Y_DC2;
  // Double32_t Residu_Y_TOF;

//   TMatrixT<double> state;
//   TMatrixT<double> cov;
  
  THyphiTrack();
  THyphiTrack(const THyphiTrack& H);
  ~THyphiTrack();

  virtual void Clear(Option_t *option ="");
  // void SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec);
  // void SetPara(const TString& name,Int_t MC,Double_t chi2,Double_t chi2_x, Double_t chi2_y,Double_t mass, Int_t pdg,const TVector3& vec,Int_t charge,Int_t barid,Double_t beta,const TVector3& hit_tr1,const TVector3& hit_tr2, const TVector3 hit_dc2, const TVector3 num_dc2,Double_t pval2,Int_t tofsbar );

//   inline void setState(const TMatrixT<double>& aState) {
//     state = aState;
//   }
//   inline void setCov(const TMatrixT<double>& aCov) {
//     cov = aCov;
//   }
  ClassDef(THyphiTrack,5)
    
};

#endif
