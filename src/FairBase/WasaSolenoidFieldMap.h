/** HypConstField.h
 ** @author M.Al-Turany <m.al/turany@gsi.de>
 ** @since 30.01.2007
 ** @version1.0
 ** A constant  magnetic field
 **/
/** modified for HypHI **/


#ifndef WASASOLENOIDFIELDMAP_H
#define WASASOLENOIDFIELDMAP_H 1


#include "FairField.h"
#include "TString.h"
#include "MField.h"

class WasaSolenoidFieldMap : public FairField
{

 public:    

  /** Default constructor **/
  WasaSolenoidFieldMap();

  WasaSolenoidFieldMap(const char* name, const char* title, const TString& namefile_field, double maxF, double signD);

  /** Destructor **/
  virtual ~WasaSolenoidFieldMap();

  void SetMaxField(double maxField);
  void Init();

  void SetPositionFromGeoManager(const TString& name_node);

  Bool_t IsInside(Double_t x, Double_t y, Double_t z);
  /** Get components of field at a given point 
   ** @param x,y,z   Point coordinates [cm]
   **/
  virtual Double_t GetBx(Double_t x, Double_t y, Double_t z);
  virtual Double_t GetBy(Double_t x, Double_t y, Double_t z);
  virtual Double_t GetBz(Double_t x, Double_t y, Double_t z);

  void GetBxyz(const Double_t point[3], Double_t* bField);


  /** Accessors to field region **/
  Double_t GetPositionX() const { return fPosX; }
  Double_t GetPositionY() const { return fPosY; }
  Double_t GetPositionZ() const { return fPosZ; }

  /** Screen output **/
  virtual void Print(Option_t *option="") const;
  
  ClassDef(WasaSolenoidFieldMap, 1);

 private:

  /** Limits of the field region **/
  Double_t fRmin;   
  Double_t fRmax;
  Double_t fZmin;
  Double_t fZmax;

  TString nameField;
  Double_t maxField;
  Double_t signDir;

  /** Field centre position in global coordinates  **/
  Double_t fPosX, fPosY, fPosZ; 

  std::unique_ptr<MField> FieldMap;

  Double_t fieldAtOrig[3];

};


#endif
