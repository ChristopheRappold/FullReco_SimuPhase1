/** HypConstField.h
 ** @author M.Al-Turany <m.al/turany@gsi.de>
 ** @since 30.01.2007
 ** @version1.0
 ** A constant  magnetic field
 **/
/** modified for HypHI **/


#ifndef HYPSOLENOIDFIELD_H
#define HYPSOLENOIDFIELD_H 1


#include "FairField.h"
#include "TString.h"

class FrsSolenoidHypField : public FairField
{

 public:    

  /** Default constructor **/
  FrsSolenoidHypField();


  /** Standard constructor 
   ** @param name   Object name
   ** @param xMin,xMax   x region of field (global coordinates)
   ** @param yMin,yMax   y region of field (global coordinates)
   ** @param zMin,zMax   z region of field (global coordinates)
   ** @param bX,bY,bZ    Field values [kG]
   **/
  FrsSolenoidHypField(const char* name, 
		      Double_t rMin, Double_t rMax,
		      Double_t zMin, Double_t zMax, 
		      Double_t pos_x,Double_t pos_y,Double_t pos_z,
		      Double_t bX, Double_t bY, Double_t bZ);

  /** Destructor **/
  virtual ~FrsSolenoidHypField();
	
	
  /** Set the field region
   ** @param xMin,xMax   x region of field (global coordinates)
   ** @param yMin,yMax   y region of field (global coordinates)
   ** @param zMin,zMax   z region of field (global coordinates)
   **/
  void SetFieldRegion(Double_t rMin, Double_t rMax, Double_t zMin, Double_t zMax,
                      Double_t pos_x,Double_t pos_y,Double_t pos_z);


  void SetPositionFromGeoManager(const TString& name_geo_mat, const TString& name_geo_field);
  void SetPositionFromGeoManager(const TString& name_node);
  /** Set the field values
   ** @param bX,bY,bZ    Field values [kG]
   **/
  void SetField(Double_t bX, Double_t bY, Double_t bZ);
  

  Bool_t IsInside(Double_t x, Double_t y, Double_t z);
  /** Get components of field at a given point 
   ** @param x,y,z   Point coordinates [cm]
   **/
  virtual Double_t GetBx(Double_t x, Double_t y, Double_t z);
  virtual Double_t GetBy(Double_t x, Double_t y, Double_t z);
  virtual Double_t GetBz(Double_t x, Double_t y, Double_t z);

  void GetBxyz(const Double_t point[3], Double_t* bField);


  /** Accessors to field region **/
  Double_t GetRmin() const { return fRmin; }
  Double_t GetRmax() const { return fRmax; }
  Double_t GetZmin() const { return fZmin; }
  Double_t GetZmax() const { return fZmax; }
  Double_t GetPositionX() const { return fPosX; }
  Double_t GetPositionY() const { return fPosY; }
  Double_t GetPositionZ() const { return fPosZ; }
  
  /** Accessors to field values **/
  Double_t GetBx() const { return fBx; }
  Double_t GetBy() const { return fBy; }
  Double_t GetBz() const { return fBz; }


  /** Screen output **/
  virtual void Print(Option_t *option="") const;
  
  ClassDef(FrsSolenoidHypField, 1);

 private:

  /** Limits of the field region **/
  Double_t fRmin;   
  Double_t fRmax;
  Double_t fZmin;
  Double_t fZmax;
  
  /** Field centre position in global coordinates  **/
  Double_t fPosX, fPosY, fPosZ; 

  
  /** Field components inside the field region **/
  Double_t fBx;
  Double_t fBy;
  Double_t fBz;

  

};


#endif
