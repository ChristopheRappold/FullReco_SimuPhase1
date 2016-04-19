// -------------------------------------------------------------------------
// -----                    HypConstField source file                  -----
// -----                Created 30/01/07  by M. Al/Turany              -----
// -------------------------------------------------------------------------

#include "HypConstField.h"
#include "TMath.h"

#include <iomanip>
#include <iostream>


using namespace std;

// -----   Default constructor   -------------------------------------------
HypConstField::HypConstField() 
: fXmin(0),   
  fXmax(0),
  fYmin(0),
  fYmax(0),
  fZmin(0),
  fZmax(0),
  fPosX(0),
  fPosY(0),
  fPosZ(0),
  fAngleRot(0),
  fBx(0),
  fBy(0),
  fBz(0)
{
  fType = 0;
}
// -------------------------------------------------------------------------



// -----   Standard constructor   ------------------------------------------
HypConstField::HypConstField(const char* name, 
                             Double_t xMin, Double_t xMax, 
                             Double_t yMin, Double_t yMax, 
                             Double_t zMin, Double_t zMax, 
                             Double_t pos_x,Double_t pos_y,Double_t pos_z,Double_t angle,
                             Double_t bX, Double_t bY, Double_t bZ) 
    : FairField(name),
      fXmin(xMin),fXmax(xMax),
      fYmin(yMin),fYmax(yMax),
      fZmin(zMin),fZmax(zMax),
      fPosX(pos_x),fPosY(pos_y),fPosZ(pos_z),fAngleRot(angle),
      fBx(bX),fBy(bY),fBz(bZ)

{
  // convert T to kG : 1 T = 10 kG
  fBx*=10.;
  fBy*=10.;  
  fBz*=10.;

  fType = 0;
}
// -------------------------------------------------------------------------


// -----   Destructor   ----------------------------------------------------
HypConstField::~HypConstField() { }
// -------------------------------------------------------------------------



// -----   Set field region   ----------------------------------------------
void HypConstField::SetFieldRegion(Double_t xMin, Double_t xMax, 
                                   Double_t yMin, Double_t yMax, 
                                   Double_t zMin, Double_t zMax,
                                   Double_t pos_x,Double_t pos_y,Double_t pos_z,Double_t angle) 
{
  fXmin = xMin;
  fXmax = xMax;
  fYmin = yMin;
  fYmax = yMax;
  fZmin = zMin;
  fZmax = zMax;

  fPosX = pos_x;
  fPosY = pos_y;
  fPosZ = pos_z;
  
  fAngleRot = angle;

}
// -------------------------------------------------------------------------



// -----   Set field values   ----------------------------------------------
void HypConstField::SetField(Double_t bX, Double_t bY, Double_t bZ) {
  fBx   = bX;
  fBy   = bY;
  fBz   = bZ;

  // convert T to kG : 1 T = 10 kG
  fBx*=10.;
  fBy*=10.;  
  fBz*=10.;
}
// -------------------------------------------------------------------------



// -----   Get x component of field   --------------------------------------
Double_t HypConstField::GetBx(Double_t x, Double_t y, Double_t z) {
  
  
  if ( ! IsInside(x,y,z) ) 
    return 0.;
  else
    return fBx;
}
// -------------------------------------------------------------------------



// -----   Get y component of field   --------------------------------------
Double_t HypConstField::GetBy(Double_t x, Double_t y, Double_t z) {
  if ( !IsInside(x,y,z) ) 
    return 0.;
  else
    return fBy;
}
// -------------------------------------------------------------------------



// -----   Get z component of field   --------------------------------------
Double_t HypConstField::GetBz(Double_t x, Double_t y, Double_t z) {
  if ( !IsInside(x,y,z) ) 
    return 0.;
  else
    return fBz;
}
// -------------------------------------------------------------------------

Bool_t HypConstField::IsInside(Double_t x, Double_t y, Double_t z) 
{
  
  Double_t xl = x - fPosX;
  Double_t yl = y - fPosY;
  Double_t zl = z - fPosZ;
  
  Double_t xl2 = xl*TMath::Cos(fAngleRot)-zl*TMath::Sin(fAngleRot);
  Double_t yl2 = yl;
  Double_t zl2 = zl*TMath::Cos(fAngleRot)+xl*TMath::Sin(fAngleRot);
  
  // ---  Check for being outside the map range
  if ( ! ( xl2 >= fXmin && xl2 <= fXmax && yl2 >= fYmin && yl2 <= fYmax &&
          zl2 >= fZmin && zl2 <= fZmax ) ) 
  {
    return kFALSE;
  }
  else  
    return kTRUE;
  
}


// -----   Screen output   -------------------------------------------------
void HypConstField::Print(Option_t *option) const{
  cout << "======================================================" << endl;
  cout << "----  " << fTitle << " : " << fName << endl;
  cout << "----" << endl;
  cout << "----  Field type    : constant" << endl;
  cout << "----" << endl;
  cout << "----  Field regions : " << endl;
  cout << "----        x = " << setw(4) << fXmin << " to " << setw(4) 
       << fXmax << " cm" << endl;
  cout << "----        y = " << setw(4) << fYmin << " to " << setw(4) 
       << fYmax << " cm" << endl;
  cout << "----        z = " << setw(4) << fZmin << " to " << setw(4)
       << fZmax << " cm" << endl;
  cout << "---- Position = " << setw(4) << fPosX << " " << setw(4) << fPosY << " "<< setw(4) <<fPosZ<< " cm" << endl;
  cout << "---- Angle = " << setw(4) << fAngleRot << endl;
  cout.precision(4);
  cout << "----  B = ( " << fBx << ", " << fBy << ", " << fBz << " ) kG"
       << endl;
  cout << "======================================================" << endl;
}



ClassImp(HypConstField)
