// -------------------------------------------------------------------------
// -----                    FrsSolenoidHypField source file                  -----
// -----                Created 30/01/07  by M. Al/Turany              -----
// -------------------------------------------------------------------------

#include "FrsSolenoidHypField.h"
#include "TMath.h"
#include "TGeoTube.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"


#include <iomanip>
#include <iostream>


using namespace std;

// -----   Default constructor   -------------------------------------------
FrsSolenoidHypField::FrsSolenoidHypField() 
  : FairField("HypSolenoidField"),
    fRmin(0),   
    fRmax(0),
    fZmin(0),
    fZmax(0),
    fPosX(0),
    fPosY(0),
    fPosZ(0),
    fBx(0),
    fBy(0),
    fBz(0)
{
  fType = 0;
}
// -------------------------------------------------------------------------



// -----   Standard constructor   ------------------------------------------
FrsSolenoidHypField::FrsSolenoidHypField(const char* name, 
					 Double_t rMin, Double_t rMax, 
					 Double_t zMin, Double_t zMax, 
					 Double_t pos_x,Double_t pos_y,Double_t pos_z,
					 Double_t bX, Double_t bY, Double_t bZ) 
  : FairField(name),
    fRmin(rMin),fRmax(rMax),
    fZmin(zMin),fZmax(zMax),
    fPosX(pos_x),fPosY(pos_y),fPosZ(pos_z),
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
FrsSolenoidHypField::~FrsSolenoidHypField() { }
// -------------------------------------------------------------------------



// -----   Set field region   ----------------------------------------------
void FrsSolenoidHypField::SetFieldRegion(Double_t rMin, Double_t rMax, 
					 Double_t zMin, Double_t zMax,
					 Double_t pos_x,Double_t pos_y,Double_t pos_z) 
{
  fRmin = rMin;
  fRmax = rMax;
  fZmin = zMin;
  fZmax = zMax;

  fPosX = pos_x;
  fPosY = pos_y;
  fPosZ = pos_z;
  
}
// -------------------------------------------------------------------------



// -----   Set field values   ----------------------------------------------
void FrsSolenoidHypField::SetField(Double_t bX, Double_t bY, Double_t bZ) {
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
Double_t FrsSolenoidHypField::GetBx(Double_t x, Double_t y, Double_t z) {
  
  
  if ( ! IsInside(x,y,z) ) 
    return 0.;
  else
    return fBx;
}
// -------------------------------------------------------------------------



// -----   Get y component of field   --------------------------------------
Double_t FrsSolenoidHypField::GetBy(Double_t x, Double_t y, Double_t z) {
  if ( !IsInside(x,y,z) ) 
    return 0.;
  else
    return fBy;
}
// -------------------------------------------------------------------------



// -----   Get z component of field   --------------------------------------
Double_t FrsSolenoidHypField::GetBz(Double_t x, Double_t y, Double_t z) {
  if ( !IsInside(x,y,z) ) 
    return 0.;
  else
    return fBz;
}

void FrsSolenoidHypField::GetBxyz(const Double_t point[3], Double_t* bField)
{
  bField[0] = 0.;
  bField[1] = 0.;
  bField[2] = 0.;
  
  Bool_t FirstM = IsInside(point[0],point[1],point[2]) ;
   

  if ( FirstM ) 
    {
      bField[0] = fBx;
      bField[1] = fBy;
      bField[2] = fBz;
    }

    
}
// -------------------------------------------------------------------------


// -------------------------------------------------------------------------

Bool_t FrsSolenoidHypField::IsInside(Double_t x, Double_t y, Double_t z) 
{
  
  Double_t xl = x - fPosX;
  Double_t yl = y - fPosY;
  Double_t zl = z - fPosZ;
  
  Double_t rl = TMath::Sqrt(xl*xl+yl*yl);

  
  // ---  Check for being outside the map range
  if ( ! ( rl >= fRmin && rl <= fRmax && zl >= fZmin && zl <= fZmax ) ) 
    {
      return kFALSE;
    }
  else  
    return kTRUE;
  
}


// -----   Screen output   -------------------------------------------------
void FrsSolenoidHypField::Print(Option_t *option) const{
  cout << "======================================================" << endl;
  cout << "----  " << fTitle << " : " << fName << endl;
  cout << "----" << endl;
  cout << "----  Field type    : constant" << endl;
  cout << "----" << endl;
  cout << "----  Field regions : " << endl;
  cout << "----        r = " << setw(4) << fRmin << " to " << setw(4) 
       << fRmax << " cm" << endl;
  cout << "----        z = " << setw(4) << fZmin << " to " << setw(4)
       << fZmax << " cm" << endl;
  cout << "---- Position = " << setw(4) << fPosX << " " << setw(4) << fPosY << " "<< setw(4) <<fPosZ<< " cm" << endl;
  cout.precision(4);
  cout << "----  B = ( " << fBx << ", " << fBy << ", " << fBz << " ) kG"
       << endl;
  cout << "======================================================" << endl;
}


void FrsSolenoidHypField::SetPositionFromGeoManager(const TString& name_geo_mat, const TString& name_geo_field) {

  if(0==gGeoManager)
    {
      cout<<" E - No gGeoManager initialized ! "<<endl;
      Fatal("SetPosition","No gGeoManager");
    }
  
  TObjArray* listMatrices = gGeoManager->GetListOfMatrices();
  TGeoMatrix* MagneticField = (TGeoMatrix*) listMatrices->FindObject(name_geo_mat);

  TObjArray* listOfShape = gGeoManager->GetListOfShapes();
  TGeoTube* tempMag = (TGeoTube*) listOfShape->FindObject(name_geo_field);
  fRmin = tempMag->GetRmin();
  fRmax = tempMag->GetRmax();
  fZmin = -tempMag->GetDz();
  fZmax = -fZmin;

  Double_t lloc[]={0.,0.,0.};
  Double_t lmaster[]={0.,0.,0.};
  MagneticField->LocalToMaster(lloc,lmaster);
  fPosX = lmaster[0];
  fPosY = lmaster[1];
  fPosZ = lmaster[2];

  
}

void FrsSolenoidHypField::SetPositionFromGeoManager(const TString& name_node)
{

  if(gGeoManager==nullptr)
    {
      cout<<" E - No gGeoManager initialized ! "<<endl;
      Fatal("SetPositionFromGeoManager","No gGeoManager");
      return;
    }
  
  TGeoNode* mother = gGeoManager->GetTopNode();
  TObjArray* ListNodes = mother->GetNodes();
  TGeoNodeMatrix* Solenoid_node = dynamic_cast<TGeoNodeMatrix*>(ListNodes->FindObject(name_node));
  TGeoMatrix* tempMatrix = nullptr;
  if(Solenoid_node == nullptr && ListNodes->GetEntries()==1)
    {
      TGeoNodeMatrix* tempMother = dynamic_cast<TGeoNodeMatrix*>(ListNodes->At(0));
      tempMatrix = (TGeoMatrix*) tempMother->GetMatrix();
      
      ListNodes = tempMother->GetNodes();
      Solenoid_node = dynamic_cast<TGeoNodeMatrix*>(ListNodes->FindObject(name_node));
    }
  
  if(Solenoid_node == nullptr)
    {
      cout<<" E - No Solenoid node found ! "<<endl;
      Fatal("SetPositionFromGeoManager","No Solenoid_node found");
      return;
    }
    
  TGeoTube* tempMag = (TGeoTube*) Solenoid_node->GetVolume()->GetShape();
  TGeoMatrix* MagneticField = (TGeoMatrix*) Solenoid_node->GetMatrix();
  
  fRmin = tempMag->GetRmin();
  fRmax = tempMag->GetRmax();
  fZmin = -tempMag->GetDz();
  fZmax = -fZmin;

  Double_t lloc[]={0.,0.,0.};
  Double_t lmaster[]={0.,0.,0.};
  if(tempMatrix==nullptr)
    MagneticField->LocalToMaster(lloc,lmaster);
  else
    {
      TGeoHMatrix M1(*MagneticField), T1(*tempMatrix);
      TGeoHMatrix h = T1*M1;
      h.LocalToMaster(lloc,lmaster);
    }

  fPosX = lmaster[0];
  fPosY = lmaster[1];
  fPosZ = lmaster[2];

  
}




ClassImp(FrsSolenoidHypField)
