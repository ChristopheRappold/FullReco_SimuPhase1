// -------------------------------------------------------------------------
// -----                    WasaSolenoidFieldMap source file                  -----
// -----                Created 30/01/07  by M. Al/Turany              -----
// -------------------------------------------------------------------------

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TMath.h"
#include "WasaSolenoidFieldMap.h"

#include <iomanip>
#include <iostream>

using namespace std;

// -----   Default constructor   -------------------------------------------
WasaSolenoidFieldMap::WasaSolenoidFieldMap()
  : FairField("WasaSolenoidFieldMap"), fRmin(0), fRmax(0), fZmin(0), fZmax(0),nameField(""),maxField(0.),signDir(1.0),fPosX(0), fPosY(0), fPosZ(0)

{
  fType = 2;
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
WasaSolenoidFieldMap::WasaSolenoidFieldMap(const char* name, const char* title, const TString& nameF, double maxF, double signD)
  : FairField(name,title), fRmin(0), fRmax(0), fZmin(0), fZmax(0) , nameField(nameF), maxField(10.*maxF), signDir(signD), fPosX(0), fPosY(0), fPosZ(0)

{
//  convert T to kG : 1 T = 10 kG : maxField = 10.*maxF because maxF is in Tesla;
  fType = 2;
  FieldMap = std::make_unique<MField>(name);
  FieldMap->ReadParameter(nameField);

}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
WasaSolenoidFieldMap::~WasaSolenoidFieldMap() {}
// -------------------------------------------------------------------------

void WasaSolenoidFieldMap::SetMaxField(double maxF)
{
//   // convert T to kG : 1 T = 10 kG
  maxField = 10.*maxF; // maxF -> in Tesla;
}

void WasaSolenoidFieldMap::Init()
{
  FieldMap->SetScale(std::fabs(maxField));
  FieldMap->InitializeParameter();

  fRmin = FieldMap->MinMax_R[0];
  fRmax = FieldMap->MinMax_R[1];
  fZmin = FieldMap->MinMax_Z[0];
  fZmax = FieldMap->MinMax_Z[1];

  //Double_t originField[3] = {0.,0.,0.};
  Double_t originField[3] = {fPosX,fPosY,fPosZ};
  GetFieldValue(originField, fieldAtOrig);
}

// -----   Get x component of field   --------------------------------------
Double_t WasaSolenoidFieldMap::GetBx(Double_t x, Double_t y, Double_t z)
{
  Double_t p[3] = {x,y,z};
  Double_t f[3] = {0.,0.,0.};
  GetBxyz(p,f);
  return f[0];
}
// -------------------------------------------------------------------------

// -----   Get y component of field   --------------------------------------
Double_t WasaSolenoidFieldMap::GetBy(Double_t x, Double_t y, Double_t z)
{
  Double_t p[3] = {x,y,z};
  Double_t f[3] = {0.,0.,0.};
  GetBxyz(p,f);
  return f[1];
}
// -------------------------------------------------------------------------

// -----   Get z component of field   --------------------------------------
Double_t WasaSolenoidFieldMap::GetBz(Double_t x, Double_t y, Double_t z)
{
  Double_t p[3] = {x,y,z};
  Double_t f[3] = {0.,0.,0.};
  GetBxyz(p,f);
  return f[2];
}

void WasaSolenoidFieldMap::GetBxyz(const Double_t point[3], Double_t* bField)
{
  bField[0] = 0.;
  bField[1] = 0.;
  bField[2] = 0.;

  Double_t pos[3];

  pos[0] = signDir*(point[0] - fPosX);
  pos[1] = (point[1] - fPosY);
  pos[2] = signDir*(point[2] - fPosZ);

  Double_t R = TMath::Hypot(pos[0],pos[1]);

  // ---  Check for being outside the map range
  if(!(R >= fRmin && R <= fRmax && pos[2] >= fZmin && pos[2] <= fZmax))
    return;

  FieldMap->Wfld(pos,bField); // it is already in kGauss
  bField[0] *= signDir;
  bField[2] *= signDir;

}


// -----   Screen output   -------------------------------------------------
void WasaSolenoidFieldMap::Print(Option_t* option) const
{
  cout << "======================================================" << endl;
  cout << "----  " << fTitle << " : " << fName << endl;
  cout << "----" << endl;
  cout << "----  Field type    : Field Map" << endl;
  cout << "----" << endl;
  cout << "----  Field regions : " << endl;
  cout << "----  sign Direction : "<<signDir<< endl;
  cout << "----        r = " << setw(4) << fRmin << " to " << setw(4) << fRmax << " cm" << endl;
  cout << "----        z = " << setw(4) << fZmin << " to " << setw(4) << fZmax << " cm" << endl;
  cout << "---- Position = " << setw(4) << fPosX << " " << setw(4) << fPosY << " " << setw(4) << fPosZ << " cm" << endl;
  cout.precision(4);
  cout << "----  B = ( " << fieldAtOrig[0] << ", " << fieldAtOrig[1] << ", " << fieldAtOrig[2] << " ) kG" << endl;
  cout << "======================================================" << endl;
}

void WasaSolenoidFieldMap::SetPositionFromGeoManager(const TString& name_node)
{

  if(gGeoManager == nullptr)
    {
      cout << " E - No gGeoManager initialized ! " << endl;
      Fatal("SetPositionFromGeoManager", "No gGeoManager");
      return;
    }

  TGeoNode* mother              = gGeoManager->GetTopNode();
  TObjArray* ListNodes          = mother->GetNodes();
  TGeoNodeMatrix* Solenoid_node = dynamic_cast<TGeoNodeMatrix*>(ListNodes->FindObject(name_node));
  TGeoMatrix* tempMatrix        = nullptr;
  if(Solenoid_node == nullptr && ListNodes->GetEntries() == 1)
    {
      TGeoNodeMatrix* tempMother = dynamic_cast<TGeoNodeMatrix*>(ListNodes->At(0));
      tempMatrix                 = (TGeoMatrix*)tempMother->GetMatrix();

      ListNodes     = tempMother->GetNodes();
      Solenoid_node = dynamic_cast<TGeoNodeMatrix*>(ListNodes->FindObject(name_node));
    }

  if(Solenoid_node == nullptr)
    {
      cout << " E - No Solenoid node found ! " << endl;
      Fatal("SetPositionFromGeoManager", "No Solenoid_node found");
      return;
    }

  TGeoTube* tempMag         = (TGeoTube*)Solenoid_node->GetVolume()->GetShape();
  TGeoMatrix* MagneticField = (TGeoMatrix*)Solenoid_node->GetMatrix();

  Double_t lloc[]    = {0., 0., 0.};
  Double_t lmaster[] = {0., 0., 0.};
  if(tempMatrix == nullptr)
    MagneticField->LocalToMaster(lloc, lmaster);
  else
    {
      TGeoHMatrix M1(*MagneticField), T1(*tempMatrix);
      TGeoHMatrix h = T1 * M1;
      h.LocalToMaster(lloc, lmaster);
    }

  fPosX = lmaster[0];
  fPosY = lmaster[1];
  fPosZ = lmaster[2];
}

ClassImp(WasaSolenoidFieldMap)
