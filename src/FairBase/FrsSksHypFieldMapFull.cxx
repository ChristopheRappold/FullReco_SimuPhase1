// -------------------------------------------------------------------------
// -----                      FrsSksHypFieldMapFull source file                  -----
// -----         Created 12/01/04  by M. Al/Turany (FairField.cxx)      -----
// -------------------------------------------------------------------------


#include <iomanip>
#include <iostream>
#include <fstream>
#include "stdlib.h"

#include "TArrayF.h"
#include "TFile.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TVector3.h"

#include "FrsSksHypFieldMapFull.h"
#include "FrsSksHypFieldMapFullData.h"

using namespace std;

// -------------   Default constructor  ----------------------------------
FrsSksHypFieldMapFull::FrsSksHypFieldMapFull() 
  : FairField(),
    fFileName(""),
    fScale(1.0),
    funit(1.0),
    fPosX(-1.10348092e+02), fPosY(1.63900368e+02), fPosZ(0),
    fAngleRotField(-20.*TMath::DegToRad()),
    fXmin(0), fXmax(0), fXstep(0),
    fYmin(0), fYmax(0), fYstep(0),
    fZmin(0), fZmax(0), fZstep(0),
    fAngleMagnetX(0), fAngleMagnetY(0),
    fPosMagnetRefX(0), fPosMagnetRefY(0), fPosMagnetRefZ(0),
    fPosFieldRefX(0), fPosFieldRefY(0), fPosFieldRefZ(0),
    fromGeoMatrix(false),secondMagnet(false),MagneticField(NULL),SecondMagneticField(NULL),SecondMFAll(NULL),
    fXminSecond(0),fXmaxSecond(0),
    fYminSecond(0),fYmaxSecond(0),
    fZminSecond(0),fZmaxSecond(0),
    fSecondBx(0),fSecondBy(0),fSecondBz(0),
    fNx(0),fNy(0),fNz(0),
    fBx(NULL), fBy(NULL), fBz(NULL)   
    //fBy(NULL)   
{
  SetName("");
  fType = 1;
}
// ------------------------------------------------------------------------



// -------------   Standard constructor   ---------------------------------
FrsSksHypFieldMapFull::FrsSksHypFieldMapFull(const char* mapName)//, const char* fileType)
  : FairField(mapName),
    fFileName(TString("")),
    fScale(1.0),
    funit(1.0), 
    fPosX(-1.10348092e+02), fPosY(1.63900368e+02), fPosZ(0),
    fAngleRotField(-20.*TMath::DegToRad()),
    fXmin(0), fXmax(0), fXstep(0),
    fYmin(0), fYmax(0), fYstep(0),
    fZmin(0), fZmax(0), fZstep(0),
    fAngleMagnetX(0), fAngleMagnetY(0),
    fPosMagnetRefX(0), fPosMagnetRefY(0), fPosMagnetRefZ(0),
    fPosFieldRefX(0), fPosFieldRefY(0), fPosFieldRefZ(0),
    fromGeoMatrix(false),secondMagnet(false),MagneticField(NULL),SecondMagneticField(NULL),SecondMFAll(NULL),
    fSecondBx(0),fSecondBy(0),fSecondBz(0),
    fNx(0),fNy(0),fNz(0),
    fBx(NULL), fBy(NULL), fBz(NULL)   
    //fBy(NULL)   
{
  //SetName(mapName);
  //TString dir = "./";
    
  
  fFileName = mapName;
  
  TObjArray* list_name_temp = fFileName.Tokenize("/");
  TObjString* name_temp_1 = static_cast<TObjString*>(list_name_temp->Last());
  TString name_temp_2(name_temp_1->String());
  
  TObjArray* list_name_temp_2 = name_temp_2.Tokenize(".");
  if(list_name_temp_2->GetEntries()!=2)
    cout<<"The name of the FieldMapFull contains more than one dot ! change name ! "<<fFileName<<endl;
  
  TObjString* name_temp_3 = static_cast<TObjString*>(list_name_temp_2->First());
  SetName(name_temp_3->String());
  
  //if ( fileType[0] == 'R' ) fFileName += ".root";
  //else                      fFileName += ".dat";
  fType = 3;
}
// ------------------------------------------------------------------------


/*
// ------------   Constructor from FrsSksHypFieldPar   --------------------------
FrsSksHypFieldMapFull::FrsSksHypFieldMapFull(FrsSksHypFieldPar* fieldPar) 
  : FairField(),
	fFileName(TString("")),
    fScale(1.0),
    funit(1.0), // to keep in T
    fPosX(0), fPosY(0), fPosZ(0),
    fXmin(0), fXmax(0), fXstep(0),
    fYmin(0), fYmax(0), fYstep(0),
    fZmin(0), fZmax(0), fZstep(0),
    fAngleRot(0.),
    fNx(0),fNy(0),fNz(0),
    //fBx(NULL), fBy(NULL), fBz(NULL)   
    fBy(NULL)   
{
  fType = 1;
  if ( ! fieldPar ) {
    cerr << "-W- HypConstField::HypConstField: empty parameter container!"
	 << endl;
    SetName("");
	fType     = -1;
  }
  else {
    TString Name=GetName();
    fieldPar->MapName(Name);
    fPosX  = fieldPar->GetPositionX();
    fPosY  = fieldPar->GetPositionY();
    fPosZ  = fieldPar->GetPositionZ();
    fScale = fieldPar->GetScale();
    fAngleRot = fieldPar->GetAngleRot();
    TString dir = "./";
    fFileName = dir + Name + ".root";
    fType = fieldPar->GetType();
  }
}
// ------------------------------------------------------------------------
*/


// ------------   Destructor   --------------------------------------------
FrsSksHypFieldMapFull::~FrsSksHypFieldMapFull() {
  if ( fBx ) delete fBx;
  if ( fBy ) delete fBy;
  if ( fBz ) delete fBz;
  if ( SecondMFAll ) delete SecondMFAll;

}
// ------------------------------------------------------------------------



// -----------   Intialisation   ------------------------------------------
void FrsSksHypFieldMapFull::Init() {
  if      (fFileName.EndsWith(".root"))     ReadRootFile(fFileName, fName);
  else if (fFileName.EndsWith(".dat"))  ReadAsciiFile(fFileName);
  else {
    cerr << "-E- FrsSksHypFieldMapFull::Init: No proper file name defined! ("
	 << fFileName << ")" << endl;
    Fatal("Init", "No proper file name");
  }
}

void FrsSksHypFieldMapFull::InitSecond(double Bx, double By, double Bz)
{
  fSecondBx = Bx;
  fSecondBy = By;
  fSecondBz = Bz;
}
// ------------------------------------------------------------------------



// -----------   Get x component of the field   ---------------------------
Double_t FrsSksHypFieldMapFull::GetBx(Double_t x, Double_t y, Double_t z) 
{
  return 0.;
  
  Int_t ix    = 0;
  Int_t iy    = 0;
  Int_t iz    = 0;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;

  Int_t ix2    = 0;
  Int_t iy2    = 0;
  Int_t iz2    = 0;
  Double_t dx2 = 0.;
  Double_t dy2 = 0.;
  Double_t dz2 = 0.;

  Bool_t FirstM =  IsInside(x, y, z, ix, iy, iz, dx, dy, dz);
  Bool_t SecondM =  IsInsideSecond(x, y, z, ix2, iy2, iz2, dx2, dy2, dz2);

  Double_t Temp_Bx = 0.;
  
  if (FirstM ) 
    {
      
      // Get Bx field values at grid cell corners
      fHa[0][0][0] = fBx->At(ix    *fNy*fNz + iy    *fNz + iz);
      fHa[1][0][0] = fBx->At((ix+1)*fNy*fNz + iy    *fNz + iz);
      fHa[0][1][0] = fBx->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
      fHa[1][1][0] = fBx->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
      fHa[0][0][1] = fBx->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
      fHa[1][0][1] = fBx->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
      fHa[0][1][1] = fBx->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
      fHa[1][1][1] = fBx->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));
      
      // Return interpolated field value
      //return Interpolate(dx, dy, dz);
      
      Temp_Bx += Interpolate(dx, dy, dz);
    }

  if(SecondM)
    {
      Temp_Bx += fSecondBx;
    }

  return Temp_Bx;
// =======
//   if ( IsInside(x, y, z, ix, iy, iz, dx, dy, dz) ) {

//   // Get Bx field values at grid cell corners
//   fHa[0][0][0] = fBx->At(ix    *fNy*fNz + iy    *fNz + iz);
//   fHa[1][0][0] = fBx->At((ix+1)*fNy*fNz + iy    *fNz + iz);
//   fHa[0][1][0] = fBx->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
//   fHa[1][1][0] = fBx->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
//   fHa[0][0][1] = fBx->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
//   fHa[1][0][1] = fBx->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
//   fHa[0][1][1] = fBx->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
//   fHa[1][1][1] = fBx->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));

//   // Return interpolated field value
//   return Interpolate(dx, dy, dz);

//   }

//   return 0.;
// >>>>>>> start git repo from code in home/gsi
}
// ------------------------------------------------------------------------



// -----------   Get y component of the field   ---------------------------
Double_t FrsSksHypFieldMapFull::GetBy(Double_t x, Double_t y, Double_t z) {

  Int_t ix    = 0;
  Int_t iy    = 0;
  Int_t iz    = 0;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;

  Int_t ix2    = 0;
  Int_t iy2    = 0;
  Int_t iz2    = 0;
  Double_t dx2 = 0.;
  Double_t dy2 = 0.;
  Double_t dz2 = 0.;

  Bool_t FirstM =  IsInside(x, y, z, ix, iy, iz, dx, dy, dz);
  Bool_t SecondM =  IsInsideSecond(x, y, z, ix2, iy2, iz2, dx2, dy2, dz2);

  Double_t Temp_By = 0.;

  if ( FirstM ) 
    {

      // Get By field values at grid cell corners
      fHa[0][0][0] = fBy->At(ix    *fNy*fNz + iy    *fNz + iz);
      fHa[1][0][0] = fBy->At((ix+1)*fNy*fNz + iy    *fNz + iz);
      fHa[0][1][0] = fBy->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
      fHa[1][1][0] = fBy->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
      fHa[0][0][1] = fBy->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
      fHa[1][0][1] = fBy->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
      fHa[0][1][1] = fBy->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
      fHa[1][1][1] = fBy->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));
      
      // Return interpolated field value
      //return Interpolate(dx, dy, dz);
      Temp_By += Interpolate(dx, dy, dz);
    }
  if( SecondM )
    {
      Temp_By += fSecondBy;
    }

  return Temp_By;
// =======
//   if ( IsInside(x, y, z, ix, iy, iz, dx, dy, dz) ) {

//   // Get By field values at grid cell corners
//   fHa[0][0][0] = fBy->At(ix    *fNy*fNz + iy    *fNz + iz);
//   fHa[1][0][0] = fBy->At((ix+1)*fNy*fNz + iy    *fNz + iz);
//   fHa[0][1][0] = fBy->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
//   fHa[1][1][0] = fBy->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
//   fHa[0][0][1] = fBy->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
//   fHa[1][0][1] = fBy->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
//   fHa[0][1][1] = fBy->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
//   fHa[1][1][1] = fBy->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));

//   // Return interpolated field value
//   return Interpolate(dx, dy, dz);

//   }

//   return 0.;
// >>>>>>> start git repo from code in home/gsi
}
// ------------------------------------------------------------------------



// -----------   Get z component of the field   ---------------------------
Double_t FrsSksHypFieldMapFull::GetBz(Double_t x, Double_t y, Double_t z) 
{
  
  Int_t ix    = 0;
  Int_t iy    = 0;
  Int_t iz    = 0;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;

  Int_t ix2    = 0;
  Int_t iy2    = 0;
  Int_t iz2    = 0;
  Double_t dx2 = 0.;
  Double_t dy2 = 0.;
  Double_t dz2 = 0.;

  Bool_t FirstM =  IsInside(x, y, z, ix, iy, iz, dx, dy, dz);
  Bool_t SecondM =  IsInsideSecond(x, y, z, ix2, iy2, iz2, dx2, dy2, dz2);

  Double_t Temp_Bz = 0.;

  if ( FirstM )
    {
      // Get Bz field values at grid cell corners
      fHa[0][0][0] = fBz->At(ix    *fNy*fNz + iy    *fNz + iz);
      fHa[1][0][0] = fBz->At((ix+1)*fNy*fNz + iy    *fNz + iz);
      fHa[0][1][0] = fBz->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
      fHa[1][1][0] = fBz->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
      fHa[0][0][1] = fBz->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
      fHa[1][0][1] = fBz->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
      fHa[0][1][1] = fBz->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
      fHa[1][1][1] = fBz->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));
      
      // Return interpolated field value
      //return Interpolate(dx, dy, dz);
      Temp_Bz += Interpolate(dx, dy, dz);
    }
  if( SecondM )
    {
      Temp_Bz += fSecondBz;
    }

  return Temp_Bz;
// =======
//   if ( IsInside(x, y, z, ix, iy, iz, dx, dy, dz) ) {

//   // Get Bz field values at grid cell corners
//   fHa[0][0][0] = fBz->At(ix    *fNy*fNz + iy    *fNz + iz);
//   fHa[1][0][0] = fBz->At((ix+1)*fNy*fNz + iy    *fNz + iz);
//   fHa[0][1][0] = fBz->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
//   fHa[1][1][0] = fBz->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
//   fHa[0][0][1] = fBz->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
//   fHa[1][0][1] = fBz->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
//   fHa[0][1][1] = fBz->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
//   fHa[1][1][1] = fBz->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));

//   // Return interpolated field value
//   return Interpolate(dx, dy, dz);

//   }

//   return 0.;
// >>>>>>> start git repo from code in home/gsi
}
// ------------------------------------------------------------------------



// -----------   Check whether a point is inside the map   ----------------
Bool_t FrsSksHypFieldMapFull::IsInside(Double_t x, Double_t y, Double_t z,
			     Int_t& ix, Int_t& iy, Int_t& iz,
			     Double_t& dx, Double_t& dy, Double_t& dz) {

 
  // Transform into the coordinate system of the sksplus magnet
  //tempPos->RotateY(-fAngleMagnetY);
  Double_t local[] = {0.,0.,0.};
  if(fromGeoMatrix==false)
    {
      Double_t temp_x0 = x - fPosMagnetRefX; 
      Double_t temp_y0 = y - fPosMagnetRefY; 
      Double_t temp_z0 = z - fPosMagnetRefZ; 
      
      Double_t s = TMath::Sin(-fAngleMagnetY);
      Double_t c = TMath::Cos(-fAngleMagnetY);
      
      Double_t temp_z1 = c*temp_z0 - s*temp_x0;
      Double_t temp_x1 = s*temp_z0 + c*temp_x0;
      Double_t temp_y1 = temp_y0;
      //tempPos->RotateX(-fAngleMagnetX);
 
      s = TMath::Sin(-fAngleMagnetX);
      c = TMath::Cos(-fAngleMagnetX);
      
      Double_t temp_y2 = c*temp_y1 - s*temp_z1;
      Double_t temp_z2 = s*temp_y1 + c*temp_z1;
      Double_t temp_x2 = temp_x1;
      
      
      // --- Transform into local coordinate system
      
      s = TMath::Sin(fAngleRotField);
      c = TMath::Cos(fAngleRotField);
      
      Double_t temp_xm0 = c*temp_x2 - s*temp_y2;
      Double_t temp_ym0 = s*temp_x2 + c*temp_y2;
      Double_t temp_zm0 = temp_z2;
      
      local[0] = temp_xm0 - fPosX;
      local[1] = temp_ym0 - fPosY;
      local[2] = temp_zm0 - fPosZ;
    }
  else
    {
      Double_t master[] = {x,y,z};
      MagneticField->MasterToLocal(master,local);
      local[0] += fPosX;
      local[1] += fPosY;
      local[2] += fPosZ;
    }
  // ---  Check for being outside the map range
  if ( ! ( local[0] >= fXmin && local[0] <= fXmax && local[1] >= fYmin && local[1] <= fYmax &&
	   local[2] >= fZmin && local[2] <= fZmax ) ) {
    ix = iy = iz = 0;
    dx = dy = dz = 0.;
    return kFALSE;
  }

  // --- Determine grid cell
  ix = Int_t( (local[0]-fXmin) / fXstep );
  iy = Int_t( (local[1]-fYmin) / fYstep );
  iz = Int_t( (local[2]-fZmin) / fZstep );


  // Relative distance from grid point (in units of cell size)
  dx = (local[0]-fXmin) / fXstep - Double_t(ix);
  dy = (local[1]-fYmin) / fYstep - Double_t(iy);
  dz = (local[2]-fZmin) / fZstep - Double_t(iz);


  //cout<<" I - "<<local[0]<<" "<<local[1]<<" "<<local[2]<<" "<<ix<<" "<<iy<<" "<<iz<<" "<<fXmin<<" "<<fYmin<<" "<<fZmin<<" | "<<fBx->At(ix*fNy*fNz+iy*fNz+iz)<<" "<<fBy->At(ix*fNy*fNz+iy*fNz+iz)<<" "<<fBz->At(ix*fNy*fNz+iy*fNz+iz)<<endl;

  // --- Determine grid cell
  // ix += Int_t(fNx/2);
  // iy += Int_t(fNy/2);
  // iz += Int_t(fNz/2);
  
  
  return kTRUE;

}
// ------------------------------------------------------------------------



// -----------   Check whether a point is inside the map   ----------------
Bool_t FrsSksHypFieldMapFull::IsInsideSecond(Double_t x, Double_t y, Double_t z,
					     Int_t& ix, Int_t& iy, Int_t& iz,
					     Double_t& dx, Double_t& dy, Double_t& dz) {

 
  if(0==SecondMFAll)
    {
      SecondMFAll = new TGeoHMatrix(*SecondMagneticField);
      SecondMFAll->Multiply(MagneticField);
    }


  Double_t local[] = {0.,0.,0.};
  if(secondMagnet==false)
    {
      return kFALSE;
    }
  else
    {
      Double_t master[] = {x,y,z};
      //SecondMagneticField->MasterToLocal(master,local);
      SecondMFAll->MasterToLocal(master,local);
      
      // local[0] += fPosX;
      // local[1] += fPosY;
      // local[2] += fPosZ;
    }
  // ---  Check for being outside the map range
  if ( ! ( local[0] >= fXminSecond && local[0] <= fXmaxSecond && local[1] >= fYminSecond && local[1] <= fYmaxSecond &&
	   local[2] >= fZminSecond && local[2] <= fZmaxSecond ) ) {
    ix = iy = iz = 0;
    dx = dy = dz = 0.;
    return kFALSE;
  }

  // --- Determine grid cell
  ix = Int_t( (local[0]-fXminSecond) / fXstep );
  iy = Int_t( (local[1]-fYminSecond) / fYstep );
  iz = Int_t( (local[2]-fZminSecond) / fZstep );


  // Relative distance from grid point (in units of cell size)
  dx = (local[0]-fXminSecond) / fXstep - Double_t(ix);
  dy = (local[1]-fYminSecond) / fYstep - Double_t(iy);
  dz = (local[2]-fZminSecond) / fZstep - Double_t(iz);


  //cout<<" I - "<<local[0]<<" "<<local[1]<<" "<<local[2]<<" "<<ix<<" "<<iy<<" "<<iz<<" "<<fXmin<<" "<<fYmin<<" "<<fZmin<<" | "<<fBx->At(ix*fNy*fNz+iy*fNz+iz)<<" "<<fBy->At(ix*fNy*fNz+iy*fNz+iz)<<" "<<fBz->At(ix*fNy*fNz+iy*fNz+iz)<<endl;

  // --- Determine grid cell
  // ix += Int_t(fNx/2);
  // iy += Int_t(fNy/2);
  // iz += Int_t(fNz/2);
  
  
  return kTRUE;

}
// ------------------------------------------------------------------------

// ----------   Write the map to an ASCII file   --------------------------
void FrsSksHypFieldMapFull::WriteAsciiFile(const char* fileName) {

  // Open file
  cout << "-I- FrsSksHypFieldMapFull: Writing field map to ASCII file " 
       << fileName << endl;
  ofstream mapFile(fileName);
  if ( ! mapFile.is_open() ) {
    cerr << "-E- FrsSksHypFieldMapFull:ReadAsciiFile: Could not open file! " << endl;
    return;
  }

  // Write field map grid parameters
  mapFile.precision(4);
  mapFile << showpoint;
  if ( fType == 1 ) mapFile << "nosym" << endl;
  if ( fType == 2 ) mapFile << "Solenoid" << endl;
  if ( fType == 3 ) mapFile << "Dipole" << endl;
  if ( fType == 4 ) mapFile << "Trans" << endl;
  if ( funit == 10.0  ) mapFile << "T" << endl;
  if ( funit == 0.001 ) mapFile << "G" << endl;
  if ( funit == 1.0   ) mapFile << "kG" << endl;

  mapFile << fXmin << " " << fXmax << " " << fNx << endl;
  mapFile << fYmin << " " << fYmax << " " << fNy << endl;
  mapFile << fZmin << " " << fZmax << " " << fNz << endl;

  // Write field values
  Double_t factor = funit * fScale;  // Takes out scaling 
  cout << right;
  Int_t nTot = fNx * fNy * fNz;
  cout << "-I- FrsSksHypFieldMapFull: " << fNx*fNy*fNz << " entries to write... " 
       << setw(3) << 0 << " % ";
  Int_t index=0;
  div_t modul;
  Int_t iDiv = TMath::Nint(nTot/100.);
  for(Int_t ix=0; ix<fNx; ix++) {
    for(Int_t iy=0; iy<fNy; iy++) {
      for(Int_t iz=0; iz<fNz; iz++) {
	index =ix*fNy*fNz + iy*fNz + iz;
	modul = div(index,iDiv);
	if ( modul.rem == 0 ) {
	  Double_t perc = TMath::Nint(100.*index/nTot);
	  cout << "\b\b\b\b\b\b" << setw(3) << perc << " % " << flush;
	}
	mapFile << fBx->At(index)/factor << " " <<fBy->At(index)/factor 
		<< " " << fBz->At(index)/factor << endl;
      } // z-Loop
    }   // y-Loop
  }     // x-Loop
  cout << "   " << index+1 << " written" << endl;
  mapFile.close();		

}	
// ------------------------------------------------------------------------



// -------   Write field map to a ROOT file   -----------------------------
void FrsSksHypFieldMapFull::WriteRootFile(const char* fileName,
				const char* mapName) {

  FrsSksHypFieldMapFullData* data = new FrsSksHypFieldMapFullData(mapName, *this);
  TFile* oldFile = gFile;
  TFile* file = new TFile(fileName, "RECREATE");
  data->Write();
  file->Close();
  if(oldFile) oldFile->cd();

}
// ------------------------------------------------------------------------



// -----  Set the position of the field centre in global coordinates  -----
void FrsSksHypFieldMapFull::SetPosition(Double_t x, Double_t y, Double_t z,Double_t angle_x,Double_t angle_y) {
  fAngleMagnetX = angle_x;
  fAngleMagnetY = angle_y;
  fPosMagnetRefX = x;
  fPosMagnetRefY = y;
  fPosMagnetRefZ = z;

}
// ------------------------------------------------------------------------

// -----  Set the position of the field centre in global coordinates  -----
void FrsSksHypFieldMapFull::SetPositionFromGeoManager(const TString& name_geo_field) {

  if(0==gGeoManager)
    {
      cout<<" E - No gGeoManager initialized ! "<<endl;
      Fatal("SetPosition","No gGeoManager");
    }
  
  TObjArray* listMatrices = gGeoManager->GetListOfMatrices();
  MagneticField = (TGeoMatrix*) listMatrices->FindObject(name_geo_field);
  fromGeoMatrix = true;
  Double_t lloc[]={0.,0.,0.};
  Double_t lmaster[]={0.,0.,0.};
  MagneticField->LocalToMaster(lloc,lmaster);
  fPosFieldRefX = lmaster[0];
  fPosFieldRefY = lmaster[1];
  fPosFieldRefZ = lmaster[2];
  Double_t lloc2[]={176.7,18.4,0.};
  
  MagneticField->LocalToMaster(lloc2,lmaster);
  fPosMagnetRefX = lmaster[0];
  fPosMagnetRefY = lmaster[1];
  fPosMagnetRefZ = lmaster[2];

  
  
}

// -----  Set the position of the field centre in global coordinates  -----
void FrsSksHypFieldMapFull::SetPositionSecondFromGeoManager(const TString& name_geo_mat, const TString& name_geo_field) {

  if(0==gGeoManager)
    {
      cout<<" E - No gGeoManager initialized ! "<<endl;
      Fatal("SetPosition","No gGeoManager");
    }
  
  TObjArray* listMatrices = gGeoManager->GetListOfMatrices();
  SecondMagneticField = (TGeoMatrix*) listMatrices->FindObject(name_geo_mat);

  secondMagnet = true;

  TObjArray* listOfShape = gGeoManager->GetListOfShapes();
  TGeoBBox* tempSecondMag = (TGeoBBox*) listOfShape->FindObject(name_geo_field);
  fXminSecond = -tempSecondMag->GetDX();
  fYminSecond = -tempSecondMag->GetDY();
  fZminSecond = -tempSecondMag->GetDZ();

  fXmaxSecond = -fXminSecond;
  fYmaxSecond = -fYminSecond;
  fZmaxSecond = -fZminSecond;

  if(0!=MagneticField)
    {
      SecondMFAll = new TGeoHMatrix(*SecondMagneticField);
      //cout<<" SetPos :"<<endl;
      //SecondMFAll->Print();
      SecondMFAll->MultiplyLeft(MagneticField);
      //cout<<" SetPos after Multi :"<<endl;
      //SecondMFAll->Print();
    }
  
}

// ------------------------------------------------------------------------



// ---------   Screen output   --------------------------------------------
void FrsSksHypFieldMapFull::Print(Option_t *option) const
{
  TString type = "MapFull";
  if ( fType == 2 ) type = "Soleniod Map ";
  if ( fType == 3 ) type = "Dipole Map ";
  if ( fType == 4 ) type = "Trans Map ";
  cout << "======================================================" << endl;
  cout.precision(4);
  cout << showpoint;
  cout << "----  " << fTitle << " : " << fName << endl;
  cout << "----" << endl;
  cout << "----  Field type     : " << type << endl;
  cout << "----" << endl;
  cout << "----  Field map grid : " << endl;
  cout << "----  x = " << setw(4) << fXmin << " to " << setw(4) << fXmax 
       << " cm, " << fNx << " grid points, dx = " << fXstep << " cm" << endl;
  cout << "----  y = " << setw(4) << fYmin << " to " << setw(4) << fYmax 
       << " cm, " << fNy << " grid points, dy = " << fYstep << " cm" << endl;
  cout << "----  z = " << setw(4) << fZmin << " to " << setw(4) << fZmax 
       << " cm, " << fNz << " grid points, dz = " << fZstep << " cm" << endl;
  cout << endl;
  if(fromGeoMatrix)
    cout << "----  Field centre position: ( " << setw(6) << fPosFieldRefX << ", " << setw(6) << fPosFieldRefY << ", " << setw(6) << fPosFieldRefZ << ") cm" << endl << "----  Field centre position2: ( " << setw(6) << fPosX << ", " << setw(6) << fPosY << ", " << setw(6) << fPosZ << ") cm" << endl;
  else
    cout << "----  Field centre position: ( " << setw(6) << fPosX << ", " << setw(6) << fPosY << ", " << setw(6) << fPosZ << ") cm" << endl;

  cout << "----  Magnet Ref position: ( " << setw(6) << fPosMagnetRefX << ", " << setw(6) << fPosMagnetRefY << ", " << setw(6) << fPosMagnetRefZ << ") cm" << endl;
  cout << "----  Field scaling factor: " << fScale << endl;

  if(fromGeoMatrix)
    MagneticField->Print();
  TVector3 temp_mag(0.,194.,0.);
  temp_mag.RotateZ(-20*TMath::DegToRad());
  Double_t fPosXmag(-1.10348092e+02), fPosYmag(1.63900368e+02);
  temp_mag-=TVector3(fPosX+fPosXmag,fPosY+fPosYmag,fPosZ);
  
  temp_mag.Print();
  Double_t toLoc[] = {temp_mag.X(),temp_mag.Y(),temp_mag.Z()};
  Double_t AtMas[] = {0.,0.,0.};

  if(fromGeoMatrix)
    MagneticField->LocalToMaster(toLoc,AtMas);

  // cout<<" MM : "<<AtMas[0]<<" "<<AtMas[1]<<" "<<AtMas[2]<<endl;
  // Double_t Bxyz[] = {0.,0.,0.};

  // GetBxyz(AtMas,Bxyz);
  
  // cout << "----" << endl;
  // cout << "----  Field at entrance magnet is ( " << setw(6) << Bxyz[0] << " "<< Bxyz[1] <<" "<< Bxyz[2] << ") kG" << endl;


  temp_mag.SetXYZ(0.,194.,0.);
  temp_mag.RotateZ(-45*TMath::DegToRad());
  
  temp_mag-=TVector3(fPosX+fPosXmag,fPosY+fPosYmag,fPosZ);
  
  temp_mag.Print();
  Double_t toLoc2[] = {temp_mag.X(),temp_mag.Y(),temp_mag.Z()};
  
  if(fromGeoMatrix)
    MagneticField->LocalToMaster(toLoc2,AtMas);

  // cout<<" MM : "<<AtMas[0]<<" "<<AtMas[1]<<" "<<AtMas[2]<<endl;

  // GetBxyz(AtMas,Bxyz);

  // cout << "----" << endl;
  // cout << "----  Field at exit magnet is ( " << setw(6) << Bxyz[0] << " "<< Bxyz[1] <<" "<< Bxyz[2] << ") kG" << endl;

  if(secondMagnet)
    {
      cout<<" Second Field Size : ["<<fXminSecond<<" "<<fXmaxSecond<<"] ["<<fYminSecond<<" "<<fYmaxSecond<<"] ["<<fZminSecond<<" "<<fZmaxSecond<<"]"<<endl; 
      SecondMagneticField->Print();

    }
  if(SecondMFAll)
    {
      cout<< " Second All Matrix :"<<endl;
      SecondMFAll->Print();
    }
  
  //cout << "----  Field at origin is ( " << setw(6) << by << ") T" << endl;
 cout << "======================================================" << endl;
}
// ------------------------------------------------------------------------  



// ---------    Reset parameters and data (private)  ----------------------
void FrsSksHypFieldMapFull::Reset() {
  fPosX = fPosY = fPosZ = 0.;
  fXmin = fYmin = fZmin = 0.;
  fXmax = fYmax = fZmax = 0.;
  fXstep = fYstep = fZstep = 0.;
  fNx = fNy = fNz = 0;
  fScale = 1.;
  funit = 10.0;
  if ( fBx ) { delete fBx; fBx = NULL; }
  if ( fBy ) { delete fBy; fBy = NULL; }
  if ( fBz ) { delete fBz; fBz = NULL; }
}
// ------------------------------------------------------------------------  



// -----   Read field map from ASCII file (private)   ---------------------
void FrsSksHypFieldMapFull::ReadAsciiFile(const char* fileName) {

  Double_t bx=0., by=0., bz=0.;
  Double_t xx, yy, zz;
  // Open file
  cout << "-I- FrsSksHypFieldMapFull: Reading field map from ASCII file " 
       << fileName << endl;
  ifstream mapFile(fileName);
  if ( ! mapFile.is_open() ) {
    cerr << "-E- FrsSksHypFieldMapFull:ReadAsciiFile: Could not open file! " << endl;
    Fatal("ReadAsciiFile","Could not open file");
  }

  // Read map type
  TString type;
  mapFile >> type;

  Int_t iType = 0;
  if ( type == "nosym" ) iType = 1;
  if ( type == "Solenoid") iType = 2;
  if ( type == "Dipole"  ) iType = 3;
  if ( type == "Trans"  ) iType = 4;
  if ( fType != iType ) {
    cout << "-E- FrsSksHypFieldMapFull::ReadAsciiFile: Incompatible map types!"
	 << endl;
    cout << "    Field map is of type " << fType 
	 << " but map on file is of type " << iType << endl;
    Fatal("ReadAsciiFile","Incompatible map types");
  }
  // Read Units
  TString unit;
  mapFile >> unit;
  if ( unit == "G" ) funit = 0.001;
  else if ( unit == "T"  ) funit = 10.0;
  else if ( unit == "kG"  ) funit=1.0;
  else {
    cout << "-E- FieldMapFull::ReadAsciiFile: No units!"
	 << endl;
        Fatal("ReadAsciiFile","No units defined");
  }


  // Read grid parameters
  mapFile >> fNx >> fNy >> fNz ;
  mapFile >> fXmin >> fYmin >> fZmin ;
  mapFile >> fXstep >> fYstep >> fZstep;

  // mapFile >>fXmin >> fXmax >> fNx;
  // mapFile >>fYmin >> fYmax >> fNy;
  // mapFile >>fZmin >> fZmax >> fNz;
  // fXstep = ( fXmax - fXmin ) / Double_t( fNx - 1 );
  // fYstep = ( fYmax - fYmin ) / Double_t( fNy - 1 );
  // fZstep = ( fZmax - fZmin ) / Double_t( fNz - 1 );

  fXmax = (Double_t(fNx -1))* fXstep + fXmin ;
  fYmax = (Double_t(fNy -1))* fYstep + fYmin ;
  fZmax = (Double_t(fNz -1))* fZstep + fZmin ;

  fPosX = 0.5*(fXmin+fXmax);
  fPosY = 0.5*(fYmin+fYmax);
  fPosZ = 0.5*(fZmin+fZmax);

  // Create field arrays
  fBx = new TArrayF(fNx * fNy * fNz);
  fBy = new TArrayF(fNx * fNy * fNz);
  fBz = new TArrayF(fNx * fNy * fNz);
  // Read the field values
  Double_t factor = fScale * funit;   // Factor 1/1000 for G -> kG
  cout << right;
  Int_t nTot = fNx * fNy * fNz;
  cout << "-I- FrsSksHypFieldMapFull: " << nTot << " entries to read... " 
       << setw(3) << 0 << " % ";
  Int_t index = 0;
  div_t modul;
  Int_t iDiv = TMath::Nint(nTot/100.);


   for (Int_t ix=0; ix<fNx; ix++) 
     for (Int_t iy = 0; iy<fNy; iy++) 
       for (Int_t iz = 0; iz<fNz; iz++) 
	 {
	   index = ix*fNy*fNz + iy*fNz + iz;
	   fBx->AddAt(0.0, index);
	   fBy->AddAt(0.0, index);
	   fBz->AddAt(0.0, index);
	 }


  // for (Int_t ix=0; ix<fNx; ix++) {
  //   for (Int_t iy = 0; iy<fNy; iy++) {
  //     for (Int_t iz = 0; iz<fNz; iz++) {
  while(mapFile >> xx >> yy >> zz >> bx >> by >> bz)
    {

      
      int ix = (xx-fXmin)/fXstep;
      int iy = (yy-fYmin)/fYstep;
      int iz = (zz-fZmin)/fZstep;
      
      if (! mapFile.good()) 
	cerr << "-E- FrsSksHypFieldMapFull::ReadAsciiFile: " << "I/O Error at " << ix << " " << iy << " " << iz << endl;

      index = ix*fNy*fNz + iy*fNz + iz;
      modul = div(index,iDiv);
      if ( modul.rem == 0 ) 
	{
	  Double_t perc = TMath::Nint(100.*index/nTot);
	  cout << "\b\b\b\b\b\b" << setw(3) << perc << " % " << flush;
	}
      //mapFile >> xx>>yy>>zz>>  bx >> by >> bz ;
      //mapFile >>  bx >> by >> bz ;
      //cout  << " x= " <<xx <<" y= " << yy<<" z= " << zz<<" b= " <<by<< endl;
      if(TMath::Abs(xx-176.)<1e-3 && TMath::Abs(yy-18.)<1e-3 && TMath::Abs(zz)<1e-3)
	cout  << " x= " <<xx <<" y= " << yy<<" z= " << zz<<" bx= " <<  bx <<" by= " <<by <<" bz= " << bz<< endl;
      fBx->AddAt(factor*bx, index);
      fBy->AddAt(factor*by, index);
      fBz->AddAt(factor*bz, index);
      if ( mapFile.eof() ) 
	{
	  cerr << endl << "-E- FrsSksHypFieldMapFull::ReadAsciiFile: EOF"
	       << " reached at " << ix << " " << iy << " " << iz << endl;
	  mapFile.close();
	  break;
	}
    }
  cout << "   " << index+1 << " read" << endl;

  mapFile.close();

}
// ------------------------------------------------------------------------



// -------------   Read field map from ROOT file (private)  ---------------
void FrsSksHypFieldMapFull::ReadRootFile(const char* fileName, 
			       const char* mapName) {

  // Store gFile pointer
  TFile* oldFile = gFile;

  // Open root file
  cout << "-I- FrsSksHypFieldMapFull: Reading field map from ROOT file " 
       << fileName << endl; 
  TFile* file = new TFile(fileName, "READ");		
  if (file->IsZombie()) {
    cerr << "-E- FrsSksHypFieldMapFull::ReadRootfile: Cannot read from file! " 
	 << endl;
    Fatal("ReadRootFile","Cannot read from file");
  }

  // Get the field data object
  FrsSksHypFieldMapFullData* data = NULL;
  file->GetObject(mapName, data);
  if ( ! data ) {
    cout << "-E- FrsSksHypFieldMapFull::ReadRootFile: data object " << mapName
	 << " not found in file! " << endl;
    exit(-1);
  }

  // Get the field parameters
  SetField(data);

  // Close the root file and delete the data object
  file->Close();
  delete data;
  if ( oldFile ) oldFile->cd();
  delete file;
}
// ------------------------------------------------------------------------



// ------------   Set field parameters and data (private)  ----------------
void FrsSksHypFieldMapFull::SetField(const FrsSksHypFieldMapFullData* data) {

  // Check compatibility
  if ( data->GetType() != fType ) {
    cout << "-E- FrsSksHypFieldMapFull::SetField: Incompatible map types!"
	 << endl;
    cout << "    Field map is of type " << fType 
	 << " but map on file is of type " << data->GetType() << endl;
    Fatal("SetField","Incompatible map types");
  }
  
  
  fXmin = data->GetXmin();
  fYmin = data->GetYmin();
  fZmin = data->GetZmin();
  fXmax = data->GetXmax();
  fYmax = data->GetYmax();
  fZmax = data->GetZmax();

  fPosX = 0.5*(fXmin+fXmax);
  fPosY = 0.5*(fYmin+fYmax);
  fPosZ = 0.5*(fZmin+fZmax);

  fNx = data->GetNx();
  fNy = data->GetNy();
  fNz = data->GetNz();
  fXstep = ( fXmax - fXmin ) / Double_t( fNx - 1 );
  fYstep = ( fYmax - fYmin ) / Double_t( fNy - 1 );
  fZstep = ( fZmax - fZmin ) / Double_t( fNz - 1 );
  if ( fBx ) delete fBx;
  if ( fBy ) delete fBy;
  if ( fBz ) delete fBz;
  fBx = new TArrayF(*(data->GetBx()));
  fBy = new TArrayF(*(data->GetBy()));
  fBz = new TArrayF(*(data->GetBz()));

  std::cout<<"Loading FieldMapFull :"<<fXmin<<" "<<fXmax<<" "<<fYmin<<" "<<fYmax<<" "<<fZmin<<" "<<fZmax<<std::endl;

  // Scale and convert from G(or T) to kG
  Double_t factor = fScale * funit;
  Int_t index = 0;
  for (Int_t ix=0; ix<fNx; ix++) {
    for (Int_t iy=0; iy<fNy; iy++) {
      for (Int_t iz=0; iz<fNz; iz++) {
	index = ix*fNy*fNz + iy*fNz + iz;
	if ( fBx ) (*fBx)[index] = (*fBx)[index] * factor;
	if ( fBy ) (*fBy)[index] = (*fBy)[index] * factor;
	if ( fBz ) (*fBz)[index] = (*fBz)[index] * factor;
      }
    }
  }

}
// ------------------------------------------------------------------------  



// ------------   Interpolation in a grid cell (private)  -----------------
Double_t FrsSksHypFieldMapFull::Interpolate(Double_t dx, Double_t dy, Double_t dz) {

  // Interpolate in x coordinate
  fHb[0][0] = fHa[0][0][0] + ( fHa[1][0][0]-fHa[0][0][0] ) * dx;
  fHb[1][0] = fHa[0][1][0] + ( fHa[1][1][0]-fHa[0][1][0] ) * dx;
  fHb[0][1] = fHa[0][0][1] + ( fHa[1][0][1]-fHa[0][0][1] ) * dx;
  fHb[1][1] = fHa[0][1][1] + ( fHa[1][1][1]-fHa[0][1][1] ) * dx;

  // Interpolate in y coordinate
  fHc[0] = fHb[0][0] + ( fHb[1][0] - fHb[0][0] ) * dy;
  fHc[1] = fHb[0][1] + ( fHb[1][1] - fHb[0][1] ) * dy;

  // Interpolate in z coordinate
  return fHc[0] + ( fHc[1] - fHc[0] ) * dz;

}
// ------------------------------------------------------------------------

void FrsSksHypFieldMapFull::GetBxyz(const Double_t point[3], Double_t* bField)
{
    Int_t ix    = 0;
    Int_t iy    = 0;
    Int_t iz    = 0;
    Double_t dx = 0.;
    Double_t dy = 0.;
    Double_t dz = 0.;

    Int_t ix2    = 0;
    Int_t iy2    = 0;
    Int_t iz2    = 0;
    Double_t dx2 = 0.;
    Double_t dy2 = 0.;
    Double_t dz2 = 0.;

    Bool_t FirstM = IsInside(point[0],point[1],point[2], ix, iy, iz, dx, dy, dz) ;
    Bool_t SecondM = IsInsideSecond(point[0],point[1],point[2], ix2, iy2, iz2, dx2, dy2, dz2) ;
    
    bField[0] = 0.;
    bField[1] = 0.;
    bField[2] = 0.;

    if ( FirstM )
      //if ( IsInside(point[0],point[1],point[2], ix, iy, iz, dx, dy, dz) ) 
      {
	// Get Bx field values at grid cell corners
	fHa[0][0][0] = fBx->At(ix    *fNy*fNz + iy    *fNz + iz);
	fHa[1][0][0] = fBx->At((ix+1)*fNy*fNz + iy    *fNz + iz);
	fHa[0][1][0] = fBx->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
	fHa[1][1][0] = fBx->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
	fHa[0][0][1] = fBx->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
	fHa[1][0][1] = fBx->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
	fHa[0][1][1] = fBx->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
	fHa[1][1][1] = fBx->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));
	
	// Return interpolated field value
	bField[0] += Interpolate(dx, dy, dz);
	//bField[0] = Interpolate(dx, dy, dz);
	
	// Get By field values at grid cell corners
	fHa[0][0][0] = fBy->At(ix    *fNy*fNz + iy    *fNz + iz);
	fHa[1][0][0] = fBy->At((ix+1)*fNy*fNz + iy    *fNz + iz);
	fHa[0][1][0] = fBy->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
	fHa[1][1][0] = fBy->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
	fHa[0][0][1] = fBy->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
	fHa[1][0][1] = fBy->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
	fHa[0][1][1] = fBy->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
	fHa[1][1][1] = fBy->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));
	
	// Return interpolated field value
	bField[1] += Interpolate(dx, dy, dz);
	//bField[1] = Interpolate(dx, dy, dz);

	  // Get Bz field values at grid cell corners
	fHa[0][0][0] = fBz->At(ix    *fNy*fNz + iy    *fNz + iz);
	fHa[1][0][0] = fBz->At((ix+1)*fNy*fNz + iy    *fNz + iz);
	fHa[0][1][0] = fBz->At(ix    *fNy*fNz + (iy+1)*fNz + iz);
	fHa[1][1][0] = fBz->At((ix+1)*fNy*fNz + (iy+1)*fNz + iz);
	fHa[0][0][1] = fBz->At(ix    *fNy*fNz + iy    *fNz + (iz+1));
	fHa[1][0][1] = fBz->At((ix+1)*fNy*fNz + iy    *fNz + (iz+1));
	fHa[0][1][1] = fBz->At(ix    *fNy*fNz + (iy+1)*fNz + (iz+1));
	fHa[1][1][1] = fBz->At((ix+1)*fNy*fNz + (iy+1)*fNz + (iz+1));
	
	// Return interpolated field value
	bField[2] += Interpolate(dx, dy, dz);
      }

    if(SecondM)
      {
	bField[0] += fSecondBx;
	bField[1] += fSecondBy;
	bField[2] += fSecondBz;

	//cout<<"Inside SecondM : ("<<point[0]<<", "<<point[1]<<", "<<point[2]<<") :"<<bField[1]<<" "<<bField[0]<<" "<<bField[2]<<endl;
// =======
// 	bField[2] = Interpolate(dx, dy, dz);


//       }
//     else
//       {
// 	bField[0] = 0.;
// 	bField[1] = 0.;
// 	bField[2] = 0.;
// >>>>>>> start git repo from code in home/gsi
      }

    if(fromGeoMatrix)
      {
	Double_t newbField[] = {0.,0.,0.};
	MagneticField->LocalToMasterVect(bField,newbField);
	bField[0] = newbField[0]; 
	bField[1] = newbField[1];
	bField[2] = newbField[2];
      }

}


ClassImp(FrsSksHypFieldMapFull)
