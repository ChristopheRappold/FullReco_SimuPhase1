// -------------------------------------------------------------------------
// -----                      HypFieldMap source file                  -----
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


#include "HypFieldMap.h"
#include "HypFieldMapData.h"

using namespace std;

// -------------   Default constructor  ----------------------------------
HypFieldMap::HypFieldMap() 
  : FairField(),
    fFileName(""),
    fScale(1.0),
    funit(1.0),
    fPosX(0), fPosY(0), fPosZ(0),
    fXmin(0), fXmax(0), fXstep(0),
    fYmin(0), fYmax(0), fYstep(0),
    fZmin(0), fZmax(0), fZstep(0),
    fAngleRot(0.),
    fNx(0),fNy(0),fNz(0),
    //fBx(NULL), fBy(NULL), fBz(NULL)   
    fBy(NULL)   
{
  SetName("");
  fType = 1;
}
// ------------------------------------------------------------------------



// -------------   Standard constructor   ---------------------------------
HypFieldMap::HypFieldMap(const char* mapName)//, const char* fileType)
  : FairField(mapName),
    fFileName(TString("")),
    fScale(1.0),
    funit(10.0), 
    fPosX(0), fPosY(0), fPosZ(0),
    fXmin(0), fXmax(0), fXstep(0),
    fYmin(0), fYmax(0), fYstep(0),
    fZmin(0), fZmax(0), fZstep(0),
    fAngleRot(0.),
    fNx(0),fNy(0),fNz(0),
    //fBx(NULL), fBy(NULL), fBz(NULL)   
    fBy(NULL)   
{
  //SetName(mapName);
  //TString dir = "./";
    
  
  fFileName = mapName;
  
  TObjArray* list_name_temp = fFileName.Tokenize("/");
  TObjString* name_temp_1 = static_cast<TObjString*>(list_name_temp->Last());
  TString name_temp_2(name_temp_1->String());
  
  TObjArray* list_name_temp_2 = name_temp_2.Tokenize(".");
  if(list_name_temp_2->GetEntries()!=2)
    cout<<"The name of the FieldMap contains more than one dot ! change name ! "<<fFileName<<endl;
  
  TObjString* name_temp_3 = static_cast<TObjString*>(list_name_temp_2->First());
  SetName(name_temp_3->String());
  
  //if ( fileType[0] == 'R' ) fFileName += ".root";
  //else                      fFileName += ".dat";
  fType = 1;
}
// ------------------------------------------------------------------------


/*
// ------------   Constructor from HypFieldPar   --------------------------
HypFieldMap::HypFieldMap(HypFieldPar* fieldPar) 
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
HypFieldMap::~HypFieldMap() {
  //if ( fBx ) delete fBx;
  if ( fBy ) delete fBy;
  //if ( fBz ) delete fBz;
}
// ------------------------------------------------------------------------



// -----------   Intialisation   ------------------------------------------
void HypFieldMap::Init() {
  if      (fFileName.EndsWith(".root"))     ReadRootFile(fFileName, fName);
  else if (fFileName.EndsWith(".dat"))  ReadAsciiFile(fFileName);
  else {
    cerr << "-E- HypFieldMap::Init: No proper file name defined! ("
	 << fFileName << ")" << endl;
    Fatal("Init", "No proper file name");
  }
}
// ------------------------------------------------------------------------



// -----------   Get x component of the field   ---------------------------
Double_t HypFieldMap::GetBx(Double_t x, Double_t y, Double_t z) 
{
  return 0.;
}  /*
  Int_t ix    = 0;
  Int_t iy    = 0;
  Int_t iz    = 0;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;

  if ( IsInside(x, y, z, ix, iy, iz, dx, dy, dz) ) {

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
  return Interpolate(dx, dy, dz);

  }

  return 0.;
}
// ------------------------------------------------------------------------
*/


// -----------   Get y component of the field   ---------------------------
Double_t HypFieldMap::GetBy(Double_t x, Double_t y, Double_t z) {

  Int_t ix    = 0;
  Int_t iy    = 0;
  Int_t iz    = 0;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;

  if ( IsInside(x, y, z, ix, iy, iz, dx, dy, dz) ) {

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
  return Interpolate(dx, dy, dz);

  }

  return 0.;
}
// ------------------------------------------------------------------------



// -----------   Get z component of the field   ---------------------------
Double_t HypFieldMap::GetBz(Double_t x, Double_t y, Double_t z) 
{
  return 0.;
}
  /*  Int_t ix    = 0;
  Int_t iy    = 0;
  Int_t iz    = 0;
  Double_t dx = 0.;
  Double_t dy = 0.;
  Double_t dz = 0.;

  if ( IsInside(x, y, z, ix, iy, iz, dx, dy, dz) ) {

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
  return Interpolate(dx, dy, dz);

  }

  return 0.;
}
// ------------------------------------------------------------------------
*/


// -----------   Check whether a point is inside the map   ----------------
Bool_t HypFieldMap::IsInside(Double_t x, Double_t y, Double_t z,
			     Int_t& ix, Int_t& iy, Int_t& iz,
			     Double_t& dx, Double_t& dy, Double_t& dz) {

  // --- Transform into local coordinate system
  Double_t xl = x - fPosX;
  Double_t yl = y - fPosY;
  Double_t zl = z - fPosZ;
  
  Double_t xl2 = xl*TMath::Cos(fAngleRot)-zl*TMath::Sin(fAngleRot);
  Double_t yl2 = yl;
  Double_t zl2 = zl*TMath::Cos(fAngleRot)+xl*TMath::Sin(fAngleRot);
  
  // ---  Check for being outside the map range
  if ( ! ( xl2 >= fXmin && xl2 <= fXmax && yl2 >= fYmin && yl2 <= fYmax &&
	   zl2 >= fZmin && zl2 <= fZmax ) ) {
    ix = iy = iz = 0;
    dx = dy = dz = 0.;
    return kFALSE;
  }
 
  // --- Determine grid cell
  ix = Int_t( xl2 / fXstep );
  iy = Int_t( yl2 / fYstep );
  iz = Int_t( zl2 / fZstep );


  // Relative distance from grid point (in units of cell size)
  dx = xl2 / fXstep - Double_t(ix);
  dy = yl2 / fYstep - Double_t(iy);
  dz = zl2 / fZstep - Double_t(iz);

  // --- Determine grid cell
  ix += Int_t(fNx/2);
  iy += Int_t(fNy/2);
  iz += Int_t(fNz/2);
  
  
  return kTRUE;

}
// ------------------------------------------------------------------------



// ----------   Write the map to an ASCII file   --------------------------
void HypFieldMap::WriteAsciiFile(const char* fileName) {

  // Open file
  cout << "-I- HypFieldMap: Writing field map to ASCII file " 
       << fileName << endl;
  ofstream mapFile(fileName);
  if ( ! mapFile.is_open() ) {
    cerr << "-E- HypFieldMap:ReadAsciiFile: Could not open file! " << endl;
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
  cout << "-I- HypFieldMap: " << fNx*fNy*fNz << " entries to write... " 
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
	mapFile << /*fBx->At(index)/factor << " " <<*/ fBy->At(index)/factor 
		/*<< " " << fBz->At(index)/factor*/ << endl;
      } // z-Loop
    }   // y-Loop
  }     // x-Loop
  cout << "   " << index+1 << " written" << endl;
  mapFile.close();		

}	
// ------------------------------------------------------------------------



// -------   Write field map to a ROOT file   -----------------------------
void HypFieldMap::WriteRootFile(const char* fileName,
				const char* mapName) {

  HypFieldMapData* data = new HypFieldMapData(mapName, *this);
  TFile* oldFile = gFile;
  TFile* file = new TFile(fileName, "RECREATE");
  data->Write();
  file->Close();
  if(oldFile) oldFile->cd();

}
// ------------------------------------------------------------------------



// -----  Set the position of the field centre in global coordinates  -----
void HypFieldMap::SetPosition(Double_t x, Double_t y, Double_t z,Double_t angle) {
  fPosX = x;
  fPosY = y;
  fPosZ = z;
  fAngleRot = angle;
}
// ------------------------------------------------------------------------



// ---------   Screen output   --------------------------------------------
void HypFieldMap::Print(Option_t *option) const {
  TString type = "Map";
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
  cout << "----  Field centre position: ( " << setw(6) << fPosX << ", "
       << setw(6) << fPosY << ", " << setw(6) << fPosZ << ") cm" << endl;
  cout << "----  Field scaling factor: " << fScale << endl;
  //Double_t bx = GetBx(0.,0.,0.);
  //Double_t by = GetBy(0.,0.,0.);
  //Double_t bz = GetBz(0.,0.,0.);
  //cout << "----" << endl;
  //cout << "----  Field at origin is ( " << setw(6) << by << ") kG" << endl;
  //cout << "----  Field at origin is ( " << setw(6) << by << ") T" << endl;
 cout << "======================================================" << endl;
}
// ------------------------------------------------------------------------  



// ---------    Reset parameters and data (private)  ----------------------
void HypFieldMap::Reset() {
  fPosX = fPosY = fPosZ = 0.;
  fXmin = fYmin = fZmin = 0.;
  fXmax = fYmax = fZmax = 0.;
  fXstep = fYstep = fZstep = 0.;
  fNx = fNy = fNz = 0;
  fScale = 1.;
  funit = 10.0;
  //if ( fBx ) { delete fBx; fBx = NULL; }
  if ( fBy ) { delete fBy; fBy = NULL; }
  //if ( fBz ) { delete fBz; fBz = NULL; }
}
// ------------------------------------------------------------------------  



// -----   Read field map from ASCII file (private)   ---------------------
void HypFieldMap::ReadAsciiFile(const char* fileName) {

  Double_t by=0.;//bx=0., by=0., bz=0.;
  Double_t  xx, yy, zz;
  // Open file
  cout << "-I- HypFieldMap: Reading field map from ASCII file " 
       << fileName << endl;
  ifstream mapFile(fileName);
  if ( ! mapFile.is_open() ) {
    cerr << "-E- HypFieldMap:ReadAsciiFile: Could not open file! " << endl;
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
    cout << "-E- HypFieldMap::ReadAsciiFile: Incompatible map types!"
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
    cout << "-E- FieldMap::ReadAsciiFile: No units!"
	 << endl;
        Fatal("ReadAsciiFile","No units defined");
  }


  // Read grid parameters
 
  mapFile >>fXmin >> fXmax >> fNx;
  mapFile >>fYmin >> fYmax >> fNy;
  mapFile >>fZmin >> fZmax >> fNz;
  fXstep = ( fXmax - fXmin ) / Double_t( fNx - 1 );
  fYstep = ( fYmax - fYmin ) / Double_t( fNy - 1 );
  fZstep = ( fZmax - fZmin ) / Double_t( fNz - 1 );
  
  // Create field arrays
  //fBx = new TArrayF(fNx * fNy * fNz);
  fBy = new TArrayF(fNx * fNy * fNz);
  //fBz = new TArrayF(fNx * fNy * fNz);

  // Read the field values
  Double_t factor = fScale * funit;   // Factor 1/1000 for G -> kG
  cout << right;
  Int_t nTot = fNx * fNy * fNz;
  cout << "-I- HypFieldMap: " << nTot << " entries to read... " 
       << setw(3) << 0 << " % ";
  Int_t index = 0;
  div_t modul;
  Int_t iDiv = TMath::Nint(nTot/100.);
  for (Int_t ix=0; ix<fNx; ix++) {
    for (Int_t iy = 0; iy<fNy; iy++) {
      for (Int_t iz = 0; iz<fNz; iz++) {
	if (! mapFile.good()) cerr << "-E- HypFieldMap::ReadAsciiFile: "
				   << "I/O Error at " << ix << " "
				   << iy << " " << iz << endl;
	index = ix*fNy*fNz + iy*fNz + iz;
	modul = div(index,iDiv);
	if ( modul.rem == 0 ) {
	  Double_t perc = TMath::Nint(100.*index/nTot);
	  cout << "\b\b\b\b\b\b" << setw(3) << perc << " % " << flush;
	}
	mapFile >> xx>>yy>>zz>>  /*bx >> */ by /* >> bz */;
	//mapFile >>  bx >> by >> bz ;
	//cout  << " x= " <<xx <<" y= " << yy<<" z= " << zz<<" bx= " <<  bx <<" by= " <<by <<" bz= " << bz<< endl;
  if(xx==0 && yy==0 && zz==0)
    cout  << " x= " <<xx <<" y= " << yy<<" z= " << zz<<" by= " <<by<< endl;
	//fBx->AddAt(factor*bx, index);
	fBy->AddAt(factor*by, index);
	//fBz->AddAt(factor*bz, index);
	if ( mapFile.eof() ) {
	  cerr << endl << "-E- HypFieldMap::ReadAsciiFile: EOF"
	       << " reached at " << ix << " " << iy << " " << iz << endl;
	  mapFile.close();
	  break;
	}
      }   // z-Loop
    }     // y-Loop0)
  }       // x-Loop

  cout << "   " << index+1 << " read" << endl;

  mapFile.close();

}
// ------------------------------------------------------------------------



// -------------   Read field map from ROOT file (private)  ---------------
void HypFieldMap::ReadRootFile(const char* fileName, 
			       const char* mapName) {

  // Store gFile pointer
  TFile* oldFile = gFile;

  // Open root file
  cout << "-I- HypFieldMap: Reading field map from ROOT file " 
       << fileName << endl; 
  TFile* file = new TFile(fileName, "READ");		
  if (file->IsZombie()) {
    cerr << "-E- HypFieldMap::ReadRootfile: Cannot read from file! " 
	 << endl;
    Fatal("ReadRootFile","Cannot read from file");
  }

  // Get the field data object
  HypFieldMapData* data = NULL;
  file->GetObject(mapName, data);
  if ( ! data ) {
    cout << "-E- HypFieldMap::ReadRootFile: data object " << mapName
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
void HypFieldMap::SetField(const HypFieldMapData* data) {

  // Check compatibility
  if ( data->GetType() != fType ) {
    cout << "-E- HypFieldMap::SetField: Incompatible map types!"
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
  fNx = data->GetNx();
  fNy = data->GetNy();
  fNz = data->GetNz();
  fXstep = ( fXmax - fXmin ) / Double_t( fNx - 1 );
  fYstep = ( fYmax - fYmin ) / Double_t( fNy - 1 );
  fZstep = ( fZmax - fZmin ) / Double_t( fNz - 1 );
  //if ( fBx ) delete fBx;
  if ( fBy ) delete fBy;
  //if ( fBz ) delete fBz;
  //fBx = new TArrayF(*(data->GetBx()));
  fBy = new TArrayF(*(data->GetBy()));
  //fBz = new TArrayF(*(data->GetBz()));

  std::cout<<"Loading FieldMap :"<<fXmin<<" "<<fXmax<<" "<<fYmin<<" "<<fYmax<<" "<<fZmin<<" "<<fZmax<<std::endl;

  // Scale and convert from G(or T) to kG
  Double_t factor = fScale * funit;
  Int_t index = 0;
  for (Int_t ix=0; ix<fNx; ix++) {
    for (Int_t iy=0; iy<fNy; iy++) {
      for (Int_t iz=0; iz<fNz; iz++) {
	index = ix*fNy*fNz + iy*fNz + iz;
	//if ( fBx ) (*fBx)[index] = (*fBx)[index] * factor;
	if ( fBy ) (*fBy)[index] = (*fBy)[index] * factor;
	//if ( fBz ) (*fBz)[index] = (*fBz)[index] * factor;
      }
    }
  }

}
// ------------------------------------------------------------------------  



// ------------   Interpolation in a grid cell (private)  -----------------
Double_t HypFieldMap::Interpolate(Double_t dx, Double_t dy, Double_t dz) {

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



ClassImp(HypFieldMap)
