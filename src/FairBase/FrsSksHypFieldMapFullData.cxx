#include "TArrayF.h"
#include "FrsSksHypFieldMapFull.h"
#include "FrsSksHypFieldMapFullData.h"



// -------------   Default constructor  ----------------------------------
FrsSksHypFieldMapFullData::FrsSksHypFieldMapFullData()
 :fType(1),
  fXmin(0), 
  fXmax(0),
  fYmin(0), 
  fYmax(0),
  fZmin(0), 
  fZmax(0),
  fUnit(0),
  fNx(0), 
  fNy(0), 
  fNz(0),
  fBx(0), 
  fBy(0),
  fBz(0)
{

}
// ------------------------------------------------------------------------



// -------------   Standard constructor   ---------------------------------
FrsSksHypFieldMapFullData::FrsSksHypFieldMapFullData(const char* mapName)
  :TNamed(mapName, "Hyp Field MapFull Data"),
   fType(1),
   fXmin(0), 
   fXmax(0),
   fYmin(0), 
   fYmax(0),
   fZmin(0), 
   fZmax(0),
   fUnit(0),
   fNx(0), 
   fNy(0), 
   fNz(0),
   fBx(0), 
   fBy(0), 
   fBz(0)
{
}
// ------------------------------------------------------------------------



// -----   Constructor from FrsSksHypFieldMapFull   ------------------------------
FrsSksHypFieldMapFullData::FrsSksHypFieldMapFullData(const char* name,
				 const FrsSksHypFieldMapFull& map) 
  :TNamed(name, "Hyp Field MapFull Data"),
   fType( map.GetType()),
   fXmin( map.GetXmin()),
   fXmax( map.GetXmax()),
   fYmin( map.GetYmin()),
   fYmax( map.GetYmax()),
   fZmin( map.GetZmin()),
   fZmax( map.GetZmax()),
   fUnit( map.GetUnit()), 
   fNx(   map.GetNx()),
   fNy(   map.GetNy()),
   fNz(   map.GetNz()),
   fBx( new TArrayF(*(map.GetBx()))),
   fBy( new TArrayF(*(map.GetBy()))),
   fBz( new TArrayF(*(map.GetBz())))
 
{

  // Take out scaling factor and convert from kG to T
  Double_t factor = map.GetScale() ;//* 10.; 
  Int_t index = 0;
  for (Int_t ix=0; ix<fNx; ix++) {
    for (Int_t iy=0; iy<fNy; iy++) {
      for (Int_t iz=0; iz<fNz; iz++) {
	index = ix*fNy*fNz + iy*fNz + iz;
	if ( fBx ) (*fBx)[index] = (*fBx)[index] / factor;
	if ( fBy ) (*fBy)[index] = (*fBy)[index] / factor;
	if ( fBz ) (*fBz)[index] = (*fBz)[index] / factor;
      }    // z loop
    }      // y loop
  }        // x loop

}
// ------------------------------------------------------------------------



// ------------   Destructor   --------------------------------------------
FrsSksHypFieldMapFullData::~FrsSksHypFieldMapFullData() {
  if ( fBx ) delete fBx;
  if ( fBy ) delete fBy;
  if ( fBz ) delete fBz;
}
// ------------------------------------------------------------------------


ClassImp(FrsSksHypFieldMapFullData)
