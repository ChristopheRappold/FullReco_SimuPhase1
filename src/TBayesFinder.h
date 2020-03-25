#ifndef TBAYESFINDER 
#define TBAYESFINDER 

#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"

#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TRandom3.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "THyphiAttributes.h"

typedef TDataProcess<FullRecoEvent,MCAnaEventG4Sol> TDataProcessInterface;

namespace BayesFind {

  struct DataLayer {
    double cenX = {0.};
    double cenY = {0.};
    double cenZ = {0.};
    		      
    double minX = {0.};
    double minY = {0.};
    double minZ = {0.};
    		      
    double maxX = {0.};
    double maxY = {0.};
    double maxZ = {0.};

    double radius = {0.};
  };

  struct DataTrack {
    TVector3 hit = TVector3(-999.,-999.,-999.);
    TVector3 mom = TVector3(-999.,-999.,-999.);
    int Id;
  };

};

class TBayesFinder :  public TDataProcessInterface
{
  public :
  const THyphiAttributes& att;


  TBayesFinder(const THyphiAttributes& attr);
  ~TBayesFinder();

  //int Init(Ana_Hist* h);
  int operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) override;
  int SoftExit(int) override;
  int FinderTrack(FullRecoEvent& RecoEvent);

  std::vector<double> radiusCDC;
  //TVector3 Plane_time;

  std::vector<TGeoHMatrix> MatMD;
  std::vector< std::vector<TGeoNodeMatrix*> > ME;
  std::vector< std::vector< BayesFind::DataLayer >> LayerGeo;
  TRandom3 rand;
};


#endif
