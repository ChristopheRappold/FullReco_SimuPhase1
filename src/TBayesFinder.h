#ifndef TBAYESFINDER 
#define TBAYESFINDER 

#include "ReturnRes.hh"
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

class TBayesFinder final :  public TDataProcessInterface
{
  public :
  const THyphiAttributes& att;


  TBayesFinder(const THyphiAttributes& attr);
  ~TBayesFinder();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) final;
 private:
  int Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree) final;
  ReturnRes::InfoM SoftExit(int) final;
  void SelectHists() final;

  int FinderTrack(FullRecoEvent& RecoEvent);

  std::vector<double> radiusCDC;
  //TVector3 Plane_time;

  std::vector<TGeoHMatrix> MatMD;
  std::vector< std::vector<TGeoNodeMatrix*> > ME;
  std::vector< std::vector< BayesFind::DataLayer >> LayerGeo;
  TRandom3 rand;

  struct LocalHists
  {
    TH2F* h_xy;
    TH2F* h_PxPy;
    TH2F* h_xy_extrap;
    TH2F* h_PxPy_extrap;
    TH2F* h_TrackFindingStat;
    TH2F* h_SolenoidGeo[3];
  };
  LocalHists LocalHisto;


};


#endif
