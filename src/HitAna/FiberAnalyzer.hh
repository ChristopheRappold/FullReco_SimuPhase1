#ifndef FIBER_ANALYZER_HH
#define FIBER_ANALYZER_HH
#include "FiberHitAna.hh"
#include "FiberTrackAna.hh"
#include "FiberHitXUV.hh"
#include "TrackHit.hh"
#include "ParaManager.hh"
#include "TVector3.h"
#include "AbsMeasurement.h"

class FiberAnalyzer
{
  public:
    FiberAnalyzer();
    ~FiberAnalyzer();
    std::vector< std::vector< std::vector< FiberHitAna* > > > Clusterize(std::vector< std::vector< std::vector< FiberHitAna* > > > &cont);
    std::vector< std::vector< FiberHitXUV* > > FindHit( std::vector< std::vector< std::vector< FiberHitAna* > > > &cont, ParaManager *par);
    std::vector< FiberTrackAna* > DeleteDup(               std::vector< FiberTrackAna* > &cont);
    std::vector< FiberTrackAna* > DeleteDupCombi(          std::vector< FiberTrackAna* > &cont);
    std::vector< FiberTrackAna* > DeleteSame(              std::vector< FiberTrackAna* > &cont);
    std::vector< FiberTrackAna* > DeleteInclusive(         std::vector< FiberTrackAna* > &cont);
    std::vector< TrackHit* >      DeleteDupTrackHit(       std::vector< TrackHit*      > &cont);
    std::vector< TrackHit* >      DeleteDupTrackHitMDC(    std::vector< TrackHit*      > &cont);
    std::vector< TrackHit* >      DeleteInclusiveTrackHit( std::vector< TrackHit*      > &cont);
    std::vector< FiberHitXUV* >   DeleteDupXUV(            std::vector< FiberHitXUV*   > &cont);
    TVector3 GetVertexPoint( const TVector3 & Xin, const TVector3 & Xout, const TVector3 & Pin, const TVector3 & Pout, double & dist );
    double CalcPhiDif(const double phi1, const double phi2);
    double GetPSB_R(  int _seg);
    double GetPSB_Phi(int _seg);
    std::map< std::string, std::vector<FiberTrackAna*> > FiberTracking(std::vector< std::vector< std::vector< FiberHitAna* > > >& FiberHitClCont,
                                                          ParaManager *par, std::vector<std::unique_ptr<genfit::AbsMeasurement> >& ListHits_PSB);

};

#endif
