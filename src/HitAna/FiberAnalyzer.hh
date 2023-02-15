#ifndef FIBER_ANALYZER_HH
#define FIBER_ANALYZER_HH
#include "FiberHitAna.hh"
#include "FiberTrackAna.hh"
#include "FiberHitXUV.hh"
#include "TVector3.h"

class FiberAnalyzer
{
  public:
    FiberAnalyzer();
    ~FiberAnalyzer();
    std::vector< std::vector< std::vector< FiberHitAna* > > > Clusterize(std::vector< std::vector< std::vector< FiberHitAna* > > > &cont);
    std::vector< std::vector< FiberHitXUV* > > FindHit( std::vector< std::vector< std::vector< FiberHitAna* > > > &cont, ParaManager *par);
    std::vector< FiberTrackAna* > DeleteDup( std::vector< FiberTrackAna* > &cont);
    std::vector< FiberTrackAna* > DeleteDupCombi(       std::vector< FiberTrackAna* > &cont);
    std::vector< FiberHitXUV* >   DeleteDupXUV(         std::vector< FiberHitXUV*   > &cont);
    TVector3 GetVertexPoint( const TVector3 & Xin, const TVector3 & Xout, const TVector3 & Pin, const TVector3 & Pout, double & dist );

};


#endif
