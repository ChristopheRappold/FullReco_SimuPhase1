#ifndef HypMeasurement_h
#define HypMeasurement_h

#include "AbsMeasurement.h"
#include "AbsHMatrix.h"


namespace genfit {

class HypMeasurement : public PlanarMeasurement {

 public:
  HypMeasurement(int nDim = 2);
  HypMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, double angle, int detId, int hitId, TrackPoint* trackPoint);

  virtual ~HypMeasurement() {;}

  virtual AbsMeasurement* clone() const {return new HypMeasurement(*this);}

  

  //virtual SharedPlanePtr constructPlane(const StateOnPlane& state) const;

  //virtual std::vector<MeasurementOnPlane*> constructMeasurementsOnPlane(const StateOnPlane& state) const;

  //virtual void setPlane(const SharedPlanePtr& physicalPlane, int planeId = -1) {physicalPlane_ = physicalPlane; planeId_ = planeId;}

  ClassDef(SpacepointMeasurement,1)
};

} /* End of namespace genfit */
/** @} */

#endif // HypMeasurement_h
