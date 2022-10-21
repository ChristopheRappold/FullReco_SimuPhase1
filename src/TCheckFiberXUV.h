#ifndef TCHECKFIBERXUV
#define TCHECKFIBERXUV

#include "ReturnRes.hh"
#include "TDataProcess.h"
#include "Debug.hh"
#include "FullRecoEvent.hh"
#include "THyphiAttributes.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"

#include "TFile.h"
//#include "TRandom3.h"
#include "TTree.h"
//#include "TH1I.h"
//#include "TH2I.h"
#include "TVector2.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>
#include <cmath>

template<class Out>
using TDataProcessInterface = TDataProcess<FullRecoEvent,Out>;

template<class Out>
class TCheckFiberXUV final :  public TDataProcessInterface<Out>
{
  public :
  const THyphiAttributes& att;

  TCheckFiberXUV(const THyphiAttributes& attr);
  ~TCheckFiberXUV();

  //int Init(Ana_Hist* h);
  void InitMT() final;
  ReturnRes::InfoM operator() (FullRecoEvent& RecoEvent,Out* OutTree) override;
 private:
  int Exec(FullRecoEvent& RecoEvent,Out* OutTree) override;

  ReturnRes::InfoM SoftExit(int) override;
  void SelectHists() final;
  int CheckHitXUV(const FullRecoEvent& RecoEvent);

  void FindHitXUV(genfit::AbsMeasurement* hitx, genfit::AbsMeasurement* hitu, genfit::AbsMeasurement* hitv, int id_det, TVector2* HitXY);


  std::vector<std::vector<TVector2> >     vect_HitXY = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<TVector2> > vect_realHitXY = {{}, {}, {}, {}, {}, {}, {}};
  std::vector<std::vector<double> >  vect_realHitAng = {{}, {}, {}, {}, {}, {}, {}};

  std::vector<double> cut_d = {};

  double ang[7][3] = {
    {  0.,  30., -30.}, /* UFT1 */
    {  0.,  30., -30.}, /* UFT2 */
    {  0.,  30., -30.}, /* UFT3 */
    {  0., -60.,  60.}, /* MFT1 */
    {  0.,  60., -60.}, /* MFT2 */
    { 30., -30.,   0.}, /* DFT1 */
    {  0.,  30., -30.}}; /* DFT2 */

  double EnergyThreshold     = 0.; // in MeV


  struct LocalHists
  {
    TH1F* h_ResidualFiberHitX[7];
    TH1F* h_ResidualFiberHitY[7];
    TH2F* h_ResidualFiberHitXY[7];
    TH2F* h_ResidualFiberHitX_Angle[7];
    TH2F* h_ResidualFiberHitY_Angle[7];
  };
  LocalHists LocalHisto;
};


#endif
