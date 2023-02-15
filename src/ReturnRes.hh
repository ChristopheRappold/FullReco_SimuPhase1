#ifndef TRETURNRES
#define TRETURNRES

#include "msgpack.hpp"

namespace ReturnRes
{
  enum InfoM : int
    {
      Fine = 0,
      MultiS2_Start,
      StartTimingCut,
      ChamberHitLimit,
      TriggerCut,
      NoBeam,
      BuildError,
      PrimVtxError,
      KalmanError,
      DecayVtxError,
      EndRun,
      CloseRun,
      SIZEOF_INFOM
    };

};

MSGPACK_ADD_ENUM(ReturnRes::InfoM);

#endif
