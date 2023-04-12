#ifndef FIBER_HIT_XUV_HH
#define FIBER_HIT_XUV_HH
#include "FiberHitAna.hh"

class FiberHitXUV
{
  public:
    FiberHitXUV(double _pos_x, double _pos_y, double _d, FiberHitAna *hit_x, FiberHitAna *hit_u, FiberHitAna *hit_v, int _id = 0);
    ~FiberHitXUV();
    void Print(void);
    double GetPosX(){ return pos_x;};
    double GetPosY(){ return pos_y;};
    double GetD(){    return d;    };
    int    GetID(){   return id;   };
    FiberHitAna* GetHitX(){return p_x;}
    FiberHitAna* GetHitU(){return p_u;}
    FiberHitAna* GetHitV(){return p_v;}
    FiberHitAna* GetHit0();
    FiberHitAna* GetHit1();
    FiberHitAna* GetHit2();
    FiberHitAna* GetHit(int _i);

  private:
    double pos_x = -9999.;
    double pos_y = -9999.;
    double d     = -9999.;
    int    id    = -9999;
    FiberHitAna *p_x;
    FiberHitAna *p_u;
    FiberHitAna *p_v;
};


#endif
