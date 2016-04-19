#include "Ana_EventNew_v14.hh"
#include "THypernucleus.hh"
//#include "THyphiSimpleHit.hh"
#include "THyphiHitDet.hh"
//#include "TTrackSimple_v2.hh"
#include "THyphiTrack_v2.hh"

using namespace std;

//#define DEBUG_DAISUKE
ClassImp(Ana_Event)

TClonesArray *Ana_Event::gMC_Particle = 0;
//TClonesArray *Ana_Event::gTR0 = 0;
TClonesArray *Ana_Event::gfTrack = 0;
TClonesArray *Ana_Event::gfHyp = 0;

TClonesArray *Ana_Event::gTR0x = 0;
TClonesArray *Ana_Event::gTR0y = 0;
TClonesArray *Ana_Event::gTR1x = 0;
TClonesArray *Ana_Event::gTR1y = 0;
TClonesArray *Ana_Event::gTR2x = 0;
TClonesArray *Ana_Event::gTR2y = 0;

TClonesArray *Ana_Event::gDC1x = 0;
TClonesArray *Ana_Event::gDC1xp= 0;
TClonesArray *Ana_Event::gDC1u = 0;
TClonesArray *Ana_Event::gDC1up= 0;
TClonesArray *Ana_Event::gDC1v = 0;
TClonesArray *Ana_Event::gDC1vp= 0;

TClonesArray *Ana_Event::gDC2x = 0;
TClonesArray *Ana_Event::gDC2xp= 0;
TClonesArray *Ana_Event::gDC2y = 0;
TClonesArray *Ana_Event::gDC2yp= 0;
TClonesArray *Ana_Event::gDC2u = 0;

TClonesArray *Ana_Event::gTOFs = 0;
TClonesArray *Ana_Event::gTOFp = 0;
TClonesArray *Ana_Event::gBg = 0;

Ana_Event::Ana_Event()
{
   if(!gMC_Particle) gMC_Particle = new TClonesArray("TMcParticle",30);
   //if(!gTR0) gTR0 = new TClonesArray("THyphiSimpleHit",40);
   //if(!gTR0) gTR0 = new TClonesArray("THyphiHitDet",40);
   if(!gTR0x) gTR0x = new TClonesArray("THyphiHitDet",40);
   if(!gTR0y) gTR0y = new TClonesArray("THyphiHitDet",40);
   if(!gTR1x) gTR1x = new TClonesArray("THyphiHitDet",40);
   if(!gTR1y) gTR1y = new TClonesArray("THyphiHitDet",40);
   if(!gTR2x) gTR2x = new TClonesArray("THyphiHitDet",40);
   if(!gTR2y) gTR2y = new TClonesArray("THyphiHitDet",40);

   if(!gDC1x) gDC1x = new TClonesArray("THyphiHitDet",40);
   if(!gDC1xp) gDC1xp = new TClonesArray("THyphiHitDet",40);
   if(!gDC1u) gDC1u = new TClonesArray("THyphiHitDet",40);
   if(!gDC1up) gDC1up = new TClonesArray("THyphiHitDet",40);
   if(!gDC1v) gDC1v = new TClonesArray("THyphiHitDet",40);
   if(!gDC1vp) gDC1vp = new TClonesArray("THyphiHitDet",40);

   if(!gDC2x) gDC2x = new TClonesArray("THyphiHitDet",40);
   if(!gDC2xp) gDC2xp = new TClonesArray("THyphiHitDet",40);
   if(!gDC2y) gDC2y = new TClonesArray("THyphiHitDet",40);
   if(!gDC2yp) gDC2yp = new TClonesArray("THyphiHitDet",40);
   if(!gDC2u) gDC2u = new TClonesArray("THyphiHitDet",40);

   if(!gTOFs) gTOFs = new TClonesArray("THyphiHitDet",40);
   if(!gTOFp) gTOFp = new TClonesArray("THyphiHitDet",40);
   if(!gBg) gBg = new TClonesArray("THyphiHitDet",40);

   //if(!gfTrack) gfTrack = new TClonesArray("TTrackSimple",30);
   if(!gfTrack) gfTrack = new TClonesArray("THyphiTrack",30);
   if(!gfHyp) gfHyp = new TClonesArray("THypernucleus",20);


   //   TR0 = gTR0;
   TR0x = gTR0x;
   TR0y = gTR0y;
   TR1x = gTR1x;
   TR1y = gTR1y;
   TR2x = gTR2x;
   TR2y = gTR2y;
   DC1x = gDC1x;
   DC1xp = gDC1xp;
   DC1u = gDC1u;
   DC1up = gDC1up;
   DC1v = gDC1v;
   DC1vp = gDC1vp;
   DC2x = gDC2x;
   DC2xp = gDC2xp;
   DC2y = gDC2y;
   DC2yp = gDC2yp;
   DC2u = gDC2u;

   TOFs = gTOFs;
   TOFp = gTOFp;
   Bg = gBg;

   fMC_Particle = gMC_Particle;
   fTrack = gfTrack;
   fHyp = gfHyp;

   Nmc=0;
   //   Nhit_tr0=0;
   Nhit_tr0x=0;
   Nhit_tr0y=0;
   Nhit_tr1x=0;
   Nhit_tr1y=0;
   Nhit_tr2x=0;
   Nhit_tr2y=0;
   Nhit_dc1x=0;
   Nhit_dc1xp=0;
   Nhit_dc1u=0;
   Nhit_dc1up=0;
   Nhit_dc1v=0;
   Nhit_dc1vp=0;

   Nhit_dc2x=0;
   Nhit_dc2xp=0;
   Nhit_dc2y=0;
   Nhit_dc2yp=0;
   Nhit_dc2u=0;

   Nhit_tofs=0;
   Nhit_tofp=0;
   Nhit_bg=0;

   Ntracks=0;
   Nhyp=0;
   trigger = 0;
   Setup();
}

Ana_Event::~Ana_Event()
{
  Clear();
  Reset();
}

void Ana_Event::Clear(Option_t *option)
{
#ifdef DEBUG_DAISUKE
  std::cout<<"Ana_Event Clear start "<<std::endl;
#endif
   fMC_Particle->Clear("C");
   //   TR0->Clear("C");
#ifdef DEBUG_DAISUKE
   std::cout<<"Ana_Event Clear line "<<__LINE__<<std::endl;
#endif

   TR0x->Clear("C");
   TR0y->Clear("C");
   TR1x->Clear("C");
   TR1y->Clear("C");
   TR2x->Clear("C");
   TR2y->Clear("C");
#ifdef DEBUG_DAISUKE
   std::cout<<"Ana_Event Clear line "<<__LINE__<<std::endl;
#endif

   DC1x->Clear("C");
   DC1xp->Clear("C");
   DC1u->Clear("C");
   DC1up->Clear("C");
   DC1v->Clear("C");
   DC1vp->Clear("C");
#ifdef DEBUG_DAISUKE
   std::cout<<"Ana_Event Clear line "<<__LINE__<<std::endl;
#endif

   DC2x->Clear("C");
   DC2xp->Clear("C");
   DC2y->Clear("C");
   DC2yp->Clear("C");
   DC2u->Clear("C");

   TOFs->Clear("C");
   TOFp->Clear("C");
   Bg->Clear("C");

   fTrack->Clear("C");
   fHyp->Clear("C");
#ifdef DEBUG_DAISUKE
   std::cout<<"Ana_Event Clear line "<<__LINE__<<std::endl;
#endif

   Setup();
#ifdef DEBUG_DAISUKE
  std::cout<<"Ana_Event Clear finished "<<std::endl;
#endif
}

void Ana_Event::Reset()
{
#ifdef DEBUG_DAISUKE
  std::cout<<"Ana_Event Reset start "<<gfTrack<<" "<<std::endl;
#endif
  delete gMC_Particle; gMC_Particle = 0;
  delete gfTrack; gfTrack = 0;
  //  delete gTR0; gTR0 = 0;
  delete gTR0x; gTR0x = 0;
  delete gTR0y; gTR0y = 0;
  delete gTR1x; gTR1x = 0;
  delete gTR1y; gTR1y = 0;
  delete gTR2x; gTR2x = 0;
  delete gTR2y; gTR2y = 0;
  delete gDC1x; gDC1x = 0;
  delete gDC1xp; gDC1xp = 0;
  delete gDC1u; gDC1u = 0;
  delete gDC1up; gDC1up = 0;
  delete gDC1v; gDC1v = 0;
  delete gDC1vp; gDC1vp = 0;
  delete gDC2x; gDC2x = 0;
  delete gDC2xp; gDC2xp = 0;
  delete gDC2y; gDC2y = 0;
  delete gDC2yp; gDC2yp = 0;
  delete gDC2u; gDC2u = 0;

  delete gTOFs; gTOFs = 0;
  delete gTOFp; gTOFp = 0;
  delete gBg; gBg = 0;

  delete gfHyp; gfHyp = 0;
#ifdef DEBUG_DAISUKE
  std::cout<<"Ana_Event Reset finished"<<std::endl;
#endif

}


int Ana_Event::Setup()
{
#ifdef DEBUG_DAISUKE
  std::cout<<"Ana_Event Setup start"<<std::endl;
#endif

  Nmc=0;
  //  Nhit_tr0=0;
   Nhit_tr0x=0;
   Nhit_tr0y=0;
   Nhit_tr1x=0;
   Nhit_tr1y=0;
   Nhit_tr2x=0;
   Nhit_tr2y=0;
   Nhit_dc1x=0;
   Nhit_dc1xp=0;
   Nhit_dc1u=0;
   Nhit_dc1up=0;
   Nhit_dc1v=0;
   Nhit_dc1vp=0;

   Nhit_dc2x=0;
   Nhit_dc2xp=0;
   Nhit_dc2y=0;
   Nhit_dc2yp=0;
   Nhit_dc2u=0;

   Nhit_tofs=0;
   Nhit_tofp=0;
   Nhit_bg=0;

  Ntracks=0;
  Nhyp=0;
  trigger = 0;

#ifdef DEBUG_DAISUKE
  std::cout<<"Ana_Event Setup finished"<<std::endl;
#endif

  return 0;
}

int Ana_Event::Add_MC(const TMcParticle& M)
{
   TClonesArray &MC_Particle_Ref = *fMC_Particle;
   
   new(MC_Particle_Ref[Nmc]) TMcParticle(M);
   Nmc++;

   int status=0;
   
   return status;
   
}
//int Ana_Event::Add_HitTr0(const THyphiSimpleHit& H)
// int Ana_Event::Add_HitTr0(const THyphiHitDet& H)
// {
//    TClonesArray &HitTr0_Ref = *TR0;
   
//    new(HitTr0_Ref[Nhit_tr0]) THyphiHitDet(H);
//    Nhit_tr0++;
   
//    int status=0;
   
//    return status;
// }

int Ana_Event::Add_Hit(const THyphiHitDet& H,TString detector)
{

  std::string funcname ="Add_Hit(const THyphiHitDet& H,TString detector) :  ";
  std::string name_detector(detector.Data());
  funcname += name_detector;
  if(name_detector=="TR0x")
    {
      TClonesArray &HitTr0x_Ref = *TR0x;
      
      new(HitTr0x_Ref[Nhit_tr0x]) THyphiHitDet(H);
      Nhit_tr0x++;
    }
  else if(name_detector=="TR0y")
    {
      TClonesArray &HitTr0y_Ref = *TR0y;
      
      new(HitTr0y_Ref[Nhit_tr0y]) THyphiHitDet(H);
      Nhit_tr0y++;
    }
  else if(name_detector=="TR1x")
    {
      TClonesArray &HitTr1x_Ref = *TR1x;
      
      new(HitTr1x_Ref[Nhit_tr1x]) THyphiHitDet(H);
      Nhit_tr1x++;
    }
  else if(name_detector=="TR1y")
    {
      TClonesArray &HitTr1y_Ref = *TR1y;
      
      new(HitTr1y_Ref[Nhit_tr1y]) THyphiHitDet(H);
      Nhit_tr1y++;
    }
  else if(name_detector=="TR2x")
    {
      TClonesArray &HitTr2x_Ref = *TR2x;
      
      new(HitTr2x_Ref[Nhit_tr2x]) THyphiHitDet(H);
      Nhit_tr2x++;
    }
  else if(name_detector=="TR2y")
    {
      TClonesArray &HitTr2y_Ref = *TR2y;
      
      new(HitTr2y_Ref[Nhit_tr2y]) THyphiHitDet(H);
      Nhit_tr2y++;
    }
   
  else if(name_detector=="DC1x")
    {
      TClonesArray &HitDC1x_Ref = *DC1x;
      
      new(HitDC1x_Ref[Nhit_dc1x]) THyphiHitDet(H);
      Nhit_dc1x++;
    }
  else if(name_detector=="DC1xp")
    {
      TClonesArray &HitDC1xp_Ref = *DC1xp;
      
      new(HitDC1xp_Ref[Nhit_dc1xp]) THyphiHitDet(H);
      Nhit_dc1xp++;
    }
  else if(name_detector=="DC1u")
    {
      TClonesArray &HitDC1u_Ref = *DC1u;
      
      new(HitDC1u_Ref[Nhit_dc1u]) THyphiHitDet(H);
      Nhit_dc1u++;
    }
  else if(name_detector=="DC1up")
    {
      TClonesArray &HitDC1up_Ref = *DC1up;
      
      new(HitDC1up_Ref[Nhit_dc1up]) THyphiHitDet(H);
      Nhit_dc1up++;
    }
  else if(name_detector=="DC1v")
    {
      TClonesArray &HitDC1v_Ref = *DC1v;
      
      new(HitDC1v_Ref[Nhit_dc1v]) THyphiHitDet(H);
      Nhit_dc1v++;
    }
  else if(name_detector=="DC1vp")
    {
      TClonesArray &HitDC1vp_Ref = *DC1vp;
      
      new(HitDC1vp_Ref[Nhit_dc1vp]) THyphiHitDet(H);
      Nhit_dc1vp++;
    }
  else if(name_detector=="DC2x")
    {
      TClonesArray &HitDC2x_Ref = *DC2x;
      
      new(HitDC2x_Ref[Nhit_dc2x]) THyphiHitDet(H);
      Nhit_dc2x++;
    }
  else if(name_detector=="DC2xp")
    {
      TClonesArray &HitDC2xp_Ref = *DC2xp;
      
      new(HitDC2xp_Ref[Nhit_dc2xp]) THyphiHitDet(H);
      Nhit_dc2xp++;
    }
  else if(name_detector=="DC2y")
    {
      TClonesArray &HitDC2y_Ref = *DC2y;
      
      new(HitDC2y_Ref[Nhit_dc2y]) THyphiHitDet(H);
      Nhit_dc2y++;
    }
  else if(name_detector=="DC2yp")
    {
      TClonesArray &HitDC2yp_Ref = *DC2yp;
      
      new(HitDC2yp_Ref[Nhit_dc2yp]) THyphiHitDet(H);
      Nhit_dc2yp++;
    }
  else if(name_detector=="DC2u")
    {
      TClonesArray &HitDC2u_Ref = *DC2u;
      
      new(HitDC2u_Ref[Nhit_dc2u]) THyphiHitDet(H);
      Nhit_dc2u++;
    }
  else if(name_detector=="TOFs")
    {
      TClonesArray &HitTOFs_Ref = *TOFs;
      
      new(HitTOFs_Ref[Nhit_tofs]) THyphiHitDet(H);
      Nhit_tofs++;
    }
  else if(name_detector=="TOFp")
    {
      TClonesArray &HitTOFp_Ref = *TOFp;
      
      new(HitTOFp_Ref[Nhit_tofp]) THyphiHitDet(H);
      Nhit_tofp++;
    }
  else if(name_detector=="Bg")
    {
      TClonesArray &HitBg_Ref = *Bg;
      
      new(HitBg_Ref[Nhit_bg]) THyphiHitDet(H);
      Nhit_bg++;
    }
  else
    {
      std::cout<<funcname<<" wrond detector"<<std::endl;
      
    }

   int status=0;
   
   return status;
}


//int Ana_Event::Add_Track(const TTrackSimple& T)
int Ana_Event::Add_Track(const THyphiTrack& T)
{
   TClonesArray &Track_Ref = *fTrack;
   
   //new(Track_Ref[Ntracks]) TTrackSimple(T);
   new(Track_Ref[Ntracks]) THyphiTrack(T);
   Ntracks++;
   
   int status=0;
   
   return status;
}

int Ana_Event::Add_Hyp(const THypernucleus& H)
{
   TClonesArray &Hyp_Ref = *fHyp;
   
   
   new(Hyp_Ref[Nhyp]) THypernucleus(H);
   Nhyp++;
   
   int status=0;
   
   return status;
}

