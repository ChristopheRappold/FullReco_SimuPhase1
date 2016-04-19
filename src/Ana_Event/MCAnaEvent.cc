#include "MCAnaEvent.hh"
//#include "THypernucleus.hh"
//#include "THyphiSimpleHit.hh"
//#include "TMcHit.hh"
//#include "TTrackSimple_v2.hh"
//#include "THyphiTrack_v3.hh"

using namespace std;

//#define DEBUG_DAISUKE
ClassImp(MCAnaEvent)

TClonesArray *MCAnaEvent::gMC_Particle = 0;
//TClonesArray *MCAnaEvent::gTR0 = 0;
TClonesArray *MCAnaEvent::gfTrack = 0;
//TClonesArray *MCAnaEvent::gfHyp = 0;

TClonesArray *MCAnaEvent::gTR0 = 0;
TClonesArray *MCAnaEvent::gTR1 = 0;
TClonesArray *MCAnaEvent::gTR2 = 0;

TClonesArray *MCAnaEvent::gDC1 = 0;

TClonesArray *MCAnaEvent::gDC2 = 0;
TClonesArray *MCAnaEvent::gDC3 = 0;

TClonesArray *MCAnaEvent::gDC2stop = 0;

TClonesArray *MCAnaEvent::gTOFp = 0;
TClonesArray *MCAnaEvent::gSTOP = 0;
TClonesArray *MCAnaEvent::gSTOP2 = 0;

MCAnaEvent::MCAnaEvent():TObject()
{
   if(!gMC_Particle) gMC_Particle = new TClonesArray("TMcParticle",10);

   if(!gTR0) gTR0 = new TClonesArray("TMcHit",20);
   if(!gTR1) gTR1 = new TClonesArray("TMcHit",20);
   if(!gTR2) gTR2 = new TClonesArray("TMcHit",20);

   if(!gDC1) gDC1 = new TClonesArray("TMcHit",20);

   if(!gDC2) gDC2 = new TClonesArray("TMcHit",20);
   if(!gDC3) gDC3 = new TClonesArray("TMcHit",20);
   if(!gDC2stop) gDC2stop = new TClonesArray("TMcHit",20);

   if(!gTOFp) gTOFp = new TClonesArray("TMcHit",20);
   if(!gSTOP) gSTOP = new TClonesArray("TMcHit",20);
   if(!gSTOP2) gSTOP2 = new TClonesArray("TMcHit",20);

   if(!gfTrack) gfTrack = new TClonesArray("THyphiTrack",20);
   //if(!gfHyp) gfHyp = new TClonesArray("THypernucleus",20);


   TR0 = gTR0;
   TR1 = gTR1;
   TR2 = gTR2;
   DC1 = gDC1;
   DC2 = gDC2;
   DC3 = gDC3;
   DC2stop = gDC2stop;

   TOFp = gTOFp;
   STOP = gSTOP;
   STOP2 = gSTOP2;

   fMC_Particle = gMC_Particle;
    fTrack = gfTrack;
   // fHyp = gfHyp;

   Nmc=0;
 
   Nhit_tr0=0;
   Nhit_tr1=0;
   Nhit_tr2=0;
   Nhit_dc1=0;
   
   Nhit_dc2=0;
   Nhit_dc3=0;
   Nhit_dc2stop=0;
   
   Nhit_tofp=0;
   Nhit_stop=0;
   Nhit_stop2=0;

   Nsystematic=0;
   SysX_shift = 0;
   SysZ_shift = 0;
   SysY_angle = 0;

   Field = 0;
   SysX_shift2 = 0;
   SysY_shift2 = 0;
   SysZ_angle2 = 0;
   Field2 = 0;

   Ntracks=0;
   // Nhyp=0;
   trigger = 0;
   Setup();
}

MCAnaEvent::~MCAnaEvent()
{
  Clear();
  Reset();
}

void MCAnaEvent::Clear(Option_t *option)
{
#ifdef DEBUG_DAISUKE
  std::cout<<"MCAnaEvent Clear start "<<std::endl;
#endif
   fMC_Particle->Clear("C");

#ifdef DEBUG_DAISUKE
   std::cout<<"MCAnaEvent Clear line "<<__LINE__<<std::endl;
#endif

   TR0->Clear("C");
   TR1->Clear("C");
   TR2->Clear("C");
#ifdef DEBUG_DAISUKE
   std::cout<<"MCAnaEvent Clear line "<<__LINE__<<std::endl;
#endif

   DC1->Clear("C");
#ifdef DEBUG_DAISUKE
   std::cout<<"MCAnaEvent Clear line "<<__LINE__<<std::endl;
#endif

   DC2->Clear("C");
   DC3->Clear("C");
   DC2stop->Clear("C");

   TOFp->Clear("C");
   STOP->Clear("C");
   STOP2->Clear("C");

    fTrack->Clear("C");
   // fHyp->Clear("C");
#ifdef DEBUG_DAISUKE
   std::cout<<"MCAnaEvent Clear line "<<__LINE__<<std::endl;
#endif

   Setup();
#ifdef DEBUG_DAISUKE
  std::cout<<"MCAnaEvent Clear finished "<<std::endl;
#endif
}

void MCAnaEvent::Reset()
{
#ifdef DEBUG_DAISUKE
  std::cout<<"MCAnaEvent Reset start "<<gfTrack<<" "<<std::endl;
#endif
  delete gMC_Particle; gMC_Particle = 0;
  delete gfTrack; gfTrack = 0;
  //  delete gTR0; gTR0 = 0;
  delete gTR0; gTR0 = 0;
  delete gTR1; gTR1 = 0;
  delete gTR2; gTR2 = 0;
  delete gDC1; gDC1 = 0;
  delete gDC2; gDC2 = 0;
  delete gDC3; gDC3 = 0;
  delete gDC2stop; gDC2stop = 0;

  delete gTOFp; gTOFp = 0;
  delete gSTOP; gSTOP = 0;
  delete gSTOP2; gSTOP2 = 0;

  //delete gfHyp; gfHyp = 0;
#ifdef DEBUG_DAISUKE
  std::cout<<"MCAnaEvent Reset finished"<<std::endl;
#endif

}


int MCAnaEvent::Setup()
{
#ifdef DEBUG_DAISUKE
  std::cout<<"MCAnaEvent Setup start"<<std::endl;
#endif
   Nmc=0;
 
   Nhit_tr0=0;
   Nhit_tr1=0;
   Nhit_tr2=0;
   Nhit_dc1=0;
   
   Nhit_dc2=0;
   Nhit_dc3=0;
   Nhit_dc2stop=0;
   
   Nhit_tofp=0;
   Nhit_stop=0;
   Nhit_stop2=0;

   Nsystematic=0;
   SysX_shift = 0;
   SysZ_shift = 0;
   SysY_angle = 0;

   Field = 0;
   SysX_shift2 = 0;
   SysY_shift2 = 0;
   SysZ_angle2 = 0;
   Field2 = 0;
   
   Ntracks=0;
   // Nhyp=0;
   trigger = 0;

#ifdef DEBUG_DAISUKE
  std::cout<<"MCAnaEvent Setup finished"<<std::endl;
#endif

  return 0;
}

int MCAnaEvent::Add_MC(const TMcParticle& M)
{
  //std::cout<<" Add_MC "<<M.type<<std::endl;

  TClonesArray &MC_Particle_Ref = *fMC_Particle;
  
  new(MC_Particle_Ref[Nmc]) TMcParticle(M);
  Nmc++;

  int status=0;
   
  return status;
   
}

int MCAnaEvent::Add_Hit(const TMcHit& H,TString detector)
{

  std::string funcname ="Add_Hit(const TMcHit& H,TString detector) :  ";
  std::string name_detector(detector.Data());
  funcname += name_detector;
  if(name_detector=="TR0")
    {
      TClonesArray &HitTr0_Ref = *TR0;
      
      new(HitTr0_Ref[Nhit_tr0]) TMcHit(H);
      Nhit_tr0++;
    }
  else if(name_detector=="TR1")
    {
      TClonesArray &HitTr1_Ref = *TR1;
      
      new(HitTr1_Ref[Nhit_tr1]) TMcHit(H);
      Nhit_tr1++;
    }
  else if(name_detector=="TR2")
    {
      TClonesArray &HitTr2_Ref = *TR2;
      
      new(HitTr2_Ref[Nhit_tr2]) TMcHit(H);
      Nhit_tr2++;
    }
  else if(name_detector=="DC1")
    {
      TClonesArray &HitDC1_Ref = *DC1;
      
      new(HitDC1_Ref[Nhit_dc1]) TMcHit(H);
      Nhit_dc1++;
    }
  else if(name_detector=="DC2")
    {
      TClonesArray &HitDC2_Ref = *DC2;
      
      new(HitDC2_Ref[Nhit_dc2]) TMcHit(H);
      Nhit_dc2++;
    }
  else if(name_detector=="DC3")
    {
      TClonesArray &HitDC3_Ref = *DC3;
      
      new(HitDC3_Ref[Nhit_dc3]) TMcHit(H);
      Nhit_dc3++;
    }
  else if(name_detector=="DC2stop")
    {
      TClonesArray &HitDC2stop_Ref = *DC2stop;
      
      new(HitDC2stop_Ref[Nhit_dc2stop]) TMcHit(H);
      Nhit_dc2stop++;
    }
  else if(name_detector=="TOFp")
    {
      TClonesArray &HitTOFp_Ref = *TOFp;
      
      new(HitTOFp_Ref[Nhit_tofp]) TMcHit(H);
      Nhit_tofp++;
    }
  else if(name_detector=="STOP")
    {
      TClonesArray &HitSTOP_Ref = *STOP;
      
      new(HitSTOP_Ref[Nhit_stop]) TMcHit(H);
      Nhit_stop++;
    }
  else if(name_detector=="STOP2")
    {
      TClonesArray &HitSTOP2_Ref = *STOP2;
      
      new(HitSTOP2_Ref[Nhit_stop2]) TMcHit(H);
      Nhit_stop2++;
    }
  else
    {
      std::cout<<funcname<<" wrond detector"<<std::endl;
      
    }

   int status=0;
   
   return status;
}


//int MCAnaEvent::Add_Track(const TTrackSimple& T)
int MCAnaEvent::Add_Track(const THyphiTrack& T)
{
   TClonesArray &Track_Ref = *fTrack;

   //new(Track_Ref[Ntracks]) TTrackSimple(T);
   new(Track_Ref[Ntracks]) THyphiTrack(T);
   Ntracks++;

   int status=0;

   return status;
}

// int MCAnaEvent::Add_Hyp(const THypernucleus& H)
// {
//    TClonesArray &Hyp_Ref = *fHyp;
   
   
//    new(Hyp_Ref[Nhyp]) THypernucleus(H);
//    Nhyp++;
   
//    int status=0;
   
//    return status;
// }

