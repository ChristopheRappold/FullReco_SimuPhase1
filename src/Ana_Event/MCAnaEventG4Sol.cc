#include "MCAnaEventG4Sol.hh"

using namespace std;

ClassImp(MCAnaEventG4Sol)

TClonesArray *MCAnaEventG4Sol::gMC_Particle = 0;

TClonesArray *MCAnaEventG4Sol::gInSi = 0;
TClonesArray *MCAnaEventG4Sol::gTR = 0;
TClonesArray *MCAnaEventG4Sol::gFiber = 0;
TClonesArray *MCAnaEventG4Sol::gCDC = 0;
TClonesArray *MCAnaEventG4Sol::gCDH = 0;
TClonesArray *MCAnaEventG4Sol::gFwdTracker = 0;
TClonesArray *MCAnaEventG4Sol::gRPC = 0;
TClonesArray *MCAnaEventG4Sol::gFMF2 = 0;

TClonesArray *MCAnaEventG4Sol::gPSFE = 0;
TClonesArray *MCAnaEventG4Sol::gPSBE = 0;
TClonesArray *MCAnaEventG4Sol::gPSCE = 0;

TClonesArray *MCAnaEventG4Sol::gfTrack = 0;

MCAnaEventG4Sol::MCAnaEventG4Sol()//:TObject()
{
  if(!gMC_Particle) gMC_Particle = new TClonesArray("TMcParticle",10);

  if(!gInSi) gInSi = new TClonesArray("TMcHit",20);
  if(!gTR) gTR = new TClonesArray("TMcHit",20);
  if(!gFiber) gFiber = new TClonesArray("TMcHit",20);
  if(!gCDC) gCDC = new TClonesArray("TMcHit",20);
  if(!gCDH) gCDH = new TClonesArray("TMcHit",20);
  if(!gFwdTracker) gFwdTracker = new TClonesArray("TMcHit",20);
  if(!gRPC) gRPC = new TClonesArray("TMcHit",20);
  if(!gFMF2) gFMF2 = new TClonesArray("TMcHit",20);

  if(!gPSBE) gPSBE = new TClonesArray("TMcHit",20);
  if(!gPSFE) gPSFE = new TClonesArray("TMcHit",20);
  if(!gPSCE) gPSCE = new TClonesArray("TMcHit",20);

  if(!gfTrack) gfTrack = new TClonesArray("THyphiTrack",20);
  //if(!gfHyp) gfHyp = new TClonesArray("THypernucleus",20);


  InSi = gInSi;
  TR = gTR;
  Fiber = gFiber;
  CDC = gCDC;
  CDH = gCDH;
  FwdTracker = gFwdTracker;
  RPC = gRPC;
  FMF2 = gFMF2;

  PSFE = gPSFE; 
  PSBE = gPSBE; 
  PSCE = gPSCE;
  
  fMC_Particle = gMC_Particle;
  fTrack = gfTrack;
  // fHyp = gfHyp;

  Nmc=0;
 
  NInSi=0;
  NTr=0;
  NFiber=0;
  NCdc=0;
  NCdh=0;
  NFwdtracker=0;
  NRpc=0;
  NFmf2=0;
  NPsbe=0;
  NPsfe=0;
  NPsce=0;
  
  Field = 0.;
  ReducFactor = 1.;
  NtrackCand = 0;
  Ntrack = 0;
  trigger = 0;

  Setup();
}

MCAnaEventG4Sol::~MCAnaEventG4Sol()
{
  Clear();
  Reset();
}

void MCAnaEventG4Sol::Clear(Option_t *option)
{

  fMC_Particle->Clear("C");

  InSi->Clear("C");
  TR->Clear("C");
  Fiber->Clear("C");
  CDC->Clear("C");
  CDH->Clear("C");
  FwdTracker->Clear("C");
  RPC->Clear("C");
  FMF2->Clear("C");

  PSBE->Clear("C");
  PSFE->Clear("C");
  PSCE->Clear("C"); 
  
  TrackCand->Clear("C");

  fTrack->Clear("C");

  Setup();
}

void MCAnaEventG4Sol::Reset()
{

  delete gMC_Particle; gMC_Particle = 0;
  delete gfTrack; gfTrack = 0;

  delete gInSi; gInSi = 0;
  delete gTR; gTR = 0;
  delete gFiber; gFiber = 0;
  delete gCDC; gCDC = 0;
  delete gCDH; gCDH = 0;
  delete gFwdTracker; gFwdTracker = 0;
  delete gRPC; gRPC = 0;
  delete gFMF2; gFMF2 = 0;
  delete gPSFE; gPSFE = 0;
  delete gPSBE; gPSBE = 0;
  delete gPSCE; gPSCE = 0;

}


int MCAnaEventG4Sol::Setup()
{
  Nmc=0;
 
  NInSi=0;
  NTr=0;
  NFiber=0;
  NCdc=0;
  NCdh=0;
  NFwdtracker=0;
   
  NRpc=0;
  NFmf2=0;

  NPsbe=0;
  NPsfe=0;
  NPsce=0;

  Field = 0;
  ReducFactor =1.;
  NtrackCand =0;
  Ntrack =0;
  trigger = 0;


  return 0;
}

// int MCAnaEventG4Sol::Add_MC(const TMcParticle& M)
// {
//   //std::cout<<" Add_MC "<<M.type<<std::endl;

//   TClonesArray &MC_Particle_Ref = *fMC_Particle;
  
//   new(MC_Particle_Ref[Nmc]) TMcParticle(M);
//   Nmc++;

//   int status=0;
   
//   return status;
   
// }

// int MCAnaEventG4Sol::Add_Hit(const TMcHit& H,TString detector)
// {

//   std::string funcname ="Add_Hit(const TMcHit& H,TString detector) :  ";
//   std::string name_detector(detector.Data());
//   funcname += name_detector;
//   if(name_detector=="InSi")
//     {
//       TClonesArray &HitTr0_Ref = *InSi;
      
//       new(HitTr0_Ref[Nhit_tr0]) TMcHit(H);
//       Nhit_tr0++;
//     }
//   else if(name_detector=="CDC")
//     {
//       TClonesArray &HitTr1_Ref = *CDC;
      
//       new(HitTr1_Ref[Nhit_tr1]) TMcHit(H);
//       Nhit_tr1++;
//     }
//   else if(name_detector=="CDH")
//     {
//       TClonesArray &HitTr2_Ref = *CDH;
      
//       new(HitTr2_Ref[Nhit_tr2]) TMcHit(H);
//       Nhit_tr2++;
//     }
//   else if(name_detector=="FwdTracker")
//     {
//       TClonesArray &HitFwdTracker_Ref = *FwdTracker;
      
//       new(HitFwdTracker_Ref[Nhit_dc1]) TMcHit(H);
//       Nhit_dc1++;
//     }
//   else if(name_detector=="RPC")
//     {
//       TClonesArray &HitRPC_Ref = *RPC;
      
//       new(HitRPC_Ref[Nhit_dc2]) TMcHit(H);
//       Nhit_dc2++;
//     }
//   else if(name_detector=="FMF2")
//     {
//       TClonesArray &HitFMF2_Ref = *FMF2;
      
//       new(HitFMF2_Ref[Nhit_dc3]) TMcHit(H);
//       Nhit_dc3++;
//     }
//   else if(name_detector=="RPCstop")
//     {
//       TClonesArray &HitRPCstop_Ref = *RPCstop;
      
//       new(HitRPCstop_Ref[Nhit_dc2stop]) TMcHit(H);
//       Nhit_dc2stop++;
//     }
//   else if(name_detector=="TOFp")
//     {
//       TClonesArray &HitTOFp_Ref = *TOFp;
      
//       new(HitTOFp_Ref[Nhit_tofp]) TMcHit(H);
//       Nhit_tofp++;
//     }
//   else if(name_detector=="STOP")
//     {
//       TClonesArray &HitSTOP_Ref = *STOP;
      
//       new(HitSTOP_Ref[Nhit_stop]) TMcHit(H);
//       Nhit_stop++;
//     }
//   else if(name_detector=="STOP2")
//     {
//       TClonesArray &HitSTOP2_Ref = *STOP2;
      
//       new(HitSTOP2_Ref[Nhit_stop2]) TMcHit(H);
//       Nhit_stop2++;
//     }
//   else
//     {
//       std::cout<<funcname<<" wrond detector"<<std::endl;
      
//     }

//   int status=0;
   
//   return status;
// }


// //int MCAnaEventG4Sol::Add_Track(const TTrackSimple& T)
// int MCAnaEventG4Sol::Add_Track(const THyphiTrack& T)
// {
//   TClonesArray &Track_Ref = *fTrack;

//   //new(Track_Ref[Ntrack ]) TTrackSimple(T);
//   new(Track_Ref[Ntrack ]) THyphiTrack(T);
//   Ntrack ++;

//   int status=0;

//   return status;
// }

// int MCAnaEventG4Sol::Add_Hyp(const THypernucleus& H)
// {
//    TClonesArray &Hyp_Ref = *fHyp;
   
   
//    new(Hyp_Ref[Nhyp]) THypernucleus(H);
//    Nhyp++;
   
//    int status=0;
   
//    return status;
// }

