#include "Ana_Event.hh"
#include "THypernucleus.hh"

using namespace std;

ClassImp(Ana_Event)

TClonesArray *Ana_Event::gMC_Mom = 0;
TClonesArray *Ana_Event::gMC_Vtx = 0;
//TClonesArray *Ana_Event::gfSeedMom = 0;
TClonesArray *Ana_Event::gfTrackMom = 0;
//TClonesArray *Ana_Event::gfTrackPos = 0;
TClonesArray *Ana_Event::gfHyp = 0;
//TClonesArray *Ana_Event::gfHypMom = 0;
//TClonesArray *Ana_Event::gfHypVtx = 0;
//TClonesArray *Ana_Event::gfHypD1Mom = 0;
//TClonesArray *Ana_Event::gfHypD2Mom = 0;
//TClonesArray *Ana_Event::gfHypD3Mom = 0;


// THypernucleus::THypernucleus():MomMass(0.,0.,0.,0.),Vtx(0.,0.,0.,0.),MomMassD1(0.,0.,0.,0.),MomMassD2(0.,0.,0.,0.),MomMassD3(0.,0.,0.,0.)
// {
// }

// THypernucleus::THypernucleus(const THypernucleus& H)
// {
//    type=H.type;
//    Ndecay=H.Ndecay;
//    //   Chi2=H.Chi2;
//    Pvalue=H.Pvalue;
//    InvMass=H.InvMass;
//    Dist=H.Dist;
//    MomMass=H.MomMass;
//    Vtx=H.Vtx;
//    MomMassD1=H.MomMassD1;
//    MomMassD2=H.MomMassD2;
//    MomMassD3=H.MomMassD3;
   
// }


// //THypernucleus& THypernucleus::operator=(const THypernucleus& H)
// //{
// //   THypernucleus temp(H);
// //   return temp;
// //}


// THypernucleus::~THypernucleus()
// {
// }

// void THypernucleus::Clear(Option_t *option)
// {

//    MomMass.SetXYZT(0.,0.,0.,0.);
//    Vtx.SetXYZT(0.,0.,0.,0.);
//    MomMassD1.SetXYZT(0.,0.,0.,0.);
//    MomMassD2.SetXYZT(0.,0.,0.,0.);
//    MomMassD3.SetXYZT(0.,0.,0.,0.);

// }


Ana_Event::Ana_Event()
{
   if(!gMC_Mom) gMC_Mom = new TClonesArray("TVector3",30);
   if(!gMC_Vtx) gMC_Vtx = new TClonesArray("TVector3",30);
   //if(!gfSeedMom) gfSeedMom = new TClonesArray("TVector3",30);
   if(!gfTrackMom) gfTrackMom = new TClonesArray("TVector3",30);
   //if(!gfTrackPos) gfTrackPos = new TClonesArray("TVector3",30);
   if(!gfHyp) gfHyp = new TClonesArray("THypernucleus",10);
   //if(!gfHypMom) gfHypMom = new TClonesArray("TVector3",30);
   //if(!gfHypVtx) gfHypVtx = new TClonesArray("TVector3",30);
   //if(!gfHypD1Mom) gfHypD1Mom = new TClonesArray("TVector3",30);
   //if(!gfHypD2Mom) gfHypD2Mom = new TClonesArray("TVector3",30);
   //if(!gfHypD3Mom) gfHypD3Mom = new TClonesArray("TVector3",30);
   
   fMC_Mom = gMC_Mom;
   fMC_Vtx = gMC_Vtx;
   //fSeedMom = gfSeedMom;
   fTrackMom = gfTrackMom;
   //fTrackPos = gfTrackPos;
   fHyp = gfHyp;
   /*fHypMom = gfHypMom;
   fHypVtx = gfHypVtx;
   fHypD1Mom = gfHypD1Mom;
   fHypD2Mom = gfHypD2Mom;
   fHypD3Mom = gfHypD3Mom;
*/


   Nmc=0;
   //Nseed=0;
   NId=0;
   Ntracks=0;
   Nhyp=0;

   Setup();
}

Ana_Event::~Ana_Event()
{

   Clear();
   /*   MC_Mom=0;
   MC_Vtx=0;

   fSeedMom =0;
   fTrackMom =0;
   fTrackPos =0;

   fHypMom =0;
   fHypVtx =0;
   */
}

void Ana_Event::Clear(Option_t *option)
{
   fMC_Mom->Clear("C");
   fMC_Vtx->Clear("C");
   //fSeedMom->Clear("C");
   fTrackMom->Clear("C");
   //fTrackPos->Clear("C");
   fHyp->Clear("C");
   //fHypMom->Clear("C");
   //fHypVtx->Clear("C");
   //fHypD1Mom->Clear("C");
   //fHypD2Mom->Clear("C");
   //fHypD3Mom->Clear("C");
   Setup();
}

void Ana_Event::Reset()
{

delete gMC_Mom; gMC_Mom = 0;
delete gMC_Vtx; gMC_Vtx = 0;
//delete gfSeedMom; gfSeedMom = 0;
//delete gfTrackMom; gfTrackMom = 0;
//delete gfTrackPos; gfTrackPos = 0;
delete gfHyp; gfHyp = 0;
//delete gfHypMom; gfHypMom = 0;
//delete gfHypVtx; gfHypVtx = 0;
//delete gfHypD1Mom; gfHypD1Mom = 0;
//delete gfHypD2Mom; gfHypD2Mom = 0;
//delete gfHypD3Mom; gfHypD3Mom = 0;

}


int Ana_Event::Setup()
{
  fMC_name.clear();
  fId.clear();
  fIdName.clear();

  fTR0id.clear();
  //fTR0_X.clear();
  //fTR0_Y.clear();
  //fTR0_Z.clear();
  //fTR0_t.clear();
  fTR0_E.clear();

  //fTR1id.clear();
  //fTR1_X.clear();
  //fTR1_Y.clear();
  //fTR1_Z.clear();
  //fTR1_t.clear();
  //fTR1_E.clear();

  //fTR2id.clear();
  //fTR2_X.clear();
  //fTR2_Y.clear();
  //fTR2_Z.clear();
  //fTR2_t.clear();
  //fTR2_E.clear();

  //fTOFid.clear();
  //fTOF_X.clear();
  //fTOF_Y.clear();
  //fTOF_Z.clear();
  //fTOF_t.clear();
  //fTOF_E.clear();

  //fDC1id.clear();
  //fDC1_X.clear();
  //fDC1_Y.clear();
  //fDC1_Z.clear();
  //fDC1_t.clear();
  //fDC1_E.clear();

  //fDC2id.clear();
  //fDC2_X.clear();
  //fDC2_Y.clear();
  //fDC2_Z.clear();
  //fDC2_t.clear();
  //fDC2_E.clear();

  //fId_tof.clear();
  //fId_bfm.clear();
  //fId_dc.clear();

  //fTrackCharge.clear();
  fTrackMass.clear();
  fTrackChi2.clear();
  //fTrackPV.clear();

  //fHypName.clear();
  //fHypD1Name.clear();
  //fHypD2Name.clear();
  //fHypD3Name.clear();
  //fHypPV.clear();
  //fHypDist.clear();
  //fHypInvMass.clear();
  
  return 0;
}

int Ana_Event::Add_MC(TVector3 Mom,TVector3 Vtx)
{
   TClonesArray &MC_mom_Ref = *fMC_Mom;
   TClonesArray &MC_vtx_Ref = *fMC_Vtx;
   
   TVector3* testMom = new(MC_mom_Ref[Nmc]) TVector3(Mom);
   TVector3* testVtx = new(MC_vtx_Ref[Nmc]) TVector3(Vtx);
   Nmc++;

   int status=0;
   if(testMom==0 || testVtx==0)
   status++;
   if(testMom==0 && testVtx==0)
   status++;
   
   return status;
   
}
/*
int Ana_Event::Add_Seed(TVector3 Mom)
{
  TClonesArray &Seed_Ref = *fSeedMom;
   
   TVector3* testMom = new(Seed_Ref[Nseed++]) TVector3(Mom);   

   int status=0;
   if(testMom==0)
   status++;
   
   return status;
  
}
*/
int Ana_Event::Add_Track(TVector3 Mom,TVector3 Pos)
{
   TClonesArray &Track_Mom_Ref = *fTrackMom;
   //TClonesArray &Track_Pos_Ref = *fTrackPos;
   
   TVector3 *testMom = new(Track_Mom_Ref[Ntracks]) TVector3(Mom);
   //TVector3 *testPos = new(Track_Pos_Ref[Ntracks]) TVector3(Pos);
   Ntracks++;
   
   int status=0;
   //if(testMom==0 || testPos==0)
   status++;
   if(testMom==0) //&& testPos==0)
   status++;
   
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

/*int Ana_Event::Add_HypAndD(TVector3 Mom,TVector3 Vtx,TVector3 D1mom, TVector3 D2mom)
{
   TClonesArray &Hyp_mom_Ref = *fHypMom;
   TClonesArray &Hyp_vtx_Ref = *fHypVtx;
   TClonesArray &Hyp_D1mom_Ref = *fHypD1Mom;
   TClonesArray &Hyp_D2mom_Ref = *fHypD2Mom;

   
   TVector3* testMom = new(Hyp_mom_Ref[Nhyp]) TVector3(Mom);
   TVector3* testPos = new(Hyp_vtx_Ref[Nhyp]) TVector3(Vtx);
   TVector3* testD1Mom = new(Hyp_D1mom_Ref[Nhyp]) TVector3(D1mom);
   TVector3* testD2Mom = new(Hyp_D2mom_Ref[Nhyp]) TVector3(D2mom);
   Nhyp++;
   
   int status=0;
   if(testMom==0 || testPos==0 || testD1Mom==0 || testD2Mom==0)
   status++;
   if(testMom==0 && testPos==0 && testD1Mom==0 && testD2Mom==0)
   status++;
   
   return status;
}

int Ana_Event::Add_HypAndD(TVector3 Mom,TVector3 Vtx,TVector3 D1mom,TVector3 D2mom,TVector3 D3mom)
{
   TClonesArray &Hyp_mom_Ref = *fHypMom;
   TClonesArray &Hyp_vtx_Ref = *fHypVtx;
   TClonesArray &Hyp_D1mom_Ref = *fHypD1Mom;
   TClonesArray &Hyp_D2mom_Ref = *fHypD2Mom;
   TClonesArray &Hyp_D3mom_Ref = *fHypD3Mom;

   
   TVector3* testMom = new(Hyp_mom_Ref[Nhyp]) TVector3(Mom);
   TVector3* testPos = new(Hyp_vtx_Ref[Nhyp]) TVector3(Vtx);
   TVector3* testD1Mom = new(Hyp_D1mom_Ref[Nhyp]) TVector3(D1mom);
   TVector3* testD2Mom = new(Hyp_D2mom_Ref[Nhyp]) TVector3(D2mom);
   TVector3* testD3Mom = new(Hyp_D3mom_Ref[Nhyp]) TVector3(D3mom);
   Nhyp++;
   
   int status=0;
   if(testMom==0 || testPos==0 || testD1Mom==0 || testD2Mom==0 || testD3Mom==0)
   status++;
   if(testMom==0 && testPos==0 && testD1Mom==0 && testD2Mom==0 && testD3Mom==0)
   status++;
   
   return status;
}
*/
