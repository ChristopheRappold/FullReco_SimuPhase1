#include "TMergerOutput_ZMQ.h"

#include "Debug.hh"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD

using namespace std;
using namespace ZMQ;

TMergerOutput_ZMQ::TMergerOutput_ZMQ(const THyphiAttributes& attribut) : TDataMerger("merger_out"), att(attribut)
{
  att._logger->info("TMergerOutput_ZMQ::TMergerOutput_ZMQ");
}

TMergerOutput_ZMQ::~TMergerOutput_ZMQ() {}

ReturnRes::InfoM TMergerOutput_ZMQ::operator()(DataFitterOut& In, MCAnaEventG4Sol* OutTree)
{
  int result = Exec(In, OutTree);

  return SoftExit(result);
}

void TMergerOutput_ZMQ::SelectHists()
{
  LocalHisto.h_stats = AnaHisto->CloneAndRegister(AnaHisto->h_stats);
}

ReturnRes::InfoM TMergerOutput_ZMQ::SoftExit(int return_build) { return ReturnRes::Fine; }

int TMergerOutput_ZMQ::Exec(DataFitterOut& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  OutTree->Field = att.Field_Strength;

  switch(RecoEvent.status)
    {
    case ReturnRes::Fine:
      PushToTreeTracks(RecoEvent, OutTree);
    default:
      PushToTreeParticles(RecoEvent, OutTree);
      PushToTreeHits(RecoEvent, OutTree);
    }

  return 0;
}

void TMergerOutput_ZMQ::PushToTreeParticles(DataFitterOut& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  att._logger->info("Merger> PushToTreeParticles : size {}", RecoEvent.previousStep.DumpParticles.size());

  for(auto& OParticle : RecoEvent.previousStep.DumpParticles)
    {
      TMcParticle* OutParticle =
          dynamic_cast<TMcParticle*>(OutTree->fMC_Particle->ConstructedAt(OutTree->fMC_Particle->GetEntries()));

      OutParticle->type      = OParticle.type;
      OutParticle->Mc_id     = OParticle.Mc_id;
      OutParticle->Mother_id = OParticle.Mother_id;
      OutParticle->Pdg       = OParticle.Pdg;
      OutParticle->Charge    = OParticle.Charge;
      OutParticle->MomMass.SetXYZM(OParticle.MomMass[0], OParticle.MomMass[1], OParticle.MomMass[2],
                                   OParticle.MomMass[3]);
      OutParticle->Vtx.SetXYZT(OParticle.Vtx[0], OParticle.Vtx[1], OParticle.Vtx[2], OParticle.Vtx[3]);
      OutParticle->Weigth = OParticle.Weigth;
      OutParticle->GeoAcc = OParticle.GeoAcc;
    }

  OutTree->Nmc = OutTree->fMC_Particle->GetEntries();
}

void TMergerOutput_ZMQ::PushToTreeHits(DataFitterOut& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  att._logger->info("Merger> PushToTreeHits : size {}", RecoEvent.previousStep.DumpHits.size());

  auto fillOutHit = [](TClonesArray* out, const std::vector<OutHit>& AllHits) {
    for(const auto& hit : AllHits)
      {
        TMcHit* OutHit  = dynamic_cast<TMcHit*>(out->ConstructedAt(out->GetEntries()));
        OutHit->name    = hit.name;
        OutHit->LayerID = hit.LayerID;
        OutHit->HitID   = hit.HitID;
        OutHit->MCHit.SetXYZ(hit.MCHit[0], hit.MCHit[1], hit.MCHit[2]);
        OutHit->Hit.SetXYZ(hit.Hit[0], hit.Hit[1], hit.Hit[2]);
        OutHit->MC_id  = hit.MC_id;
        OutHit->Pdg    = hit.Pdg;
        OutHit->Charge = hit.Charge;
        OutHit->MCparticle.SetXYZM(hit.MCparticle[0], hit.MCparticle[1], hit.MCparticle[2], hit.MCparticle[3]);
        OutHit->Brho              = hit.Brho;
        OutHit->MagnetInteraction = hit.MagnetInteraction;
        // std::cout<<" Out> LayerID:"<<LayerID<<" "<<HitID<<std::endl;
      }
  };

  for(int TypeDet = G4Sol::InSi0, id_det_last = G4Sol::FMF2Stop2; TypeDet <= id_det_last; ++TypeDet)
    {
      if(TypeDet >= G4Sol::InSi0 && TypeDet <= G4Sol::InSi3)
        fillOutHit(OutTree->InSi, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet >= G4Sol::TR1 && TypeDet <= G4Sol::TR2)
        fillOutHit(OutTree->TR, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet >= G4Sol::CDC_layer0 && TypeDet <= G4Sol::CDC_layer14)
        fillOutHit(OutTree->CDC, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet >= G4Sol::MG01 && TypeDet <= G4Sol::MG17)
        fillOutHit(OutTree->CDC, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet == G4Sol::CDHBar)
        fillOutHit(OutTree->CDH, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet == G4Sol::RPC_l || TypeDet == G4Sol::RPC_h)
        fillOutHit(OutTree->RPC, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet == G4Sol::PSFE)
        fillOutHit(OutTree->PSFE, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet == G4Sol::PSCE)
        fillOutHit(OutTree->PSCE, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet == G4Sol::PSBE)
        fillOutHit(OutTree->PSBE, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet >= G4Sol::TrFwd0 && TypeDet <= G4Sol::TrFwd2)
        fillOutHit(OutTree->FwdTracker, RecoEvent.previousStep.DumpHits[TypeDet]);

      if(TypeDet >= G4Sol::FMF2Stop0 && TypeDet <= G4Sol::FMF2Stop2)
        fillOutHit(OutTree->FMF2 , RecoEvent.previousStep.DumpHits[TypeDet]);
    }

  OutTree->NInSi       = OutTree->InSi->GetEntries();
  OutTree->NTr         = OutTree->TR->GetEntries();
  OutTree->NCdc        = OutTree->CDC->GetEntries();
  OutTree->NCdh        = OutTree->CDH->GetEntries();
  OutTree->NRpc        = OutTree->RPC->GetEntries();
  OutTree->NFwdtracker = OutTree->FwdTracker->GetEntries();
  OutTree->NFmf2       = OutTree->FMF2->GetEntries();
  OutTree->NPsbe       = OutTree->PSBE->GetEntries();
  OutTree->NPsfe       = OutTree->PSFE->GetEntries();
  OutTree->NPsce       = OutTree->PSCE->GetEntries();
}

void TMergerOutput_ZMQ::PushToTreeTracks(DataFitterOut& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  att._logger->info("Merger> PushToTreeTracks : size {}", RecoEvent.DumpTracks.size());

  for(auto& TrackRes : RecoEvent.DumpTracks)
    {

      THyphiTrack* OutTrack = dynamic_cast<THyphiTrack*>(OutTree->fTrack->ConstructedAt(OutTree->fTrack->GetEntries()));

      OutTrack->type        = TrackRes.type;
      OutTrack->MC_status   = TrackRes.MC_status;
      OutTrack->Chi2        = TrackRes.Chi2;
      OutTrack->Chi2_X      = TrackRes.Chi2_X;
      OutTrack->Chi2_Y      = TrackRes.Chi2_Y; // TrackRes.firstHit;
      OutTrack->Mass        = TrackRes.Mass;
      OutTrack->pdgcode     = TrackRes.pdgcode;
      OutTrack->MomMass.SetXYZM(TrackRes.MomMass[0],TrackRes.MomMass[1],TrackRes.MomMass[2],TrackRes.MomMass[3]);
      OutTrack->Mom.SetXYZ(TrackRes.Mom[0],TrackRes.Mom[1],TrackRes.Mom[2]);

      OutTrack->BarId  = TrackRes.BarId;
      OutTrack->Charge = TrackRes.Charge;
      OutTrack->dE     = TrackRes.dE;
      OutTrack->Beta   = TrackRes.Beta;
      OutTrack->RefPoint.SetXYZ(TrackRes.RefPoint[0],TrackRes.RefPoint[1],TrackRes.RefPoint[2]);
      OutTrack->Pval2      = TrackRes.Pval2;
      OutTrack->PathLength = TrackRes.PathLength;
      OutTrack->TOF        = TrackRes.TOF;

      OutTrack->MomIni.SetXYZ(TrackRes.MomIni[0], TrackRes.MomIni[1], TrackRes.MomIni[2]);
      OutTrack->BetaIni       = TrackRes.BetaIni;
      OutTrack->MassIni       = TrackRes.MassIni;
      OutTrack->TOFIni        = TrackRes.TOFIni;
      OutTrack->PathLengthIni = TrackRes.PathLengthIni;
      OutTrack->RChiIni       = TrackRes.RChiIni;

      OutTrack->Sim2Vtx.SetXYZT(TrackRes.Sim2Vtx[0],TrackRes.Sim2Vtx[1],TrackRes.Sim2Vtx[2],TrackRes.Sim2Vtx[3]);

      for(int i=0;i<6;++i)
	OutTrack->State[i] = TrackRes.State[i];
      
      for(int row = 0; row < 6; row++)
        for(int col = 0; col < 6; col++)
          OutTrack->Cov[row][col] = TrackRes.Cov[row][col];
    }

  OutTree->Ntrack = OutTree->fTrack->GetEntries();
}
