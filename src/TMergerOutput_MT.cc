#include "TMergerOutput_MT.h"

#include "Debug.hh"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD

using namespace std;

TMergerOutput_MT::TMergerOutput_MT(const THyphiAttributes& attribut) : TDataMerger("merger_out"), att(attribut)
{
  att._logger->info("TMergerOutput_MT::TMergerOutput_MT");
}

TMergerOutput_MT::~TMergerOutput_MT() {}

ReturnRes::InfoM TMergerOutput_MT::operator()(FullRecoEvent& RecoEvent, ReturnRes::InfoM Status, MCAnaEventG4Sol* OutTree)
{
  int result = Exec(RecoEvent, Status, OutTree);

  return SoftExit(result);
}

void TMergerOutput_MT::SelectHists()
{
  LocalHisto.h_stats = AnaHisto->CloneAndRegister(AnaHisto->h_stats);
}

ReturnRes::InfoM TMergerOutput_MT::SoftExit(int ) { return ReturnRes::Fine; }

int TMergerOutput_MT::Exec(FullRecoEvent& RecoEvent, ReturnRes::InfoM Status, MCAnaEventG4Sol* OutTree)
{
  OutTree->Field = att.Field_Strength;

  switch(Status)
    {
    case ReturnRes::Fine:
      PushToTreeTracks(RecoEvent, OutTree);
    default:
      PushToTreeParticles(RecoEvent, OutTree);
      PushToTreeHits(RecoEvent, OutTree);
    }

  return 0;
}

void TMergerOutput_MT::PushToTreeParticles(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{
  for(auto OParticle : RecoEvent.ToDumpParticles)
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

void TMergerOutput_MT::PushToTreeHits(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{

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
        fillOutHit(OutTree->InSi, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet >= G4Sol::TR1 && TypeDet <= G4Sol::TR2)
        fillOutHit(OutTree->TR, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet >= G4Sol::CDC_layer0 && TypeDet <= G4Sol::CDC_layer14)
        fillOutHit(OutTree->CDC, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet >= G4Sol::MG01 && TypeDet <= G4Sol::MG17)
        fillOutHit(OutTree->CDC, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet == G4Sol::CDHBar)
        fillOutHit(OutTree->CDH, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet == G4Sol::RPC_l || TypeDet == G4Sol::RPC_h)
        fillOutHit(OutTree->RPC, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet == G4Sol::PSFE)
        fillOutHit(OutTree->PSFE, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet == G4Sol::PSCE)
        fillOutHit(OutTree->PSCE, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet == G4Sol::PSBE)
        fillOutHit(OutTree->PSBE, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet >= G4Sol::TrFwd0 && TypeDet <= G4Sol::TrFwd2)
        fillOutHit(OutTree->FwdTracker, RecoEvent.ToDumpHits[TypeDet]);

      if(TypeDet >= G4Sol::FMF2Stop0 && TypeDet <= G4Sol::FMF2Stop2)
        fillOutHit(OutTree->FMF2 , RecoEvent.ToDumpHits[TypeDet]);
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

void TMergerOutput_MT::PushToTreeTracks(FullRecoEvent& RecoEvent, MCAnaEventG4Sol* OutTree)
{

  for(auto TrackRes : RecoEvent.DAF_results)
    {
      int TrackID             = TrackRes.first;
      const ResSolDAF& FitRes = TrackRes.second;
      auto TInfo              = RecoEvent.TrackInfo.find(TrackID);
      auto IsDecay            = RecoEvent.TrackMother.find(TrackID);
      int Decay               = IsDecay == RecoEvent.TrackMother.end() ? 0 : 1;

      THyphiTrack* OutTrack = dynamic_cast<THyphiTrack*>(OutTree->fTrack->ConstructedAt(OutTree->fTrack->GetEntries()));
      auto PDG_particle     = TDatabasePDG::Instance()->GetParticle(FitRes.pdg_guess);
      OutTrack->type        = PDG_particle->GetName();
      OutTrack->MC_status   = TrackID + 10000 * Decay;
      OutTrack->Chi2        = FitRes.chi2;
      OutTrack->Chi2_X      = FitRes.ndf;
      OutTrack->Chi2_Y      = FitRes.Ncentral; // FitRes.firstHit;
      OutTrack->Mass        = FitRes.mass;
      OutTrack->pdgcode     = FitRes.pdg_guess;
      OutTrack->MomMass.SetXYZM(FitRes.momX, FitRes.momY, FitRes.momZ, FitRes.mass);
      OutTrack->Mom.SetXYZ(FitRes.momX, FitRes.momY, FitRes.momZ);

      OutTrack->BarId  = FitRes.lastHit;
      OutTrack->Charge = FitRes.charge;
      OutTrack->dE     = TInfo->second[FitRes.lastHit].Eloss;
      OutTrack->Beta   = FitRes.beta;
      OutTrack->RefPoint.SetXYZ(FitRes.posX, FitRes.posY, FitRes.posZ);
      OutTrack->Pval2      = FitRes.pvalue;
      OutTrack->PathLength = FitRes.path_length;
      OutTrack->TOF        = FitRes.tof;

      OutTrack->MomIni.SetXYZ(FitRes.momX_init, FitRes.momY_init, FitRes.momZ_init);
      OutTrack->BetaIni       = FitRes.beta2;
      OutTrack->MassIni       = FitRes.mass2;
      OutTrack->TOFIni        = FitRes.tof2;
      OutTrack->PathLengthIni = TInfo->second[FitRes.lastHit].tracklength; // FitRes.path_length2;
      OutTrack->RChiIni       = FitRes.fitter;
      if(Decay == 1)
        OutTrack->Sim2Vtx.SetXYZT(std::get<1>(IsDecay->second), std::get<2>(IsDecay->second),
                                  std::get<3>(IsDecay->second), std::get<4>(IsDecay->second));

      OutTrack->State[0] = FitRes.posX;
      OutTrack->State[1] = FitRes.posY;
      OutTrack->State[2] = FitRes.posZ;
      OutTrack->State[3] = FitRes.momX;
      OutTrack->State[4] = FitRes.momY;
      OutTrack->State[5] = FitRes.momZ;
      for(int row = 0; row < 6; row++)
        for(int col = 0; col < 6; col++)
          OutTrack->Cov[row][col] = FitRes.cov_matrix[row][col];
    }

  OutTree->Ntrack = OutTree->fTrack->GetEntries();
}
