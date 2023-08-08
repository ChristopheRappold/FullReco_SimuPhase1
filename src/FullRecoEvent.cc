#include "FullRecoEvent.hh"

#include <iostream>

const std::string PDG_fromName::ElName2[] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu",
    "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds"};

FullRecoEvent::FullRecoEvent(unsigned int idTh) : idThread(idTh)
{
  std::cout << " *** > FullRecoEvent object created" << std::endl;
  Clear();
}

FullRecoEvent::~FullRecoEvent()
{
  std::cout << " *** > FullRecoEvent object deleted" << std::endl;
  Clear();
}

void FullRecoEvent::Clear(int toclean)
{
  idEvent = -1;
  ToDumpParticles.clear();
  ToDumpHits.clear();
  // ToDumpTracks.clear();
  // Mom_Particle.clear();
  DAF_results.clear();
  // for(auto& det : ListHits)
  //   for(auto& hit : det)
  //     {
  // 	if(hit!=nullptr)
  // 	  delete hit;
  //     	hit = nullptr;
  //     }
  ListHits.clear();
  OldListHits.clear();
  ListHitsInfo.clear();
  
  SegmentHit1Ds.clear();

  ListHitsToTracks.clear();

  TrackDAF.clear();
  TrackDAFSim.clear();
  TrackDAFInit.clear();
  TrackDAFInitSim.clear();
  TrackInfo.clear();

  TracksFound.clear();
  IdHitsToMeasurement.clear();

  TrackMother.clear();
  DaughtersTrackDAFInit.clear();

  MWDCTracks.clear();

  BeamTracks.clear();
  PrimaryTracks.clear();

  FragmentTracks.clear();
  FragmentPID = -999;

  for(auto x : FiberTrackCont)
    for(size_t i=0; i<x.second.size(); ++i)
      delete x.second[i];

  for(size_t i = 0; i < FiberHitClCont.size(); ++i)
    {
      for(size_t j = 0; j < FiberHitClCont[i].size(); ++j)
        {
          size_t num = FiberHitClCont[i][j].size();
          for(size_t k = 0; k < num; ++k)
            delete FiberHitClCont[i][j][k];
        }
    }

  FiberTrackCont.clear();
  FiberHitClCont.clear();

  for(size_t i = 0; i < MDCHitCont.size(); ++i)
    for(size_t j = 0; j < MDCHitCont[i].size(); ++j)
      delete MDCHitCont[i][j];
  MDCHitCont.clear();

  for(size_t i = 0; i <PSBHitCont.size(); ++i)
    delete PSBHitCont[i];
  PSBHitCont.clear();

  for(size_t i = 0; i <PSFEHitCont.size(); ++i)
    delete PSFEHitCont[i];
  PSFEHitCont.clear();


  //Added when merging with master
  //FragmentTracks.clear();
  //PionTracks.clear();


  Si_HitsEnergyLayer.clear();
/*
  Hits_Si1.clear();
  Hits_Si2.clear();
  Hits_Si3.clear();
  Hits_Si4.clear();
  
  HitsX_Si1.clear();
  HitsY_Si1.clear();
  HitsX_Si2.clear();
  HitsY_Si2.clear();
*/

  Mother_MomE.SetXYZM(0.,0.,0.,0.);
  InteractionPoint.fill(0.);
  DecayVertex.fill(0.);
  Hyp_LifeTime = -1.;

  PrimVtxRecons.SetXYZ(0., 0., 0.);
  CovMatrix_IP.fill(0.);

  DecayVtxRecons.SetXYZ(0., 0., 0.);
  CovMatrix_SV.fill(0.);

  Hyp_Vect.clear();
}



double PrimaryVtxTrack::GetTheta()
{
  if(a > -990. && b > -990.)
    return TMath::ATan2(std::sqrt( a*a + b*b ), 1.) * 180. / M_PI;
  else
    return -999.;
}

double PrimaryVtxTrack::GetPhi()
{
  if(a > -990. && b > -990.)
    return TMath::ATan2(b, a) * 180. / M_PI;
  else
    return -999.;
}
