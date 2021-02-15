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

  TrackDAF.clear();
  TrackDAFSim.clear();
  TrackDAFInit.clear();
  TrackInfo.clear();
  TrackMother.clear();
  
  Si_HitsEnergyLayer.clear();
  InteractionPoint = {0.,0.,0.};
}
