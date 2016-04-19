#include "FullRecoEvent.hh"
#include <iostream>


FullRecoEvent::FullRecoEvent()
{
  std::cout<<" *** > FullRecoEvent object created"<<std::endl; 
  Clear();
}

FullRecoEvent::~FullRecoEvent()
{
  std::cout<<" *** > FullRecoEvent object deleted"<<std::endl; 
  Clear();

}

void FullRecoEvent::Clear(int toclean)
{
  //Mom_Particle.clear();
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
  TrackInfo.clear();
}

