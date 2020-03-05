#include "CheckField.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"

#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include <iostream>
#include <tuple>

#include "FieldManager.h"

CheckField::CheckField(const THyphiAttributes& attribut):TDataProcessInterface("check_field"),att(attribut),done(0)
{ }

void CheckField::InitMT() {att._logger->error("E> Not supposed to be multithreaded !"); }

ReturnRes::InfoM CheckField::operator() (FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree)
{

  if(done==0)
    {
      Exec(RecoEvent,OutTree);
      done=1;
    }
  return SoftExit(0);
}

int CheckField::Exec(FullRecoEvent& RecoEvent,MCAnaEventG4Sol* OutTree)
{
  return Check();
}

ReturnRes::InfoM CheckField::SoftExit(int result_full)
{
  return ReturnRes::Fine;    
}

void CheckField::SelectHists()
{
  for(size_t i=0;i<3;++i)
    {
      LocalHisto.FieldXY[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldXY[i]); 
      LocalHisto.FieldXZ[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldXZ[i]); 
      LocalHisto.FieldYZ[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldYZ[i]); 
      LocalHisto.FieldXY_n[i] = AnaHisto->CloneAndRegister(AnaHisto->FieldXY_n[i]); 
      LocalHisto.FieldXZ_n[i] = AnaHisto->CloneAndRegister(AnaHisto->FieldXZ_n[i]); 
      LocalHisto.FieldYZ_n[i] = AnaHisto->CloneAndRegister(AnaHisto->FieldYZ_n[i]); 
    }
}

int CheckField::Check()
{
  att._logger->info("Start Field check :");

  TGeoVolume* v = gGeoManager->GetVolume(0);
  TGeoShape* s = v->GetShape();
  s->ComputeBBox();
  TGeoBBox* sb = dynamic_cast<TGeoBBox*>(s);
  if(sb==nullptr)
    {
      att._logger->error("E> CheckField: dynamic_cast failed !");
      v->Print();
      return 0;
    }

  //const Double_t* origin = sb->GetOrigin();
  double Dx = sb->GetDX();
  double Dy = sb->GetDY();
  double Dz = sb->GetDZ();
  att._logger->info(" between : {} {} {}",Dx, Dy, Dz);
  auto GenfitMag = genfit::FieldManager::getInstance();
  // const auto max_i = 10000000;
  // auto i = max_i;
  
  // while(i>0)
  //   {
  //     if(i%500000==0)
  // 	std::cout<<"Processing Check#"<<i<<"\n";

  //     double x = origin[0]+gRandom->Uniform(-Dx,Dx);
  //     double y = origin[1]+gRandom->Uniform(-Dy,Dy);
  //     double z = origin[2]+gRandom->Uniform(-Dz,Dz);

  size_t i = 0;
  double y=0.;
  for(double x = -Dx; x<Dx ; x += 1.) 
    //for(double y = -Dy; y<Dy ; y += 1.) 
      for(double z = -Dz; z<Dz ; z += 1.) 
	{
	  if(i%10000000==0)
	    att._logger->info("Processing Check# {}",i);
	  ++i;
	  double b[3]={0.,0.,0.};
	  
	  GenfitMag->getFieldVal(x,y,z,b[0],b[1],b[2]);
	  for(int k=0;k<3;++k)
	    if(TMath::Abs(b[k])>1e-10)
	      {
		LocalHisto.FieldXY[k]->Fill(x,y,b[k]);
		LocalHisto.FieldXZ[k]->Fill(x,z,b[k]);
		LocalHisto.FieldYZ[k]->Fill(y,z,b[k]);

		LocalHisto.FieldXY_n[k]->Fill(x,y);
		LocalHisto.FieldXZ_n[k]->Fill(x,z);
		LocalHisto.FieldYZ_n[k]->Fill(y,z);
	      }
	}
  
  return 0;
}
