#include "CheckField.h"

#include "Ana_Event/MCAnaEventG4Sol.hh"
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

template<class Out>
CheckField<Out>::CheckField(const THyphiAttributes& attribut):TDataProcessInterface<Out>("check_field"),att(attribut),done(0)
{ }

template<class Out>
void CheckField<Out>::InitMT() {att._logger->error("E> Not supposed to be multithreaded !"); }

template<class Out>
ReturnRes::InfoM CheckField<Out>::operator() (FullRecoEvent& RecoEvent,Out* OutTree)
{

  if(done==0)
    {
      Exec(RecoEvent,OutTree);
      done=1;
    }
  return SoftExit(0);
}

template<class Out>
int CheckField<Out>::Exec(FullRecoEvent& RecoEvent,Out* OutTree)
{
  return Check();
}

template<class Out>
ReturnRes::InfoM CheckField<Out>::SoftExit(int result_full)
{
  return ReturnRes::Fine;    
}

template<class Out>
void CheckField<Out>::SelectHists()
{
  for(size_t i=0;i<3;++i)
    {
      LocalHisto.FieldXY[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldXY[i]); 
      LocalHisto.FieldXZ[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldXZ[i]); 
      LocalHisto.FieldYZ[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldYZ[i]); 
      LocalHisto.FieldXYmax[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldXYmax[i]);
      LocalHisto.FieldXZmax[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldXZmax[i]);
      LocalHisto.FieldYZmax[i]   = AnaHisto->CloneAndRegister(AnaHisto->FieldYZmax[i]);
      LocalHisto.FieldXY_n[i] = AnaHisto->CloneAndRegister(AnaHisto->FieldXY_n[i]); 
      LocalHisto.FieldXZ_n[i] = AnaHisto->CloneAndRegister(AnaHisto->FieldXZ_n[i]); 
      LocalHisto.FieldYZ_n[i] = AnaHisto->CloneAndRegister(AnaHisto->FieldYZ_n[i]); 

      LocalHisto.FieldXY[i]->Sumw2()  ;
      LocalHisto.FieldXZ[i]->Sumw2()  ;
      LocalHisto.FieldYZ[i]->Sumw2()  ;
      LocalHisto.FieldXY_n[i]->Sumw2();
      LocalHisto.FieldXZ_n[i]->Sumw2();
      LocalHisto.FieldYZ_n[i]->Sumw2();
    }
}

template<class Out>
int CheckField<Out>::Check()
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
  auto GenfitMag = genfit::FieldManager::getInstance();

  double minX = att.Config.IsAvailable("CheckField_minX") ? att.Config.Get<double>("CheckField_minX") : -Dx;
  double maxX = att.Config.IsAvailable("CheckField_maxX") ? att.Config.Get<double>("CheckField_maxX") : Dx;
  double minY = att.Config.IsAvailable("CheckField_minY") ? att.Config.Get<double>("CheckField_minY") : -Dy;
  double maxY = att.Config.IsAvailable("CheckField_maxY") ? att.Config.Get<double>("CheckField_maxY") : Dy;
  double minZ = att.Config.IsAvailable("CheckField_minZ") ? att.Config.Get<double>("CheckField_minZ") : -Dz;
  double maxZ = att.Config.IsAvailable("CheckField_maxZ") ? att.Config.Get<double>("CheckField_maxZ") : Dz;

  att._logger->info(" between : [{}, {}] [{}, {}] [{}, {}]",minX, maxX, minY, maxY, minZ, maxZ);

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
  //double y = 0.;

  att._logger->debug("binXYZ {} {} {}",LocalHisto.FieldXY[0]->GetXaxis()->GetNbins(),LocalHisto.FieldXY[0]->GetYaxis()->GetNbins(),LocalHisto.FieldXZ[0]->GetYaxis()->GetNbins());

  for(int iX = 1; iX <= LocalHisto.FieldXY[0]->GetNbinsX();++iX)
    for(int iY = 1; iY <= LocalHisto.FieldXY[0]->GetNbinsY();++iY)
      for(int iZ = 1; iZ <= LocalHisto.FieldXZ[0]->GetNbinsY();++iZ)
  // for(double x = minX; x<maxX ; x += 1.)
  //   //for(double y = minY; y<maxY ; y += 1.)
  //     for(double z = minZ; z<maxZ ; z += 1.)
	{
	  double x = LocalHisto.FieldXY[0]->GetXaxis()->GetBinCenter(iX);
	  double y = LocalHisto.FieldXY[0]->GetYaxis()->GetBinCenter(iY);
	  double z = LocalHisto.FieldXZ[0]->GetYaxis()->GetBinCenter(iZ);
	  if(x<minX || x>maxX || y<minY || y>maxY || z<minZ || z>maxZ)
	    continue;

	  if(i%10000000==0)
	    att._logger->info("Processing Check# {}",i);
	  ++i;
	  double b[3]={0.,0.,0.};
	  
	  GenfitMag->getFieldVal(x,y,z,b[0],b[1],b[2]);
	  for(int k=0;k<3;++k)
	    if(TMath::Abs(b[k])>1e-5)
	      {
		// LocalHisto.FieldXY[k]->Fill(x,y,b[k]);
		// LocalHisto.FieldXZ[k]->Fill(x,z,b[k]);
		// LocalHisto.FieldYZ[k]->Fill(y,z,b[k]);

		// LocalHisto.FieldXY_n[k]->Fill(x,y);
		// LocalHisto.FieldXZ_n[k]->Fill(x,z);
		// LocalHisto.FieldYZ_n[k]->Fill(y,z);

		LocalHisto.FieldXY[k]->AddBinContent(LocalHisto.FieldXY[k]->GetBin(iX,iY),b[k]);
		LocalHisto.FieldXZ[k]->AddBinContent(LocalHisto.FieldXZ[k]->GetBin(iX,iZ),b[k]);
		LocalHisto.FieldYZ[k]->AddBinContent(LocalHisto.FieldYZ[k]->GetBin(iY,iZ),b[k]);

		LocalHisto.FieldXY_n[k]->AddBinContent(LocalHisto.FieldXY_n[k]->GetBin(iX,iY),1.);
		LocalHisto.FieldXZ_n[k]->AddBinContent(LocalHisto.FieldXZ_n[k]->GetBin(iX,iZ),1.);
		LocalHisto.FieldYZ_n[k]->AddBinContent(LocalHisto.FieldYZ_n[k]->GetBin(iY,iZ),1.);

		double tempB1 = LocalHisto.FieldXYmax[k]->GetBinContent(iX,iY);
		if(tempB1<=b[k])
		  LocalHisto.FieldXYmax[k]->SetBinContent(iX,iY,b[k]);

		double tempB2 = LocalHisto.FieldXZmax[k]->GetBinContent(iX,iZ);
		if(tempB2<=b[k])
		  LocalHisto.FieldXZmax[k]->SetBinContent(iX,iZ,b[k]);

		double tempB3 = LocalHisto.FieldYZmax[k]->GetBinContent(iY,iZ);
		if(tempB3<=b[k])
		  LocalHisto.FieldYZmax[k]->SetBinContent(iY,iZ,b[k]);


	      }
	}

  auto f_norm = [](TH2F* h1, TH2F* h_n)
  {
    for(int iX = 1; iX<h1->GetNbinsX();++iX)
      for(int iY = 1; iY<h1->GetNbinsY();++iY)
	{
	  double tempB = h1->GetBinContent(iX,iY);
	  double tempN = h_n->GetBinContent(iX,iY);

	  if(TMath::Abs(tempN)>0.1)
	    h1->SetBinContent(iX,iY,tempB/tempN);
	  else
	    h1->SetBinContent(iX,iY,0.);
	}
  };

  for(int k=0;k<3;++k)
    {
      f_norm(LocalHisto.FieldXY[k], LocalHisto.FieldXY_n[k]);
      f_norm(LocalHisto.FieldXZ[k], LocalHisto.FieldXZ_n[k]);
      f_norm(LocalHisto.FieldYZ[k], LocalHisto.FieldYZ_n[k]);
    }


  return 0;
}


template class CheckField<MCAnaEventG4Sol>;
