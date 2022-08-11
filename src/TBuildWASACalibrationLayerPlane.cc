#include "TBuildWASACalibrationLayerPlane.h"

#include "Debug.hh"
#include "TGeoManager.h"

#include <list>
#include <map>
#include <set>
#include <vector>

//#define DEBUG_BUILD
//#define DEBUG_BUILD3

using namespace std;

TBuildWASACalibrationLayerPlane::TBuildWASACalibrationLayerPlane(const THyphiAttributes& attribut) : TDataBuilder("build_det"), att(attribut)
{
  att._logger->info("TBuildWASACalibrationLayerPlane::TBuildWASACalibrationLayerPlane");


}

TBuildWASACalibrationLayerPlane::~TBuildWASACalibrationLayerPlane() {}

#ifdef ROOT6
ReturnRes::InfoM TBuildWASACalibrationLayerPlane::operator()(const EventWASAUnpack& event,
                                                         FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree)
{
  int result = Exec(event, RecoEvent, OutTree);

  return SoftExit(result);
}
#else
ReturnRes::InfoM TBuildWASACalibrationLayerPlane::operator()(const EventWASAUnpack& event,
                                                         FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree)
{
  int result = Exec(event, RecoEvent, OutTree);

  return SoftExit(result);
}

#endif
void TBuildWASACalibrationLayerPlane::SelectHists() { LocalHisto.h_stats = AnaHisto->CloneAndRegister(AnaHisto->h_stats); }

ReturnRes::InfoM TBuildWASACalibrationLayerPlane::SoftExit(int return_build)
{
  if(return_build == -1)
    {
      att._logger->warn("!> Multiplicity > 2 on Start : event rejected");
      LocalHisto.h_stats->Fill("start M>2", 1);
      return ReturnRes::MultiS2_Start;
    }
  else if(return_build == -2)
    {
      att._logger->warn("!> TDC Timing Start cut : event rejected");
      LocalHisto.h_stats->Fill("start Timing cut", 1);
      return ReturnRes::StartTimingCut;
    }
  else if(return_build == -3)
    {
      att._logger->warn("!> Chamber Hit > 1000 : event rejected");
      LocalHisto.h_stats->Fill("chamber hit>1000", 1);
      return ReturnRes::ChamberHitLimit;
    }
  else if(return_build == -9)
    {
      att._logger->warn("!> No Beam : event rejected");
      LocalHisto.h_stats->Fill("No Beam", 1);
      return ReturnRes::NoBeam;
    }
  else if(return_build != 0)
    {
      att._logger->warn("Error in Build Detector !");
      LocalHisto.h_stats->Fill("Error", 1);
      return ReturnRes::BuildError;
    }
  LocalHisto.h_stats->Fill("start Ok", 1);

  return ReturnRes::Fine;
}
#ifdef ROOT6
int TBuildWASACalibrationLayerPlane::Exec(const EventWASAUnpack& event,
                                      FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree)
#else
int TBuildWASACalibrationLayerPlane::Exec(const EventWASAUnpack& event,
                                      FullRecoEvent& RecoEvent, Ana_WasaEvent* OutTree)
#endif
{
  att._logger->info("Start Build Calibration of Wasa Detector !");
  // Sub classes of calibrations
  // Fiber : Input event.s2fiber 

  // MDC : Input event.s2mdc
  
  // PSB - PSBE - PSFE - T0 : Input event.s2tq1 event.s2tq2
  
  // CsI
  
  // S4 - FRSTPC


  // After calibration : 
  // creation of the genfit measurements
  // collection of additional information that does not go into genfit measurements (i.e ToF, dE, ion optics ...)

  // Save calibration into Ana_WasaEvent for possible restart after this point.











// #ifdef DEBUG_BUILD2
//                   std::cout << "fiber" << std::endl;
//                   std::string tmpName = orderDetName.find(TypeDet)->second;
//                   std::cout << "name : " << tmpName << std::endl;
//                   std::cout << "LayerID : " << LayerID << std::endl;
//                   std::cout << "HitPosX : " << hit.HitPosX << std::endl;
//                   std::cout << "HitPosY : " << hit.HitPosY << std::endl;
//                   std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
//                   gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->Print();
//                   gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix()->Print();
//                   gGeoManager->GetVolume("MFLD")
//                       ->GetNode(motherName.c_str())
//                       ->GetVolume()
//                       ->GetNode((volumeName + "_0").c_str())
//                       ->Print();
//                   gGeoManager->GetVolume("MFLD")
//                       ->GetNode(motherName.c_str())
//                       ->GetVolume()
//                       ->GetNode((volumeName + "_0").c_str())
//                       ->GetMatrix()
//                       ->Print();
//                   gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->Print();
//                   gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix()->Print();
// #endif
//                   TGeoMatrix* g1 =
//                       gGeoManager->GetVolume(volumeName.c_str())->GetNode(LayerID * 2 + 1)->GetMatrix(); // fiber core
//                   TGeoMatrix* g2 = gGeoManager->GetVolume("MFLD")
//                                        ->GetNode(motherName.c_str())
//                                        ->GetVolume()
//                                        ->GetNode((volumeName + "_0").c_str())
//                                        ->GetMatrix(); // fiber layer
//                   TGeoMatrix* g3 =
//                       gGeoManager->GetVolume("MFLD")->GetNode(motherName.c_str())->GetMatrix(); // fiber station
//                   TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();     // MFLD
//                   TGeoHMatrix H1(*g1), H2(*g2), H3(*g3), H4(*g4);
//                   TGeoHMatrix H = H2 * H1;
//                   H             = H3 * H;
//                   H             = H4 * H;
//                   TGeoHMatrix w1("w1");
//                   TGeoHMatrix w2("w2");
//                   w1.SetDz(-10);
//                   w2.SetDz(10);
//                   TGeoHMatrix Hw1 = H * w1;
//                   TGeoHMatrix Hw2 = H * w2;
// #ifdef DEBUG_BUILD2
//                   H.Print();
//                   Hw1.Print();
//                   Hw2.Print();
// #endif
//                   double* edge1 = Hw1.GetTranslation();
//                   double* edge2 = Hw2.GetTranslation();
//                   double* shift = H.GetTranslation();
//                   TVector3 o(0., 0., shift[2]), zdir(0., 0., 1.);
//                   TVector3 fiber_dir(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
//                   fiber_dir  = fiber_dir.Unit();
//                   TVector3 u = fiber_dir.Cross(zdir);
//                   TVector3 v = fiber_dir;
//                   genfit::SharedPlanePtr plane(new genfit::DetPlane(o, u, v));

//                   TVectorD hitCoords(1);
//                   hitCoords(0) = u.Dot(TVector3(shift[0], shift[1], 0));
//                   TMatrixDSym hitCov(1);
//                   hitCov(0, 0) = resolution_fiber * resolution_fiber;
//                   measurement =
//                       std::make_unique<genfit::PlanarMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
//                   dynamic_cast<genfit::PlanarMeasurement*>(measurement.get())->setPlane(plane);

//                   hitCoordsTree(0) = hit.HitPosX;
//                   hitCoordsTree(1) = hit.HitPosY;
//                   hitCoordsTree(2) = hit.HitPosZ;
//                 }

//               else if(IsWire(TypeDet))
//                 {
// #ifdef DEBUG_BUILD2
//                   std::cout << "wire" << std::endl;
//                   std::string tmpName = orderDetName.find(TypeDet)->second;
//                   std::cout << "name : " << tmpName << std::endl;
//                   std::cout << "LayerID : " << LayerID << std::endl;
//                   std::cout << "HitPosX : " << hit.HitPosX << std::endl;
//                   std::cout << "HitPosY : " << hit.HitPosY << std::endl;
//                   std::cout << "HitPosZ : " << hit.HitPosZ << std::endl;
//                   gGeoManager->GetVolume("INNER")
//                       ->GetNode(TypeDet - G4Sol::MG01 + 1)
//                       ->GetVolume()
//                       ->GetNode(LayerID - 1)
//                       ->Print();
//                   gGeoManager->GetVolume("INNER")
//                       ->GetNode(TypeDet - G4Sol::MG01 + 1)
//                       ->GetVolume()
//                       ->GetNode(LayerID - 1)
//                       ->GetMatrix()
//                       ->Print();
//                   gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->Print();
//                   gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetMatrix()->Print();
//                   gGeoManager->GetVolume("MFLD")->GetNode(0)->Print();
//                   gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix()->Print();
//                   gGeoManager->GetVolume("WASA")->GetNode(0)->Print();
//                   gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix()->Print();
// #endif
//                   TGeoMatrix* g1 = gGeoManager->GetVolume("INNER")
//                                        ->GetNode(TypeDet - G4Sol::MG01 + 1)
//                                        ->GetVolume()
//                                        ->GetNode(LayerID - 1)
//                                        ->GetMatrix(); // ME, MG
//                   TGeoShape* tempShape =
//                       gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetVolume()->GetShape();
//                   TGeoMatrix* g2 =
//                       gGeoManager->GetVolume("INNER")->GetNode(TypeDet - G4Sol::MG01 + 1)->GetMatrix(); // MD
//                   TGeoMatrix* g3 = gGeoManager->GetVolume("MFLD")->GetNode(0)->GetMatrix();             // INNER
//                   TGeoMatrix* g4 = gGeoManager->GetVolume("WASA")->GetNode(0)->GetMatrix();             // MFLD
//                   TGeoHMatrix H1(*g1), H2(*g2), H3(*g3), H4(*g4);
//                   TGeoHMatrix H = H2 * H1;
//                   H             = H3 * H;
//                   H             = H4 * H;
//                   double* shift = H.GetTranslation();
//                   TGeoHMatrix w1("w1");
//                   TGeoHMatrix w2("w2");
//                   Double_t minZ, maxZ;
//                   tempShape->GetAxisRange(3, minZ, maxZ);
//                   w1.SetDz(minZ);
//                   w2.SetDz(maxZ);
//                   TGeoHMatrix Hw1 = H * w1;
//                   TGeoHMatrix Hw2 = H * w2;
// #ifdef DEBUG_BUILD2
//                   H.Print();
//                   Hw1.Print();
//                   Hw2.Print();
// #endif
//                   double* edge1 = Hw1.GetTranslation();
//                   double* edge2 = Hw2.GetTranslation();

//                   TVector3 x1(shift[0], shift[1], shift[2]);
//                   TVector3 p1(edge2[0] - edge1[0], edge2[1] - edge1[1], edge2[2] - edge1[2]);
//                   TVector3 x2(hit.HitPosX, hit.HitPosY, hit.HitPosZ);
//                   TVector3 p2(hit.MomX, hit.MomY, hit.MomZ);
//                   double dl = CloseDist(x1, x2, p1, p2);

//                   double dlmax = 0;
//                   switch(TypeDet - G4Sol::MG01 + 1)
//                     {
//                     case 1:
//                     case 2:
//                     case 3:
//                     case 4:
//                     case 5:
//                       dlmax = 0.2;
//                       break;
//                     case 6:
//                     case 7:
//                     case 8:
//                     case 9:
//                     case 10:
//                     case 11:
//                       dlmax = 0.3;
//                       break;
//                     case 12:
//                     case 13:
//                     case 14:
//                     case 15:
//                     case 16:
//                     case 17:
//                       dlmax = 0.4;
//                       break;
//                     default:
//                       att._logger->warn("Error in WireMeasurement !");
//                       break;
//                     }
//                   double temp_dl = gRandom->Gaus(dl, resolution_dl);
//                   bool doneRand  = false;
//                   while(doneRand)
//                     {
//                       if(temp_dl < 0 || temp_dl > dlmax)
//                         temp_dl = gRandom->Gaus(dl, resolution_dl);
//                       else
//                         doneRand = true;
//                     }
//                   // if(temp_dl<0)     dl = 0;
//                   // if(temp_dl>dlmax) dl = dlmax;

//                   TVectorD hitCoords(7);
//                   hitCoords(0) = edge1[0];
//                   hitCoords(1) = edge1[1];
//                   hitCoords(2) = edge1[2];
//                   hitCoords(3) = edge2[0];
//                   hitCoords(4) = edge2[1];
//                   hitCoords(5) = edge2[2];
//                   hitCoords(6) = temp_dl;
//                   if(edge1[2] > edge2[2])
//                     {
//                       hitCoords(0) = edge2[0];
//                       hitCoords(1) = edge2[1];
//                       hitCoords(2) = edge2[2];
//                       hitCoords(3) = edge1[0];
//                       hitCoords(4) = edge1[1];
//                       hitCoords(5) = edge1[2];
//                     }
//                   TMatrixDSym hitCov(7);
//                   hitCov(6, 6) = resolution_dl * resolution_dl;
//                   measurement =
//                       std::make_unique<genfit::WireMeasurement>(hitCoords, hitCov, int(TypeDet), LayerID, nullptr);
//                   dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setLeftRightResolution(0);
//                   dynamic_cast<genfit::WireMeasurement*>(measurement.get())->setMaxDistance(dlmax);

//                   hitCoordsTree(0) = hit.HitPosX;
//                   hitCoordsTree(1) = hit.HitPosY;
//                   hitCoordsTree(2) = hit.HitPosZ;
//                 }
 

  return 0;
}

