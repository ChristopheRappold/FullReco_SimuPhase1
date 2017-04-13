#ifndef THYPHIATTRIBUTE
#define THYPHIATTRIBUTE
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>
//#include "TRandom3.h"
#include "FairBase/FrsSksHypFieldMapFull.h"
#include "FairBase/HypConstField.h"
#include "FairBase/HypFieldMap.h"
#include "FairBase/HypFieldMapFull.h"

#include "FairBase/FrsSolenoidHypField.h"

struct DataSim
{
  std::vector<std::string>* nameDet;
  std::map<std::string, double>* simParameters;
};

class THyphiAttributes
{

  public:
  int Nb_CPU;
  int Nb_Fraction;

  bool G4_simu;
  bool G4_TimeResolution;
  bool G4_GeoResolution;

  bool back_tracking;

  double Target_PositionX;
  double Target_PositionY;
  double Target_PositionZ;

  double Target_Size;
  double Field_Strength;

  std::vector<std::string> name_GeoVolumes;

  bool beam_only;
  bool Debug_DAF;
  bool DoNoMaterial;

  FairField* Field;

  const DataSim& InputPar;

  THyphiAttributes();
  THyphiAttributes(const std::list<std::string>& type, const std::list<std::string>& option, double FS, const DataSim& InputParameters);
  ~THyphiAttributes();

  int Init_Para();
};

#endif
