#include "FullRecoConfig.hh"

#include <boost/foreach.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <fstream>
#include <getopt.h>
#include <sstream>
#include <cmath>



static struct option optlong[] = {{"help", 0, NULL, 'h'},  {"cpu", 1, NULL, 'c'},   {"num", 1, NULL, 'n'},
                                  {"event", 1, NULL, 'e'}, {"start", 1, NULL, 's'}, {"geo", 1, NULL, 'g'},
                                  {"log", 1, NULL, 'l'},   {"config", 1, NULL, 'f'},{"gnn", 1, NULL, 'a'}};

FullRecoConfig::FullRecoConfig() { SetDefault(); }

FullRecoConfig::FullRecoConfig(int argc, char** argv)
{
  SetDefault();
  status = ParseCmd(argc, argv);
}

void FullRecoConfig::ParseLine(std::stringstream& lineStream, std::vector<std::string>& res)
{
  std::string lastStr;
  lineStream >> lastStr;
  if(lastStr.empty())
    return;
  else
    {
      res.emplace_back(lastStr);
      return ParseLine(lineStream, res);
    }
}

void FullRecoConfig::ParseConfig(const std::string& namefile)
{
  std::cout << "start reading" << std::endl;
  std::ifstream ifs(namefile.c_str());
  if(ifs.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs, temp_line))
        {
          std::stringstream stream(temp_line);
          std::string testComment(stream.str());
          auto it_comment = testComment.find(CommentSymbol);
          if(it_comment != std::string::npos)
            {
              // std::cout<<"!> Skip "<<test<<std::endl;
              continue;
            }

          std::vector<std::string> out_list;
          ParseLine(stream, out_list);

          switch(out_list.size())
            {
            case 2:
              tree.put(out_list[0], out_list[1]);
              break;
            case 3:
              {
                std::string key1(out_list[0] + ".unit");
                tree.put(out_list[0], out_list[1]);
                tree.put(key1, out_list[2]);
                std::string key2(key1 + "." + out_list[2]);
                double valUnit = GetDimension(out_list[2]);
                tree.put(key2, valUnit);
                break;
              }
            default:
              std::cout << "!> stream of line from config file with unknown parsing " << out_list.size() << " \n";
              break;
            }
        }
    }
}

int FullRecoConfig::ParseCmd(int argc, char** argv)
{
  auto print_help = [&argv](const std::string& specific_msg) {
    std::cout << specific_msg << " \n";
    std::cout << "!> Example of use:\n";
    std::cout << "!> " << argv[0];
    std::cout << "[-g Geofile] [--geo Geofile] [-c nb_cpu] [--cpu nb_cpu] [-n fraction] [--num fraction] [-s start_ev] "
                 "[--start start_ev] [-e nb_event] [--event nb_event] [-l lvllog] [--log lvllog] [-a gnntext] [--gnn gnntext] [-h]  OutputFile "
                 "RootInputFile \n";
    std::cout << "\n";
  };

  if(argc < 3)
    {
      print_help("E> Wrong number of parameters!");
      return -1;
    }

  int option_char;
  int Nb_CPU      = 1;
  int Nb_event    = -1;
  int Nb_fraction = 1;
  int Start       = 0;
  int Log_lvl     = 1;
  std::string nameGeo("./geo/GeoSolenoid.root");
  std::string nameConf("testconfig.cfg");
  std::string nameGnn("test.txt");

  while((option_char = getopt_long(argc, argv, "+hc:n:e:s:g:l:f:a:", optlong, NULL)) != EOF)
    switch(option_char)
      {
      case 'h':
        print_help("!> Help:");
        return -1;
        break;
      case 'c':
        Nb_CPU = std::atoi(optarg);
        tree.put("Nb_CPU", Nb_CPU);
        break;
      case 'n':
        Nb_fraction = std::atoi(optarg);
        tree.put("Nb_Fraction", Nb_fraction);
        break;
      case 'e':
        Nb_event = std::atoi(optarg);
        tree.put("Nb_Event_Cfg", Nb_event);
        break;
      case 's':
        Start = std::atoi(optarg);
        tree.put("Start", Start);
        break;
      case 'g':
        nameGeo = std::string(optarg);
        tree.put("Geo", nameGeo);
        break;
      case 'l':
        Log_lvl = std::atoi(optarg);
        tree.put("Log_Lvl", Log_lvl);
        break;
      case 'f':
        nameConf = std::string(optarg);
        tree.put("Config", nameConf);
        break;
      case 'a':
        nameGnn = std::string(optarg);
        tree.put("GNN_Node", nameGnn);
        break;
      case '?':
      default:
        print_help("E> Wrong number of parameters!");
        return -1;
      }

  std::string name_in, name_out;

  if(optind == argc)
    {
      std::stringstream ss1;
      ss1 << "E> Input and output Rootfile are missing !" << optind << " " << argc;
      print_help(ss1.str());
      return -1;
    }
  else
    {
      name_out = argv[optind];
      name_in  = argv[optind + 1];

      tree.put("Output_Namefile", name_out);
      tree.put("Input_Namefile", name_in);
    }

  ParseConfig(nameConf);

  return 0;
}

std::string FullRecoConfig::display(const int depth, const pt::ptree& t)
{
  std::stringstream ss1;
  BOOST_FOREACH(pt::ptree::value_type const& v, t.get_child(""))
    {
      pt::ptree subtree   = v.second;
      std::string nodestr = t.get<std::string>(v.first);

      // print current node
      ss1 << std::string("").assign(depth * 2, ' ') << "* ";
      ss1 << v.first;
      if(nodestr.length() > 0)
        ss1 << " -> \"" << t.get<std::string>(v.first) << "\"";
      ss1 << std::endl;

      // recursive go down the hierarchy
      ss1 << display(depth + 1, subtree);
    }
  return ss1.str();
}

std::string FullRecoConfig::CheckConfig() { return display(0, tree); }

void FullRecoConfig::Save(std::string& Sjson) const
{
  std::stringstream outs;
  pt::write_json(outs, tree,false);
  Sjson = outs.str();
}

int FullRecoConfig::Reload(const std::vector<char>& confjson) const
{
  std::stringstream ins;
  for(auto c : confjson)
    ins << c;

  pt::ptree tempTree;
  pt::read_json(ins,tempTree);

  std::string geo1 = Get<std::string>("Geo");
  std::string geo2 = tempTree.get<std::string>("Geo");

  int same = 0;
  if(geo1 != geo2)
    ++same;

  return same;
}

double FullRecoConfig::GetDimension(const std::string& dimension)
{
  // Time
  if(dimension == "ns")
    return 1.e-6;
  if(dimension == "ms")
    return 1.e-3;
  if(dimension == "s")
    return 1.;
  // Length
  if(dimension == "m")
    return 0.01;
  if(dimension == "cm")
    return 1.;
  if(dimension == "mm")
    return 0.1;
  // Fields
  if(dimension == "tesla")
    return 1.;
  if(dimension == "kG")
    return 0.1;
  // Energy
  if(dimension == "eV")
    return 1.e-6;
  if(dimension == "keV")
    return 1.e-3;
  if(dimension == "MeV")
    return 1.;
  if(dimension == "GeV")
    return 1.e3;
  // Angle
  if(dimension == "degree")
    return M_PI / 180.;
  if(dimension == "rad")
    return 1.;
  //
  std::cerr << "!> Unknown dimension " << dimension << "\n";
  std::cerr << "!> Exiting !!!\n";
  std::exit(1);
}

void FullRecoConfig::SetDefault()
{
  tree.put("Nb_CPU", 1);
  tree.put("Nb_Event_Cfg", -1);
  tree.put("Nb_Fraction", 1);
  tree.put("Start", 0);
  tree.put("Log_Lvl", 1);

  // tree.put("TargetRegionCut.unit","mm");
  // tree.put("TargetRegionCut.unit.mm",1.*mm);
}
