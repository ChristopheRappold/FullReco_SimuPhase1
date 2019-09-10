#ifndef FULLRECOCONFIG_h
#define FULLRECOCONFIG_h

#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace pt = boost::property_tree;

class FullRecoConfig
{
public:
  FullRecoConfig();
  FullRecoConfig(int argc, char** argv);
  ~FullRecoConfig() = default;
  void ParseConfig(const std::string& namefile);
  int ProperConf() { return status; }

  template <typename T>
  T Get(const std::string& key) const
  {
    T temp = tree.get<T>(key);
    if(boost::optional<std::string> unit = tree.get_optional<std::string>(key + ".unit"))
      {
        std::string unitName(*unit);
        T unitVal = tree.get<T>(key + ".unit." + unitName);
        return temp * unitVal;
      }
    else
      return temp;
  }

  template <typename T>
  void Add(const std::string& key, const T& val)
  {
    tree.put(key, val);
  }

  bool IsAvailable(const std::string& key) const
  {
    size_t count = tree.count(key);
    if(count > 0)
      return true;
    else
      return false;
  }
  std::string CheckConfig();

private:
  pt::ptree tree;
  int status;

  int ParseCmd(int argc, char** argv);
  void ParseLine(std::stringstream& lineStream, std::vector<std::string>& res);
  std::string display(const int depth, const pt::ptree& t);
  double GetDimension(const std::string&);
  void SetDefault();
};

template <>
inline std::string FullRecoConfig::Get(const std::string& key) const
{
  return tree.get<std::string>(key);
}
template <>
inline boost::optional<std::string> FullRecoConfig::Get(const std::string& key) const
{
  return tree.get_optional<std::string>(key);
}

#endif
