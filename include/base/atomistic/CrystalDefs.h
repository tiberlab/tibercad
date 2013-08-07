#include "Specie.h"
#include <set>
#include <map>

class CrystalDefs {

  private:


  public:

    static bool is_anion(const std::string material, const Specie sp);

    static bool is_cation(const std::string material, const Specie sp);

    static std::map<std::string, Specie::Type >
     anion;
    
    static std::map<std::string, Specie::Type >
     cation;
};

