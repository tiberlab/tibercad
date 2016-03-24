#include "Specie.h"
#include <set>
#include <map>

class CrystalDefs {

  private:

  // std::map<std::string, std::map<Specie::Type, unsigned int> > ion_map;

  public:

    static bool is_anion(const std::string material, const Specie sp);

    static bool is_cation(const std::string material, const Specie sp);

    static std::map<std::string, Specie::Type > anion;
    
    static std::map<std::string, Specie::Type > cation;

  //  void add_material_types(const std::string material, const Specie sp, unsigned int label);

  //static bool is_specie(const std::string& material, const Specie sp, unsigned int label);  

};

