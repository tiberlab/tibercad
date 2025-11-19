#include "tibercad/atomistic/Specie.h"

#include <set>
#include <vector>
#include <map>

//! This class contains static definitions and utilities for Crystals
class CrystalDefs {

  public:

    static bool is_anion(const std::string material, const Specie sp);

    static bool is_cation(const std::string material, const Specie sp);

    //! Convert any type of symmetry name to its International Symbol
    static std::string convert_to_international_symbol(const std::string& symmetry);

    //! Convert Schoenflies notation to International Symbol
    static std::string schoenflies_to_IS(const std::string& schoenflies);

    //! Convert Space group Symbol or number to International Symbol
    static std::string spacegroup_to_IS(const std::string& sp_grp);

    //! Convert International Symbol to crystal class
    static std::string IS_to_crystal_class(const std::string& int_sym);

    //! Convert crystal class to International Symbol
    static std::string crystal_class_to_IS(const std::string& cclass);

    //! Convert International Symbol to crystal system
    static std::string IS_to_crystal_system(const std::string& int_sym);

    //! Convert given lattice type string to a standard form
    /*!
     * Standard form is one of the 14 Bravais lattices
     */
    static std::string get_bravais_lattice(const std::string& name);

    static std::string bravais_short_to_long_name(const std::string& in);


    static std::map<std::string, Specie::Type > anion;
    
    static std::map<std::string, Specie::Type > cation;



  private:
    
    static unsigned int spacegroup_to_IS_id(const std::string& sp_grp);

    //! International Symbols
    static const std::vector<std::string> _international_sym;

    //! Schoenflies Symbols, ordered as International Symbols
    static const std::vector<std::string> _schoenflies;

    //! Crystal classes, order as international symbols
    static const std::vector<std::string> _crystal_classes;

    //! Crystal class index to crystal system
    static const std::vector<unsigned int> _class_to_system;

    //! The space groups
    static const std::vector<std::string> _space_groups;

};

