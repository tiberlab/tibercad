// $Id$

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "PhysicalObject.h"
#include "TypeDefs.h"
#include "Specie.h"

// LibMesh includes
#include "point.h"

// C++ includes
#include <string>
#include <set>
#include <vector>
#include <cassert>

// forward declarations
class Dopant;
class RotatedCrystal;


//! Contains all needed data for a material
/*!
 * For any material/structure combination that is used for a simulation
 * a \c Material object is built which contains all needed physical
 * models. Every simulation will have exactly one model in the models list
 * which it can use for its calculations. Additionally, a Material object
 * contains also the list of donors and acceptors (see Dopant) and a
 * RotatedCrystal object.
 */
class Material : public PhysicalObject
{

  public:


    //! An iterator to iterate over all dopants
    typedef std::set<Dopant*>::iterator dopant_iterator;

    typedef std::set<Specie>::iterator crystal_species_iterator;
  
    //! A const iterator to iterate over all dopants
    //typedef std::set<Dopant*>::const_iterator const_dopant_iterator;


    //! Destructor
    /*!
     * Deletes all \c PhysicalModel objects
     */
    virtual ~Material(void);


    //! Create a bulk material with name \c name and options
    static Material* create(const std::string& name,
        const ModelOptions& options);


    //! Tells whether it is an alloy
    bool is_alloy(void) const;


    //! Do some preparatory work at creation time
    void preinit(void);


    //! Add a dopant
    void add_dopant(Dopant* dopant);


    //! Set the structure
    void set_structure(const std::string& structure);


    //! Get the crystal structure
    const std::string& get_structure(void) const;


    //! Get a reference to the RotatedCrystal
    const RotatedCrystal& get_rotated_crystal(void) const;


    //! Get the total n-doping
    double get_total_donor_density(void) const;


    //! Get the total p-doping
    double get_total_acceptor_density(void) const;


    //! Get the total doping density
    /*!
     * The return value is \f$N_d + N_a\f$
     */
    double get_total_doping_density(void) const;


    //! Get the total net doping density
    /*!
     * The return value is \f$N_d - N_a\f$
     */
    double get_net_doping_density(void) const;


    //! Get the first iterator for the donors
    dopant_iterator donors_begin(void) const;

    //! Get the past-the-end iterator for the donors
    dopant_iterator donors_end(void) const;


    //! Get the first iterator for the acceptors
    dopant_iterator acceptors_begin(void) const;

    //! Get the past-the-end iterator for the acceptors
    dopant_iterator acceptors_end(void) const;


    //! Print some information
    void info(void) const;
 
    //! Tells if a specie belongs to the material
    bool has_specie(Specie) const;

    //! tells if a specie occupies a defined position in lattice
    bool is_specie(Specie, unsigned int) const; 
   
    //! gives label (atom position in unit cell) from specie
    unsigned int get_label(Specie) const;

    //! count all atomic species of a given label 
    unsigned int count_species(unsigned int) const;

    //! returns all possible atomic labels  
    unsigned int count_labels() const;

    //! iterators over all species in one specific position
    crystal_species_iterator species_begin(unsigned int label) const;
  
    crystal_species_iterator species_end(unsigned int label) const;
    
    //! Print the list of atomic species in material
    void print_species(void) const;

  protected:

    //! Construct a material
    /*!
     *
     * \param name the name of the Material
     */
    Material(const std::string& name, const ModelOptions& options,
        bool alloy = false);


    //! The real preinit function
    /*!
     * Prepares the rotated crystal object
     */
    virtual void do_preinit(void);


    //! The real init function
    virtual void do_init(void);


    //! Setup the doping
    void setup_doping(void);


    //! Can be used to print specific information
    virtual void do_info(void) const {};


    //! Get a writable pointer to the RotatedCrystal
    RotatedCrystal* get_crystal(void);


    //! Set the RotatedCrystal
    void set_crystal(RotatedCrystal* crystal);

    //! Fill list of species, to be used during initialization
    virtual void fill_species(void);

    
     //! List of all species in material (irrespective of their label) 
    std::set<Specie> _species;

    /*! Map atomic labels into a set of species. 
     *  e.g. InGaAs has atoms with label 1 -> (In, Ga)
     *              and atoms with label 2 -> As
     */ 
    std::map<unsigned int, std::set<Specie>> _crystal_type_map;

  private:


    //! The crystal structure
    /*!
     * The crystal structure as wz, zb etc
     */
    std::string _structure;


    //! True if this is an alloy
    bool _is_alloy;


    //! The RotatedCrystal object
    RotatedCrystal* _rotated_crystal;


    //! The list of donors
    std::set<Dopant*> _donors;

    //! The list acceptors
    std::set<Dopant*> _acceptors;


    //! Clear all doping
    void clear_doping(void) TBDLLOCAL;


};


//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------

inline
bool
Material::is_alloy(void) const
{
  return _is_alloy;
}



inline
const std::string&
Material::get_structure(void) const
{
  return  _structure;
}


inline
const RotatedCrystal&
Material::get_rotated_crystal(void) const
{
  return *_rotated_crystal;
}


inline
RotatedCrystal*
Material::get_crystal(void)
{
  return _rotated_crystal;
}


inline
void
Material::set_structure(const std::string& structure)
{
  _structure = structure;
}




inline
double
Material::get_total_doping_density(void) const
{
  return (get_total_donor_density() + get_total_acceptor_density());
}


inline
double
Material::get_net_doping_density(void) const
{
  return (get_total_donor_density() - get_total_acceptor_density());
}


inline
Material::dopant_iterator
Material::donors_begin(void) const
{
  return _donors.begin();
}


inline
Material::dopant_iterator
Material::donors_end(void) const
{
  return _donors.end();
}


inline
Material::dopant_iterator
Material::acceptors_begin(void) const
{
  return _acceptors.begin();
}


inline
Material::dopant_iterator
Material::acceptors_end(void) const
{
  return _acceptors.end();
}

  
inline
Material::crystal_species_iterator 
Material::species_begin(unsigned int label) const
{
  if (_crystal_type_map.count(label)>0)
    return (_crystal_type_map.find(label)->second).begin();
  else
    return (_crystal_type_map.find(1)->second).end();
}

inline
Material::crystal_species_iterator 
Material::species_end(unsigned int label) const
{
  if (_crystal_type_map.count(label)>0)
    return (_crystal_type_map.find(label)->second).end();
  else
    return (_crystal_type_map.find(1)->second).end();
}

inline
unsigned int
Material::count_species(unsigned int label) const
{
  if (_crystal_type_map.count(label)>0)
    return (_crystal_type_map.find(label)->second).size();
  else
    return 0;
}

inline
unsigned int
Material::count_labels(void) const
{
  return _crystal_type_map.size();
}


#endif // _MATERIAL_H_



    /*
    OLD ITERATOR:

    class crystal_species_iterator
    {
      public:
      crystal_species_iterator(const std::map<unsigned int, std::set<Specie>>& map, 
                               unsigned int label,
                               bool end = false) :
        _map(map),
        _iter(map.begin()),
        _label(label)
      {
        if (!end && _iter != _map.end())
        {
          while ( (_iter != _map.end()) && (label != _iter->second ) )
            ++_iter;
        }
        else
          _iter = _map.end();       
      }
      
      //! Copy constructor
      crystal_species_iterator(const crystal_species_iterator& it) :
        _map(it._map),
        _iter(it._iter),
        _label(it._label)
      { };
        
      //! Prefix increment
      crystal_species_iterator& operator++(void)
      {
        do{
          ++_iter;
        }
        while ( (_iter != _map.end()) && (_iter->second != _label) );
        return *this;
      }
      
      // Assignment
      //crystal_species_iterator& operator=(const  crystal_species_iterator& rhs)
      //{
      //  _map = rhs._map;
      //  _iter = rhs._iter;
      //  _label = rhs._label;
      //  return *this;
      //}

      //! Comparison
      bool operator==(const crystal_species_iterator& rhs)
      {
        return (_iter == rhs._iter);
      }
      
      //! Comparison
      bool operator!=(const crystal_species_iterator& rhs)
      {
        return (_iter != rhs._iter);
      }
      
      //! Dereference the iterator to get the Specie
      Specie operator*(void)
      {
        return _iter->first;
      }
      
    private:
      const std::map<unsigned int, std::set<Specie>>& _map;
      unsigned int _label;
      std::set<Specie>::const_iterator _iter;
    };
    */
