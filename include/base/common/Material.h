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
    void fill_species(void);


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

    //! List of all species 
    std::set<Specie> _species;

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
bool
Material::has_specie(Specie sp) const 
{
  if (_species.find(sp) == _species.end())
    return true;
  else return false; 
}


#endif // _MATERIAL_H_
