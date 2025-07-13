// $Id$

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "PhysicalObject.h"
#include "TypeDefs.h"
#include "Specie.h"

// LibMesh includes
#include "libmesh/point.h"
#include "libmesh/tensor_value.h"

// C++ includes
#include <string>
#include <set>
#include <vector>
#include <cassert>

// forward declarations
class Dopant;
class BulkCrystal;


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


    //! Get the bulk crystal
    const BulkCrystal* get_bulk_crystal(void) const;

    //! Get the Bloch basis functions defining the CB
    const std::vector<std::string>& get_cb_bloch_functions(void) const; 

    //! Get the Bloch basis functions defining the VB
    const std::vector<std::string>& get_vb_bloch_functions(void) const; 

    //! Get the atomic orbital basis functions defining the CB
    /*!
     * The first index refers to the atom in the basis
     */
    const std::vector<std::vector<std::string>>& get_cb_atomic_orbitals(void) const; 

    //! Get the atomic orbital basis functions defining the VB
    /*!
     * The first index refers to the atom in the basis
     */
    const std::vector<std::vector<std::string>>& get_vb_atomic_orbitals(void) const; 

    //! Get the rotation matrix
    /*!
     * For a crystal, this is the rotation from standard crystal
     * coordinates into the calculation system.
     */
    const libMesh::RealTensor& get_rotation_matrix(void) const;


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
    
    //! iterators over all species in one specific position
    crystal_species_iterator species_begin(void) const;

    crystal_species_iterator species_end(void) const;

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
    virtual void do_preinit(void) {};


    //! The real init function
    virtual void do_init(void);


    //! Setup the doping
    void setup_doping(void);


    //! Can be used to print specific information
    virtual void do_info(void) const {};


    //! The rotation matrix
    /*!
     * The rotation matrix is taken from BulkCrystal,
     * whenever it can be defined, otherwise it will be
     * the unity matrix.
     * In future, we might add an additional rotation due
     * to macroscopic deformation.
     */
    libMesh::RealTensor _rotation_matrix {1.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0, 1.0};

    //! Set the RotatedCrystal
    void set_crystal(BulkCrystal* crystal);

    //! Fill list of species, to be used during initialization
    virtual void fill_species(void);

    
     //! List of all species in material (irrespective of their label) 
    std::set<Specie> _species;

    /*! Map atomic labels into a set of species. 
     *  e.g. InGaAs has atoms with label 1 -> (In, Ga)
     *              and atoms with label 2 -> As
     */ 
    std::map<unsigned int, std::set<Specie>> _crystal_type_map;

    //! The Bloch basis states contributig to the CB
    std::vector<std::string> _cb_bloch_states;

    //! The Bloch basis states contributig to the VB
    std::vector<std::string> _vb_bloch_states;

    //! The atomic orbitals contributing to the CB
    /*!
     * Here the first index refers to the atom in the basis
     */
    std::vector<std::vector<std::string>> _cb_atomic_orbitals;

    //! The atomic orbitals contributing to the VB
    /*!
     * Here the first index refers to the atom in the basis
     */
    std::vector<std::vector<std::string>> _vb_atomic_orbitals;


  private:


    //! The crystal structure
    /*!
     * The crystal structure as wz, zb etc
     */
    std::string _structure {""};


    //! True if this is an alloy
    bool _is_alloy;


    //! The bulk crystal object
    BulkCrystal* _bulk_crystal {nullptr};


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
const BulkCrystal*
Material::get_bulk_crystal(void) const
{
  return _bulk_crystal;
}


inline
const libMesh::RealTensor&
Material::get_rotation_matrix(void) const
{
  return _rotation_matrix;
}



inline
void
Material::set_structure(const std::string& structure)
{
  _structure = structure;
}

inline
const std::vector<std::string>&
Material::get_cb_bloch_functions(void) const
{
  return _cb_bloch_states;
} 

inline
const std::vector<std::string>&
Material::get_vb_bloch_functions(void) const
{
  return _vb_bloch_states;
}

inline
const std::vector<std::vector<std::string>>&
Material::get_cb_atomic_orbitals(void) const
{
  return _cb_atomic_orbitals;
}

inline
const std::vector<std::vector<std::string>>&
Material::get_vb_atomic_orbitals(void) const
{
  return _vb_atomic_orbitals;
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
Material::crystal_species_iterator
Material::species_begin(void) const
{
  return (_species.begin());
}

inline
Material::crystal_species_iterator
Material::species_end(void) const
{
  return (_species.end());
}

inline
unsigned int
Material::count_species(unsigned int label) const
{
  if (_crystal_type_map.count(label) > 0)
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

