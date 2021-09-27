// $Id$

#ifndef _ATOMISTICSTRUCTURE_H_
#define _ATOMISTICSTRUCTURE_H_

#include "tensor.h"
#include "ModelOptions.h"
#include "TypeDefs.h"
#include "Device.h"
#include "Database.h"
#include "InitFailedException.h"
#include "AtomisticBasis.h"

#include <boost/bind/bind.hpp>
#include <boost/function.hpp>

class Material;
namespace libMesh
{
  class UnstructuredMesh;
}

//! Contains all needed data for an atomic structure
/*!
 *For any Atomistic region an atomistic structure is
 *built, containing atom positions and other informations
 *needed by any atomistic physical model.
 */
class AtomisticStructure: public AtomisticBasis
{
public:

//! A class for Dftb options
  class AtomisticStructureOptions
  {
  public:

    //!Options inizialization
    AtomisticStructureOptions(void);

    //!Options destructor
    ~AtomisticStructureOptions(void);

    //!Tells if association has been already done
    bool is_associated;
  };



  //! Destructor for AtomisticStructure class object
  ~AtomisticStructure(void);

  //! Create a material with name /c name
  static AtomisticStructure* create();

  //! Create a material with name /c name
  static AtomisticStructure* create(const AtomisticStructure& as);

  //! Get the structure options
  ModelOptions& get_options(void);

  //! Get the structure options
  const ModelOptions& get_options(void) const;

  //! Get the structure name
  const std::string& get_name(void) const;


  //! Get set of regions covered by atomistic structure (IDs)
  const std::set<ID>& get_IDset(void) const;

   //! Set the device we're working with
  void set_device(const Device* const device);

  //! Get Device reference
  const Device*  get_device(void) const;

  //! Create a grid aligned on the atoms
  /*!
   * This is tested for tetrahedrically bonded materials (III-V, III-nitrides)
   *
   * \param mesh the mesh to be created
   * \param labels the atoms inside the basis to include in the mesh
   * \param keep_node_order if true, the nodes will be ordered in the same way as the atoms
   */
  void create_conformal_grid(libMesh::UnstructuredMesh& mesh,
      std::set<ID> labels = std::set<ID>(),
      bool keep_node_order = false) const;


  //! Get scale factor (from mesh_units to amstrong mesh_units/1e-10)
  const double& get_scale(void) const;

  //! Initialize a structure (to be read from input file)
  void init(const std::string& name, const Device* const device, const ModelOptions& options);

  //! Initialize a structure (to be read from file)
  void init(const std::string& filename);

  //! Combine structures
  void combine_structures(const std::string& name, const ModelOptions& options);

  //! Remove close atoms 
  void remove_bad_atoms(void );



  //! Print structure to file (format depends on extension used)
  void print_structure(const std::string& path);

  //! A tool for printing atomic charges on output
  void print_structure(const std::string& path, double const* const charges);

  //! Print structure in internal tgn format (contains bondmap and element map)
  void print_tgn(const std::string& path) const;

  //! Print upg file (etb_dataset may be changed into type options)
  void print_upg(const std::string& path, const std::string& etb_dataset, 
                                          bool band_offsets = false);


  //! extract alloy distribution statistics
  /*!
   * \param stats the map to put statistics in
   * \param regions the regions to be used
   * \param cutoff the cutoff radius in nm (default 0.5 nm)
   * \param y if not 0, control volume is a cube instead of sphere
   * \param z together with y defines the control cube
   * \param atomlist (optional) will be filled with all atoms in the control volume
   */
  void extract_statistics(std::map<Specie, std::vector<unsigned int>>& stats,
      const std::set<ID>& regions, double cutoff = 0.5,
      double y = 0, double z = 0,
      const ModelOptions& opt = ModelOptions()) const;

  //! Extract statistics in volume around a single atom
  /*!
   * \param atom center of control volume
   * \param counts the counts of the different species
   * \param regions the regions to be used
   * \param cutoff the cutoff radius in nm (default 0.5 nm)
   */
  void extract_statistics(unsigned int atom, std::map<Specie, unsigned int>& counts,
      const std::set<ID>& regions, double cutoff) const;


  //! Get number of non hydrogen atoms
  unsigned int get_N_without_H(void) const;

  //! Calculate the number of cations
  unsigned int compute_N_cations(void) const;

  //! Find the atom index nearest to a point in an element
  /*!
   * Will return -1 if no atom can be found up to a distance of \c cutoff nm
   */
  int find_nearest_atom(const Elem* elem, const Point& point, double cutoff);

  //! Set the model options
  void set_options(const ModelOptions& options);

  //! Get specie of atom i
  const Specie& get_specie(unsigned int i) const;

  //! Get atom Material
  const Material* get_material(const Atom& atom, bool parent = false) const;

  //! Get bond Material
  /*
   * If a couple of atoms belonging to a binary compound is specified,
   * the material can be assigned without ambiguity to the Anion.
   * The check is now hardcoded for a restricted class of materials (Nitrides
   * and Arsenides). With additional information in Material class, this can be made
   * more generic.
   */
  const Material* get_material(const Atom& atom1, const Atom& atom2,
      bool parent = false) const;

  //! Tells if random alloy approximation is used
  bool is_random_alloy(void);

  // ! Get the reference material
  //const Material* get_reference_material(void) const;

  //! Apply reordering to atoms
  void reorder(const std::vector<unsigned int>& P);

  //! Restrict the atomistic structure to given sub-regions
  void dorestrict(const std::set<ID>& rgn_ids);

   //! Get the atoms in a given mesh element
  const std::vector<unsigned int>& get_atoms_in_elem(const Elem* element) const;

  
  //! Used to set the virtual types. 
  void set_virtual_types(const std::set<std::string>& types);

  //! clear the atom_type set
  void clear_virtual_types(){_virtual_atom_types.clear();};

  //! Get atom type index (types are stored in _virtual_atom_types)
  unsigned int get_virtual_type_index(const std::string &  type) const;
  
  //! used by print upg (temporarily here)
  void assign_virtual_species(void);

  //! Register a callback to be executed after structure creation
  /*!
   * This is used to synchronize random alloy distributions between
   * atomistic and continuous representation
   */
  static void register_callback(std::string& name,
      boost::function<void(void)> callback);


  void interface_interactions(const Material* mat1, 
                              const Material* mat2, 
                              std::vector<std::string>& str,
                              std::vector<double>& frac,
                              const Atom& at1,
                              const Atom& at2);

  //! extract alloy statistics
  void extract_alloy_statistics(const ModelOptions& opt);

  //! plot alloy composition on a vtk file
  void plot_alloy_composition(const ModelOptions& opt);

  //! compute radial distribution
  void radial_distribution(std::string suffix = "");



private:

  //! Constructor for AtomisticStructure class object
  /*!
   * A name to keep reference of structure must
   * be specified
   *
   */
  AtomisticStructure();

  AtomisticStructure(const std::string& name);
  
  //! copy constructor
  AtomisticStructure(const AtomisticStructure& other);
  
  //! Calculate the number of atoms excluding hydrogens, usefule for passivated structures
  void compute_N_without_H(void);
  
  //!Override lattice vectors from structure generation
  void parse_lattice_vectors(void); 
  
  //! Initialize the structure using mesh infos
  void init_periodicity(void);
  
  //! Initialize the structure using mesh infos
  void init_mesh_structure(void);
  
  //!Build mesh regions infos
  void parse_regions(void);


  //!Build element to atoms association map
  void build_elem_to_atoms(void);

  //! AtomisticStructureOptions object pointer
  AtomisticStructureOptions _atomistic_structure_options;

  //! Associate at any alement atoms contained
  std::map<const libMesh::Elem*, std::vector<unsigned int> > _elem_to_atoms;

  //!Associate elements: any atom keep tracks of the elements he belongs to
  void associate_elements();

  //! Read structure from file
  void read_structure(const std::string& path);

  //! Read tgn structure
  void read_xyz(const std::string& path, const Tensor1& transl);
  
  //! Read xyz structure
  void read_gen(const std::string& path, const Tensor1& transl);

  //! Read gen structure
  void read_tgn(const std::string& path, const Tensor1& transl);

  //! set the labels (used e.g. when reading xyz)
  void set_labels(void);

  //! Options for the structure (from Atomistic Region)
  ModelOptions _options;

  //! Name of the structure (will be the same of associated Atomistic region)
  std::string _name;

  //! Scale factor (from mesh_units to amstrong mesh_units/1e-10)
  double _scale;

  //! Set of mesh regions covered by atomistic structure (numbers: mesh regions)
  std::set <ID> _IDset;

  //! Tell if the object has been already initialized
  bool _is_initialized;

  //! True if the alloy must be treated as a Random Alloy
  /*!
   * (Species are used to distinguish the parent material the atom
   * belongs to). By default VCA is used
   */
  bool _random_alloy;

  //! Contains reference to device we're working with
  const Device* _device;

  //!Number of atoms excluding hydrogens
  unsigned int _N_without_H;

  //! The map with callbacks
  static std::map<std::string,
    std::list<boost::function<void(void)>>> _callback_functions;


  //! Manage structure printing
  void print_driver(void);

  //! Set virtual atom types in structure as for VCA InGaAs (InGa), As 
  std::vector<std::string> _virtual_atom_types;

  std::map<std::string, unsigned int> _virtual_type_idx;

  // calculate radial distribution function for a given specie
  void compute_g(const Specie&, double Rc, double dr, std::vector<std::map<Specie,unsigned int>>& g);

};

//----------------------------------------------------
// Inline member functions
//----------------------------------------------------


inline
ModelOptions& AtomisticStructure::get_options(void)
{
    assert( !(_options.is_empty()) );

    return _options;
}


inline
const ModelOptions& AtomisticStructure::get_options(void) const
{
    return _options;
}

inline
bool AtomisticStructure::is_random_alloy(void)
{
  return _random_alloy;
}


inline
void AtomisticStructure::set_options(const ModelOptions& options)
{
  _options = options;
}


inline
const std::string& AtomisticStructure::get_name(void) const
{
  return _name;
}


inline
const double& AtomisticStructure::get_scale(void) const
{
  return _scale;
}






inline
const Device* AtomisticStructure::get_device(void) const
{

  assert ( (_device != NULL) );

  return _device;

}

inline
void AtomisticStructure::set_device(const Device* const device)
{
 _device = device;
}


inline
const std::set<ID>&
AtomisticStructure::get_IDset(void) const
{
  return _IDset;
}


inline
const Specie&
AtomisticStructure::get_specie(unsigned int i) const
{
  return _atoms[i].get_specie();
}

inline
unsigned int
AtomisticStructure::get_N_without_H(void) const
{
  return _N_without_H;
}


//inline 
//const Material* 
//AtomisticStructure::get_reference_material(void) const
//{
//  return _reference_material;
//}

#endif // _ATOMISTICSTRUCTURE_H_
