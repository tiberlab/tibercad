// $Id$

#include "AtomisticStructure.h"
#include "AtomisticGenerator.h"
#include "BondMap.h"
#include "Messages.h"
#include "MeshUtils.h"
#include "Utils.h"
#include "TiberCad.h"
#include "Material.h"
#include "DataOutput.h"
#include "RuntimeException.h"

#include "libmesh/mesh_tetgen_interface.h"
// unfortunately tetgen defines REAL
#undef REAL

//C++ includes
//--------------------
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <tr1/random>
//---------------------

using namespace std;
using namespace libMesh;

map<string, list<boost::function<void(void)>>>
AtomisticStructure::_callback_functions;



AtomisticStructure::AtomisticStructure()
:_device(NULL),
 _atomistic_structure_options(),
 _scale(1.0),
 _random_alloy(false),
 _N_without_H(0)
{
}


AtomisticStructure::AtomisticStructure(const std::string& name)
:_name(name),
 _atomistic_structure_options(),
 _scale(1.0),
 _device(NULL),
 _random_alloy(false),
 _N_without_H(0)
{
}

AtomisticStructure::AtomisticStructure(const AtomisticStructure& other) :
  AtomisticBasis(other),
  _options(other._options),
  _name(other._name),
  _scale(other._scale),
  _IDset(other._IDset),
  _is_initialized(other._is_initialized),
  _random_alloy(other._random_alloy),
  _device(other._device),
  _N_without_H(other._N_without_H),
  _virtual_atom_types(other._virtual_atom_types),
  _virtual_type_idx(other._virtual_type_idx)
{
}


AtomisticStructure::~AtomisticStructure(void)
{
}


AtomisticStructure*
AtomisticStructure::create()
{
  AtomisticStructure* st = new AtomisticStructure();
  return st;
}

AtomisticStructure*
AtomisticStructure::create(const AtomisticStructure& as)
{
  AtomisticStructure* st = new AtomisticStructure(as);

  return st;
}



AtomisticStructure::AtomisticStructureOptions::AtomisticStructureOptions(void)
:is_associated(false)
{}


AtomisticStructure::AtomisticStructureOptions::~AtomisticStructureOptions(void)
{}


void
AtomisticStructure::init(const std::string& name,
    const Device* const device, const ModelOptions& options)
{
  std::string random_alloy;

  _name = name;
  _device = device;
  set_options(options);

  //Setting scale factor
  _scale = ( ( _device->get_mesh_units() ) / 1e-10 );

  //There's 3 ways of building a structure:
  // by material (bulk), by regions (mesh based)
  // or by file.
  // In the input section we can specify any, here we
  // define which of the three ways is given
  //
  // If file path is defined, structure is read from file.
  // Physical regions can be specified for a faster projection
  // between mesh and atoms. If it's not, all regions are considered
  // by default. If device is not defined, no regions are defined.
  //
  // If reference_region is defined and physical_regions is
  // not defined, a bulk description of material in the
  // reference region is done. AtomisticStructure build a copy of
  // Material, so options can be override
  //
  // If reference_region is defined and physical_regions is defined
  // the atomistic generator for complex structures is invoked
  //

  Messages m;
  m.newline();
  m.info("Building Atomistic Structure " + get_name());
  m.indent();

  random_alloy = _options.get_option("random_alloy", "false");
  if (random_alloy == "true") _random_alloy = true;

  _elem_to_atoms.clear();

  if (_options.find_option("load_structure") || _options.find_option("load"))
  {
    std::string filename;

    if (_options.find_option("load_structure")) filename = _options.get_option("load_structure","none");
    if (_options.find_option("load")) filename = _options.get_option("load","none");

    //------------------------------------------------------------
    ostringstream os;
    os << "Reading structure from file " << filename <<
        ". Any further information will be neglected" << std::endl;
    Messages::info(os.str(), true);
    os.str(std::string());
    //---------------------------------------------------------------
    read_structure(filename);

    Messages::info("Parse regions");
    parse_regions();

    //NOTE: I'm calculating the bond map again anyway because otherwise the translation vectors are not
    //correctly reproduced for periodic structures. If we really want to import the bond map,
    //this needs to be changed
    if (_bondmap != NULL)
    {
      delete _bondmap;
      _bondmap = NULL;
    }

    Messages::info("Build Bond Map");
    build_bond_map();

    //Messages::info("Associate elements");
    if (!_atomistic_structure_options.is_associated)
      associate_elements();

    set_labels();

    dorestrict(_IDset);

  }
  else if (_options.find_option("combine_structures"))
  {

    // combine two structures into a single one

    vector<string> pieces;
    _options.get_option("combine_structures", pieces);

    ostringstream os;
    os << "Combining atomistic structures : ";
    for (auto&& i : pieces)
      os << i << ", ";

    Messages::info(os.str());
    Messages::warning("Currently combined structures will share "
        "their atoms with the original structures!");

    const MeshBase& mesh = get_device()->get_mesh();

    // we need to recalculate the lattice vectors
    libMesh::RealVectorValue a;
    libMesh::RealVectorValue b;
    libMesh::RealVectorValue c;

    // the sum of sizes
    Point size(0);

    // the combined bbox
    Point bb_min(0);
    Point bb_max(0);

    bool first = true;
    for (auto&& i : pieces)
    {
      const AtomisticStructure* const as = _device->get_atomistic_structure(i);

      if (first)
      {
        *this = *as;
        _name = name;
        set_options(options);
        as->get_lattice_vectors(a, b, c);
        const auto& bbox = as->get_bounding_box();
        bb_min = bbox.first;
        bb_max = bbox.second;
        size = bb_max - bb_min;

        first = false;
      }
      else
      {

        _IDset.insert(as->_IDset.begin(), as->_IDset.end());
        _virtual_atom_types.insert(_virtual_atom_types.end(),
            as->_virtual_atom_types.begin(), as->_virtual_atom_types.end());

        _atoms.insert(_atoms.end(), as->_atoms.begin(), as->_atoms.end());

        libMesh::RealVectorValue ai, bi, ci;
        as->get_lattice_vectors(ai, bi, ci);
        switch (mesh.mesh_dimension())
        {
          case 3:
            c += ci;
          case 2:
            b += bi;
          case 1:
            a += ai;
          default:
            break;
        }

        const auto& bbox = as->get_bounding_box();
        size += bbox.second - bbox.first;
        for (unsigned int i = 0; i < 3; ++i)
        {
          if (bbox.first(i) < bb_min(i))
            bb_min(i) = bbox.first(i);

          if (bbox.second(i) > bb_max(i))
            bb_max(i) = bbox.second(i);
        }

      }
    }
    cerr << size << endl;
    cerr << bb_min << endl;
    cerr << bb_max << endl;

    set_lattice_vectors(a, b, c);

    refresh();

  }
  else
  {
    // Build mesh based representation
    Utils::Timer tt;
    tt.reset();

    parse_regions();
    init_mesh_structure();

    Messages::info("Atomistic structure build time: "+tt.elapsed_string());

  }
      
   //Calculate the number of atoms excluding hydrogens 
   //(Useful for passivated semiconductors)
   compute_N_without_H();

   m.newline();
   ostringstream os;
   os << "Atomistic Structure containing " << get_N_atoms() <<
         " atoms has been built. " <<std::endl;
   os << "Size not counting passivation hydrogens: "<< get_N_without_H()<<std::endl;
   m.info(os.str());
   os.str(std::string());
 
   m.info("Supercell structure");
   m.info("Lattice vectors (A):");
   m.indent();
 
   libMesh::RealVectorValue a, b, c;
   get_lattice_vectors(a, b, c);
   os << "a1 = (";
   a.write_unformatted(os, false);
   os << ")\na2 = (";
   b.write_unformatted(os, false);
   os << ")\na3 = (";
   c.write_unformatted(os, false);
   os << ")\n";
   m.info(os.str());
   os.str(std::string());


   Messages::info("Output structure(s)");
   print_driver();


   // Device call_back function 
   map<string, list<boost::function<void(void)>>>::iterator mit =
       _callback_functions.find(get_name());
 
   if (mit != _callback_functions.end())
   {
     list<boost::function<void(void)>>::iterator it((mit->second).begin());
     for ( ; it != (mit->second).end(); ++it)
     {
       (*it)();
     }
   }
}



void
AtomisticStructure::parse_lattice_vectors(void)
{
  if (_options.find_option("lattice_vectors"))
  {
    std::vector<double> lattice_vectors;
    _options.get_option("lattice_vectors", lattice_vectors);
 
    if (lattice_vectors.size() == 3)
    {  
       //set_periodicity( true, true, true);
       _lattice_vectors[0] = lattice_vectors[0];
       _lattice_vectors[4] = lattice_vectors[1];
       _lattice_vectors[8] = lattice_vectors[2];
    }
 
    else if (lattice_vectors.size() == 9)
    {
       for (int i = 0; i < 9; i++)
       {
         //set_periodicity( true, true, true);
         _lattice_vectors[i] = lattice_vectors[i]; 
       }
    }
 
    else Messages::error("lattice_vectors must have 3 or 9 components");
  }
}


void
AtomisticStructure::init(const std::string& )
{
}


void
AtomisticStructure::print_driver(void)
{
  if (_options.find_option("print"))
  {
    std::vector<std::string> extensions;
    _options.get_option("print", extensions);
  
    for (int i = 0; i < extensions.size(); i++)
    {

       std::string name(_name + "." + extensions[i]);
       print_structure(name);
    }
    
  }
       
  if (_options.has_submodel("radial_distribution"))
  {
     this->radial_distribution(); 
  }
       
  if (_options.has_submodel("alloy_statistics"))
  {
     const ModelOptions& opt = (_options.submodels_begin("alloy_statistics"))->second;
     this->extract_alloy_statistics(opt); 

     if (opt.get_option("plot_alloy_composition",false))
        this->plot_alloy_composition(opt); 
  }

}

void
AtomisticStructure::parse_regions(void)
{
  std::string physreg = _options.get_option("regions", "all");
  _device->extract_physical_regions(physreg, _IDset);
}


void
AtomisticStructure::init_mesh_structure()
{
  std::ostringstream os;
  assert(_device != NULL);

  //---------------------------------------------------------------
  os << "Mesh_units is " << _device->get_mesh_units() << std::endl
      << "Scale factor is " << _scale << std::endl;
  Messages::debug(os.str());
  os.str(std::string());
  //--------------------------------------------------------------

  if (!get_options().find_option("periodicity"))
  {
    switch (this->get_device()->get_mesh().mesh_dimension())
    {
      case 1:
        get_options().set_option("periodicity", 0);
        break;
      case 2:
        get_options().set_option("periodicity", "(0,0)");
        break;
     case 3:
        get_options().set_option("periodicity", "(0,0,0)");
        break;
    }
  }
  
  vector<bool> periodicity;
  get_options().get_option("periodicity",periodicity);

  //if (periodicity.size() != this->get_device()->get_mesh().mesh_dimension())
  //{
  //   throw RuntimeException("periodicity must have the same dimension of mesh");
  //}

  periodicity.resize(3,true);
  this->set_periodicity(periodicity);
  
  //---------------------------------------------------------------
  // Extend mesh for contacts 
  //unsigned int num_sides;
  //Point normal = Device::get_normal();
  //---------------------------------------------------------------

  AtomisticGenerator* generator 
    = AtomisticGenerator::create(this, _device->get_mesh().mesh_dimension());

  //build();
  //associate_elements(_as->get_IDset());
  //cut the structure (only flags atoms)
  //iterates on _super_basis and assign species
  //check_periodic(); Eventually enlarge along dummy supercell directions
  //passivate();
  //remove_atoms();
  //build_random_alloy();  
  generator->do_init();


  // Copy the structure back from the generator here
  generator->finalize();

  _elem_to_atoms.clear();
  for (size_t i = 0; i < _atoms.size(); ++i)
    _elem_to_atoms[_atoms[i].get_elem()].push_back(i);
 
  parse_lattice_vectors();

  Messages::info("Build final Bond Map...");

  // extract atom types and create final bondmap
  refresh();

  delete generator;

}


void
AtomisticStructure::dorestrict(const std::set<ID>& rgn_ids) 
{
  _IDset = rgn_ids;

  AtomisticGenerator* generator 
    = AtomisticGenerator::create(this, _device->get_mesh().mesh_dimension());

  generator->dorestrict();

  generator->finalize();

  Messages::info("Build final Bond Map...");
  build_bond_map();
  
  compute_N_without_H();

  delete generator;

  std::ostringstream os;
  os << "New structure size: "<< get_N_atoms();
  Messages::info(os.str());
  os.str("");
  os << "Without H: "<<_N_without_H;
  Messages::info(os.str()); 
}



void
AtomisticStructure::associate_elements()
{
  // the tensor grid to real mesh mapper for fast association atom->Elem
  // NOTE: we pass the relevant ID set, since Quantum contacts could be present
  //       which have to be included in the atomistic structure
  MeshUtils::GridMapper& mapper =
      MeshUtils::GridMapper::get_mapper(_device->get_mesh(), _IDset);

  std::vector<Atom>::iterator atom = _atoms.begin();

  unsigned int dim =  _device->get_mesh().mesh_dimension();
  _elem_to_atoms.clear();

  size_t n_atoms = _atoms.size();
  Utils::Progress prog("Assign elements", n_atoms);
  unsigned int progress = 0;

  // NOTE: Hydrogens remains outside regions and are not associated to elements
  for (size_t atom = 0; atom < n_atoms; ++atom)
  {
    
    Point p(_atoms[atom].get_position());
    p *= 1.0 / _scale;
    
    // set unneeded dim to 0, so atoms are associated to the correct elements
    switch ( dim )
    {
    case 0:
      p(0) = 0.0;
    case 1:
      p(1) = 0.0;
    case 2:
      p(2) = 0.0;
    default:
      break;
    }
        
    const Elem* elem = mapper.get_element(p);
    
    if (elem != NULL)
    {
      _atoms[atom].set_elem(elem);
      _elem_to_atoms[elem].push_back(atom);
    }

    progress++;
    prog.progress_message(progress);

  }
  //std::cout << std::endl;

  Messages::debug("Finished associate_elements");

}

void
AtomisticStructure::register_callback(string& name,
    boost::function<void(void)> callback)
{
  if (!name.empty())
  {
    vector<string> tokens;

    Utils::tokenize(name, tokens, ".");

    _callback_functions[tokens[0]].push_back(callback);
  }
}


const std::vector<unsigned int>&
AtomisticStructure::get_atoms_in_elem(const Elem* element) const
{
  static std::vector<unsigned int> empty_int_vector(0);

  std::map<const Elem*, std::vector<unsigned int>>::const_iterator
    it(_elem_to_atoms.find(element));
  if (it != _elem_to_atoms.end())
    return(it->second);

  return(empty_int_vector);
}




void
AtomisticStructure::set_virtual_types(const std::set<std::string>& atom_types)
{
  if (_virtual_atom_types.size()>0)  _virtual_atom_types.clear();

  std::set<std::string>::iterator type = atom_types.begin();
  for ( ; type != atom_types.end(); type++)
    _virtual_atom_types.push_back( *type );
  
  for(unsigned int i=0; i<_virtual_atom_types.size(); i++)
  {
    _virtual_type_idx.insert( std::make_pair( _virtual_atom_types[i], i+1) );
  }
}

unsigned int
AtomisticStructure::get_virtual_type_index(const std::string& type) const
{
  unsigned int result = 0;
  for (unsigned int i = 0; i < get_N_types(); i++){
      if ( (type.compare( _atom_types[i] ) == 0) ) result = i + 1;
  }

  return result;

}



// Virtual species are defined in the alloy database (InGa, AsP, ...)
void
AtomisticStructure::assign_virtual_species(void)
{
  std::set<ID> reg_ids(get_IDset());
  // For each region we have a map: label -> virtual_type
  // So map< region, map<label,virtual_type> >
  // 
  // Finally we need a map: label -> type_idx (atom_t)
  // This is done with an intermediate map:
  // virtual_type::iterator -> type_idx
  //
  std::set<std::string> atom_types;
  std::map<ID, std::map<Atom::atom_t, std::string > > assign;
  std::pair<std::set<std::string>::iterator, bool> ret;

  Messages::info("Assign Virtual Species");

  if (!this->is_random_alloy())
  {
    std::set<ID>::iterator reg(reg_ids.begin());

    for ( ; reg != reg_ids.end(); ++reg)
    {
      const Material* mat = get_device()->get_material( (*reg) );

      Database db = mat->get_database();
      // by setting no-mixing db.get() takes the first component
      db.set_alloy_mixing(Database::NONE);

      db.set_section("atomistic_structure");
      unsigned int n_species = db.get("n_basis_specie", 0);

      Atom::atom_t label = 0;

      //Build up conversion map from file
      for (unsigned int i = 1; i <= n_species; i++)
      {
        std::string record;
        std::string s;
        std::stringstream out;

        out << i;
        s = out.str();
        record = "specie_" + s;
        std::string specie = db.get(record.c_str(),"none");

        ret = atom_types.insert(specie);

        record = "n_" + s;
        unsigned int n = db.get(record.c_str(), 0);

        for (unsigned int j = 1; j <= n; j++)
        {
          assign[*reg][++label] = *(ret.first);
        }

      }

      //No more reading from section atomistic_structure in database are needed
      db.set_section("");

    }
  }
  else
  {
    atom_types.insert(_atom_types.begin(), _atom_types.end());
    atom_types.erase("H");
  }
 
  // Passivation Hydrogen has label 0 (we insert Hx by default)
  atom_types.insert("Hx");

  // Setup a map between the strings and label index
  // set _virtual_atom_types needed later 
  set_virtual_types(atom_types); 
   
  // -- DEBUG ---
  //std::cout<<"type set:"<<std::endl;
  //std::set<std::string>::iterator it=atom_types.begin();
  //for( ; it!=atom_types.end(); it++)
  //  std::cout<<*it<<" : "<<static_cast<unsigned int>(_virtual_type_idx[*it])<<std::endl;

  // Cycle upon all atoms and change specie according to assign map
  unsigned int count=1;
   
  std::vector<Atom>::iterator atom = _atoms.begin();
  for ( ; atom != _atoms.end(); ++atom)
  {
    //cout<< count<< " label: "<<static_cast<unsigned int>((*atom).get_label())<<endl;

    // For now passivation atoms have label == 0
    if ((*atom).get_label() == 0)
    {
      (*atom).set_type( _virtual_type_idx["Hx"] );
    }
    else
    {
      std::string str(atom->get_specie().get_string());

      if (!this->is_random_alloy())
      {
        ID el_reg = (*atom).get_elem()->subdomain_id();
        str = assign[el_reg][(*atom).get_label()];
      }

      (*atom).set_type( _virtual_type_idx[str] );
    }
     
    count++;
  }
 
  
}


// Virtual species are defined in the alloy database (InGa, AsP, ...)
void
AtomisticStructure::set_labels(void)
{
  for (unsigned int i = 0; i < get_N_atoms(); i++)
  {
    Atom& atom = _atoms[i];
    if ( atom.get_elem() == NULL)
      atom.set_label(0);
    else
    {
      const Material* mat = get_material(atom, false);
      unsigned int lb = mat->get_label(atom.get_specie());
      if (lb==255)
      {
        std::ostringstream os;
        os<<"Cannot assign a valid label to atom "+atom.get_specie().get_string()+" number "<<i+1<<std::endl;
        os<<"The loaded structure is incompatible with material region "<<mat->get_name();
        Messages::error(os.str());
        throw RuntimeException("");
      } 
      atom.set_label(lb);
    }
  }
}



void
AtomisticStructure::read_structure(const std::string& path)
{

  std::ifstream file;
  Tensor1 transl;

  // Delete eventually existing structure
  if (!(_atoms.empty())) _atoms.clear();
  if (!(_atom_types.empty())) _atom_types.clear();

  file.open(path.c_str(), std::ifstream::in);

  if (!file)
  {
    Messages::error("Unable to open file "+path+". Cannot read Atomistic Structure.");
  }

  // Recognize type of input file and read it properly
  std::string extension = path.substr(path.size()-4);
  
  //A translation vector can be specified to modify supercell alignment
  std::vector<double> translation (3,0.0);
  if ( get_options().find_option("translation") )
  {
   get_options().get_option("translation", translation);
   transl(1) += translation[0]; 
   transl(2) += translation[1]; 
   transl(3) += translation[2];
  }


  // XYZ file
  if ( (extension.compare(".xyz") == 0) || (extension.compare(".XYZ") == 0) )
  {
    read_xyz(path, transl);
  }
  else if ( (extension.compare(".gen") == 0) || (extension.compare(".GEN") == 0) )
  {
    read_gen(path, transl);
  }
  else if ( (extension.compare(".tgn") == 0) || (extension.compare(".TGN") == 0) )
  {
    read_tgn(path, transl);
  }
  else
  {
    Messages::error("Structure file extension is not recognized. ");
  }

  file.close();

}

void
AtomisticStructure::read_xyz(const std::string& path, const Tensor1& transl)
{
  std::ifstream file;
  file.open(path.c_str(), std::ifstream::in);

  std::string line, record;
  unsigned int n_specie;
  Atom tmp_atom;
  Tensor1 pos;
  // First line is number of atoms
  
  getline(file, line);
  unsigned int N_atoms = atoi(line.c_str());
  _atoms.reserve(N_atoms);

  if (N_atoms == 0)
  {
    std::cerr << "No atoms in structure files or non valid integer in first line \n";
    exit(1);
  }

  // Skip second line
  getline(file, line);

  // Start reading  lines
  // Note: Multiple initialization of line_string are allowed as inside a command block.
  // It's a suggested solution as it allows inner loop on record to work well
  // without further line_string manipulation
  while (getline(file, line))
  {
    std::stringstream line_string(line);
    // Extract atom type and check if it's a new type
    line_string >> record;

    if ( ~(_atom_types.empty()) )
      {
        bool not_present = true;

        for (unsigned int i = 0; i < ( _atom_types.size()); i++)
          {
            if ( record.compare(_atom_types[i]) == 0) not_present = false;
          }

        if (not_present) _atom_types.push_back(record);
      }
    else
      {
        _atom_types.push_back(record);
      }

    tmp_atom.set_specie( record );

    for (unsigned int i = 1; i < 4; i++)
      {

        line_string >> record;
        pos(i) = atof(record.c_str());
        
      }

    tmp_atom.set_position( pos+transl );

    _atoms.push_back(tmp_atom);

  }

  if ( (_atoms.size()) != N_atoms )
    std::cerr << "Warning: in file xyz number of atoms is wrong \n";

  // Warning: XYZ file has no informations about structure periodicity
  parse_lattice_vectors();

}  // end read_xyz

void
AtomisticStructure::read_gen(const std::string& path, const Tensor1& transl)
{
  std::ifstream file;
  file.open(path.c_str(), std::ifstream::in);

  std::string line, record;
  unsigned int n_specie;
  Atom tmp_atom;
  Tensor1 pos;

  getline(file, line);
  std::stringstream line_string(line);

  unsigned int N_atoms = atoi(line.c_str());
  _atoms.reserve(N_atoms);

  if (N_atoms == 0)
  {
    std::cerr << "No atoms in structure files or non valid integer in first line. \n";
    exit(1);
  }

  line_string >> record;

  if ( (record.compare("S") == 0) || (record.compare("s") == 0))
    set_periodicity( true, true, true);
  else  if ( (record.compare("C") == 0) || (record.compare("c") == 0))
    set_periodicity(false, false, false);
  else
    std::cerr << "Warning (in GEN file at first line): "
              << "Cluster (C) or Supercell (S) must be specified. "
              << " By default a Cluster (no periodicity) is considered. \n";

  getline(file, line);

  //  //This line clean stringstream in a safe way
  line_string.str(std::string());
  line_string.clear(std::stringstream::goodbit);
  //---------------------------------------------

  line_string << line;

  while ( line_string >> record)
  {
    _atom_types.push_back(record);
  }

  // Cycle upon specified number of atoms (last rows are for periodicity vectors)
  for (unsigned int i = 1; i <= N_atoms; i++)
  {
    getline(file, line);
    std::stringstream line_string(line);

    // First value is ignored (just atoms enumeration)
    line_string >> record;
    line_string >> record;
    n_specie = atoi(record.c_str());
    tmp_atom.set_specie( _atom_types[n_specie -1] );

    for (unsigned int j = 1; j <= 3; j++)
      {
        line_string >> record;
        pos(j) =  atof(record.c_str());
      }
    tmp_atom.set_position( pos+transl );
    _atoms.push_back(tmp_atom);
  }

  // An additional line is present in GEN files. It's the coordinates origin and it's
  // not needed
  getline(file, line);

  // Read periodicity vectors anyway, if system is not periodical they will be ignored
  unsigned int count = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    getline(file, line);
    std::stringstream line_string(line);
    for (unsigned int j = 0; j < 3; j++)
    {
      line_string >> record;
      _lattice_vectors[count] = atof(record.c_str());
      count++;
    }
  }
} // end read_gen


void
AtomisticStructure::read_tgn(const std::string& path, const Tensor1& transl)
{

  std::ifstream file;
  file.open(path.c_str(), std::ifstream::in);
  
  std::string line, record;
  unsigned int n_specie;
  Atom tmp_atom;
  Tensor1 pos;
  std::stringstream line_string;
  std::ostringstream os;

  Messages::debug("Reading tgn file. Loading atom coords and bond map ");

  getline(file, line);

  line_string.clear();
  line_string.str(line);

  line_string >> record;
 

  unsigned int N_atoms = atoi(record.c_str());
  _atoms.reserve(N_atoms);

  //Prepare bond map object
  if ( _bondmap == NULL) _bondmap = new BondMap(N_atoms);
  else
  {
      delete _bondmap;
      _bondmap = new BondMap(N_atoms);
  }

  BondMap& bondmap = *_bondmap;

  if (N_atoms == 0)
    {
      std::cerr << "No atoms in structure files or non valid integer in first line. \n";
      exit(1);
    }
  _atoms.resize(N_atoms);

  line_string >> record;

  if ( (record.compare("S") == 0) || (record.compare("s") == 0))
    set_periodicity( true, true, true);
  else  if ( (record.compare("C") == 0) || (record.compare("c") == 0))
    set_periodicity(false, false, false);
  else
  {
    std::cerr << "Warning (in GEN file at first line): " 
              << "Cluster (C) or Supercell (S) must be specified. " 
              << "By default a Cluster (no periodicity) is considered. \n"<<std::endl;
  }

  getline(file, line);

  //This 2 lines clean stringstream in a safe way
  line_string.clear();
  line_string.str(line);


  while ( line_string >> record)
    {
      _atom_types.push_back(record);
    }


  // Cycle upon specified number of atoms (last rows are for periodicity vectors)
  for (unsigned int i = 1; i <= N_atoms; i++)
    {
      getline(file, line);

      line_string.clear();
      line_string.str(line);

      // First value contains atom label 
      line_string >> record;
      _atoms[i - 1].set_label(atoi(record.c_str()));

      line_string >> record;
      n_specie = atoi(record.c_str());

      for (unsigned int j = 1; j <= 3; j++)
        {
          line_string >> record;
          pos(j) =  atof(record.c_str());
        }
      _atoms[i - 1].set_specie(_atom_types[n_specie -1]);
      _atoms[i - 1].set_position(pos+transl);

      //Get bond map
      line_string >> record;
      //TODO: write a set_bond_map function.
      // Change pointers in vectors to manage constness more easily
      bondmap[i-1].resize( atoi(record.c_str()) );

      for (unsigned int j = 0; j < bondmap[i - 1].size(); j++)
        {
          line_string >> record;
          bondmap[i - 1][j] = atoi(record.c_str()) - 1;
        }

      line_string >> record;
      int tmp_id = atoi(record.c_str());

      if (tmp_id == -1) _atoms[i - 1].set_elem(NULL);
      else
        {
          _atoms[i - 1].set_elem(_device->get_mesh().elem(tmp_id));
        }


    }

  // An additional line is present in GEN files. It's the coordinates origin and it's
  // not needed
  getline(file, line);

  // Read periodicity vectors anyway, if system is not periodical they will be ignored
  unsigned int count = 0;
  for (unsigned int i = 0; i < 3; i++)
    {
      getline(file, line);
      line_string.clear();
      line_string.str(line);

      for (unsigned int j = 0; j < 3; j++)
        {
          line_string >> record;
          _lattice_vectors[count] = atof(record.c_str());
          count++;
        }
    }

  file.close();
  _atomistic_structure_options.is_associated = true;

}


void
AtomisticStructure::print_tgn(const std::string& path) const
{
  std::ofstream file;

  const BondMap& bondmap = *_bondmap;
  // -------------------------------------------

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/" + path;

  //I must build a materials map
  std::map<const Material*, unsigned int> material_map;

  std::set<ID>::iterator ID_it;

  // we need to know which regions have been used actually, otherwise
  // we will write too many materials at the end.
  unsigned int id = 1;
  const Material* mat = NULL;

  for (ID_it = _IDset.begin(); ID_it != _IDset.end(); ++ID_it)
    {
      mat = _device->get_material(*ID_it);
      material_map.insert(std::pair<const Material*, unsigned int>(mat, id));
      id++;
    }

  //Standard gen section (modified with material index)
  // --------------------------------------------
  file.open(file_name.c_str());
      
      file << _atoms.size();
      if (is_periodic()) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(6) << _atom_types[i];
        }
      file << std::endl;
      
      for (unsigned int i = 0; i < _atoms.size(); i++)
        {
          unsigned int n_specie;
          for (n_specie = 0; n_specie < _atom_types.size(); n_specie++)
            {
              if ( _atom_types[n_specie] == _atoms[i].get_specie() ) break;
            }
          
          // Print atom label 
          file << std::setw(3);
          file<< static_cast<unsigned int>(_atoms[i].get_label());

          file << std::setw(5) << n_specie + 1
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2)) ;

          if (_bondmap != NULL)
            {

              file << std::setw(5) << bondmap[i].size();

              // N.B. Indexing is in Fortran notation (first atom is labelled as 1) !!!!!!!!!!!!!!!!
              for (unsigned int j = 0; j < bondmap[i].size(); j++)
                {
                  file << std::setw(10) << bondmap[i][j] + 1;
                }
              ///////////////////////////////////////////

            }

          //ID of element is saved (note: no modifications to mesh are allowed to preserve compatibility)
          file << std::setw(14);
          if (_atoms[i].get_elem() == NULL || _atoms[i].get_label() == 0) file << -1;
          else file << _atoms[i].get_elem()->id();

          //Write Material ID (useful for some external codes...)
          file << std::setw(3);
          if (_atoms[i].get_elem() ==  NULL)
          {
            file << material_map[_device->get_material(_atoms[get_bond_map()[i][0]].get_region_ID()) ];
          }
          else
          {
            file << material_map[ (_device->get_material(_atoms[i].get_region_ID())) ];
          }
         

          file << std::endl;

        }



      //A line of zeros is put here (coordinates origin)
      file <<  std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) << "\n";

      // Periodicity vectors at the bottom
      unsigned int count = 0;
      for (unsigned int i = 0; i < 3; i++)
        {
          for (unsigned int j = 0; j < 3; j++)
            { 
              file << std::setw(20) << std::setprecision(10) << std::fixed <<
                  _lattice_vectors[count];
              count++;
            }
          file << "\n";
        }
      file.close();
}

void
AtomisticStructure::print_structure(const std::string& path)
{
  std::ofstream file;


  // -------------------------------- -----------

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/" + path;

  // --------------------------------------------



  //#ifdef DEBUG 
  //  std::cerr << "AtomisticStructure::print_structure(path) begin. \n";
  //#endif


  // Recognize type of input file and print it properly
  std::string extension = path.substr(path.size()-4);


  if ( (extension.compare(".xyz") == 0) || (extension.compare(".XYZ") == 0) )
    {
     print_xyz(path);
    }

  else if ( (extension.compare(".xyb") == 0) || (extension.compare(".XYB") == 0) )
    {
     print_xyb(path);
    }

  else if ( (extension.compare(".tgn") == 0) || (extension.compare(".TGN") == 0) )
    {
     print_tgn(path);
    }


  else if ( (extension.compare(".gen") == 0) || (extension.compare(".GEN") == 0) )
    {
     print_gen(path);
    }

  else if ( (extension.compare(".upg") == 0) || (extension.compare(".UPG") == 0) )
    {
      print_upg(file_name, "", true);
    }

  else
    {
      std::cerr << "File extension " << extension <<
          " of file does not correspond to any internal format. File not print. \n";
    }

  file.close();

  //#ifdef DEBUG
  //  std::cerr << "AtomisticStructure::print_structure(path) end. \n";
  //#endif

}



/*
* Quaternary Type A
*
* Al(x)In(y)Ga(1-x-y)P = x AlP + y InP + (1-x-y) GaP
* [Al(b)Ga(1-b)P](a) + [In(b)Ga(1-b)P](1-a)
* ab AlP + a(1-b) GaP + (1-a)b InP + (1-a)(1-b) GaP
* ab AlP + (1-a)b InP + (1-b) GaP
*
* x = ab    
* y = b - x  
*
* Quaternary Type B
*
* Al(x)Ga(1-x)As(y)P(1-y) = xy AlAs + (1-x)(1-y) GaP + x(1-y) AlP + y(1-x) GaAs
* [Al(b)Ga(1-b)As](a) + [Al(b)Ga(1-b)P](1-a)
* ab AlAs + a(1-b) GaAs + (1-a)b AlP + (1-a)(1-b) GaP
*
* a=y 
* b=x
*
*/

void
AtomisticStructure::print_upg(const std::string& path, const std::string& etb_dataset,
                                                       bool band_offsets)
{
  // assign virtual species (essentially setup _virtual_atom_types)

  assign_virtual_species();

  std::ofstream file, os;

  if (_bondmap == NULL)
    Messages::error("print upg file requires a valid bondmap");

  const BondMap& bondmap = *_bondmap;

  Messages::info("Printing upg file");

  // Recognize type of input file and print it properly
  std::string extension = path.substr(path.size()-4);
  //Gen format modified for uptight input

  file.open(path.c_str());
  
  //Material -> index 
  std::map<const Material*, unsigned int> material_map;
  
  // Create Material map--------------------------------------
  std::set<ID>::iterator ID_it;
  
  unsigned int id = 1;
  const Material* mat = NULL;
  
  for (ID_it = _IDset.begin(); ID_it != _IDset.end(); ++ID_it)
  {
    mat = _device->get_material(*ID_it);
    material_map.insert(std::pair<const Material*, unsigned int>(mat, id));
    id++;
  }
  //-----------------------------------------------------------

  std::set<Utils::Couple<const Material* > > interfaces;

  //Standard gen section (modified with material index)
  file << _atoms.size();
  
  if (is_periodic()) file << std::setw(10) << "S \n";
  else file << std::setw(10) << "C \n";
  
  for (unsigned int i = 0; i < _virtual_atom_types.size(); i++)
  {
    file << std::setw(7) << _virtual_atom_types[i];
  }
  file << std::endl;
  
  for (unsigned int i = 0; i < _atoms.size(); i++)
  {
    const Material *mati = NULL;
    const Material *matj = NULL;

    file << std::setw(10);
    if (_atoms[i].get_elem() ==  NULL)
    {
      mati = _device->get_material(_atoms[bondmap[i][0]].get_region_ID()); 
      file << material_map[mati];
    }
    else
    {
      mati = _device->get_material(_atoms[i].get_region_ID()); 
      file << material_map[ mati ];
    }

    file << std::setw(5) << static_cast<unsigned int>(_atoms[i].get_type());

    file << std::setw(20) << std::setprecision(10)
    << std::fixed << double(_atoms[i].get_position(0))
    << std::setw(20) << std::setprecision(10)
    << std::fixed  << double(_atoms[i].get_position(1))
    << std::setw(20) << std::setprecision(10)
    << std::fixed  << double(_atoms[i].get_position(2));

    
    file << std::setw(5) << bondmap[i].size();
    
    // N.B. Indexing is in Fortran notation (first atom is labelled as 1) !
    for (unsigned int j = 0; j < bondmap[i].size(); j++)
    {
      file << std::setw(10) << bondmap[i][j] + 1;
   
      if(_atoms[bondmap[i][j]].get_elem() != NULL)
      {
        matj = _device->get_material(_atoms[bondmap[i][j]].get_region_ID()); 
        if (mati != matj)
          interfaces.insert( Utils::Couple<const Material*>(mati,matj) );
      }
    }
    
    file << std::endl;
    
  }
  

  // Periodicity vectors at the bottom
  if (is_periodic())
  {
    
    //A line of zeros is put here (coordinates origin)
    file <<  std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
         << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
         << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) << "\n";
    
    unsigned int count = 0;
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
      {
        file << std::setw(20) << std::setprecision(10) << std::fixed <<
          _lattice_vectors[count];
        count++;
      }
      file << std::endl;
    }
  }
  
  //Information about materials
  
  file << "#Materials " << std::endl;

  std::map<const Material*, unsigned int> n_species_map;

  // to get them ordered correctly
  std::map<unsigned int, const Material*> inv_material_map;
  {
    std::map<const Material*, unsigned int>::iterator mat_it = material_map.begin();
    for ( ; mat_it != material_map.end(); ++mat_it)
    {
      inv_material_map[mat_it->second] = mat_it->first;
    }
  }

  std::map<unsigned int, const Material*>::iterator mat_it = inv_material_map.begin();
  for ( ; mat_it != inv_material_map.end(); ++mat_it)
  {
    const Material* mat = (*mat_it).second;
    
    Database db = mat->get_database();
    db.set_section("atomistic_structure");
    db.set_alloy_mixing(Database::NONE);
    n_species_map[mat] =  db.get("n_basis_specie", 0);

    std::string alloy_type("");
    
    if      ( n_species_map[mat] == 1 ){ alloy_type = "simple";}
    else if ( n_species_map[mat] == 2 ){ alloy_type = "binary";}
    else if ( n_species_map[mat] == 3 ){ alloy_type = "ternary";}

    if (mat->is_alloy())
    {
      db.set_section("");
      alloy_type = db.get("alloy_type", "ternary");
    } 
    if (alloy_type == "")
      Messages::error("Could not define alloy_type in AtomisticStructure.C");

    std::vector<double> Ev(2,0.0);
    
    if (mat->is_alloy())
    {

      file << (*mat_it).first
           << " " << mat->get_name()
           << " " << mat->get_structure()
           << " " << alloy_type << " ";

      if (alloy_type == "binary")
        file << 2;
      else if (alloy_type == "ternary")
        file << 2;
      else if (alloy_type == "quaternary")
        file << 4;

      if (is_random_alloy())
        file << " RND";
      else
        file << " VCA";

      std::string nameA = (static_cast<const Alloy*>(mat))->get_name_A();
      std::string nameB = (static_cast<const Alloy*>(mat))->get_name_B();
      const Material* matA = (static_cast<const Alloy*>(mat))->get_component_A();	  
      const Material* matB = (static_cast<const Alloy*>(mat))->get_component_B();	
      double x =  (static_cast<const Alloy*>(mat))->get_molar_fraction();
      //  double x = mat->get_options().get_option("x",1.0);

      if (alloy_type == "binary" || alloy_type == "ternary")
      {
        file << " "  << nameA
             << " "  << nameB
             << " "  << std::setprecision(3) << x
             << " "  << std::setprecision(3) << 1.0 - x
             << " "  << nameA << etb_dataset + ".etb"
             << " "  << nameB << etb_dataset + ".etb";
      }

      if (alloy_type=="binary")
      {
        file << " " << nameA << nameB << etb_dataset + ".etb";
      }
    
      if (alloy_type == "quaternary")
      {
        db.set_section("");
        double x_A = db.get("x_A",1.0);
        double x_B = db.get("x_B",1.0);

        file << " "  << (static_cast<const Alloy*>(matA))->get_name_A()
             << " "  << (static_cast<const Alloy*>(matA))->get_name_B()
             << " "  << (static_cast<const Alloy*>(matB))->get_name_A()
             << " "  << (static_cast<const Alloy*>(matB))->get_name_B()
             << " "  <<  std::setprecision(3) <<  x*x_A
             << " "  <<  std::setprecision(3) <<  x*(1-x_A)
             << " "  <<  std::setprecision(3) <<  (1-x)*x_B
             << " "  <<  std::setprecision(3) <<  (1-x)*(1-x_B)
             << " "  << (static_cast<const Alloy*>(matA))->get_name_A() << etb_dataset + ".etb"
             << " "  << (static_cast<const Alloy*>(matA))->get_name_B() << etb_dataset + ".etb"
             << " "  << (static_cast<const Alloy*>(matB))->get_name_A() << etb_dataset + ".etb"
             << " "  << (static_cast<const Alloy*>(matB))->get_name_B() << etb_dataset + ".etb";
      }

      if (band_offsets)
      { 
        const Database& dbA = matA->get_database();
        const Database& dbB = matB->get_database();
        dbA.set_section("valenceband");
        dbB.set_section("valenceband");
    
        file           << " " << dbA.get("E_v",0.0)
                       << " " << dbB.get("E_v",0.0);
      }

      if (alloy_type == "quaternary")
        file                          << " 0.0  0.0";

      file                            << " 0.0  0.0"
                                      << " 0.0  0.0" << std::endl;
    }
    else
    {
      db.set_section("valenceband");
      file << (*mat_it).first
           << " " << mat->get_name()
           << " " << mat->get_structure()
           << " " << alloy_type
           << " " << 1
           << " " << "CRY"
           << " " << mat->get_name()
           << " " << std::setprecision(3) << 1.0
           << " " << mat->get_name() << etb_dataset + ".etb";
      if (band_offsets)
      {               
        file << " " << db.get("E_v",0.0);
      } 
      file << " 0.0  0.0"   << std::endl;
    }

  }
 
  //std::cout<<"Print interfaces: "<<interfaces.size()<<std::endl;
  

  file << "#Interfaces " << std::endl;

  //std::set<Utils::Couple<const Material*>>::iterator int_it = interfaces.begin();
  //for (;int_it != interfaces.end(); ++int_it)

  std::multimap<Utils::Couple<ID>, Utils::Couple<Specie>> processed;


  for (size_t i = 0; i < get_N_atoms(); ++i)
  {
    //const Material* mat1 = (*int_it).first;
    //const Material* mat2 = (*int_it).second;
    ID reg_i = _atoms[i].get_region_ID();

    if (reg_i == INVALID_ID)
      continue;

    for (unsigned int j = 0; j < bondmap[i].size(); ++j)
    {
      ID reg_j = _atoms[bondmap[i][j]].get_region_ID();

      if ((reg_j == INVALID_ID) ||
          (reg_i == reg_j))
        continue;

      Utils::Couple<ID> couple(reg_i, reg_j);
      Utils::Couple<Specie> sp_couple(_atoms[i].get_specie(),
          _atoms[bondmap[i][j]].get_specie());
      if(processed.count(couple))
      {
        auto range = processed.equal_range(couple);
        //auto it = find(range.first, range.second,
        //    make_pair(couple, sp_couple));
        bool found = false;
        for (auto& it = range.first; it != range.second; ++it)
        {
          if (it->second == sp_couple)
            found = true;
        }

        if (found)
          continue;
      }

      processed.insert(make_pair(couple, sp_couple));

      const Material* mat1 = _device->get_material(reg_i);
      const Material* mat2 = _device->get_material(reg_j);

      // Put interfaces in increasing order
      if (material_map[mat1] > material_map[mat2])
      {
        swap(mat1, mat2);
      }

      file << std::setw(3) << material_map[mat1]<<"  "<< material_map[mat2]<<"  ";

      if (mat1->is_alloy() || mat2->is_alloy())
      {
        if (is_random_alloy())
          file <<  " RND";
        else
          file <<  " VCA";
      }
      else
      {
        file <<  " CRY";
      }

      std::vector<std::string> str;
      std::vector<double> frac;
      str.reserve(8); frac.reserve(8);

      interface_interactions(mat1, mat2, str, frac,
          _atoms[i], _atoms[bondmap[i][j]]);

      unsigned int size1 = str.size();
      file << " " << size1 << "  ";

      interface_interactions(mat2, mat1, str, frac,
          _atoms[bondmap[i][j]], _atoms[i]);

      file << " " << str.size()-size1 << "  ";

      for (unsigned int i=0; i < str.size(); i++)
        file <<"  " << str[i];

      for (unsigned int i=0; i < str.size(); i++)
        file <<"  " << frac[i];

      for (unsigned int i=0; i<str.size(); i++)
        file <<"  " << str[i]<<etb_dataset<<".etb";
      //if (check_data_file(Database::get_search_path()+"/"+str[i]+etb_dataset+".etb") )
      //    file <<"  " << str[i]<<etb_dataset<<".etb";
      //else
      //  Messages::error("Database file does not exist"+Database::get_search_path()+"/"+str[i]+etb_dataset+".etb");

      file<<std::endl;
    }
  }

  file.close();


  Messages::debug("upg file printed");

}
 
void
AtomisticStructure::interface_interactions(const Material* mat1, 
                                           const Material* mat2, 
                                           std::set<std::string>& str)
{
  str.clear();
  stringstream ss;

     
  {
    // sets all possible cation then all anions, e.g. AlGa-AsP = AlAs, AlP, GaAs, GaP
    // If a material has no anions (e.g. Si, Ge) species are inverted (not that robust!)
    // e.g.: SiGe-Si = SiSi SiGe (GeSi inverted)
    //       SiGe-SiGe = SiSi SiGe 
    Material::crystal_species_iterator sp_it1 = mat1->species_begin(1);
    Material::crystal_species_iterator sp_end1 = mat1->species_end(1);
    for (;sp_it1 != sp_end1; ++sp_it1)
    {
      Material::crystal_species_iterator sp_it2 = mat2->species_begin(2);
      Material::crystal_species_iterator sp_end2 = mat2->species_end(2);
      if ( sp_it2 != sp_end2 )
      {
        for (;sp_it2 != sp_end2; ++sp_it2)
        {
          ss << (*sp_it1)<<(*sp_it2);
          str.insert(ss.str());
          ss.str("");
        }
      }
      else
      {
        sp_it2 = mat2->species_begin(1);
        ss << (*sp_it2)<<(*sp_it1);
        str.insert(ss.str());
        ss.str("");
      }
    }
  }
  
  {
    Material::crystal_species_iterator sp_it1 = mat2->species_begin(1);
    Material::crystal_species_iterator sp_end1 = mat2->species_end(1);
    for (;sp_it1 != sp_end1; ++sp_it1)
    {
      Material::crystal_species_iterator sp_it2 = mat1->species_begin(2);       
      Material::crystal_species_iterator sp_end2 = mat1->species_end(2);
      if ( sp_it2 != sp_end2 )
      {
        for (;sp_it2 != sp_end2; ++sp_it2)
        {
          ss << (*sp_it1)<<(*sp_it2);
          str.insert(ss.str());
          ss.str("");
        }
      }
      else
      {
        sp_it2 = mat1->species_begin(1);
        ss << (*sp_it2)<<(*sp_it1);
        str.insert(ss.str());
        ss.str("");
      }
    }
  }

}

void
AtomisticStructure::interface_interactions(const Material* mat1, 
                                           const Material* mat2, 
                                           std::vector<std::string>& str,
                                           std::vector<double>& frac,
                                           const Atom& at1,
                                           const Atom& at2)
{

  str.clear();
  stringstream ss;
  double conc1, conc2;

  Specie sp1(at1.get_specie());
  Specie sp2(at2.get_specie());
  
  // Number of interctions (ca1)*(an2)
  // sets all possible cation then all anions, e.g. AlGa-AsP = AlAs, AlP, GaAs, GaP
  // repetitions are possible!
  // e.g.: SiGe-SiGe = SiSi SiGe GeSi GeGe
  unsigned int label1 = at1.get_label();

  Material::crystal_species_iterator sp_it1 = mat1->species_begin(label1);
  Material::crystal_species_iterator sp_end1 = mat1->species_end(label1);
  for (;sp_it1 != sp_end1; ++sp_it1)
  {
    if (*sp_it1 == sp1)
    {
      conc1 = 1.0;
      if (mat1->is_alloy())
        conc1=(static_cast<const Alloy*>(mat1))->get_molar_fraction(label1, sp1);

      unsigned int label2 = at2.get_label();

      Material::crystal_species_iterator sp_it2 = mat2->species_begin(label2);
      Material::crystal_species_iterator sp_end2 = mat2->species_end(label2);


      for (;sp_it2 != sp_end2; ++sp_it2)
      {
        if (*sp_it2 == sp2)
        {

          if (label2 < label1)
            swap(sp1, sp2);

          ss << (sp1)<<(sp2);
          str.push_back(ss.str());
          ss.str("");

          conc2 = 1.0;
          if (mat2->is_alloy())
            conc2=(static_cast<const Alloy*>(mat2))->get_molar_fraction(label2, sp2);

          frac.push_back(conc1*conc2);
        }
      }
    }


  }

 
}


void
AtomisticStructure::print_structure(const std::string& path, double const* const charges)
{
  std::ofstream file;

  //#ifdef DEBUG
  //  std::cerr << "AtomisticStructure::print_structure(path) begin. \n";
  //#endif

  file.open(path.c_str());

  // Recognize type of input file and print it properly
  std::string extension = path.substr(path.size()-4);

  if ( (extension.compare(".xyz") == 0) || (extension.compare(".XYZ") == 0) )
    {
      file << _atoms.size() << std::endl << std::endl;

      for (unsigned int i = 0; i < _atoms.size(); i++)
        {
          file << std::setw(2) << _atoms[i].get_specie()
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(charges[i]) << "\n";
        }

    }

  else if ( (extension.compare(".gen") == 0) || (extension.compare(".GEN") == 0) )
    {
      file << _atoms.size();

      if (is_periodic()) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(4) << _atom_types[i];
        }
      file << std::endl;

      for (unsigned int i = 0; i < _atoms.size(); i++)
        {
          unsigned int n_specie;
          for (n_specie = 0; n_specie < _atom_types.size(); n_specie++)
            {
              if ( _atom_types[n_specie] == _atoms[i].get_specie() ) break;
            }
          file << std::setw(10) << i + 1 << std::setw(5) << n_specie + 1
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(charges[i]) << "\n";
        }

      //A line of zeros is put here (coordinates origin)
      file <<  std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) << "\n";

      // Periodicity vectors at the bottom
      if (is_periodic())
        {
          unsigned int count = 0;
          for (unsigned int i = 0; i < 3; i++)
            {
              for (unsigned int j = 0; j < 3; j++)
                {
                  file << std::setw(20) << std::setprecision(10) << std::fixed <<
                      _lattice_vectors[count];
                  count++;
                }
              file << "\n";
            }
        }

    }

  else
    {
      std::cerr << "File extension does not correspond to any internal format. File not print. \n";
    }

  file.close();

  //#ifdef DEBUG
  //  std::cerr << "AtomisticStructure::print_structure(path) end. \n";
  //#endif

}



void
AtomisticStructure::compute_N_without_H(void)
{
  unsigned int N = 0;
  const unsigned int size = _atoms.size();
  for (unsigned int i = 0; i < size; i++)
  {
    if (_atoms[i].get_label() != 0)
    {
      N++;
    }
  }
  _N_without_H = N;
}
 
unsigned int
AtomisticStructure::compute_N_cations(void) const
{
  unsigned int N = 0;
  for (unsigned int i = 0; i < _atoms.size(); i++)
  {
    if (_atoms[i].is_cation()) N++;
  }
 
  return N;
}



//! Get atom Material
const Material*
AtomisticStructure::get_material(const Atom& atom, bool parent) const
{
 const Material* mat = get_device()->get_material(atom.get_region_ID());

 if (parent && mat->is_alloy())
   {
     const Material* matA = (static_cast<const Alloy*>(mat))->get_component_A();	  
     const Material* matB = (static_cast<const Alloy*>(mat))->get_component_B();	
     if (matA->has_specie(atom.get_specie()) &&
         (!matB->has_specie(atom.get_specie())))
       mat = matA;
     else if (matB->has_specie(atom.get_specie()) &&
         (!matA->has_specie(atom.get_specie())))       
       mat = matB;
     else
     {
       Messages::error("Ambiguity for alloy component assignation"
           " in AtomisticStructure::get_material(Atom&, bool)");
     }
   }

 return mat;

}

const Material*
AtomisticStructure::get_material(const Atom& atom1, const Atom& atom2,
    bool ) const
{
  const Material* mat = nullptr;

  // first look in the materials
  const Material* mat1 = get_device()->get_material(atom1.get_region_ID());
  const Material* mat2 = get_device()->get_material(atom2.get_region_ID());

  if (mat1->has_specie(atom1.get_specie()) && mat1->has_specie(atom2.get_specie()))
  {
    // can take it from here
    mat = mat1;
  }

  if (mat == nullptr)
  {
    // try here
    if (mat2->has_specie(atom1.get_specie()) && mat2->has_specie(atom2.get_specie()))
    {
      // can take it from here
      mat = mat2;
    }
  }

  // if mat exists and is an alloy, go into components
  if ((mat != nullptr) && (mat->is_alloy()))
  {
    const Alloy* alloy = static_cast<const Alloy*>(mat);
    mat = alloy->get_parent(make_pair(atom1.get_specie(), atom2.get_specie()));
  }

  // if we didn't find it, we check interface materials (which is created if needed)

 //If not, we need to decide based on some other criteria. Up to now we're able to
 //decide only for III-V or II-VI alloys with different cations (eg. Ga-As belong to GaAs)
 //if (mat1->is_cation(atom1.get_specie()))
 //  return get_material(atom1, parent);
 //else if (mat2->is_cation(atom2.get_specie()))
 //  return get_material(atom2, parent);
 //else
 //(Alex) Test using the label on atoms rather than hardcoded CrystalDefs. 
  /*
 if (atom1.is_cation())
   return get_material(atom1, parent);
 else if (atom2.is_cation())
   return get_material(atom2, parent);
 else 
 {
   //If no value was already returned, throw an exception
   Messages::error("Material for couple of atoms is decided "
       "depending on the cation species. I cannot find a valid cation ");
   throw RuntimeException("Cannot find valid cation");
 }
 */
 return(mat);
}

void
AtomisticStructure::reorder(const std::vector<unsigned int>& P)
{
  if (P.size() != _atoms.size()) 
    Messages::error("Reorder size is wrong");


  // Inverse Permutation vector for reordering bondmap
  std::vector<int> invP(_atoms.size());
  for(unsigned int i=0; i<_atoms.size(); i++)
  {
    invP[P[i]]=i;
  }

  // Define a new vector of atoms and Bondamp
  std::vector<Atom> new_atoms(_atoms.size()); 

  BondMap new_bondmap(_atoms.size());
  BondMap& bondmap = *_bondmap;

  for (unsigned int i=0; i< _atoms.size(); i++)
  {
    new_atoms[i] = _atoms[P[i]];
    //std::cout<<i<<" "<<new_atoms[i].get_position(0)<<" "<<_atoms[P[i]].get_position(0)
    //            <<" | "<<new_atoms[i].get_elem()<<" "<<_atoms[P[i]].get_elem()
    //            <<" | "<<new_atoms[i].get_specie()<<" "<<_atoms[P[i]].get_specie()<<"  ";

    new_bondmap[i].resize( bondmap[P[i]].size() );

    for (unsigned int j = 0; j < new_bondmap[i].size(); j++)
    {
       new_bondmap[i][j] = invP[bondmap[P[i]][j]];
       //std::cout<<" "<<new_bondmap[i][j];
    }
    //std::cout<<std::endl;
  }

  // Copy the old atoms on the new atoms
  for (unsigned int i=0; i< _atoms.size(); i++)
  {
    _atoms[i] = new_atoms[i]; 
    
    bondmap[i].resize(new_bondmap[i].size());

    for (unsigned int j = 0; j < new_bondmap[i].size(); j++)
    {
      bondmap[i][j] = new_bondmap[i][j];
    }
  }

}


void
AtomisticStructure::create_conformal_grid(libMesh::UnstructuredMesh& mesh,
    set<ID> labels, bool keep_node_order) const
{
#ifdef LIBMESH_HAVE_TETGEN
  Messages::info("Creating FEM mesh conforming with the atomistic structure."
      " This may take some time... ", false);
  mesh.set_mesh_dimension(3);
  const std::vector<Atom>& structure = get_structure_atoms();

  mesh.reserve_nodes(structure.size());
  //mesh.reserve_elem(_structure_atoms.size());

  const BondMap& bm = *_bondmap;

  map<ID, ID> nodes;
  unsigned int ctr = 0;
  unsigned int elem_ctr = 0;

  // WARNING this works only for tetrahedric cells
  for (unsigned int i = 0; i < structure.size(); i++)
  {
    const Atom& atom = get_structure_atom(i);

    if ((atom.get_specie() != Specie::H) &&
        (labels.empty() || labels.count(atom.get_label())))
    {
      Point p(atom.get_position());
      p *= 0.1;
      mesh.add_point(p);
      /*
      const vector<unsigned int>& neigh = bm[i];
      if (neigh.size() == 4)
      {
        bool gonext = false;
        for (unsigned int j = 0; j < neigh.size(); ++j)
        {
          if(get_structure_atom(neigh[j]).get_specie() ==  Specie::H)
            gonext = true;
        }
        if (gonext)
          continue;

        Elem* elem = Elem::build(TET4).release();
        elem->set_id(elem_ctr);
        mesh.add_elem(elem);

        vector<unsigned int> local_nodes(4);
        for (unsigned int j = 0; j < neigh.size(); ++j)
        {
          map<ID,ID>::iterator it(nodes.find(neigh[j]));
          if (it == nodes.end())
          {
            it = (nodes.insert(make_pair(neigh[j], ctr))).first;

            const Atom& atomj = get_structure_atom(neigh[j]);
            Node* node = mesh.add_point(atomj.get_position(), ctr);
            ctr++;
          }
          local_nodes[j] = it->second;
          elem->set_node(j) = mesh.node_ptr(local_nodes[j]);
        }

        elem_ctr++;
      }
      */
    }
  }

  TetGenMeshInterface tetgenif(mesh);
  tetgenif.triangulate_pointset();

  MeshBase::element_iterator it(mesh.elements_begin());
  while (it != mesh.elements_end())
  {
    Elem* el = *it;
//    if (el->volume() < 0)
//    {
//      unsigned int n1 = el->node(0);
//      el->set_node(0) = mesh.node_ptr(el->node(1));
//      el->set_node(1) = mesh.node_ptr(n1);
//    }

    ++it;

    if (el->volume() <= 1e-12)
    {
      // eliminate all degenerate elements
      mesh.delete_elem(el);
    }
    else
    {
      Point centroid(el->centroid());
      const Elem* dev_el = MeshUtils::search_element(
          &(get_device()->get_mesh()), centroid);
      ID id = INVALID_ID;
      if (dev_el != NULL)
        id = dev_el->subdomain_id();

      el->subdomain_id() = id;

      // eliminate all elements that seem to lie outside of the structure
      if (id == INVALID_ID)
        mesh.delete_elem(el);
    }
  }

  mesh.prepare_for_use(keep_node_order);

  Messages::info("done");
#endif
}





void
AtomisticStructure::extract_statistics(map<Specie, vector<unsigned int>>& stats,
    const set<ID>& regions, double cutoff, double y, double z, const ModelOptions& opt) const
{

  //cutoff = _options.get_option("control_volume_radius", cutoff);

  int ref_atom = _options.get_option("reference_atom", -1);

  int bin_size = opt.get_option("bin_size", 0);
  set<unsigned int> binning_atom;
  opt.get_option("binning_atom", binning_atom);

  Messages m;
  m.info("Extracting alloy statistics for structure " + this->get_name() + ": ");
  m.indent();
  ostringstream os;
  if (bin_size > 0)
  {
    os << "bin size: " << bin_size << ", counting atoms with label ";
    for (auto&& ba : binning_atom)
      os << ba << " ";
  }
  else
    os << "control sphere radius: " << cutoff << " nm";

  m.info(os.str());

  if (ref_atom >= 0)
  {
    os.str("");
    os << "control volumes centered at atoms with label (= reference atom): " << ref_atom;
    m.info(os.str());
  }


  // A random starting seed is needed to actually have different sequences
  // we try to use something that is different also if launching simulations
  // at the same time
  int seed = opt.get_option("random_generator_seed",
      static_cast<int>(time(NULL) * std::tr1::random_device()()));
  {
    std::ostringstream os;
    os << "Initializing  MT19937 random generator with seed " << seed;
        //std::ios::hex << seed;
    Messages::info(os.str());
  }
  std::tr1::mt19937 generator(seed);


  // we can induce an artificial detection efficiency, to model e.g. APT
  map<Specie, double> detection_efficiency;
  {
    vector<string> det_eff;
    opt.get_option("detection_efficiency", det_eff);
    if (det_eff.size() % 2 != 0)
    {
      Messages::error("detection_efficiency for alloy statistics is ill defined: "
          "need pairs (species, efficiency)");
    }
    else
    {
      for (int i = 0; i < det_eff.size(); i += 2)
      {
        double eff = Utils::convert<double>(det_eff[i+1]);
        if ((eff > 1.0) || (eff < 0))
          Messages::error("Wrong detection efficiency for " + det_eff[i] + ": "
              + det_eff[i+1]);
        else
          detection_efficiency[det_eff[i]] = eff;
      }
    }
  }

  std::vector<unsigned char> detect;
  if (!detection_efficiency.empty())
  {
    Messages m;
    m.indent();
    m.info("Using detection efficiencies:");
    m.indent();
    for (auto&& it : detection_efficiency)
    {
      std::ostringstream os;
      os << it.first << " : " << it.second << Messages::endl;
      Messages::info(os.str());
    }

    detect.resize(get_N_atoms(), 0);

    for (unsigned int i = 0; i < get_N_atoms(); i++)
    {
      const Atom& atm = get_structure_atom(i);

      if (atm.get_label() == 0)
        continue;

      ID subdomain = atm.get_region_ID();

      if (regions.count(subdomain))
      {
        if (!stats.count(atm.get_specie()))
        {
          if (!detection_efficiency[atm.get_specie()])
            detection_efficiency[atm.get_specie()] = 1;

          double rnd = static_cast<double>(generator()) / generator.max();
          if (rnd < detection_efficiency[atm.get_specie()])
            detect[i] = 1;

        }
      }
    }
  }


  for (unsigned int i = 0; i < get_N_atoms(); i++)
  {
    const Atom& atm = get_structure_atom(i);

    // we skip atoms with no "real" label (such as passivation H)
    if (atm.get_label() == 0)
      continue;


    ID subdomain = atm.get_region_ID();

    if (regions.count(subdomain))
    {
      if (!stats.count(atm.get_specie()))
      {

        // this species is found for the first time
        size_t len = 1;
        if (stats.begin() != stats.end())
          len = (stats.begin()->second).size();

        stats[atm.get_specie()].resize(len, 0);
      }

      unsigned int counter = 0;

      if (detect.empty() || detect[i])
      {
        stats[atm.get_specie()][0] += 1;
      }

      if ((ref_atom < 0) || (atm.get_label() == ref_atom))
      {

        map<Specie, vector<unsigned int>>::iterator mit(stats.begin());
        map<Specie, vector<unsigned int>>::iterator mend(stats.end());
        for ( ; mit != mend; ++mit)
          (mit->second).push_back(0);

        if (detect.empty() || detect[i])
        {
          stats[atm.get_specie()].back() += 1;
          if (!binning_atom.empty() && binning_atom.count(atm.get_label()))
            counter++;
        }

        neighbor_iterator it(neighbors_begin(i, 10 * cutoff, 10 * y, 10 * z));
        neighbor_iterator end(neighbors_end(i));
        for ( ; (it != end) && (counter < bin_size); ++it)
        {
          const Atom& neigh = *(*it);

          if (!regions.count(neigh.get_region_ID()))
            continue;

          if (detect.empty() || detect[it.atom_index()])
          {
            if (!binning_atom.empty() && binning_atom.count(neigh.get_label()))
              counter++;

            if (!stats.count(neigh.get_specie()))
            {
              stats[neigh.get_specie()].resize(stats[atm.get_specie()].size(), 0);
              stats[neigh.get_specie()][0] = 1;
            }

            stats[neigh.get_specie()].back() += 1;
          }
        }
      }
    }
  }
}


void
AtomisticStructure::extract_statistics(unsigned int atom,
    map<Specie, unsigned int>& counts, const set<ID>& regions, double cutoff) const
{
  const Atom& at = this->get_structure_atom(atom);
  counts[at.get_specie()] += 1;

  neighbor_iterator it(neighbors_begin(atom, 10 * cutoff));
  neighbor_iterator end(neighbors_end(atom));
  for ( ; it != end; ++it)
  {
    const Atom& neigh = *(*it);

    if (!regions.count(neigh.get_region_ID()))
      continue;

    counts[neigh.get_specie()] += 1;
  }
}


int
AtomisticStructure::find_nearest_atom(const Elem* elem, const Point& point, double )
{
  ID id = elem->subdomain_id();

  // atomic coordinates are in Angstrom
  Point center = _scale * point;


  const vector<unsigned int>& atoms = get_atoms_in_elem(elem);
  int atom = -1;

  //
  // it may not even contain atoms
  if (!atoms.empty())
  {
    // ok, there is at least one atom inside
    double min_dist = 1e9;

    // look for the atom nearest to the center
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      if (get_structure_atom(atoms[i]).get_specie() == Specie::H)
        continue;

      double dist = Point(center -
          get_structure_atom(atoms[i]).get_position()).size();
      if (dist < min_dist)
      {
        min_dist = dist;
        atom = atoms[i];
      }
    }
    //cerr << " found in elem: " << get_structure_atom(atom).get_position();
  }
  else
  {
    // we try to find some nearby atom, proceeding in shells
    set<const Elem*> processed_elems;
    set<const Elem*> to_process;
    set<const Elem*> next_to_process;
    to_process.insert(elem);

    double min_dist = 1e9;

    while (!to_process.empty())
    {
      set<const Elem*>::iterator it(to_process.begin());
      const Elem* next_el = *it;

      const vector<unsigned int>& atoms = get_atoms_in_elem(next_el);
      if (!atoms.empty())
      {
        // ok, there is at least one atom inside
        // atomic coordinates are in Angstrom

        // look for the atom nearest to the center
        for (unsigned int i = 0; i < atoms.size(); ++i)
        {
          if (get_structure_atom(atoms[i]).get_specie() == Specie::H)
            continue;

          double dist = Point(center -
              get_structure_atom(atoms[i]).get_position()).size();

          if (dist < min_dist)
          {
            min_dist = dist;
            atom = atoms[i];
          }
        }

        // and now we can jump out and look at the neighbors
        break;
      }

      for (int s = 0; s < next_el->n_sides(); s++)
      {

        const Elem* neigh = next_el->neighbor(s);

        if ((neigh != NULL) && //!processed_elems.count(neigh))
            !processed_elems.count(neigh) &&
            !to_process.count(neigh) &&
            (Point(point - neigh->centroid()).size() <
                min_dist / _scale))
        {
          next_to_process.insert(neigh);
        }
      }

      processed_elems.insert(next_el);
      to_process.erase(it);

      if (to_process.empty())
      {
        to_process = next_to_process;
        next_to_process.clear();
      }
    }
  }

  return(atom);
}
   
//
// COMPUTE STATISTICS ON ATOMISTIC STRUCTURE
//
void
AtomisticStructure::extract_alloy_statistics(const ModelOptions& opt)
{

  // we can induce an artificial detection efficiency, to model e.g. APT
  map<Specie, double> detection_efficiency;
  vector<string> det_eff;
  opt.get_option("detection_efficiency", det_eff);
  if (det_eff.size() % 2 != 0)
  {
    Messages::error("detection_efficiency for alloy statistics is ill defined: "
        "need pairs (species, efficiency)");
  }
  else
  {
    for (int i = 0; i < det_eff.size(); i += 2)
    {
      double eff = Utils::convert<double>(det_eff[i+1]);
      if ((eff > 1.0) || (eff < 0))
        Messages::error("Wrong detection efficiency for " + det_eff[i] + ": "
            + det_eff[i+1]);
      else
        detection_efficiency[det_eff[i]] = eff;
    }
  }

  // A random starting seed is needed to actually have different sequences
  // we try to use something that is different also if launching simulations
  // at the same time
  int seed = opt.get_option("random_generator_seed",
      static_cast<int>(time(NULL) * std::tr1::random_device()()));
  {
    std::ostringstream os;
    os << "Initializing  MT19937 random generator with seed " << seed;
        //std::ios::hex << seed;
    Messages::info(os.str());
  }
  std::tr1::mt19937 generator(seed);


  if (opt.get_option("control_volume", string("sphere")) == "column")
  {

    libMesh::RealVectorValue a, b, c;
    get_lattice_vectors(a, b, c);

    double x = a.size() / 10;
    double y = b.size() / 10;
    double z = c.size() / 10;

    string dir = opt.get_option("orientation", "x");

    Point center(get_structure_atom(0).get_position());
    center *= 0.1;

    // first find atoms on a plane orthogonal to the column orientation
    switch (dir[0])
    {

      case 'y':
        center(1) += y / 2.0;
        y = 0.2;
        break;

      case 'z':
        center(2) += z / 2.0;
        z = 0.2;
        break;

      default: // 'x'
        center(0) += x / 2.0;
        x = 0.2;
        break;
    }

    ostringstream os;
    os << "(" << x << "," << y << "," << z << ")";
    ModelOptions p_opts;
    p_opts["shape"] = "box";
    p_opts["sides"] = os.str();
    os.str("");
    os << "(" << center(0) << "," << center(1) << "," << center(2) << ")";
    p_opts["center"] = os.str();


    vector<unsigned int> atomlist;
    get_subset(atomlist, p_opts);

    /*
    {
      ofstream of("test.xyz");
      of << atomlist.size() << "\n\n";
      for (unsigned int i = 0; i < atomlist.size(); ++i)
      {
        const Atom& atom = get_structure_atom(atomlist[i]);
        of << atom.get_specie() << "  " <<
            atom.get_position()(0) << "  " <<
            atom.get_position()(1) << "  " <<
            atom.get_position()(2) << "\n";
      }
    }
    */

    map<Specie, vector<unsigned int>> stats;

    ModelOptions c_opts(opt);
    c_opts["shape"] = "column";

    for (unsigned int i = 0; i < atomlist.size(); ++i)
    {
      Point pos = get_structure_atom(atomlist[i]).get_position();
      pos *= 0.1;

      ostringstream os;
      os << "(" << pos(0) << "," << pos(1) << "," << pos(2) << ")";
      c_opts["center"] = os.str();

      vector<unsigned int> column_atoms;
      get_subset(column_atoms, c_opts);

      /*
      {
        ostringstream fname;
        fname << "test" << i << ".xyz";
        ofstream ff(fname.str());
        ff << column_atoms.size() << "\n\n";
        for (unsigned int j = 0; j < column_atoms.size(); ++j)
        {
          const Atom& atom = get_structure_atom(column_atoms[j]);
          ff << atom.get_specie() << "  " <<
              atom.get_position()(0) << "  " <<
              atom.get_position()(1) << "  " <<
              atom.get_position()(2) << "\n";
        }
      }
      */


      for (unsigned int j = 0; j < column_atoms.size(); ++j)
      {
        const Atom& atom = get_structure_atom(column_atoms[j]);

        //if (!regions.count(neigh.get_region_ID()))
        //  continue;

        double rnd = static_cast<double>(generator()) / generator.max();

        if (rnd < detection_efficiency[atom.get_specie()])
        {
          if (!stats.count(atom.get_specie()))
          {
            stats[atom.get_specie()].resize(atomlist.size(), 0);
            stats[atom.get_specie()][i] = 1;
          }

          stats[atom.get_specie()][i] += 1;
        }
      }
    }

    ofstream of(TiberCad::get_output_dir() + "/" + get_name() +
        "_columnwise_statistics.dat");


    of << "% alloy statistics for structure " << get_name() <<
        ", columnwise along " << dir << "\n";
    //of << "% (first row gives total numbers)\n";
    of << "% ";

    int NN = 0;
    map<Specie, vector<unsigned int>>::iterator mit(stats.begin());
    const map<Specie, vector<unsigned int>>::iterator mend(stats.end());
    for ( ; mit != mend; ++mit)
    {
      of << mit->first << "  ";
      NN = (mit->second).size();
    }
    of << "\n";

    //of << "% ";
    //for (mit = stats.begin(); mit != mend; ++mit)
    //  of << (mit->second)[0] << " ";
    //of << "\n";

    vector<unsigned int> sums(NN, 0);
    unsigned int max = 0;
    for (int i = 0; i < NN; ++i)
    {
      for (mit = stats.begin(); mit != mend; ++mit)
        sums[i] += (mit->second)[i];

      if (sums[i] > max)
        max = sums[i];
    }

    for (int i = 0; i < NN; ++i)
    {
      //if (sums[i] >= (max - 1))
      {
        for (mit = stats.begin(); mit != mend; ++mit)
          of << (mit->second)[i] << " ";
        of << "\n";
      }
    }

  }
  else
  {


    double cutoff = opt.get_option("control_volume_radius", 0.5);
    string cutoff_str = opt.get_option("control_volume_radius", "0.5");

    vector<string> reg_names(1, "all");
    _options.get_option("regions", reg_names);

    for (unsigned int i = 0; i < reg_names.size(); ++i)
    {

      IDSet reg_ids;
      _device->extract_physical_regions(reg_names[i], reg_ids);
      IDSet::iterator id_it(reg_ids.begin());
      while (id_it != reg_ids.end())
      {
        const Material* mat = _device->get_material(*id_it);

        // we extract statistics only for alloys
        IDSet::iterator to_be_deleted(id_it);
        ++id_it;

        if (!mat->is_alloy())
          reg_ids.erase(to_be_deleted);
      }

      // if there are no ids, we have no alloy
      if (reg_ids.empty())
        continue;

      string reg_name = reg_names[i];

      ofstream of(TiberCad::get_output_dir() + "/" + get_name() +
          "_" + reg_name + "_statistics_R" + cutoff_str + ".dat");

      // define the control volume
      string shape = opt.get_option("control_volume", "sphere");
      double x = cutoff;
      double y = 0;
      double z = 0;

      if (shape == "cube")
      {
        x = opt.get_option("side_length", cutoff);
        y = z = x;
      }
      else if (shape == "column")
      {
        libMesh::RealVectorValue a, b, c;
        get_lattice_vectors(a, b, c);
        string dir = opt.get_option("orientation", "x");
        double side = opt.get_option("side_length", cutoff);
        switch (dir[0])
        {

          case 'y':
            y = b.size() / 10.0;
            x = z = side;
            break;

          case 'z':
            z = c.size() / 10.0;
            y = x = side;
            break;

          default: // 'x'
            x = a.size() / 10.0;
            y = z = side;
            break;
        }
      }

      map<Specie, vector<unsigned int>> stats;
      extract_statistics(stats, reg_ids, x, y, z, opt);

      of << "% alloy statistics for structure " << get_name() <<
          ", region " << reg_name << "\n";
      of << "% (first row gives total numbers)\n";
      of << "% ";

      int NN = 0;
      map<Specie, vector<unsigned int>>::iterator mit(stats.begin());
      const map<Specie, vector<unsigned int>>::iterator mend(stats.end());
      for ( ; mit != mend; ++mit)
      {
        of << mit->first << "  ";
        NN = (mit->second).size();
      }
      of << "\n";

      of << "% ";
      for (mit = stats.begin(); mit != mend; ++mit)
        of << (mit->second)[0] << " ";
      of << "\n";

      vector<unsigned int> sums(NN, 0);
      unsigned int max = 0;
      for (int i = 1; i < NN; ++i)
      {
        for (mit = stats.begin(); mit != mend; ++mit)
          sums[i] += (mit->second)[i];

        if (sums[i] > max)
          max = sums[i];
      }

      for (int i = 1; i < NN; ++i)
      {
        //if (sums[i] >= (max - 1))
        {
          for (mit = stats.begin(); mit != mend; ++mit)
            of << (mit->second)[i] << " ";
          of << "\n";
        }
      }

    }
  }
   
}

void
AtomisticStructure::plot_alloy_composition(const ModelOptions& opt)
{
  // CREATE A VTK FILE PROJECTING LOCAL ALLOY COMPOSITION 
  double cutoff = opt.get_option("control_volume_radius", 0.5);
  string cutoff_str = opt.get_option("control_volume_radius", "0.5");
  libMesh::UniquePtr<libMesh::UnstructuredMesh> mesh(new libMesh::Mesh(TiberCad::get_mpi_comm(), 3));

  int ref_atom = _options.get_option("reference_atom", -1);
  IDSet refatoms;
  if (ref_atom >= 0)
    refatoms.insert(ref_atom);
  create_conformal_grid(*mesh, refatoms, true);

  libMesh::UniquePtr<DataOutput> writer(DataOutput::create(
      _options.get_option("meshdata_format", "vtk")));

  writer->set_mesh(*mesh);
  writer->set_output_directory(TiberCad::get_output_dir());
  writer->set_filename(get_name() + "_alloycomposition_R" + cutoff_str);
  //writer->write(1);

  map<ID, map<SolutionDescriptor, vector<double>>> solmap;
  map<Specie, SolutionDescriptor> species_to_descr;

  // setup a map Specie->SolutionDescriptor
  const vector<string>& atom_types = get_atom_types();
  unsigned int ctr = 0;
  for (unsigned int i = 0; i < atom_types.size(); ++i)
  {
    Specie sp(atom_types[i]);
    if (!species_to_descr.count(sp) && !(sp == Specie::H))
    {
      SolutionDescriptor desc(atom_types[i], ctr,
        SolutionDescriptor::REAL, SolutionDescriptor::NODES);

      species_to_descr[sp] = desc;
      ++ctr;
    }
  }

  IDSet reg_ids;
  IDSet::iterator id_it(_IDset.begin());
  const IDSet::iterator id_end(_IDset.end());
  for ( ; id_it != id_end; ++id_it)
  {

    reg_ids.insert(*id_it);

    // we extract statistics only for alloys
    const Material* mat = _device->get_material(*id_it);
    if (!mat->is_alloy())
      continue;
    map<Specie, SolutionDescriptor>::iterator s_it(species_to_descr.begin());
    const map<Specie, SolutionDescriptor>::iterator s_end(species_to_descr.end());
    for ( ; s_it != s_end; ++s_it)
      solmap[*id_it][s_it->second].resize(0);
  }

  // define the control volume
  string shape = opt.get_option("control_volume", "sphere");
  double x = cutoff;
  double y = 0;
  double z = 0;

  if (shape == "cube")
  {
    x = opt.get_option("side_length", cutoff);
    y = z = x;
  }
  else if (shape == "column")
  {
    libMesh::RealVectorValue a, b, c;
    get_lattice_vectors(a, b, c);
    string dir = opt.get_option("orientation", "x");
    double side = opt.get_option("side_length", cutoff);
    switch (dir[0])
    {

      case 'y':
        y = b.size() / 10.0;
        x = z = side;
        break;

      case 'z':
        z = c.size() / 10.0;
        y = x = side;
        break;

      default: // 'x'
        x = a.size() / 10.0;
        y = z = side;
        break;
    }
  }

  map<Specie, vector<unsigned int>> stats;
  extract_statistics(stats, reg_ids, x, y, z, opt);

  for (id_it = _IDset.begin(); id_it != id_end; ++id_it)
  {
    ID domain = *id_it;

    // we extract statistics only for alloys
    const Material* mat = _device->get_material(*id_it);
    if (!mat->is_alloy())
      continue;

    // to keep track of already used nodes
    set<unsigned int> used_nodes;

    map<SolutionDescriptor, vector<double>>::iterator solit(solmap[domain].begin());
    const map<SolutionDescriptor, vector<double>>::iterator solend(solmap[domain].end());
    for ( ; solit != solend; ++solit)
    {
      (solit->second).clear();
      (solit->second).reserve(mesh->n_nodes());
    }


    MeshBase::element_iterator elit(mesh->active_local_elements_begin());
    const MeshBase::element_iterator elend(mesh->active_local_elements_end());
    for ( ; elit != elend; ++elit)
    {
      const Elem* elem = *elit;

      if (elem->subdomain_id() == domain)
      {
        for (unsigned int n = 0; n < elem->n_nodes(); ++n)
        {
          if (!used_nodes.count(elem->node(n)))
          {
            used_nodes.insert(elem->node(n));

            map<Specie, vector<unsigned int>>::iterator it(stats.begin());
            const map<Specie, vector<unsigned int>>::iterator end(stats.end());
            for ( ; it != end; ++it)
            {
              auto desc(species_to_descr.find(it->first));
              if (desc != species_to_descr.end())
              {
                SolutionDescriptor& descr = species_to_descr[it->first];
                // +1 because the first values are the totals !
                solmap[domain][descr].push_back((it->second)[elem->node(n) + 1]);
              }
            }
          }
        }
      }
    }
    writer->set_data(solmap[domain], domain);
  }

  writer->write();
  
}

void
AtomisticStructure::radial_distribution(std::string suffix)
{
  // COMPUTE AUTOCORRELATION FUNCTION
  const ModelOptions& opt = (_options.submodels_begin("radial_distribution"))->second;
  double Rc = opt.get_option("cutoff_radius", 1.0) * 10.0;
  string Rc_str = opt.get_option("cutoff_radius", "1.0");
  double dr = opt.get_option("resolution",0.001) * 10.0; 
  vector<string> species;
  opt.get_option("species",species);

  Messages::info("output g(r)");
  std::set<ID> reg_ids(get_IDset());
  // get total volume
  MeshBase::element_iterator el(_device->get_mesh().elements_begin());
  const MeshBase::element_iterator el_end(_device->get_mesh().elements_end());
  double Vol = 0.0;
  for ( ; el != el_end; ++el)
  {
    Elem* elem = *el;
    if (reg_ids.count(elem->subdomain_id()))
         Vol += elem->volume();
  } 
  // transform mesh_units^3 to Ang^3
  Vol *= _scale * _scale * _scale; 
  
  unsigned int N_atoms = get_N_atoms();

  for (unsigned int i=0; i<species.size(); i++) 
  {
    ofstream of(TiberCad::get_output_dir() + "/" + get_name() +
          "_" + species[i] + "_radial" + suffix + ".dat");

    Specie sp(species[i]);
    vector<map<Specie,unsigned int>> g;
    compute_g(sp, Rc, dr, g);

    of << "# r "; 
    map<Specie,unsigned int>::iterator it = g[1].begin();
    map<Specie,unsigned int>::iterator itend = g[1].end();
    for ( ; it != itend; ++it)
       of << (it->first).get_string() << "   ";

    of << "\n";      

    unsigned int Ng = floor(Rc/dr);

    for (unsigned int i=1; i<Ng; i++)
    {
      double r = i * dr; 
      of << r << "  "; 
      it = g[i].begin();
      itend = g[i].end();
      for ( ; it != itend; ++it)
      {
        of << it->second * Vol/(N_atoms*N_atoms*4.0*M_PI*r*r*dr)  << "  "; 
      }
      of << "\n";
    }
  }

}


void
AtomisticStructure::compute_g(const Specie& sp, double Rc, double dr, vector<map<Specie,unsigned int>>& g)
{
  // cutoff: Rc
  // spacing: dr
  // vector<int> g  
  const BondMap& bm = *_bondmap;
  unsigned int Nb = floor(Rc/dr);
  g.resize(Nb); 
  // a set with all species:
  set<Specie> spset;

  for (unsigned int atom=0; atom < get_N_atoms(); atom++)
  { 
    spset.insert(_atoms[atom].get_specie());
    if (_atoms[atom].get_specie() == sp)
    {
      Point p1(_atoms[atom].get_position());
      neighbor_iterator it = neighbors_begin(atom, Rc);
      neighbor_iterator itend = neighbors_end(atom, Rc);

      for ( ; it != itend; ++it)
      {
         Specie sp2( (*it)->get_specie() );
         Point p2( (*it)->get_position() + it.atom_translation() );
         double d = sqrt((p1 - p2)*(p1 - p2));
         unsigned int bb = floor(d/dr);
         if (bb < g.size()) g[bb][sp2] += 1;
      }     
    }
  }

  // add zeros to all points not added in the map by previous loop
  for (unsigned int i=0; i< g.size(); i++)
  {
     set<Specie>::iterator spit = spset.begin();
     set<Specie>::iterator spend = spset.end();
     for ( ; spit != spend; ++spit)
     {
       if (g[i].count(*spit) == 0){ g[i][*spit] = 0; }
     }
  }

}


