// $Id$

#include "AtomisticStructure.h"
#include "AtomisticGenerator.h"
#include "AtomisticGenerator1D.h"
#include "AtomisticGenerator2D.h"
#include "AtomisticGenerator3D.h"
#include "BondMap.h"
#include "Messages.h"
#include "MeshUtils.h"
#include "TiberCad.h"
#include "Material.h"
#include "RuntimeException.h"

//C++ includes
//--------------------
#include <vector>
#include <set>
#include<iostream>
#include<fstream>
#include<sstream>
#include <map>
//---------------------



AtomisticStructure::AtomisticStructure()
:_bondmap(NULL),
 _device(NULL),
 _random_alloy(false),
 _clustering(false),
  _atom_types()
{
  N_atoms = 0;
  for (unsigned int i = 0; i < 9; i++)
    {
      _periodicity_vectors[i] = 0.0;
    }
}


AtomisticStructure::AtomisticStructure(const std::string& name)
:_name(name),
 _bondmap(NULL),
 _atomistic_structure_options(),
 _scale(1.0),
 _device(NULL),
 _random_alloy(false),
 _atom_types(),
 _N_without_H(0)
{
  // Default initializations
  N_atoms = 0;
  for (unsigned int i = 0; i < 9; i++)
    {
      _periodicity_vectors[i] = 0.0;
    }

}


AtomisticStructure::~AtomisticStructure(void)
{
  if (_bondmap != NULL) delete _bondmap;
}


AtomisticStructure*
AtomisticStructure::create()
{
  AtomisticStructure* st =  NULL;
  st = new AtomisticStructure();
  return st;
}


AtomisticStructure::AtomisticStructureOptions::AtomisticStructureOptions(void)
:is_periodical(false),
 is_associated(false)
{}


AtomisticStructure::AtomisticStructureOptions::~AtomisticStructureOptions(void)
{}


void
AtomisticStructure::init(const std::string& name,
    const Device* const device, const ModelOptions& options)
{
  std::ostringstream os;
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
  // Read from file

  random_alloy = _options.get_option("random_alloy", "false");
  if (random_alloy == "true") _random_alloy = true;

  _clustering = _options.get_option("clustering", false);


  if (_options.find_option("load_structure")||_options.find_option("load"))
    {
      std::string filename;

      //------------------------------------------------------------
      os << "Reading structure from file " << filename <<
          ". Any further information will be neglected" << std::endl;
      Messages::info(os.str(), true);
      os.str(std::string());
      //---------------------------------------------------------------

      if (_options.find_option("load_structure")) filename = _options.get_option("load_structure","none");
      if (_options.find_option("load")) filename = _options.get_option("load","none");

      init(filename);

      parse_regions();

      if (_atomistic_structure_options.is_associated == false) associate_elements();

      //TODO: I'm calculating the bond map again anyway because otherwise the translation vectors are not
      //correctly reproduced for periodic structures. If we really want to import the bond map,
      //this needs to be changed
      if (_bondmap != NULL)
      {
        delete _bondmap;
        _bondmap = NULL;
      }
      if (_bondmap == NULL) build_bond_map();


      print_driver();

    }
  // Build material representation
  //else if ( (_options.find_option("reference_region"))
  //    && (! _options.find_option("regions")) )
  //{

  //}

  // Build mesh based representation
  else if ( (_options.find_option("reference_region"))
      && ( _options.find_option("regions")) )
    {
      init_mesh_structure();
    }

  else
    {
      //------------------------------------------------------------
      os << "Error in Atomistic structure " << _name <<
          "input block. You need to define reference_region and regions." << std::endl;
      Messages::info(os.str(), true);
      os.str(std::string());
      //---------------------------------------------------------------
    }

  //Calculate the number of atoms excluding hydrogens 
  //(Useful for passivated semiconductors)
  compute_N_without_H();

}

void
AtomisticStructure::parse_lattice_vectors(void)
//Overrides supercell options
{
  if (_options.find_option("lattice_vectors"))
  {
    std::vector<double> lattice_vectors;
    _options.get_option("lattice_vectors", lattice_vectors);
 
    if (lattice_vectors.size() == 3)
    {  
       _atomistic_structure_options.is_periodical = true;
       _periodicity_vectors[0] = lattice_vectors[0];
       _periodicity_vectors[4] = lattice_vectors[1];
       _periodicity_vectors[8] = lattice_vectors[2];
    }
 
    else if (lattice_vectors.size() == 9)
    {
       for (int i = 0; i < 9; i++)
       {
         _atomistic_structure_options.is_periodical = true;
         _periodicity_vectors[i] = lattice_vectors[i]; 
       }
    }
 
    else Messages::error("lattice_vectors must have 3 or 9 components");
  }

} 

void
AtomisticStructure::init(const std::string& filename)
{
  std::ostringstream os;

  read_structure(filename);

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
  else
    {
      print_structure(_name + ".xyz");
    }
}

void
AtomisticStructure::parse_regions(void)
{
  std::string physreg = _options.get_option("regions", "all");
  _device->extract_physical_regions(physreg, _IDset);

  /*

  //If no physical regions are specified, "all" is the default
  if ( (get_options().is_empty()) || (!get_options().find_option("physical_regions")) )
    _regionset.insert("all");

  //Put physical regions specified in input file in _regions
  //A vector is needed as temporary container
  std::vector<std::string> region_string;
  _options.get_option("physical_regions", region_string);

  for (unsigned int i = 0; i < region_string.size(); i++)
  {
    _regionset.insert(region_string[i]);
  }
  region_string.clear();

  //If all regions are specified (value = "all", must fill with real names of all regions)
  if ( _regionset.count("all") == 1)
  {
    _regionset.clear();
    std::set < ID >::iterator region_ID_iterator = get_device()->get_active_region_ids().begin();

    for (unsigned int i = 0; i < get_device()->get_active_region_ids().size(); i++)
    {
      _regionset.insert( get_device()->get_region_name(*region_ID_iterator) );
      region_ID_iterator ++;
    }

  }

  //Build an array of physical regions ID
  for (std::set<std::string>::iterator i = _regionset.begin(); i != _regionset.end(); i++)
  {
    std::set<ID> tmp_ID;
    get_device()->get_active_region_ids( (*i), tmp_ID);
    _IDset.insert(tmp_ID.begin(), tmp_ID.end());
  }
   */
}


void
AtomisticStructure::build_bond_map(void)
{
  if (_bondmap != NULL)
    {
      delete _bondmap;
    }
  _bondmap = new BondMap;
  _bondmap->do_init(_structure_atoms.size());
  Tensor2Gen period;
  for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
        {
          period(j + 1, i + 1) = _periodicity_vectors[i*3 + j];
        }
    }
  _bondmap->do_solve(_structure_atoms, period);

}


void
AtomisticStructure::init_mesh_structure()
{
  std::ostringstream os;
  assert(_device != NULL);

  // Read structure from file
  std::string path;

  //---------------------------------------------------------------
  os << "Mesh_units is " << _device->get_mesh_units() << std::endl
      << "Scale factor is " << _scale << std::endl;
  Messages::debug(os.str());
  os.str(std::string());
  //--------------------------------------------------------------

  parse_regions();

  //--------------------------------------------------------------
  os << "Atomistic structure builder started " << path << std::endl;
  Messages::info(os.str(), true);
  os.str(std::string());
  //-----------------------------------------------------------

  _structure_atoms.clear();

  //---------------------------------------------------------------
  // Extend mesh for contacts
  //unsigned int num_sides;
  //Point normal = Device::get_normal();


  //---------------------------------------------------------------

  AtomisticGenerator* generate;

  if ( _device->get_mesh().mesh_dimension() == 1 ) generate = static_cast<AtomisticGenerator1D*> ( AtomisticGenerator::create(this, 1 ) );
  if ( _device->get_mesh().mesh_dimension() == 2 ) generate = static_cast<AtomisticGenerator2D*> ( AtomisticGenerator::create(this, 2 ) );
  if ( _device->get_mesh().mesh_dimension() == 3 ) generate = static_cast<AtomisticGenerator3D*> ( AtomisticGenerator::create(this, 3 ) );

  generate->do_init();
  parse_lattice_vectors();
  build_bond_map();

  //if (_atomistic_structure_options.is_associated == false) associate_elements();

  //Refresh some information after structure building
  N_atoms = _structure_atoms.size();

  print_driver();

}


void
AtomisticStructure::associate_elements()
{
  //Associate atoms with NULL element pointer to right mesh elements
  //TODO: it's almost a O(n^2) algorithm, but with no smarter solution
  //to manage mesh and atoms together it's the only way!!!

  Messages::debug("Starting associate_elements");

  bool set = false;
  Point p;
  unsigned int dim = get_device()->get_mesh().mesh_dimension();

  //Get iterators to all elements
  MeshBase::element_iterator  el_start = get_device()->get_mesh().elements_begin();
  MeshBase::element_iterator  el_end = get_device()->get_mesh().elements_end();
  MeshBase::element_iterator  it = el_start;

  for (unsigned int i = 0; i < get_structure_atoms().size(); i++)
    {
      //Up to now avoid association for hydrogens (as they're always (UP TO NOW) passivation atoms
      // and their associated element is not important)
      if (get_structure_atoms()[i].get_elem() == NULL)
        {
          set = false;
          p(0) = 0.0; p(1) = 0.0; p(2) = 0.0;
          p(0) = get_structure_atoms()[i].get_position(0) / get_scale();
          if ( (dim == 2) || (dim == 3) )   p(1) = get_structure_atoms()[i].get_position(1) / get_scale();
          if ( (dim == 3) )  p(2) = get_structure_atoms()[i].get_position(2) / get_scale();

          for (it = el_start; it != el_end; it++)
            {
              Elem* elem = *it;
              if (MeshUtils::may_belong_to_element(elem,p))
                {
                  if ( (elem->contains_point(p) ) )
                    {
                      _structure_atoms[i].set_elem(elem);
                    }

                  set = true;
                }
            }

          //TODO: UNCOMMENT THIS LINE, COMMENTED ONLY FOR DIRTY WORKS PURPOSE
          //(some passivation atoms stand out of mesh)
          //if (!set) Messages::warning("An atom has NULL element pointer: it stands outside mesh! ");

        }
    }

  Messages::debug("Finished associate_elements");

}


void
AtomisticStructure::read_structure(const std::string& path)
{

  std::ifstream file;
  std::string line, record;
  unsigned int n_specie;
  Atom tmp_atom;
  Tensor1 pos;

  Messages::debug("Reading structure from file");

  // Delete eventually existing structure
  if (!(_structure_atoms.empty())) _structure_atoms.clear();
  if (!(_atom_types.empty())) _atom_types.clear();

  file.open(path.c_str(), std::ifstream::in);

  if (!file)
    {
      std::cerr << "Unable to open file " << path << ". Cannot read Atomistic Structure. \n";
      exit(1);   // call system to stop
    }

  // Recognize type of input file and read it properly
  std::string extension = path.substr(path.size()-4);

  // XYZ file
  if ( (extension.compare(".xyz") == 0) || (extension.compare(".XYZ") == 0) )
    {
      // First line is number of atoms
      getline(file, line);
      N_atoms = atoi(line.c_str());

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

          tmp_atom.set_position( pos );

          _structure_atoms.push_back(tmp_atom);

        }

      if ( (_structure_atoms.size()) != N_atoms )
        std::cerr << "Warning: in file xyz number of atoms is wrong \n";


      N_types = _atom_types.size();

      // Warning: XYZ file has no informations about structure periodicity

    }

  // GEN file
  else if ( (extension.compare(".gen") == 0) || (extension.compare(".GEN") == 0) )

    {
      getline(file, line);
      std::stringstream line_string(line);

      line_string >> record;

      N_atoms = atoi(line.c_str());

      //#ifdef DEBUG
      //      std::cerr << "N_atoms is " << N_atoms << std::endl;
      //#endif

      if (N_atoms == 0)
        {
          std::cerr << "No atoms in structure files or non valid integer in first line. \n";
          exit(1);
        }

      line_string >> record;

      if ( (record.compare("S") == 0) || (record.compare("s") == 0))
        _atomistic_structure_options.is_periodical = true;
      else  if ( (record.compare("C") == 0) || (record.compare("c") == 0))
        _atomistic_structure_options.is_periodical = false;
      else
        std::cerr << "Warning (in GEN file at first line): Cluster (C) or Supercell (S) must be specified. By default a Cluster (no periodicity) is considered. \n";

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

      N_types = _atom_types.size();

      // Cycle upon specified number of atoms (last rows are for periodicity vectors)
      for (unsigned int i = 1; i <= N_atoms; i++)
        {
          getline(file, line);
          std::stringstream line_string(line);

          // First value is ignored (just atoms enumeration)
          line_string >> record;
          line_string >> record;
          n_specie = atoi(record.c_str());
          tmp_atom.set_specie ( _atom_types[n_specie -1] );

          for (unsigned int j = 1; j <= 3; j++)
            {
              line_string >> record;
              pos(j) =  atof(record.c_str());
            }
          tmp_atom.set_position( pos );
          _structure_atoms.push_back(tmp_atom);
        }

      // An additional line is present in GEN files. It's the coordinates origin and it's
      // not needed
      getline(file, line);

      // Read periodicity vectors anyway, if system is not periodical they will be ignored
      //    if (_atomistic_structure_options.is_periodical)
      //    {
      unsigned int count = 0;
      for (unsigned int i = 0; i < 3; i++)
        {
          getline(file, line);
          std::stringstream line_string(line);
          for (unsigned int j = 0; j < 3; j++)
            {
              line_string >> record;
              _periodicity_vectors[count] = atof(record.c_str());
              count++;
            }
        }
      //    }
    }

  else if ( (extension.compare(".tgn") == 0) || (extension.compare(".TGN") == 0) )
    {
      read_tgn(path);
    }

  else
    {
      std::cerr << "Structure file extension is not recognized. \n";
      exit(1);
    }

  file.close();

  //#ifdef DEBUG
  //  std::cerr << "AtomisticStructure::read_structure(path) end. \n";
  //#endif

}


void
AtomisticStructure::read_tgn(const std::string& path)
{

  std::ifstream file;
  std::string line, record;
  unsigned int n_specie;
  Atom tmp_atom;
  Tensor1 pos;
  std::stringstream line_string;
  std::ostringstream os;

  Messages::debug("Reading tgn file. Loading atom coords and bond map ");

  file.open(path.c_str());

  if ( !(file.is_open()) )
    {
      //-------------------------------------------------------------------------------
      os << "Unable to open file " << path << ". Cannot read Atomistic Structure. \n";
      Messages::error(os.str());
      //--------------------------------------------------------------------------------

      exit(1);   // call system to stop
    }

  getline(file, line);

  line_string.clear();
  line_string.str(line);

  line_string >> record;

  N_atoms = atoi(record.c_str());
  _structure_atoms.reserve(N_atoms);

  //Prepare bond map object
  if ( _bondmap == NULL) _bondmap = new BondMap;
  else
    {
      delete _bondmap;
      _bondmap = new BondMap;
    }
  _bondmap->do_init(N_atoms);

  if (N_atoms == 0)
    {
      std::cerr << "No atoms in structure files or non valid integer in first line. \n";
      exit(1);
    }
  _structure_atoms.resize(N_atoms);

  line_string >> record;

  if ( (record.compare("S") == 0) || (record.compare("s") == 0))
    _atomistic_structure_options.is_periodical = true;
  else  if ( (record.compare("C") == 0) || (record.compare("c") == 0))
    _atomistic_structure_options.is_periodical = false;
  else
    std::cerr << "Warning (in GEN file at first line): Cluster (C) or Supercell (S) must be specified. By default a Cluster (no periodicity) is considered. \n";

  getline(file, line);

  //This 2 lines clean stringstream in a safe way
  line_string.clear();
  line_string.str(line);


  while ( line_string >> record)
    {
      _atom_types.push_back(record);
    }


  N_types = _atom_types.size();


  // Cycle upon specified number of atoms (last rows are for periodicity vectors)
  for (unsigned int i = 1; i <= N_atoms; i++)
    {
      getline(file, line);

      line_string.clear();
      line_string.str(line);

      // First value is ignored (just atoms enumeration)
      line_string >> record;
      line_string >> record;
      n_specie = atoi(record.c_str());

      for (unsigned int j = 1; j <= 3; j++)
        {
          line_string >> record;
          pos(j) =  atof(record.c_str());
        }
      _structure_atoms[i - 1].set_specie(_atom_types[n_specie -1]);
      _structure_atoms[i - 1].set_position(pos);

      //Get bond map
      line_string >> record;
      //TODO: write a set_bond_map function. Change pointers in vectors to manage constness more easily
      _bondmap->get_bond_map()[i - 1].resize(atoi(record.c_str()));

      for (unsigned int j = 0; j < _bondmap->get_bond_map()[i - 1].size(); j++)
        {
          line_string >> record;
          _bondmap->get_bond_map()[i - 1][j] = atoi(record.c_str()) - 1;
        }

      line_string >> record;
      int tmp_id = atoi(record.c_str());

      if (tmp_id == -1) _structure_atoms[i - 1].set_elem(NULL);
      else
        {
          _structure_atoms[i - 1].set_elem(_device->get_mesh().elem(tmp_id));
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
          _periodicity_vectors[count] = atof(record.c_str());
          count++;
        }
    }

  file.close();
  _atomistic_structure_options.is_associated = true;

}


void
AtomisticStructure::print_structure(const std::string& path)
{
  std::ofstream file;


  // -------------------------------------------

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

      file.open(file_name.c_str());

      file << _structure_atoms.size() << std::endl << std::endl;

      for (unsigned int i = 0; i < _structure_atoms.size(); i++)
        {
          file << std::setw(2) << _structure_atoms[i].get_specie()
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(0))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(1))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(2)) << "\n";

        }

    }

  else if ( (extension.compare(".xyb") == 0) || (extension.compare(".XYB") == 0) )
    {

      file.open(file_name.c_str());

      file << _structure_atoms.size() << std::endl << std::endl;

      for (unsigned int i = 0; i < _structure_atoms.size(); i++)
        {
          file << std::setw(2) << _structure_atoms[i].get_specie()
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(0))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(1))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(2));

          if (_bondmap != NULL)
            {

              file << std::setw(5) << _bondmap->get_bond_map()[i].size();

              // N.B. Indexing is in Fortran notation (first atom is labelled as 1) !!!!!!!!!!!!!!!!
              for (unsigned int j = 0; j < _bondmap->get_bond_map()[i].size(); j++)
                {
                  file << std::setw(10) << _bondmap->get_bond_map()[i][j] + 1;
                }
              ///////////////////////////////////////////

            }

          file << std::endl;
        }


    }

  else if ( (extension.compare(".tgn") == 0) || (extension.compare(".TGN") == 0) )
    {

      file.open(file_name.c_str());

      file << _structure_atoms.size();

      if (_atomistic_structure_options.is_periodical) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(6) << _atom_types[i];
        }
      file << std::endl;

      for (unsigned int i = 0; i < _structure_atoms.size(); i++)
        {
          unsigned int n_specie;
          for (n_specie = 0; n_specie < _atom_types.size(); n_specie++)
            {
              if ( _atom_types[n_specie] == _structure_atoms[i].get_specie() ) break;
            }
          file << std::setw(10) << i + 1 << std::setw(5) << n_specie + 1
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(0))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(1))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(2)) ;

          if (_bondmap != NULL)
            {

              file << std::setw(5) << _bondmap->get_bond_map()[i].size();

              // N.B. Indexing is in Fortran notation (first atom is labelled as 1) !!!!!!!!!!!!!!!!
              for (unsigned int j = 0; j < _bondmap->get_bond_map()[i].size(); j++)
                {
                  file << std::setw(10) << _bondmap->get_bond_map()[i][j] + 1;
                }
              ///////////////////////////////////////////

            }

          //ID of element is saved (note: no modifications to mesh are allowed to preserve compatibility)
          file << std::setw(14);
          if (_structure_atoms[i].get_elem() == NULL) file << -1;
          else file << _structure_atoms[i].get_elem()->id();

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
                  _periodicity_vectors[count];
              count++;
            }
          file << "\n";
        }

    }


  else if ( (extension.compare(".gen") == 0) || (extension.compare(".GEN") == 0) )
    {

      file.open(file_name.c_str());

      file << _structure_atoms.size();

      if (_atomistic_structure_options.is_periodical) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(6) << _atom_types[i];
        }
      file << std::endl;

      for (unsigned int i = 0; i < _structure_atoms.size(); i++)
        {
          unsigned int n_specie;
          for (n_specie = 0; n_specie < _atom_types.size(); n_specie++)
            {
              if (_atom_types[n_specie] == _structure_atoms[i].get_specie() ) break;
            }
          file << std::setw(10) << i + 1 << std::setw(5) << n_specie + 1
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(0))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(1))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(2)) << "\n";
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
                  _periodicity_vectors[count];
              count++;
            }
          file << "\n";
        }

    }

  else if ( (extension.compare(".upg") == 0) || (extension.compare(".UPG") == 0) )
    {
      print_upg(file_name, "VOID");
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



void
AtomisticStructure::print_upg(const std::string& path, const std::string& etb_dataset)
{

  std::ofstream file, os;

  Messages::debug("Printing upg file for Uptight");

  // Recognize type of input file and print it properly
  std::string extension = path.substr(path.size()-4);
  //Gen format modified for uptight input
  if ( (extension.compare(".upg") == 0) || (extension.compare(".UPG") == 0) )
    {

      file.open(path.c_str());

      //I must build a materials map
      std::map<const Material*, unsigned int> material_map;

      std::set<ID>::iterator ID_it;

      unsigned int id = 1;
      const Material* mat = NULL;

      for (ID_it = _IDset.begin(); ID_it != _IDset.end(); ID_it ++)
        {
          mat = _device->get_material(*ID_it);
          material_map.insert(std::pair<const Material*, unsigned int>(mat, id));
          id++;
        }

      //Standard gen section (modified with material index)
      file << _structure_atoms.size();

      if (_atomistic_structure_options.is_periodical) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(6) << _atom_types[i];
        }
      file << std::endl;

      for (unsigned int i = 0; i < _structure_atoms.size(); i++)
        {
          unsigned int n_specie;
          for (n_specie = 0; n_specie < _atom_types.size(); n_specie++)
            {
              if (_atom_types[n_specie] == _structure_atoms[i].get_specie() ) break;
            }
          file << std::setw(10);
          if (_structure_atoms[i].get_specie() ==  Specie::H)
            {
              file << material_map[_device->get_material(_structure_atoms[get_bond_map()[i][0]].get_region_ID()) ];
            }
          else
            {
              file << material_map[ (_device->get_material(_structure_atoms[i].get_region_ID())) ];
            }

          file << std::setw(5) << n_specie + 1
              << std::setw(20) << std::setprecision(10)
          << std::fixed << double(_structure_atoms[i].get_position(0))
          << std::setw(20) << std::setprecision(10)
          << std::fixed  << double(_structure_atoms[i].get_position(1))
          << std::setw(20) << std::setprecision(10)
          << std::fixed  << double(_structure_atoms[i].get_position(2));


          if (_bondmap != NULL)
            {
              file << std::setw(5) << _bondmap->get_bond_map()[i].size();

              // N.B. Indexing is in Fortran notation (first atom is labelled as 1) !!!!!!!!!!!
              for (unsigned int j = 0; j < _bondmap->get_bond_map()[i].size(); j++)
                {
                  file << std::setw(10) << _bondmap->get_bond_map()[i][j] + 1;
                }
              ///////////////////////////////////////////

            }
          file << std::endl;

        }

      // Periodicity vectors at the bottom
      if (_atomistic_structure_options.is_periodical)
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
                      _periodicity_vectors[count];
                  count++;
                }
              file << std::endl;
            }
        }

      //Information about materials

      file << "#Materials " << std::endl;

      for (std::map<const Material*, unsigned int>::iterator mat_it = material_map.begin(); mat_it != material_map.end(); mat_it++)
        {
          //Note: material_map has material in ascending order, starting from one, despite any general enumeration
          //in TiberCAD. Check some lines above how it's built. In file they need to be stored in ascending order,
          //that's why I'm not using safer iterators and I work in this way

          //  std::map<Material*, unsigned int>::iterator mat_it = material_map.begin();
          ////Select material with this id in material_map
          //  for (mat_it = material_map.begin(); mat_it != material_map.end(); mat_it++)
          //  {
          //    if ((*mat_it).second == i) break;
          //  }
          const Material* mat = (*mat_it).first;
          Database db = mat->get_database();
          db.set_section("atomistic_structure");

          std::string alloy_type;

          //TODO: IT NEEDS TO BE EXTENDED WORKING WITH OTHER ALLOYS (QUATERNARY...)
          if (mat->is_alloy())  {alloy_type = "ternary";}
          else if (db.get("n_basis_specie", 0) == 1) {alloy_type = "simple";}
          else if (db.get("n_basis_specie", 0) == 2 ) {alloy_type = "binary";}
          else Messages::error("Could not define alloy_type variable in AtomisticStructure.C");

          //Mancano da inserire i file con i dati per Uptight
          std::string path = "./ ";
          std::string structure = "unknown";

          if (! mat->is_alloy())
            {
              file << std::setw(3)  << (*mat_it).second
                   << std::setw(12) << mat->get_name()
                   << std::setw(6)  << mat->get_structure()
                   << std::setw(12) << alloy_type
                   << std::setw(4)  << 1
                   << std::setw(5)  << "CRY"
                   << std::setw(8)  << mat->get_name()
                   << std::setw(10) << std::setprecision(3) << 1.0
                   << std::setw(10) << mat->get_name() << etb_dataset + ".etb"
                   << " 0.0  0.0"   << std::endl;
            };
          if (mat->is_alloy() && !(is_random_alloy()))
            {
              file << std::setw(3)  << (*mat_it).second
                   << std::setw(12) << mat->get_name()
                   << std::setw(6)  << mat->get_structure()
                   << std::setw(12) << alloy_type
                   << std::setw(4)  << 2
                   << std::setw(5)  << "VCA"
                   << std::setw(8)  << (static_cast<const Alloy*>(mat))->get_name_A()
                   << std::setw(8)  << (static_cast<const Alloy*>(mat))->get_name_B()
                   << std::setw(10) <<  std::setprecision(3)  << mat->get_options().get_option("x",1.0)
                                    << std::setw(10) << std::setprecision(3)
                                    <<  ( 1.0 - mat->get_options().get_option("x",1.0) )
                   << std::setw(10) << (static_cast<const Alloy*>(mat))->get_name_A()
                                    << etb_dataset + ".etb"
                   << std::setw(10) << (static_cast<const Alloy*>(mat))->get_name_B()
                                    << etb_dataset + ".etb"
                                    << " 0.0  0.0"
                                    << " 0.0  0.0" << std::endl;
            };
          if (mat->is_alloy() && (is_random_alloy()))
                      {
              file << std::setw(3)  << (*mat_it).second
                                 << std::setw(12) << mat->get_name()
                                 << std::setw(6)  << mat->get_structure()
                                 << std::setw(12) << alloy_type
                                 << std::setw(4)  << 2
                                 << std::setw(5)  << "RND"
                                 << std::setw(8)  << (static_cast<const Alloy*>(mat))->get_name_A()
                                 << std::setw(8)  << (static_cast<const Alloy*>(mat))->get_name_B()
                                 << std::setw(10) <<  std::setprecision(3)  << mat->get_options().get_option("x",1.0)
                                                  << std::setw(10) << std::setprecision(3)
                                                  <<  ( 1.0 - mat->get_options().get_option("x",1.0) )
                                 << std::setw(10) << (static_cast<const Alloy*>(mat))->get_name_A()
                                                  << etb_dataset + ".etb"
                                 << std::setw(10) << (static_cast<const Alloy*>(mat))->get_name_B()
                                                  << etb_dataset + ".etb"
                                                  << " 0.0  0.0"
                                                  << " 0.0  0.0" << std::endl;
                      };

        }

      file.close();
    }
  else
    {
      Messages::warning("File extension does not correspond to any internal format. File not print.");
    }

  Messages::debug("upg file printed");

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
      file << _structure_atoms.size() << std::endl << std::endl;

      for (unsigned int i = 0; i < _structure_atoms.size(); i++)
        {
          file << std::setw(2) << _structure_atoms[i].get_specie()
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(0))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(1))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(2))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(charges[i]) << "\n";
        }

    }

  else if ( (extension.compare(".gen") == 0) || (extension.compare(".GEN") == 0) )
    {
      file << _structure_atoms.size();

      if (_atomistic_structure_options.is_periodical) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(4) << _atom_types[i];
        }
      file << std::endl;

      for (unsigned int i = 0; i < _structure_atoms.size(); i++)
        {
          unsigned int n_specie;
          for (n_specie = 0; n_specie < _atom_types.size(); n_specie++)
            {
              if ( _atom_types[n_specie] == _structure_atoms[i].get_specie() ) break;
            }
          file << std::setw(10) << i + 1 << std::setw(5) << n_specie + 1
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(0))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(1))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].get_position(2))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(charges[i]) << "\n";
        }

      //A line of zeros is put here (coordinates origin)
      file <<  std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) << "\n";

      // Periodicity vectors at the bottom
      if (_atomistic_structure_options.is_periodical)
        {
          unsigned int count = 0;
          for (unsigned int i = 0; i < 3; i++)
            {
              for (unsigned int j = 0; j < 3; j++)
                {
                  file << std::setw(20) << std::setprecision(10) << std::fixed <<
                      _periodicity_vectors[count];
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
AtomisticStructure::set_periodicity_vectors(const Tensor2Gen& T)
{

  unsigned int count = 0;

  for (int i = 0; i < 3 ; i++)
    {
      for (int j = 0; j < 3 ; j++)
        {
          _periodicity_vectors[count] = T(j+1,i+1);
          count++;
        }
    }
}

void
AtomisticStructure::set_atom_types(const std::set<std::string>& atom_types)
{

  for (std::set<std::string>::iterator types = atom_types.begin(); types != atom_types.end(); types++)
    {
      _atom_types.push_back( *types );
    }

}


int
AtomisticStructure::get_type_index(const std::string& type)
{
  int result = 0;
  for (int i = 0; i < N_types; i++){
      if ( (type.compare( _atom_types[i] ) == 0) ) result = i + 1;
  }

  return result;

}

//TODO: not allocating arrays could be too slow, find a way to
//implement some memory reservation
void
AtomisticStructure::build_elem_to_atoms(void)
{
  //Get information from Atom objects
  for (unsigned int i = 0; i < _structure_atoms.size(); i++)
    {
      if (_structure_atoms[i].get_elem() != NULL)
        {
          _elem_to_atoms[_structure_atoms[i].get_elem()].push_back(i);
        }
    }
}


void
AtomisticStructure::compute_N_without_H(void)
{
  unsigned int N = 0;
  for (unsigned int i = 0; i < _structure_atoms.size(); i++)
    {
      if (_structure_atoms[i].get_specie() != Specie::H)
        {
          N++;
        }
    }
  _N_without_H = N;
}


//! Get atom Material
const Material*
AtomisticStructure::get_material(const Atom& atom, bool parent) const
{
 const Material* mat = get_device()->get_material(atom.get_region_ID());

 if (!(mat->is_alloy()))
    {
     return mat;
    }

 if (!(parent))
    {
     return mat;
    }

 if (parent)
   {
     const Alloy* alloy = dynamic_cast<const Alloy*>(mat);
     //TODO: At first we only consider InGaX, AlGaX alloys with this hardcoded trick.
     //To be generalized.
     //if ((atom.get_specie() == Specie::Ga))
     //  {
     //    return alloy->get_component_B();
     //  }
     //else return alloy->get_component_A();
     //
     if (alloy->get_component_A()->has_specie(atom.get_specie()) &&
         (!alloy->get_component_B()->has_specie(atom.get_specie())))
       return alloy->get_component_A();
     else if (alloy->get_component_B()->has_specie(atom.get_specie()) &&
         (!alloy->get_component_A()->has_specie(atom.get_specie())))       
       return alloy->get_component_B();
     else
       throw RuntimeException("Ambiguity for alloy component assignation"
           "in AtomisticStructure::get_material(Atom&, bool)");

   }


 return NULL;

}

//! Get atom Material
const Material*
AtomisticStructure::get_material(const Atom& atom1, const Atom& atom2,
    bool parent) const
{
 const Material* mat1 = get_device()->get_material(atom1.get_region_ID());
 const Material* mat2 = get_device()->get_material(atom2.get_region_ID());

 if (atom1.get_specie() == Specie::In || atom1.get_specie() == Specie::Al
     || atom1.get_specie() == Specie::Ga)
   return get_material(atom1, parent);
 else if (atom2.get_specie() == Specie::In || atom2.get_specie() == Specie::Al
     || atom2.get_specie() == Specie::Ga)
   return get_material(atom2, parent);
 else
   throw RuntimeException("WARNING: material for couple of atoms is decided "
       "depending on the anion specie. I cannot find a valid anion "
       "(only In, Al, Ga supported)");

}
