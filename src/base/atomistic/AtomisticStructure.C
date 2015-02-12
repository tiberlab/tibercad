// $Id$

#include "AtomisticStructure.h"
#include "AtomisticGenerator.h"
#include "AtomisticGenerator1D.h"
#include "AtomisticGenerator2D.h"
#include "AtomisticGenerator3D.h"
#include "BondMap.h"
#include "Messages.h"
#include "MeshUtils.h"
#include "Utils.h"
#include "TiberCad.h"
#include "Material.h"
#include "DataOutput.h"
#include "RuntimeException.h"

#include "mesh_tetgen_support.h"
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
//---------------------

using namespace std;

map<string, list<boost::function<void(void)>>>
AtomisticStructure::_callback_functions;



AtomisticStructure::AtomisticStructure()
:_device(NULL),
 _random_alloy(false),
 _clustering(false)
{
  N_atoms = 0;
}


AtomisticStructure::AtomisticStructure(const std::string& name)
:_name(name),
 _atomistic_structure_options(),
 _scale(1.0),
 _device(NULL),
 _random_alloy(false),
 _N_without_H(0)
{
  // Default initializations
  N_atoms = 0;
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

AtomisticStructure*
AtomisticStructure::create(const AtomisticStructure& as)
{
  AtomisticStructure* st = NULL;
  st = new AtomisticStructure();

  *st = as;

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


  if (_options.find_option("load_structure") || _options.find_option("load"))
    {
      std::string filename;

      if (_options.find_option("load_structure")) filename = _options.get_option("load_structure","none");
      if (_options.find_option("load")) filename = _options.get_option("load","none");

      //------------------------------------------------------------
      os << "Reading structure from file " << filename <<
          ". Any further information will be neglected" << std::endl;
      Messages::info(os.str(), true);
      os.str(std::string());
      //---------------------------------------------------------------
      Messages::info("Reading "+filename);
      init(filename);

      Messages::info("Parse regions");
      parse_regions();

      //TODO: I'm calculating the bond map again anyway because otherwise the translation vectors are not
      //correctly reproduced for periodic structures. If we really want to import the bond map,
      //this needs to be changed
      if (_bondmap != NULL)
      {
        delete _bondmap;
        _bondmap = NULL;
      }

      Messages::info("Build Bond Map");
      build_bond_map();

      Messages::info("Associate elements");
      if (_atomistic_structure_options.is_associated == false) associate_elements();

      Messages::info("Output structure(s)");
      print_driver();

    }
  else if (_options.find_option("regions"))
    {
    // Build mesh based representation
      Utils::Timer tt;
      tt.reset();
      Messages m;
      m.indent();
      init_mesh_structure();
      Messages::info("Output structure(s)");
      print_driver();
      Messages::newline();
      Messages::info("Atomistic structure build time: "+tt.elapsed_string());


      double cutoff = _options.get_option("control_volume_radius", 0.5);
      string cutoff_str = _options.get_option("control_volume_radius", "0.5");

      if (_options.get_option("extract_alloy_statistics", false))
      {
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

          map<Specie, vector<unsigned int>> stats;
          extract_statistics(stats, reg_ids, cutoff);

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

          for (int i = 0; i < NN; ++i)
          {
//            if (sums[i] >= ceil(0.95*max))
            if (sums[i] >= (max - 1))
            {
              for (mit = stats.begin(); mit != mend; ++mit)
                of << (mit->second)[i] << " ";
              of << "\n";
            }
          }

        }
      }

      if (_options.get_option("plot_alloy_composition", false))
      {

        AutoPtr<UnstructuredMesh> mesh(new Mesh(3));

        int ref_atom = _options.get_option("reference_atom", -1);
        IDSet refatoms;
        if (ref_atom >= 0)
          refatoms.insert(ref_atom);
        create_conformal_grid(*mesh, refatoms, true);

        AutoPtr<DataOutput> writer(DataOutput::create(
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
            cerr << "Atom : " << atom_types[i] << endl;

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




        map<Specie, vector<unsigned int>> stats;
        extract_statistics(stats, reg_ids, cutoff);

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

    }

  else
    {
      //------------------------------------------------------------
      os << "Error in Atomistic structure " << _name <<
          "input block. You need to define regions or load from file." << std::endl;
      Messages::info(os.str(), true);
      os.str(std::string());
      //---------------------------------------------------------------
    }

  //Calculate the number of atoms excluding hydrogens 
  //(Useful for passivated semiconductors)
  compute_N_without_H();

  Messages m;
  m.newline();
  m.indent();
  os << "Atomistic Structure containing " << N_atoms << 
        " atoms has been built. " <<std::endl;
  os << "Size not counting passivation hydrogens: "<< get_N_without_H()<<std::endl;
  m.info(os.str());
  os.str(std::string());

  m.info("Supercell structure");
  m.info("Lattice vectors (A):");
  m.indent();

  RealVectorValue a, b, c;
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
       _lattice_vectors[0] = lattice_vectors[0];
       _lattice_vectors[4] = lattice_vectors[1];
       _lattice_vectors[8] = lattice_vectors[2];
    }
 
    else if (lattice_vectors.size() == 9)
    {
       for (int i = 0; i < 9; i++)
       {
         _lattice_vectors[i] = lattice_vectors[i]; 
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


  // set a default for the periodicity
  if (!this->get_options().find_option("periodicity"))
  {
    switch (this->get_device()->get_mesh().mesh_dimension())
    {
      case 1:
        this->get_options().set_option("periodicity", "(0, 1, 1)");
        break;

      case 2:
        this->get_options().set_option("periodicity", "(0, 0, 1)");
        break;

      default:
        break;
    }
  }

  vector<bool> periodicity;
  this->get_options().get_option("periodicity", periodicity);
  this->set_periodicity(periodicity);


  parse_regions();

  //--------------------------------------------------------------
  Messages::info("Building Atomistic Structure " + get_name());
  //-----------------------------------------------------------

  _atoms.clear();
  _elem_to_atoms.clear();

  //---------------------------------------------------------------
  // Extend mesh for contacts 
  //unsigned int num_sides;
  //Point normal = Device::get_normal();
  //---------------------------------------------------------------

  AtomisticGenerator* generator 
    = AtomisticGenerator::create(this, _device->get_mesh().mesh_dimension());

  generator->do_init();
 
  generator->finalize();

  parse_lattice_vectors();

  for (size_t i = 0; i < _atoms.size(); ++i)
  {
    _elem_to_atoms[_atoms[i].get_elem()].push_back(i);
  }
 
  Messages::info("Build final Bond Map...");
  build_bond_map();

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

  delete generator;

}


void
AtomisticStructure::restrict(const std::set<ID>& rgn_ids) 
{
  _IDset = rgn_ids;


  AtomisticGenerator* generator 
    = AtomisticGenerator::create(this, _device->get_mesh().mesh_dimension());

  generator->restrict();

  generator->finalize();

  Messages::info("Build final Bond Map...");
  build_bond_map();
  
  compute_N_without_H();

  delete generator;

  std::cout<< "(AS) size: " << N_atoms <<std::endl; 
  std::cout<< "(AS) without H: " << _N_without_H <<std::endl; 
}



void
AtomisticStructure::associate_elements()
{
  // the tensor grid to real mesh mapper for fast association atom->Elem
  // NOTE: we pass the relevant ID set, since Quantum contacts could be present
  //       which have to be included in the atomistic structure
  MeshUtils::GridMapper& mapper =
      MeshUtils::GridMapper::get_mapper(_device->get_mesh(), _IDset);


  unsigned int dim =  _device->get_mesh().mesh_dimension();

  size_t n_atoms = _atoms.size();
  Utils::Progress prog("Assign elements", n_atoms);
  unsigned int progress = 0;

  // NOTE: Hydrogens remains outside regions and are not associated to elements
  // BondMap has not been created yet

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
  std::cout << std::endl;

  /*
  if (_bondmap != NULL)
  {
    const BondMap& bondmap = *_bondmap;

    for (unsigned int i=0; i < _atoms.size(); i++)
    {
      if ( _atoms[i].get_specie() == Specie::H )
      {
        for (unsigned int j=0; j < bondmap[i].size(); j++)
        { 
          const Elem* elem = _atoms[bondmap[i][j]].get_elem();
          if (elem != NULL) _atoms[i].set_elem(elem);
        }        
      }      
    }    
  }
  */

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
AtomisticStructure::read_structure(const std::string& path)
{

  std::ifstream file;
  std::string line, record;
  unsigned int n_specie;
  Atom tmp_atom;
  Tensor1 pos;

  Messages::debug("Reading structure from file "+path);

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

  // XYZ file
  if ( (extension.compare(".xyz") == 0) || (extension.compare(".XYZ") == 0) )
    {
      // First line is number of atoms
      getline(file, line);
      N_atoms = atoi(line.c_str());
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

          tmp_atom.set_position( pos );

          _atoms.push_back(tmp_atom);

        }

      if ( (_atoms.size()) != N_atoms )
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
      _atoms.reserve(N_atoms);
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
        set_periodicity({1, 1, 1});
      else  if ( (record.compare("C") == 0) || (record.compare("c") == 0))
        set_periodicity({0, 0, 0});
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
    }

  else if ( (extension.compare(".tgn") == 0) || (extension.compare(".TGN") == 0) )
    {
      read_tgn(path);
    }

  else
    {
      Messages::error("Structure file extension is not recognized. ");
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
    set_periodicity({1, 1, 1});
  else  if ( (record.compare("C") == 0) || (record.compare("c") == 0))
    set_periodicity({0, 0, 0});
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


  N_types = _atom_types.size();


  // Cycle upon specified number of atoms (last rows are for periodicity vectors)
  for (unsigned int i = 1; i <= N_atoms; i++)
    {
      getline(file, line);

      line_string.clear();
      line_string.str(line);

      // First value is ignored (contains physical region) 
      line_string >> record;
      line_string >> record;
      n_specie = atoi(record.c_str());

      for (unsigned int j = 1; j <= 3; j++)
        {
          line_string >> record;
          pos(j) =  atof(record.c_str());
        }
      _atoms[i - 1].set_specie(_atom_types[n_specie -1]);
      _atoms[i - 1].set_position(pos);

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

      line_string >> record;
      _atoms[i - 1].set_label(atoi(record.c_str()));
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

  unsigned int id = 1;
  const Material* mat = NULL;

  for (ID_it = _IDset.begin(); ID_it != _IDset.end(); ID_it ++)
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
          
          if (_atoms[i].get_specie() ==  Specie::H)
            {
              file << material_map[_device->get_material(_atoms[get_bond_map()[i][0]].get_region_ID()) ];
            }
          else
            {
              file << material_map[ (_device->get_material(_atoms[i].get_region_ID())) ];
            }

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
          if (_atoms[i].get_elem() == NULL) file << -1;
          else file << _atoms[i].get_elem()->id();

          // Print atom label 
          file << std::setw(3);
          file<< static_cast<unsigned int>(_atoms[i].get_label());

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
      print_upg(file_name, "");
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
AtomisticStructure::print_upg(const std::string& path, const std::string& etb_dataset,
                                                       bool band_offsets)
{

  std::ofstream file, os;
  const BondMap& bondmap = *_bondmap;

  Messages::debug("Printing upg file for Uptight");

  // Recognize type of input file and print it properly
  std::string extension = path.substr(path.size()-4);
  //Gen format modified for uptight input
  if ( (extension.compare(".upg") == 0) || (extension.compare(".UPG") == 0) )
    {

      file.open(path.c_str());

      //I must build a materials map with new enumeration
      std::map<const Material*, ID> material_map;

      //std::set<ID>::iterator ID_it;

      const Material* mat = NULL;

      //for (ID_it = _IDset.begin(); ID_it != _IDset.end(); ID_it ++)
      //  {
      //    mat = _device->get_material(*ID_it);
      //    material_map.insert(std::pair<const Material*, unsigned int>(mat, id));
      //    id++;
      //  }

      //Standard gen section (modified with material index)
      file << _atoms.size();

      if (is_periodic()) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(6) << _atom_types[i];
        }
      file << std::endl;
      
      // we need to know which regions have been used actually, otherwise
      // we will write too many materials at the end.
      unsigned int id = 1;

      for (unsigned int i = 0; i < _atoms.size(); i++)
        {
          unsigned int n_specie;
          for (n_specie = 0; n_specie < _atom_types.size(); n_specie++)
            {
              if (_atom_types[n_specie] == _atoms[i].get_specie() ) break;
            }
          file << std::setw(10);

          // add to the map if not present
          if (!material_map.count(_device->get_material(_atoms[i].get_region_ID())))
          {
            material_map.insert(std::make_pair(_device->get_material(_atoms[i].get_region_ID()), id));
            ++id;
          }

          ID reg_id = material_map[(_device->get_material(_atoms[i].get_region_ID()))];
          file << reg_id;
          //if (_atoms[i].get_specie() ==  Specie::H)
          //  {
          //    file << material_map[_device->get_material(_atoms[get_bond_map()[i][0]].get_region_ID()) ];
          //  }
          //else
          //  {
          //    file << material_map[ (_device->get_material(_atoms[i].get_region_ID())) ];
          //  }

          file << std::setw(5) << n_specie + 1
              << std::setw(20) << std::setprecision(10)
          << std::fixed << double(_atoms[i].get_position(0))
          << std::setw(20) << std::setprecision(10)
          << std::fixed  << double(_atoms[i].get_position(1))
          << std::setw(20) << std::setprecision(10)
          << std::fixed  << double(_atoms[i].get_position(2));


          if (_bondmap != NULL)
            {
              file << std::setw(5) << bondmap[i].size();

              // N.B. Indexing is in Fortran notation (first atom is labelled as 1) !!!!!!!!!!!
              for (unsigned int j = 0; j < bondmap[i].size(); j++)
                {
                  file << std::setw(10) << bondmap[i][j] + 1;
                }
              ///////////////////////////////////////////

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
          std::vector<double> Ev(2,0.0);

          //TODO: IT NEEDS TO BE EXTENDED WORKING WITH OTHER ALLOYS (QUATERNARY...)
          if (mat->is_alloy())  {alloy_type = "ternary";}
          else if (db.get("n_basis_specie", 0) == 1) {alloy_type = "simple";}
          else if (db.get("n_basis_specie", 0) == 2 ) {alloy_type = "binary";}
          else Messages::error("Could not define alloy_type variable in AtomisticStructure.C");

          if (alloy_type=="ternary")
          {
  	     const Material* mat1 = (static_cast<const Alloy*>(mat))->get_component_A();	  
             Database db1 = mat1->get_database();
	     db1.set_section("valenceband");
	     Ev[0]=db1.get("E_v",0.0);
  	     mat1 = (static_cast<const Alloy*>(mat))->get_component_B();	  
             Database db2 = mat1->get_database();
	     db2.set_section("valenceband");
	     Ev[1]=db2.get("E_v",0.0);
          }
	  else
	  { 
             db.set_section("valenceband");
	     Ev[0]=db.get("E_v",0.0);
	  }

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
                   << std::setw(10) << mat->get_name() << etb_dataset + ".etb";
              if (band_offsets)
              {               
		   file << std::setw(10) << Ev[0];
              } 
              file << " 0.0  0.0"   << std::endl;
          }

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
                                    << etb_dataset + ".etb";
              if (band_offsets)
              {               
                file << std::setw(10) << Ev[0] << " "
                     << std::setw(10) << Ev[1] << " ";
              }
              file << " 0.0  0.0"
                   << " 0.0  0.0" << std::endl;
          }

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
                                                  << etb_dataset + ".etb";
              if (band_offsets)
              {               
                  file           << std::setw(10) << Ev[0] << " "
                                 << std::setw(10) << Ev[1] << " ";
              }
                  file                            << " 0.0  0.0"
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
    if (_atoms[i].get_specie() != Specie::H)
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
     const Alloy* alloy = dynamic_cast<const Alloy*>(mat);
     if (alloy->get_component_A()->has_specie(atom.get_specie()) &&
         (!alloy->get_component_B()->has_specie(atom.get_specie())))
       mat = alloy->get_component_A();
     else if (alloy->get_component_B()->has_specie(atom.get_specie()) &&
         (!alloy->get_component_A()->has_specie(atom.get_specie())))       
       mat = alloy->get_component_B();
     else
     {
       Messages::error("Ambiguity for alloy component assignation"
           "in AtomisticStructure::get_material(Atom&, bool)");
     }
   }

 return mat;

}

//! Get atom Material
const Material*
AtomisticStructure::get_material(const Atom& atom1, const Atom& atom2,
    bool parent) const
{
 //const Material* mat1 = get_device()->get_material(atom1.get_region_ID());
 //const Material* mat2 = get_device()->get_material(atom2.get_region_ID());

 //If not, we need to decide based on some other criteria. Up to now we're able to
 //decide only for III-V or II-VI alloys with different cations (eg. Ga-As belong to GaAs)
 //if (mat1->is_cation(atom1.get_specie()))
 //  return get_material(atom1, parent);
 //else if (mat2->is_cation(atom2.get_specie()))
 //  return get_material(atom2, parent);
 //else
 //(Alex) Test using the label on atoms rather than hardcoded CrystalDefs. 
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
 return NULL;
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
AtomisticStructure::create_conformal_grid(UnstructuredMesh& mesh,
    set<unsigned int> labels, bool keep_node_order) const
{

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

    //if (atom.get_specie() == Specie::In ||
    //    atom.get_specie() == Specie::Al ||
    //    atom.get_specie() == Specie::Ga) 
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

}





void
AtomisticStructure::extract_statistics(map<Specie, vector<unsigned int>>& stats,
    const set<ID>& regions, double cutoff) const
{

  cutoff = _options.get_option("control_volume_radius", cutoff);

  int ref_atom = _options.get_option("reference_atom", -1);

  Messages m;
  m.info("Extracting alloy statistics for structure " + this->get_name() + ": ");
  m.indent();
  ostringstream os;
  os << "control sphere radius: " << cutoff << " nm";
  m.info(os.str());
  if (ref_atom >= 0)
  {
    os.str("");
    os << "control volumes centered at atoms with label (= reference atom): " << ref_atom;
    m.info(os.str());
  }

  for (unsigned int i = 0; i < get_N_atoms(); i++)
  {
    const Atom& atm = get_structure_atom(i);
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

      stats[atm.get_specie()][0] += 1;

      if ((ref_atom < 0) || (atm.get_label() == ref_atom))
      {
        unsigned int counter = 0;
        unsigned int total = 1;

        map<Specie, vector<unsigned int>>::iterator mit(stats.begin());
        map<Specie, vector<unsigned int>>::iterator mend(stats.end());
        for ( ; mit != mend; ++mit)
          (mit->second).push_back(0);

        neighbor_iterator it(neighbors_begin(i, 10 * cutoff));
        neighbor_iterator end(neighbors_end(i));
        for ( ; it != end; ++it)
        {
          const Atom& neigh = *(*it);

          if (!regions.count(neigh.get_region_ID()))
            continue;

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
AtomisticStructure::find_nearest_atom(const Elem* elem, const Point& point, double cutoff)
{
  ID id = elem->subdomain_id();

  // atomic coordinates are in Angstrom
  Point center = 10 * point;

  //if (reg_ids.count(id))
  //{

  const vector<unsigned int>& atoms = get_atoms_in_elem(elem);
  int atom = -1;

  //
  // it may not even contain atoms
  if (!atoms.empty())
  {
    // ok, there is at least one atom inside
    double min_dist = 1e9;
    unsigned int nearest = 0;

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
        nearest = i;
      }
    }
    atom = atoms[nearest];
  }
  else
  {
    // we try to find some nearby atom
    set<const Elem*> processed_elems;
    set<const Elem*> to_process;
    to_process.insert(elem);

    unsigned int nearest = 0;
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
            nearest = i;
          }
        }
        atom = atoms[nearest];
      }

      for (int s = 0; s < next_el->n_sides(); s++)
      {

        const Elem* neigh = next_el->neighbor(s);

        if ((neigh != NULL) &&
            !processed_elems.count(neigh) &&
            (Point(elem->centroid() - neigh->centroid()).size() <
                _scale * 2 * cutoff / 10))
        {
          to_process.insert(neigh);
        }
      }

      processed_elems.insert(next_el);
      to_process.erase(it);
    }
  }

  return atom;
}
