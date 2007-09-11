
#include "AtomisticStructure.h"



AtomisticStructure* 
AtomisticStructure::create(const std::string& name)
{
  AtomisticStructure* st =  NULL;
  st = new AtomisticStructure(name);

  std::cout << "Created AtomisticStructure " << st->get_name() << std::endl;
 

  return st;
}



AtomisticStructure* 
AtomisticStructure::create(const std::string& name, const ModelOptions& options)
{
  AtomisticStructure* st = create(name);

  st->set_options(options);

  return st;
}


void
AtomisticStructure::init(void)
{

#ifdef DEBUG
  std::cerr << "AtomisticStructure::init() begin \n";
#endif

  assert(_device != NULL);

  // Read structure from file
  std::string path;

  if (! _options.find_option("path") )
    std::cerr << "ERROR IN ATOMISTIC REGION DEFINITION: A PATH FOR STRUCTURE FILE MUST BE SPECIFIED" << std::endl;

  path = _options.get_option("path","none");
  
  if (path.compare("none") != 0) read_structure(path);


  if ( _options.find_option("physical_regions") )
    {
   
      //Put physical regions specified in input file in _regions
      //A vector is needed as temporary container
      std::vector<std::string> region_string; 
      _options.get_option("physical_regions", region_string);

      for (int i = 0; i < region_string.size(); i++)
	{std::cerr << "Assigning in set region_string " << region_string[i] << std::endl;
	  _region.insert(region_string[i]);}
      region_string.clear();

      //If all regions are specified (value = "all", must fill with real names of all regions)
      if ( _region.count("all") == 1)
	{
	  _region.clear();
	  std::set < ID >::iterator region_ID_iterator = _device->get_region_ids().begin();

	  for (int i = 0; i < _device->get_region_ids().size(); i++)
	    {
	      _region.insert( _device->get_region_name(*region_ID_iterator) );
	      region_ID_iterator ++;
	    }

	}

    }


  for (std::set<std::string>::iterator i= _region.begin(); i !=_region.end(); i++)
    {std::cerr << "WRITING " << std::endl;
      std::cerr << "_REGION IS " << *i << std::endl;}

#ifdef DEBUG
  std::cerr << "AtomisticStructure::init() end \n";
#endif

}


const std::set<std::string>& 
AtomisticStructure::get_region(void)
{
  return _region;
}


void 
AtomisticStructure::set_device(Device* device)
{
  _device = device;
}


void 
AtomisticStructure::read_structure(const std::string& path)
{
  std::ifstream file;
  std::string line, record;
  unsigned int N_atoms, n_specie;
  atom tmp_atom;

#ifdef DEBUG
  std::cerr << "AtomisticStructure::read_structure(path) begin \n";
#endif

  // Delete eventually existing structure
  if (!(_structure_atoms.empty())) _structure_atoms.clear();
  if (!(_atom_types.empty())) _atom_types.clear();


  file.open(path.c_str(), std::ifstream::in);

  if (!file) 
    {
      std::cerr << "Unable to open file " << path << " \n";
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

      _N_atoms = N_atoms;

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
		  if ( record.compare(_atom_types[i]) == 0) not_present = true;
		}

	      if (not_present) _atom_types.push_back(record);
	    }
	  else
	    {
	      _atom_types.push_back(record);
	    }

	  tmp_atom.specie = record;

	  for (unsigned int i = 1; i == 3; i++)
	    { 

	      line_string >> record;
	      tmp_atom.position(i) = atof(record.c_str());
	    }

	  _structure_atoms.push_back(tmp_atom);

	}

      if ( (_structure_atoms.size()) != _N_atoms ) std::cerr << "Warning: in file xyz number of atoms is wrong \n";

      // Warning: XYZ file has no informations about structure periodicity

    }

  // GEN file
  else if ( (extension.compare(".gen") == 0) || (extension.compare(".GEN") == 0) )

    {
      getline(file, line);
      std::stringstream line_string(line);
      
      line_string >> record;

      N_atoms = atoi(line.c_str());

#ifdef DEBUG
      std::cerr << "N_atoms is " << N_atoms << std::endl;
#endif

      if (N_atoms == 0) 
	{
	  std::cerr << "No atoms in structure files or non valid integer in first line. \n";
	  exit(1);
	}

      _N_atoms = N_atoms;

      line_string >> record;

      if ( (record.compare("S") == 0) || (record.compare("s") == 0)) 
	_is_periodical = true;
      else  if ( (record.compare("C") == 0) && (record.compare("c") == 0))
	std::cerr << "Warning (in GEN file at first line): Cluster (C) or Supercell (S) must be specified. By default a Cluster (no periodicity) is considered. \n";

      getline(file, line);

      line_string.str("");
      line_string << line;
            
      // I assum that GEN file is correct and there's no name repetition
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
	  tmp_atom.specie = _atom_types[n_specie -1];

	  for (unsigned int j = 1; j <= 3; j++)
	    { 

	      line_string >> record;
	      tmp_atom.position(j) = atof(record.c_str());
	    }

	  _structure_atoms.push_back(tmp_atom);
	}

      // An additional line is present in GEN files. It's the coordinates origin and it's
      // not needed
      getline(file, line);

      // If GEN refers to periodical structure I expect periodicity vectors
      if (_is_periodical)
	{
	  for (unsigned int i = 0; i < 3; i++)
	    {
	      getline(file, line);
	      std::stringstream line_string(line);
	      for (unsigned int j = 0; j < 3; j++)
		{
		  line_string >> record;
		  _periodicity_vectors[i][j] = atof(record.c_str());
		}  
	    }
	}
    }

  else 
    {
      std::cerr << "Structure file extension is not recognized. \n";
      exit(1);
    }

  file.close();

#ifdef DEBUG
  std::cerr << "AtomisticStructure::read_structure(path) end. \n";
#endif

}


void 
AtomisticStructure::print_structure(const std::string& path)
{
  std::ofstream file;

#ifdef DEBUG
  std::cerr << "AtomisticStructure::print_structure(path) begin. \n";
#endif

  file.open(path.c_str());

  // Recognize type of input file and print it properly
  std::string extension = path.substr(path.size()-4);

if ( (extension.compare(".xyz") == 0) || (extension.compare(".XYZ") == 0) )
  {
    file << _structure_atoms.size() << std::endl << std::endl;

    for (unsigned int i = 0; i < _structure_atoms.size(); i++)
      {
	file << std::setw(2) << _structure_atoms[i].specie 
	   << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].position(1)) 
	   << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].position(2)) 
	   << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].position(3)) << "\n"; 

      }

  }

 else if ( (extension.compare(".gen") == 0) || (extension.compare(".GEN") == 0) )
   {
     file << _structure_atoms.size();

     if (_is_periodical) file << std::setw(10) << "S \n";
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
	     if (_atom_types[n_specie].compare(_structure_atoms[i].specie) == 0) break;
	   }
	 file << std::setw(10) << i + 1 << std::setw(5) << n_specie + 1
	     << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].position(1)) 
	   << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].position(2)) 
	   << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_structure_atoms[i].position(3)) << "\n"; 
       }

     //A line of zeros is put here (coordinates origin)
     file <<  std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) 
	   << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) 
	   << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) << "\n"; 

     // Periodicity vectors at the bottom
      if (_is_periodical)
	{
	  for (unsigned int i = 0; i < 3; i++)
	    {
	      for (unsigned int j = 0; j < 3; j++)
		{
		  file << std::setw(20) << std::setprecision(10) << std::fixed <<
		    _periodicity_vectors[i][j];
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

#ifdef DEBUG
 std::cerr << "AtomisticStructure::print_structure(path) end. \n";
#endif

}


