#include "BondMap.h"
#include "GridCells.h"
#include "Messages.h"
#include "Database.h"
#include "Messages.h"
#include <fstream>
#include "Utils.h"


BondMap::BondMap(unsigned int structure_size, unsigned int valence)
:_period(0)
{

  (*this).resize(structure_size);
  for (unsigned int i = 0; i < structure_size; i++)
  {
    (*this)[i].reserve(valence);
  }

  _translation.resize(structure_size);
  for (unsigned int i = 0; i < structure_size; i++)
  {
    _translation[i].reserve(valence);
  }

}



BondMap::~BondMap(void)
{
}


void
BondMap::clean()
{
  //_grid_cell.clear();
}

void
BondMap::set_cutoff()
{
  std::ifstream file;
  std::string database_path = Database::get_default_search_path();
  std::string filename = "cutoff.dat";
  std::string line, record;
  std::stringstream line_stream;
  Specie s;

  filename = database_path + "/" + filename;

    
  file.open(filename.c_str());
  if (!file)
  {
    Messages::error("Exception opening/reading file cutoff.dat");
    exit(1);
  }
  
  while (!file.eof())
  {
    getline(file, line);
    line_stream.str(std::string());
    line_stream.clear(std::stringstream::goodbit);
    line_stream << line;
    line_stream >> record;
    s = record;
    line_stream >> record;
    
    _cutoff[s] = atof(record.c_str());
    
  }
  
  file.close();
}



void
BondMap::do_solve(const std::vector<Atom>& basis, const Tensor2Gen& period)
{

  Messages::debug("BondMap::do_solve");

  Tensor1 edge_min, edge_max;

  //set cutoff distancies
  set_cutoff();
  
  for (unsigned int i=0; i< basis.size(); i++)
  {
    if (_cutoff[basis[i].get_specie()] == 0.0)
    {
      Messages::error("Cutoff distance for specie "
                      + basis[i].get_specie().get_string()
                      +" is not defined. Add to materials/cutoff.dat");
    }
  }
 
  _period = period;
  //define the minimum spacing of the grid. the smaller it is, the faster is bonds calculations
  //cannot be smaller than the higher bond lenght. (in amstrong)
  GridCells cells(basis, period, 8.0);

  // Loop on all cells
  Utils::Progress prog("BondMap", cells.size());

  for (unsigned int c1 = 0; c1 < cells.size(); c1++)
  {
    GridCells::NeighborIterator it = cells.begin(c1);
    GridCells::NeighborIterator end = cells.end(c1);
 
    prog.progress_message(c1+1);

    //std::cout << "\b\b\b\b\b\b\b\b" << std::setw(3) 
    //          << static_cast<int>(100*(c1+1)/cells.size())<<"% ..."<<std::flush;

    // Loop on all 27 neighboring cells (periodicity is taken care by the iterator)
    for ( ; it != end; ++it)
    {
      unsigned int c2 = (*it).first;
      const Tensor1& shift = *((*it).second);

      // Loop over all atoms in each cell
      for (unsigned int i = 0; i < cells[c1].size(); i++)
      {     
        //if (!basis[i].belong_to_structure) cycle
        for (unsigned int j = 0; j < cells[c2].size(); j++)
        {
          //if (!basis[j].belong_to_structure) cycle
          process_atoms(basis, cells[c1][i], cells[c2][j], shift);
        }  
      }
      
    }
  }
}


void
BondMap::process_atoms(const std::vector<Atom>& basis, 
                       const unsigned int i,
                       const unsigned int j, 
                       const Tensor1& period)
{

  unsigned int put_here = 0;
  bool not_already_signed;
  double cutofftmp;
  Tensor1 position1, position2;
  libMesh::Point per(period(1), period(2), period(3));

  if ((i != j) || (per.norm() > 1e-5))
  {
   
    cutofftmp = _cutoff[basis[i].get_specie()] + _cutoff[basis[j].get_specie()];
    position1 = basis[i].get_ttype_position();
    position2 = basis[j].get_ttype_position() + period;
  
    //std::cout<<"atoms "<<i<<" "<<j<<" d="<<norm(position1-position2)<<std::endl; 
    if ( norm( position1 - position2) < cutofftmp )
    {
      
      not_already_signed = true;
      put_here = (*this)[i].size();
      
      if (put_here != 0)
      {
        for (unsigned int n = 0; n < put_here; n++){
          
          if ((*this)[i][n] == j && libMesh::Point(_translation[i][n] - per).norm() < 1e-6)
          {
            not_already_signed = false;
            break;
          }
        }
      }
      else
        not_already_signed = true;

      if (not_already_signed)
      {
        (*this)[i].push_back(j);
        _translation[i].push_back(libMesh::Point(period(1), period(2), period(3)));
      }
      
    }

  }
  
}

void
BondMap::print(const std::vector<Atom>& basis)
{
  for (unsigned int i = 0; i < basis.size(); i++)
  {
    std::cout<<"BondMap["<<i<<"]=";
    for (unsigned int j = 0; j < (*this)[i].size(); j++)
      std::cout<<(*this)[i][j]<<" ";

    std::cout<<std::endl;
  }
}




