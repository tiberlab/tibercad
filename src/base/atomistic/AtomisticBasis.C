// $Id$

#include "AtomisticBasis.h"
#include "Device.h"
#include "ModelOptions.h"

#include "libmesh/vector_value.h"


//STD library includes
#include<iostream>
#include<fstream>
#include<sstream>
//-------------------


using namespace std;



AtomisticBasis::~AtomisticBasis(void)
{
  if (_bondmap != nullptr) delete _bondmap;
}

AtomisticBasis::AtomisticBasis(void)
:_bondmap(nullptr),
 _lattice_vectors(3,0.0),
 _periodicity({0, 0, 0}),
 _origin(0)
{

}


AtomisticBasis::AtomisticBasis(const AtomisticBasis& other) :
  _bondmap(nullptr),
  _atoms(other._atoms),
  _lattice_vectors(other._lattice_vectors),
  _atom_types(other._atom_types),
  _periodicity(other._periodicity),
  _origin(other._origin)
{
  _bondmap = new BondMap(_atoms.size());
  *_bondmap = *(other._bondmap);
}


void
AtomisticBasis::get_lattice_vectors(libMesh::RealVectorValue& a,
                                    libMesh::RealVectorValue& b,
                                    libMesh::RealVectorValue& c) const
{
  a = _lattice_vectors[0];
  b = _lattice_vectors[1];
  c = _lattice_vectors[2];
}


void
AtomisticBasis::set_lattice_vectors(const libMesh::RealVectorValue& a,
                                    const libMesh::RealVectorValue& b,
                                    const libMesh::RealVectorValue& c)
{
  _lattice_vectors[0] = a;
  _lattice_vectors[1] = b;
  _lattice_vectors[2] = c;

  if (_bondmap != nullptr)
  {
    _bondmap->set_periodicity(a, b, c);
  }
}



void
AtomisticBasis::set_ttype_lattice_vectors(const Tensor2Gen& T)
{
  libMesh::RealVectorValue a(T(1,1), T(2,1), T(3,1));
  libMesh::RealVectorValue b(T(1,2), T(2,2), T(3,2));
  libMesh::RealVectorValue c(T(1,3), T(2,3), T(3,3));

  set_lattice_vectors(a, b, c);
}

Tensor2Gen 
AtomisticBasis::get_ttype_lattice_vectors(void)
{
  Tensor2Gen T;
  for (int j = 0; j < 3 ; j++)
  {
    T(j+1,1) = _lattice_vectors[0](j);
    T(j+1,2) = _lattice_vectors[1](j);
    T(j+1,3) = _lattice_vectors[2](j);
  }

  return T;
}

int
AtomisticBasis::get_type_index(const std::string& type) const
{
  int result = 0;
  for (int i = 0; i < _atom_types.size(); i++){
      if ( (type.compare( _atom_types[i] ) == 0) ) result = i + 1;
  }

  return result;

}


void
AtomisticBasis::set_atom_types(const std::set<std::string>& atom_types)
{

  for (std::set<std::string>::iterator types = atom_types.begin(); types != atom_types.end(); types++)
  {
    _atom_types.push_back( *types );
  }
}


BondMap*
AtomisticBasis::build_bond_map(bool periodicity[3]) const
{
  BondMap* bm = new BondMap(_atoms.size());

  Tensor2Gen period;
  for (unsigned int i = 0; i < 3; i++)
  {
    double scale = periodicity[i] ? 1 : 10;
    for (unsigned int j = 0; j < 3; j++)
    {
      period(j + 1, i + 1) = scale * _lattice_vectors[i](j);
    }
  }
  bm->do_solve(_atoms, period, _origin);

  return bm;
}


void
AtomisticBasis::build_bond_map(void)
{
  if (_bondmap != NULL)
  {
    delete _bondmap;
  }
  _bondmap = new BondMap(_atoms.size());

  Tensor2Gen period;
  for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
        {
          period(j + 1, i + 1) = _lattice_vectors[i](j);
        }
    }
  _bondmap->do_solve(_atoms, period, _origin);

}



void
AtomisticBasis::set_bond_map(const BondMap& bondmap)
{
  if (_bondmap == nullptr)
    _bondmap = new BondMap(_atoms.size());

  *_bondmap = bondmap;
}


std::pair<Point, Point>
AtomisticBasis::get_bounding_box(void) const
{
  Point pmax(std::numeric_limits<double>::min(),
      std::numeric_limits<double>::min(),
      std::numeric_limits<double>::min());
  Point pmin(std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max());

  for (auto&& at : _atoms)
  {

    const Point& p = at.get_position();
    if (p(0) < pmin(0))
      pmin(0) = p(0);
    else if (p(0) > pmax(0))
      pmax(0) = p(0);

    if (p(1) < pmin(1))
      pmin(1) = p(1);
    else if (p(1) > pmax(1))
      pmax(1) = p(1);

    if (p(2) < pmin(2))
      pmin(2) = p(2);
    else if (p(2) > pmax(2))
      pmax(2) = p(2);

  }
  return(make_pair(pmin, pmax));
}


void
AtomisticBasis::refresh(void)
{

  std::set<std::string> types;
  for (std::vector<Atom>::iterator it = _atoms.begin(); it != _atoms.end(); ++it)
  {
    types.insert(it->get_specie().get_string());
  } 
  clear_atom_types();
  set_atom_types(types);

  //build_bond_map();
}

//-------------------------------------------------------------------
//Print utilities
//-------------------------------------------------------------------

void
AtomisticBasis::print_xyz(const std::string& path) const
{
  std::ofstream file;
  // -------------------------------------------

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/" + path;

  // --------------------------------------------
  file.open(file_name.c_str());
  file << _atoms.size() << std::endl << std::endl;
  for (unsigned int i = 0; i < _atoms.size(); i++)
  {
    file << std::setw(2) << _atoms[i].get_specie()
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2)) << "\n";
  }
  file.close();
}

void
AtomisticBasis::print_xyb(const std::string& path) const
{

  std::ofstream file;
  // -------------------------------------------

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/" + path;

  // --------------------------------------------
  file.open(file_name.c_str());
  file << _atoms.size() << std::endl << std::endl;

      for (unsigned int i = 0; i < _atoms.size(); i++)
        {
          file << std::setw(2) << _atoms[i].get_specie()
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
          << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2));

          if (_bondmap != NULL)
            {

              file << std::setw(5) << (*_bondmap)[i].size();

              // N.B. Indexing is in Fortran notation (first atom is labelled as 1) !!!!!!!!!!!!!!!!
              for (unsigned int j = 0; j < (*_bondmap)[i].size(); j++)
                {
                  file << std::setw(10) << (*_bondmap)[i][j] + 1;
                }
              ///////////////////////////////////////////

            }

          file << std::endl;
        }
  file.close();
}


void
AtomisticBasis::print_gen(const std::string& path) const
{

  std::ofstream file;
  // -------------------------------------------

  std::string outdir = TiberCad::get_output_dir();
  std::string file_name = outdir + "/" + path;

  // --------------------------------------------
  file.open(file_name.c_str());
     file << _atoms.size();

      if (this->is_periodic()) file << std::setw(10) << "S \n";
      else file << std::setw(10) << "C \n";

      for (unsigned int i = 0; i < _atom_types.size(); i++)
        {
          file << std::setw(6) << _atom_types[i];
        }
      file << std::endl;

      for (unsigned int i = 0; i < _atoms.size(); i++)
        {
          std::string strsp = _atoms[i].get_specie().get_string();
          file << std::setw(10) << i + 1 << std::setw(5) << get_type_index(strsp) 
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(0))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(1))
              << std::setw(20) << std::setprecision(10)<< std::fixed  << double(_atoms[i].get_position(2)) << "\n";
        }

      //A line of zeros is put here (coordinates origin)
      file <<  std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0)
        << std::setw(20) << std::setprecision(10)<< std::fixed  << double(0.0) << "\n";

      // Periodicity vectors at the bottom
      for (unsigned int i = 0; i < 3; i++)
        {
          for (unsigned int j = 0; j < 3; j++)
            {
              file << std::setw(20) << std::setprecision(10) << std::fixed <<
                  _lattice_vectors[i](j);
            }
          file << "\n";
        }
  file.close();

}



void
AtomisticBasis::get_subset(vector<unsigned int>& subset,
      const ModelOptions& opt)
{

  // length units in options are assumed in nm

  Point center(0);
  opt.get_option("center", center);
  center *= 10;

  // define the volume
  string shape = opt.get_option("shape", "sphere");
  double x = opt.get_option("radius", 0.5);
  double y = 0;
  double z = 0;

  if (shape == "cube")
  {
    x = opt.get_option("side_length", x);
    y = z = x;
  }
  else if (shape == "box")
  {
    Point sides(0);
    opt.get_option("sides", sides);
    x = sides(0);
    y = sides(1);
    z = sides(2);
  }
  else if (shape == "column")
  {
    libMesh::RealVectorValue a, b, c;
    get_lattice_vectors(a, b, c);
    string dir = opt.get_option("orientation", "x");
    double side = opt.get_option("side_length", x);
    switch (dir[0])
    {

      case 'y':
        y = b.size();
        x = z = side;
        break;

      case 'z':
        z = c.size();
        y = x = side;
        break;

      default: // 'x'
        x = a.size();
        y = z = side;
        break;
    }
  }

  unsigned int center_atom = 0;
  double min_dist = std::numeric_limits<double>::max();

  // look for nearest atom
  for (unsigned int i = 0; i < get_N_atoms(); i++)
  {
    const Atom& atm = get_structure_atom(i);

    // we skip atoms with no "real" label (such as passivation H)
    if (atm.get_label() == 0)
      continue;

    Point d(atm.get_position() - center);
    double d_len = d.size();
    if (d_len < min_dist)
    {
      min_dist = d_len;
      center_atom = i;
    }
  }

  neighbor_iterator it(neighbors_begin(center_atom, 10 * x, 10 * y, 10 * z));
  neighbor_iterator end(neighbors_end(center_atom));
  for ( ; it != end; ++it)
  {
    //const Atom& neigh = *(*it);
    if ((*it)->get_label() != 0)
      subset.push_back(it.atom_index());
  }

}




AtomisticBasis::neighbor_iterator::neighbor_iterator(const AtomisticBasis& structure,
    unsigned int start, double length, double height, double width, bool begin)
: _structure(structure),
  _start(start),
  _current(start),
  _length(length),
  _height(height),
  _width(width),
  _image(0)
{
  if (begin)
  {
    _visited.insert(make_pair(_start, Point(0)));
    const std::vector<unsigned int>& curr = _structure.get_bond_map()[_current];
    for (unsigned int i = 0; i < curr.size(); i++)
      _setA.insert(make_pair(curr[i],
          _structure.get_bond_map().get_translation(_current, i)));

    _itA = _setA.begin();
    //cerr << _start << endl;
  }
  else
  {
    _current = _structure.get_N_atoms();
  }
}



AtomisticBasis::neighbor_iterator::neighbor_iterator(const neighbor_iterator& rhs)
: _structure(rhs._structure),
  _start(rhs._start),
  _current(rhs._current),
  _length(rhs._length),
  _height(rhs._height),
  _width(rhs._width),
  _visited(rhs._visited),
  _setA(rhs._setA),
  _setB(rhs._setB),
  _itA(rhs._itA)
{}


AtomisticBasis::neighbor_iterator&
AtomisticBasis::neighbor_iterator::operator++(void)
{
  // we will look for neighbours of atoms slightly outside the
  // desired range, to not loose valid atoms
  bool advance = true;
  while (advance)
  {
    // if we processed all in set A, pass to their neighbours
    if (_itA == _setA.end())
    {
      advance = false;
      if (_setB.empty())
      {
        _current = _structure.get_N_atoms();
      }
      else
      {
        _setA.swap(_setB);
        _itA = _setA.begin();
        _setB.clear();

        ++(*this);
      }
    }
    else
    {
      // look for all his neighbors

      _current = _itA->first;
      _image = (_itA->second);
      _visited.insert(make_pair(_current, _structure.get_structure_atom(_start).get_position()
          - _structure.get_structure_atom(_current).get_position()
          - _image));

      const vector<unsigned int>& nn = _structure.get_bond_map()[_current];
      for (unsigned int i = 0; i < nn.size(); ++i)
      {
        auto range(_visited.equal_range(nn[i]));


        // this neighbour would be shifted by this amount, if periodic copy
        Point new_shift(
            _structure.get_bond_map().get_translation(_current, i) + _image);

        Point dist(_structure.get_structure_atom(_start).get_position()
                    - _structure.get_structure_atom(nn[i]).get_position()
                    - new_shift);

        // did we already use this (periodic) atom?
        bool visited = false;
        HMMap::iterator it(range.first);
        for ( ; it != range.second; ++it)
          if ((it->second).relative_fuzzy_equals(dist, 1e-9))
          {
            visited = true;
            break;
          }

        //if you use only length you consider a sphere with radius equal to length
        //if you define length, height and width you consider a parallelepiped
        if (!visited)
        {
          if( _height == 0 && _width == 0 )
          {
            double d = dist.norm();

            if (d < _length * 1.5)
              _setB.insert(make_pair(nn[i], new_shift));
          }
          else
          {
            double dx  = fabs(dist(0));
            double dy  = fabs(dist(1));
            double dz  = fabs(dist(2));

            if (dx < (_length/2 + _min_bond) &&
                dy < (_height/2 + _min_bond) &&
                dz < (_width/2 + _min_bond))
              _setB.insert(make_pair(nn[i], new_shift)) ;
          }
        }

      }
      ++_itA;

      // check if _current is actually inside the desired domain. If not,
      // pass to the next one
      if( _height == 0 && _width == 0 )
      {
        double d = Point(_structure.get_structure_atom(_start).get_position()
            - _structure.get_structure_atom(_current).get_position()
            - _image).norm();

        advance = (d > _length) || (d < 1e-3);
      }
      else
      {
        double dx  = fabs(_structure.get_structure_atom(_start).get_position(0)
            - _structure.get_structure_atom(_current).get_position(0) - _image(0));
        double dy  = fabs(_structure.get_structure_atom(_start).get_position(1)
            - _structure.get_structure_atom(_current).get_position(1) - _image(1));
        double dz  = fabs(_structure.get_structure_atom(_start).get_position(2)
            - _structure.get_structure_atom(_current).get_position(2) - _image(2));

        advance = (dx > (_length/2)) || (dy > (_height/2)) || (dz > (_width/2));
      }
    }
  }

  return *this;
}

AtomisticBasis::neighbor_iterator
AtomisticBasis::neighbors_begin(unsigned int index, double length, double height, double width) const
{
  return neighbor_iterator(*this, index, length, height, width);
}

AtomisticBasis::neighbor_iterator
AtomisticBasis::neighbors_end(unsigned int index, double length, double height, double width) const
{
  return neighbor_iterator(*this, index, length, height, width, false);
}




