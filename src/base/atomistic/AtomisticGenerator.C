#include "AtomisticGenerator.h"
#include "AtomisticStructure.h"
#include "AtomisticGenerator1D.h"
#include "AtomisticGenerator2D.h"
#include "AtomisticGenerator3D.h"
#include "Macrostrain.h"
#include "BondMap.h"

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>



AtomisticGenerator::AtomisticGenerator(void)
:_bondmapobject(NULL),
_reference_material(NULL)
{

}

AtomisticGenerator::~AtomisticGenerator(void)
{
  if (_bondmapobject != NULL) delete _bondmapobject;
}


const double AtomisticGenerator::tol = 1e-2;

AtomisticGenerator*
AtomisticGenerator::create(AtomisticStructure* const as, unsigned int dimension)
{
  AtomisticGenerator* ag =  NULL;
  if (dimension == 1)  ag = AtomisticGenerator1D::create(as);
  if (dimension == 2) ag = AtomisticGenerator2D::create(as);
  if (dimension == 3) ag = AtomisticGenerator3D::create(as);

  return ag;
}


void AtomisticGenerator::print_basis(std::vector<Atom> &basis, const std::string filename){

  std::ofstream output_file;

  std::vector<Atom>::iterator basis_iterator = basis.begin();

  //#ifdef DEBUG
  //  std::cerr << "Printing structure to file...";
  //#endif

  output_file.open(filename.c_str());
  output_file << basis.size() << std::endl << std::endl;

  do{

    output_file << std::setw(2) << (*basis_iterator).get_specie()
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).get_position(1))
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).get_position(2))
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).get_position(3)) << "\n";

    basis_iterator++;

  }while(basis_iterator != basis.end());

  output_file.close();

  //#ifdef DEBUG
  //  std::cerr << "done" << std::endl;
  //#endif

};


void
AtomisticGenerator::do_init()
{

  std::cout << "Building Atomistic Structure " << _as->get_name() << std::endl;

  //Set dimensional scale
  scale = _as->get_scale();

  // Set material informations
  //-----------------------------------------------------------------------------------------
  std::string structure;
  structure = "none";

  if (!(_as->get_options().find_option("reference_region"))){
    std::cerr << "No material could be set: reference_region is mandatory in Atomistic section"
    "when no structure path is specified " << std::endl;}
  std::vector<ID> ids;
  std::string ref_region;
  ref_region = _as->get_options().get_option("reference_region", "None");
  _as->get_device()->get_region_ids(ref_region, ids);
  _reference_material = _as->get_device()->get_material(ids[0]);

  structure =  _reference_material->get_structure();
  std::cout << "Parsing atomistic structure parameters... " << std::endl;
  parse_parameters(_reference_material);


  //Set Growth direction informations
  Tensor2Gen miller(1);
  miller(2,2) = 0.0;
  miller(3,3) = 0.0;

  std::vector<int> growth_direction;

  if (! ( (_reference_material->get_options().find_option("x-growth-direction")) || (_as->get_options().find_option("x-growth-direction")) ) )
    std::cerr << "Warning, no x-growth-direction is set for atomistic structure " << _as->get_name() << " Setting (1,0,0) as default " << std::endl;

  if (_reference_material != NULL){
    if (_reference_material->get_options().find_option("x-growth-direction")){
      _reference_material->get_options().get_option("x-growth-direction", growth_direction);
      miller(1,1) = growth_direction[0]; miller(2,1) = growth_direction[1]; miller(3,1) = growth_direction[2];
      if (growth_direction.size() == 4) miller(3,1) = growth_direction[3];
    }
    if (_reference_material->get_options().find_option("y-growth-direction")){
      _reference_material->get_options().get_option("y-growth-direction", growth_direction);
      miller(1,2) = growth_direction[0]; miller(2,2) = growth_direction[1]; miller(3,2) = growth_direction[2];
      if (growth_direction.size() == 4) miller(3,2) = growth_direction[3];
    }
    if (_reference_material->get_options().find_option("z-growth-direction")){
      _reference_material->get_options().get_option("z-growth-direction", growth_direction);
      miller(1,3) = growth_direction[0]; miller(2,3) = growth_direction[1]; miller(3,3) = growth_direction[2];
      if (growth_direction.size() == 4) miller(3,3) = growth_direction[3];
    }
  }

  //If Miller indexes specified in Atomistic options, take them
  if (_as->get_options().find_option("x-growth-direction")){
    _reference_material->get_options().get_option("x-growth-direction", growth_direction);
    miller(1,1) = growth_direction[0]; miller(2,1) = growth_direction[1]; miller(3,1) = growth_direction[2];
    if (growth_direction.size() == 4) miller(3,1) = growth_direction[3];
  }
  if (_as->get_options().find_option("y-growth-direction")){
    _reference_material->get_options().get_option("y-growth-direction", growth_direction);
    miller(1,2) = growth_direction[0]; miller(2,2) = growth_direction[1]; miller(3,2) = growth_direction[2];
    if (growth_direction.size() == 4) miller(3,2) = growth_direction[3];
  }
  if (_as->get_options().find_option("z-growth-direction")){
    _reference_material->get_options().get_option("z-growth-direction", growth_direction);
    miller(1,3) = growth_direction[0]; miller(2,3) = growth_direction[1]; miller(3,3) = growth_direction[2];
    if (growth_direction.size() == 4) miller(3,3) = growth_direction[3];
  }

  set_prim_miller(miller);

  //-------------------------------------------------------------------------------------------


  // Set the vector of elements covered by structure regions, useful for change specie and cut
  MeshBase::element_iterator el = _as->get_device()->get_mesh().elements_begin();
  const MeshBase::element_iterator el_end = _as->get_device()->get_mesh().elements_end();

  for ( ; el != el_end; el++)
    {           Elem* elem = *el;
    if (_as->get_IDset().find( elem->subdomain_id() ) != _as->get_IDset().end() ) _structure_elements.push_back(elem);
    }


  //Build up supercell structure with proper options
  //----------------------------------------------------------------------------------------------
  build();


  std::string preserve;
  preserve = _as->get_options().get_option("preserve", "none");
  cut_and_change_specie(preserve);

  std::string passivation;
  passivation = _as->get_options().get_option("passivation", "no");

  if (passivation.compare("yes") == 0){passivate();}

  //BondMap pointer is used for passivation, delete it and refresh bond map
  delete _bondmapobject;

  _bondmapobject = NULL;

  //Eliminate not included atoms from structure
  //Check structure to eliminate unincluded atoms (using swap in another vector)
  //-----------------------------------------------------------
  std::vector<Atom> tmp_structure;
  tmp_structure.reserve(_structure_basis.size());

  for (unsigned int i = 0; i < _structure_basis.size(); i++)
    {
      if ((_structure_basis[i].belong_to_structure))
        {
          tmp_structure.push_back(_structure_basis[i]);
        }
    }

  _structure_basis.clear();
  _structure_basis.reserve(tmp_structure.size());
  _structure_basis.swap(tmp_structure);
  //-------------------------------------------------------------

  bond_map_gen(_structure_basis);

  //----------------------------------------------------------------------------------------------




  //Pass data to AtomisticStructure
  //--------------------------------------------------------------------------------------------------

  Atom tmp_atom;
  //TODO:not safe, better swap arrays, so then we can delete AtomisticGenerator instance
  _as->set_structure_atoms(_structure_basis);
  //for (int i = 0; i < 3 ; i++){
  //  for (int j = 0; j < 3 ; j++){
  //    _as->_periodicity_vectors[i][j] = _period(i+1,j+1);
  //  }
  // }
  _as->set_periodicity_vectors(_period);

  _as->set_N_atoms( _structure_basis.size() );

  std::set<std::string> atom_types;

  for (unsigned int i = 0; i < _structure_basis.size(); i++)
    {
      atom_types.insert(_structure_basis[i].get_specie());
    }
  _as->set_N_types ( atom_types.size() );

  _as->clear_atom_types();

  //for (std::set<std::string>::iterator types = atom_types.begin(); types != atom_types.end(); types++)
  //{
  //  _as->_atom_types.push_back(*types);
  //}
  _as->set_atom_types(atom_types);

  //Passing bondmap object (deallocation will be managed by AtomisticStructure, setting local pointer NULL)
  _as->set_bondmap(_bondmapobject);
  _bondmapobject = NULL;


  //#ifdef DEBUG
  //  std::cout << "Ending AtomisticGenerator::do_init() " << std::endl;
  //#endif

};



void
AtomisticGenerator::cut_and_change_specie(std::string preserve){

  //#ifdef DEBUG
  //  std::cerr << "Cutting atoms and changing species...";
  //#endif

  std::set<ID> IDs = _as->get_IDset();
  //std::map<unsigned int , std::string> assign;
  std::map<ID, std::map<unsigned int, std::string> > assign;
  bool done;
  ID el_reg;

  _structure_basis.clear();
  assign.clear();

  for (std::set<ID>::iterator reg = _as->get_IDset().begin(); reg != _as->get_IDset().end(); reg++)
    {

      Material* mat = _as->get_device()->get_material( (*reg) );

      Database& db = mat->get_database();

      db.set_section("atomistic_structure");

      //Build up conversion map from file
      for (int i = 1; i <= db.get("n_basis_specie", 0); i++)
        {
          std::string record;
          std::string s;
          std::stringstream out;

          out << i;
          s = out.str();

          record = "specie_" + s;
          assign[*reg][i] = db.get(record.c_str(),"none");

        }

      //No more reading from section atomistic_structure in database are needed
      db.set_section("");

      //If some doping is present, a strategy must be studied and implemented here,
      // as we can choose arbitrarily species name (e.g. calling one doping Silicon Si1 and another one Si2)

    }

  //Cycle upon all atoms and change specie according to assign map
  Point p(0.0, 0.0, 0.0);

  //Different strategies if preserving conventional cell or preserving basis are needed
  if (preserve.compare("none") == 0)
    {

      //bool not_already_included is needed because a point can be contained by more than one element
      //if it falls exactly on the boundary

      for ( std::vector<Atom>::iterator atom = _super_basis.begin();
      atom != _super_basis.end(); atom++)
        {
          done = false;

          p(0) = (*atom).get_position(1) / scale;

          if ( (_dim == 2) || (_dim == 3) )   p(1) = (*atom).get_position(2) / scale;

          if ( (_dim == 3) )  p(2) = (*atom).get_position(3) / scale;

          for (std::vector<Elem*>::iterator it = _structure_elements.begin();
          it != _structure_elements.end(); it++)
            {

              Elem* elem = *it;
              el_reg = elem->subdomain_id();

              if (Macrostrain::may_belong_to_element(elem, p)){

                if ( (elem->contains_point(p) ) ) {
                  if ( assign[el_reg].find( (*atom).get_flag() ) != assign[el_reg].end() )
                    {
                      std::string tmp =  assign[el_reg][(*atom).get_flag()];
                      (*atom).set_specie(tmp);
                      (*atom).set_flag(0);
                      (*atom).set_region_ID(el_reg);
                      (*atom).belong_to_structure = true;
                      _structure_basis.push_back(*atom);
                      done = true;
                    }
                  else
                    {
                      std::cout << "Warning, atom is included but no assignament map member could be built " << std::endl;
                    }
                  break;
                }

              }


            }

          //If atom does not belong to element, push it and flag it as non internal
          if (done == false)
            {
              if ( assign[el_reg].find( (*atom).get_flag() ) != assign[el_reg].end() )
                {
                  std::string tmp =  assign[el_reg][(*atom).get_flag()];
                  (*atom).set_specie(tmp);
                  (*atom).set_flag(0);
                  (*atom).set_region_ID(el_reg);
                  (*atom).belong_to_structure = false;
                  _structure_basis.push_back(*atom);
                }
              else
                {
                  std::cout << "Warning, atom is included but no assignment map member could be built " << std::endl;
                }
            }

        }
    }


  if (preserve.compare("lattice") == 0)
    {

      Atom tmp_atom;
      for ( std::vector<Tensor1>::iterator lattice = _super_lattice.begin();
      lattice != _super_lattice.end(); lattice++)
        {
          done = false;

          p(0) = (*lattice)(1) / scale;

          if ( (_dim == 2) || (_dim == 3) )   p(1) = (*lattice)(2) / scale;

          if ( (_dim == 3) )  p(2) = (*lattice)(3) / scale;

          for (std::vector<Elem*>::iterator it = _structure_elements.begin();
          it != _structure_elements.end(); it++)
            {

              Elem* elem = *it;
              el_reg = elem->subdomain_id();

              if (Macrostrain::may_belong_to_element(elem,p))
                {

                  if ( (elem->contains_point(p) ) )
                    {
                      done = true;
                      for ( std::vector<Atom>::iterator atom = _crystal_basis.begin();
                      atom != _crystal_basis.end(); atom++)
                        {

                          tmp_atom.set_position( (*lattice)+ _rotation*_prim_vec*(*atom).get_position() );
                          tmp_atom.set_region_ID( el_reg );

                          if ( assign[el_reg].find( (*atom).get_flag() ) != assign[el_reg].end() )
                            {
                              std::string tmp =  assign[el_reg][(*atom).get_flag()];
                              tmp_atom.set_specie(tmp);
                              tmp_atom.belong_to_structure = true;
                              _structure_basis.push_back(tmp_atom);
                            }
                          else
                            {
                              std::cout << "Warning, atom is included but no assignment map member could be built " << std::endl;
                            }
                        }
                      break;
                    }

                }


            }
          if (done == false)
            {
              for ( std::vector<Atom>::iterator atom = _crystal_basis.begin();
              atom != _crystal_basis.end(); atom++)
                {

                  tmp_atom.set_position( (*lattice)+ _rotation*_prim_vec*(*atom).get_position() );
                  tmp_atom.set_region_ID( el_reg );

                  if ( assign[el_reg].find( (*atom).get_flag() ) != assign[el_reg].end() )
                    {
                      std::string tmp =  assign[el_reg][(*atom).get_flag()];
                      tmp_atom.set_specie(tmp);
                      tmp_atom.belong_to_structure = false;
                      _structure_basis.push_back(tmp_atom);
                    }
                  else
                    {
                      std::cout << "Warning, atom is included but no assignment map member could be built " << std::endl;
                    }
                }
            }
        }
    }

  if (preserve.compare("conventional") == 0)
    {

      Atom tmp_atom;
      for ( std::vector<Tensor1>::iterator conv = _super_conv.begin(); conv != _super_conv.end(); conv++)
        {

          done = false;

          for (std::vector<Elem*>::iterator it = _structure_elements.begin(); it != _structure_elements.end(); it++)
            {

              p(0) = (*conv)(1) / scale;

              if ( (_dim == 2) || (_dim == 3) )   p(1) = (*conv)(2) / scale;

              if ( (_dim == 3) )  p(2) = (*conv)(3) / scale;

              Elem* elem = *it;
              el_reg = elem->subdomain_id();


              if (Macrostrain::may_belong_to_element(elem,p))
                {

                  if ( (elem->contains_point(p) ) )
                    {
                      done = true;

                      for ( std::vector<Tensor1>::iterator conv_lattice_basis_it = _conv_lattice_basis.begin();
                      conv_lattice_basis_it != _conv_lattice_basis.end(); conv_lattice_basis_it++)
                        {
                          for ( std::vector<Atom>::iterator atom = _crystal_basis.begin();
                          atom != _crystal_basis.end(); atom++)
                            {
                              tmp_atom.set_position( (*conv) + (*conv_lattice_basis_it) +
                                  _rotation*_prim_vec*(*atom).get_position());
                              tmp_atom.set_region_ID( el_reg );
                              if ( assign[el_reg].find( (*atom).get_flag() ) != assign[el_reg].end() )
                                {
                                  std::string tmp =  assign[el_reg][(*atom).get_flag()];
                                  tmp_atom.set_specie(tmp);
                                  tmp_atom.belong_to_structure = true;
                                  _structure_basis.push_back(tmp_atom);
                                }
                              else
                                {
                                  std::cout << "Warning, atom is included but no assignment map member could be built " << std::endl;
                                }
                            }
                        }
                      break;
                    }

                }



            }

          if (done == false)
            {
              for ( std::vector<Tensor1>::iterator conv_lattice_basis_it = _conv_lattice_basis.begin();
              conv_lattice_basis_it != _conv_lattice_basis.end(); conv_lattice_basis_it++)
                {
                  for ( std::vector<Atom>::iterator atom = _crystal_basis.begin();
                  atom != _crystal_basis.end(); atom++)
                    {
                      tmp_atom.set_position( (*conv) + (*conv_lattice_basis_it) +
                          _rotation*_prim_vec*(*atom).get_position());
                      tmp_atom.set_region_ID( el_reg );
                      if ( assign[el_reg].find( (*atom).get_flag() ) != assign[el_reg].end() )
                        {
                          std::string tmp =  assign[el_reg][(*atom).get_flag()];
                          tmp_atom.set_specie(tmp);
                          tmp_atom.belong_to_structure = false;
                          _structure_basis.push_back(tmp_atom);
                        }
                      else
                        {
                          std::cout << "Warning, atom is included but no assignment map member could be built " << std::endl;
                        }
                    }
                }
            }


        }
    }

};


//Note:: make_supercell is called only with preserve_basis and preserve_conv
//This is the complete function after some modifications in Atom structure (added
// conventional cell address and ID). It's commented and the function is rewritten
// keeping only preserve_basis and preserve_conv instructions, so it's more readable
// and modifications are easier. If some change in strategy would occur we can come back
// to this one. Note that also header has been modifying, removing preserve_basis and
// preserve_conv
//
//void AtomisticGenerator::make_supercell(double l1, double l2, double l3, bool preserve_basis, bool preserve_conv){
//
//  //Build a supercell, defined by the lenght of conventional growth cell vectors
//  std::vector<Tensor1>::iterator conv_iterator;
//  std::vector<Atom>::iterator basis_iterator;
//  int i,j,l;
//  int n1,n2,n3;
//  double conv_l1, conv_l2, conv_l3;
//  Atom basis_atom;
//  Tensor1 lattice_point;
//  Tensor2Gen supercell_vect,inv_supercell_vect;
//  Tensor1 tmp_check, tmp_conv;
//  bool check_boundary, check_boundary2;
//
//  //#ifdef DEBUG
//  //  std::cerr << "Building a supercell sized " << l1 << " " << l2 << " " << l3 << " Amstrong" << std::endl;
//  //#endif
//
//  //Check values. l1,l2,l3 cannot be unwisely large (no more than (1um)^3)
//  assert((l1*l2*l3) < 1e+12);
//
//  //Find lenght of conventional cell sides
//  conv_l1 = sqrt(_conv_vect(1,1) * _conv_vect(1,1) + _conv_vect(2,1) * _conv_vect(2,1) + _conv_vect(3,1) * _conv_vect(3,1));
//  conv_l2 = sqrt(_conv_vect(1,2) * _conv_vect(1,2) + _conv_vect(2,2) * _conv_vect(2,2) + _conv_vect(3,2) * _conv_vect(3,2));
//  conv_l3 = sqrt(_conv_vect(1,3) * _conv_vect(1,3) + _conv_vect(2,3) * _conv_vect(2,3) + _conv_vect(3,3) * _conv_vect(3,3));
//
//  n1 = int(floor(l1 / conv_l1)); n2 = int(floor(l2 / conv_l2)); n3 = int(floor(l3 / conv_l3));
//
//  if (preserve_conv) {
//
//    l1 = (n1 + 1) * conv_l1; l2 = (n2 +1) * conv_l2; l3 = (n3 + 1) * conv_l3;
//
//  }
//
//
//  //Set supercell periodical vectors
//  Tensor2Gen lmat(0);
//
//
//  lmat(1,1) = (n1 + 1);
//  lmat(2,2) = (n2 + 1); lmat(3,3) = (n3 +1);
//  // Periodicity along x or y or z direction is set to a big value (double of structure lenght) (non periodic along x)
//  //according to dimensionality of the system
//
//  if (_dim == 1) lmat(1,1) = (n1 + 1) * 2;
//  if (_dim == 2) {lmat(1,1) = (n1 + 1) * 2; lmat(2,2) = (n2 + 1) * 2;}
//  if (_dim == 3) {lmat(1,1) = (n1 + 1) * 2; lmat(2,2) = (n2 + 1) * 2; lmat(3,3) = (n3 +1) * 2;}
//
//  _period = _conv_vect * lmat;
//
//  //Define vectors with same direction of conventional cell vectors, but with size specifed by l1,l2,l3
//  supercell_vect(1,1) = _conv_vect(1,1) * (l1 / conv_l1); supercell_vect(2,1) = _conv_vect(2,1) * (l1 / conv_l1); supercell_vect(3,1) = _conv_vect(3,1) * (l1 / conv_l1);
//  supercell_vect(1,2) = _conv_vect(1,2) * (l2 / conv_l2); supercell_vect(2,2) = _conv_vect(2,2) * (l2 / conv_l2); supercell_vect(3,2) = _conv_vect(3,2) * (l2 / conv_l2);
//  supercell_vect(1,3) = _conv_vect(1,3) * (l3 / conv_l3); supercell_vect(2,3) = _conv_vect(2,3) * (l3 / conv_l3); supercell_vect(3,3) = _conv_vect(3,3) * (l3 / conv_l3);
//  inv_supercell_vect = inv(supercell_vect);
//
//  for (i = -1 ; i <= n1 + 1 ; i++){
//    for (j = -1 ; j <= n2 + 1 ; j++){
//      for (l = -1 ; l <= n3 + 1 ; l++){
//
//	conv_iterator = _conv_lattice_basis.begin();
//
//	//Fill conventional edges basis (super_conv)
//	if ( (i != -1)&&(i <= n1 )&&(j != -1)&&(j <= n2)&& (l != -1)&&(l <= n3) ){
//	  tmp_conv(1) = (i * _conv_vect(1,1)) + (j * _conv_vect(1,2)) + (l * _conv_vect(1,3));
//	  tmp_conv(2) = (i * _conv_vect(2,1)) + (j * _conv_vect(2,2)) + (l * _conv_vect(2,3));
//	  tmp_conv(3) = (i * _conv_vect(3,1)) + (j * _conv_vect(3,2)) + (l * _conv_vect(3,3));
//	  _super_conv.push_back(tmp_conv + _local_origin);}
//
//	do{
//	  //Assign lattice point position
//	  lattice_point(1) = (*conv_iterator)(1) + (i * _conv_vect(1,1)) + (j * _conv_vect(1,2)) + (l * _conv_vect(1,3));
//	  lattice_point(2) = (*conv_iterator)(2) + (i * _conv_vect(2,1)) + (j * _conv_vect(2,2)) + (l * _conv_vect(2,3));
//	  lattice_point(3) = (*conv_iterator)(3) + (i * _conv_vect(3,1)) + (j * _conv_vect(3,2)) + (l * _conv_vect(3,3));
//
//	  if (preserve_basis){
//
//	    //Check if lattice point is inside bonduary
//	    if ((i >= n1) || (i <= 0) || (j >= n2) || (j<= 0) || (l >= n3) || (l <= 0)) {
//	      tmp_check = inv_supercell_vect * lattice_point;
//	      check_boundary = ((tmp_check(1) >= -tol) && (tmp_check(1) <(1.0 - tol))) &&
//		((tmp_check(2) >= -tol) && (tmp_check(2) < (1.0 - tol))) &&
//		((tmp_check(3) >= -tol) && (tmp_check(3) < (1.0 - tol)));
//	    }
//
//	    else (check_boundary = 1);
//
//	    if (check_boundary){
//	      //Put lattice point into supercell lattice points array
//	      _super_lattice.push_back(lattice_point + _local_origin);
//	      basis_iterator=_crystal_basis.begin();
//
//	      do{
//		basis_atom = (*basis_iterator);
//		basis_atom.set_position ( _local_origin + lattice_point+
//					  _rotation*_prim_vec*(*basis_iterator).get_position() );
//
//		_super_basis.push_back(basis_atom);
//		basis_iterator++;
//
//	      }while(basis_iterator != _crystal_basis.end());
//	    }
//	  }
//
//	  else{
//
//
//	    //Put lattice point into supercell lattice points array
//	    _super_lattice.push_back(lattice_point + _local_origin);
//
//	    basis_iterator=_crystal_basis.begin();
//
//	    do{
//	      basis_atom = (*basis_iterator);
//	      basis_atom.set_position ( _local_origin + lattice_point +
//					_rotation*_prim_vec*(*basis_iterator).get_position() );
//
//	      //Check if basis atom is inside bonduary when preserve_basis is off
//	      if ((i >= n1) || (i <= 0) || (j >= n2) || (j<= 0) || (l >= n3) || (l <= 0)) {
//		tmp_check = inv_supercell_vect * basis_atom.get_position();
//		check_boundary2 = ((tmp_check(1) >= -tol) && (tmp_check(1) <(1.0 + tol))) &&
//		  ((tmp_check(2) >= -tol) && (tmp_check(2) < (1.0 + tol))) &&
//		  ((tmp_check(3) >= -tol) && (tmp_check(3) < (1.0 + tol)));
//		if (check_boundary2) _super_basis.push_back(basis_atom);
//	      }
//	      else _super_basis.push_back(basis_atom);
//	      basis_iterator++;
//	    }while(basis_iterator != _crystal_basis.end());
//	  }
//
//
//	  conv_iterator++;
//	}while(conv_iterator != _conv_lattice_basis.end());
//
//      };
//    };
//  };
//};



void AtomisticGenerator::make_supercell(double l1, double l2, double l3){

  //Build a supercell, defined by the lenght of conventional growth cell vectors
  std::vector<Tensor1>::iterator conv_iterator;
  std::vector<Atom>::iterator basis_iterator;
  int i,j,l;
  int n1,n2,n3,start_i,start_j,start_l;
  double conv_l1, conv_l2, conv_l3;
  Atom basis_atom;
  Tensor1 lattice_point;
  Tensor2Gen supercell_vect,inv_supercell_vect;
  Tensor1 tmp_check, tmp_conv;


  //#ifdef DEBUG
  //  std::cerr << "Building a supercell sized " << l1 << " " << l2 << " " << l3 << " Amstrong" << std::endl;
  //#endif

  //Check values. l1,l2,l3 cannot be unwisely large (no more than (1um)^3)
  assert((l1*l2*l3) < 1e+12);

  //Find lenght of conventional cell sides
  conv_l1 = sqrt(_conv_vect(1,1) * _conv_vect(1,1) + _conv_vect(2,1) * _conv_vect(2,1) + _conv_vect(3,1) * _conv_vect(3,1));
  conv_l2 = sqrt(_conv_vect(1,2) * _conv_vect(1,2) + _conv_vect(2,2) * _conv_vect(2,2) + _conv_vect(3,2) * _conv_vect(3,2));
  conv_l3 = sqrt(_conv_vect(1,3) * _conv_vect(1,3) + _conv_vect(2,3) * _conv_vect(2,3) + _conv_vect(3,3) * _conv_vect(3,3));

  n1 = int(floor(l1 / conv_l1)); n2 = int(floor(l2 / conv_l2)); n3 = int(floor(l3 / conv_l3));

  _conv_cells_supercell_lenght[0] = n1 + 1;
  _conv_cells_supercell_lenght[1] = n2 + 1;
  _conv_cells_supercell_lenght[2] = n3 + 1;

  l1 = (n1 + 1) * conv_l1; l2 = (n2 +1) * conv_l2; l3 = (n3 + 1) * conv_l3;

  //Set supercell periodical vectors
  Tensor2Gen lmat(0);

  lmat(1,1) = (n1 + 1); lmat(2,2) = (n2 + 1); lmat(3,3) = (n3 +1);

  // Periodicity along x or y or z direction is set to a big value (ten times structure lenght) (non periodic along x)
  //according to dimensionality of the system

  if (_dim == 1) lmat(1,1) = (n1 + 1) * 10;
  if (_dim == 2) {lmat(1,1) = (n1 + 1) * 10; lmat(2,2) = (n2 + 1) * 10;}
  if (_dim == 3) {lmat(1,1) = (n1 + 1) * 10; lmat(2,2) = (n2 + 1) * 10; lmat(3,3) = (n3 +1) * 10;}

  _period = _conv_vect * lmat;

  //Define vectors with same direction of conventional cell vectors, but with size specifed by l1,l2,l3
  supercell_vect(1,1) = _conv_vect(1,1) * (l1 / conv_l1); supercell_vect(2,1) = _conv_vect(2,1) * (l1 / conv_l1); supercell_vect(3,1) = _conv_vect(3,1) * (l1 / conv_l1);
  supercell_vect(1,2) = _conv_vect(1,2) * (l2 / conv_l2); supercell_vect(2,2) = _conv_vect(2,2) * (l2 / conv_l2); supercell_vect(3,2) = _conv_vect(3,2) * (l2 / conv_l2);
  supercell_vect(1,3) = _conv_vect(1,3) * (l3 / conv_l3); supercell_vect(2,3) = _conv_vect(2,3) * (l3 / conv_l3); supercell_vect(3,3) = _conv_vect(3,3) * (l3 / conv_l3);
  inv_supercell_vect = inv(supercell_vect);

  //std::cout << "I'm bulding a supercell with " << n1 << n2 << n3 << "conventional cells" << std::endl;

  if (_dim == 1) {start_i = -1; start_j = 0; start_l = 0; n1 = n1 + 1;}
  if (_dim == 2) {start_i = -1; start_j =-1; start_l = 0; n1 = n1 + 1; n2 = n2 + 1;}
  if (_dim == 3) {start_i = -1; start_j =-1; start_l = -1; n1 = n1 + 1; n2 = n2 + 1; n3 = n3 + 1;}

  //Need to construct a redundant supercell (for passivation purposes)
  //Note that it must be redundant only in non periodic directions
  for (i = start_i; i <= n1; i++){
    for (j = start_j; j <= n2; j++){
      for (l = start_l; l <= n3; l++){

        conv_iterator = _conv_lattice_basis.begin();

        //Fill conventional edges basis (super_conv)
        //        if ( (i != -1)&&(i <= n1 )&&(j != -1)&&(j <= n2)&& (l != -1)&&(l <= n3) ){
        tmp_conv(1) = (i * _conv_vect(1,1)) + (j * _conv_vect(1,2)) + (l * _conv_vect(1,3));
        tmp_conv(2) = (i * _conv_vect(2,1)) + (j * _conv_vect(2,2)) + (l * _conv_vect(2,3));
        tmp_conv(3) = (i * _conv_vect(3,1)) + (j * _conv_vect(3,2)) + (l * _conv_vect(3,3));
        _super_conv.push_back(tmp_conv + _local_origin);
        //        }

        do{
          //Assign lattice point position
          lattice_point(1) = (*conv_iterator)(1) + (i * _conv_vect(1,1)) + (j * _conv_vect(1,2)) + (l * _conv_vect(1,3));
          lattice_point(2) = (*conv_iterator)(2) + (i * _conv_vect(2,1)) + (j * _conv_vect(2,2)) + (l * _conv_vect(2,3));
          lattice_point(3) = (*conv_iterator)(3) + (i * _conv_vect(3,1)) + (j * _conv_vect(3,2)) + (l * _conv_vect(3,3));

          //Put lattice point into supercell lattice points array
          _super_lattice.push_back(lattice_point + _local_origin);

          basis_iterator=_crystal_basis.begin();

          do{
            basis_atom = (*basis_iterator);
            basis_atom.set_position ( _local_origin + lattice_point+
                _rotation*_prim_vec*(*basis_iterator).get_position() );

            _super_basis.push_back(basis_atom);
            basis_iterator++;

          }while(basis_iterator != _crystal_basis.end());


          conv_iterator++;

        }while(conv_iterator != _conv_lattice_basis.end());

      };
    };
  };
};



void
AtomisticGenerator::set_lattice_type(const std::string lattice_name)
{

  Tensor2Gen prim_vec_dir;

  //Set the lattice type. It defines the primitive vectors
  _lattice_type=lattice_name;

  if (_lattice_type.compare("cubic") == 0) {

    assert((_lattice_constant[0] == _lattice_constant[1]) && (_lattice_constant[1] == _lattice_constant[2]));

    prim_vec_dir(1,1) = 1.0; prim_vec_dir(2,1) = 0; prim_vec_dir(3,1) = 0;
    prim_vec_dir(1,2) = 0; prim_vec_dir(2,2) = 1; prim_vec_dir(3,2) = 0;
    prim_vec_dir(1,3) = 0; prim_vec_dir(2,3) = 0; prim_vec_dir(3,3) = 1;

    _prim_vec = prim_vec_dir * _lattice_constant[0];

  }

  else if (_lattice_type.compare("bcc") == 0) {

    assert((_lattice_constant[0] == _lattice_constant[1]) && (_lattice_constant[1] == _lattice_constant[2]));

    prim_vec_dir(1,1) = -0.5; prim_vec_dir(2,1) = 0.5; prim_vec_dir(3,1) = 0.5;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = -0.5; prim_vec_dir(3,2) = 0.5;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = -0.5;

    _prim_vec = prim_vec_dir * _lattice_constant[0];

  }

  else if (_lattice_type.compare("fcc") == 0) {

    assert((_lattice_constant[0] == _lattice_constant[1]) && (_lattice_constant[1] == _lattice_constant[2]));

    prim_vec_dir(1,1) = 0.0; prim_vec_dir(2,1) = 0.5; prim_vec_dir(3,1) = 0.5;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = 0.0; prim_vec_dir(3,2) = 0.5;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = 0.0;

    _prim_vec = prim_vec_dir * _lattice_constant[0];

  }

  else if (_lattice_type.compare("hexagonal") == 0) {

    assert(_lattice_constant[0] == _lattice_constant[1]);

    prim_vec_dir(1,1) = 0.5; prim_vec_dir(2,1) = -sqrt(3.0) / 2.0; prim_vec_dir(3,1) = 0.0;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = sqrt(3.0) / 2.0; prim_vec_dir(3,2) = 0.0;
    prim_vec_dir(1,3) = 0.0; prim_vec_dir(2,3) = 0.0; prim_vec_dir(3,3) = 1.0;

    _prim_vec(1,1) = prim_vec_dir(1,1) * _lattice_constant[0]; _prim_vec(2,1) = prim_vec_dir(2,1) * _lattice_constant[0];
    _prim_vec(1,2) = prim_vec_dir(1,2) * _lattice_constant[0]; _prim_vec(2,2) = prim_vec_dir(2,2) * _lattice_constant[0];
    _prim_vec(1,3) = 0.0; _prim_vec(2,3) = 0.0; _prim_vec(3,3) = prim_vec_dir(3,3) * _lattice_constant[2];

  }

  else if (_lattice_type.compare("anatase") == 0) {

    assert(_lattice_constant[0] == _lattice_constant[1]);

    prim_vec_dir(1,1) = 1.0; prim_vec_dir(2,1) = 0.0; prim_vec_dir(3,1) = 0.0;
    prim_vec_dir(1,2) = 0.0; prim_vec_dir(2,2) = 1.0; prim_vec_dir(3,2) = 0.0;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = 0.5;

    _prim_vec(1,1) = prim_vec_dir(1,1) * _lattice_constant[0];
    _prim_vec(2,2) = prim_vec_dir(2,2) * _lattice_constant[0];
    _prim_vec(1,3) = prim_vec_dir(1,3) * _lattice_constant[0]; _prim_vec(2,3) = prim_vec_dir(2,3) * _lattice_constant[0]; _prim_vec(3,3) = prim_vec_dir(3,3) * _lattice_constant[2];

  }

  else std::cout << "Lattice type " << _lattice_type << " doesn't exist" << std::endl;

};


void AtomisticGenerator::set_prim_miller(Tensor2Gen cut_planes)
{
  //Reduce cut planes if they're not in minimal integer form

  cut_planes = reduce_vector(cut_planes);

  //Transform miller indexes in primitive reciprocal vectors basis
  //Needed when conventional cell differs from primitive cell
  Tensor2Gen prim_miller_basis = reciprocal(_prim_vec);

  if (_lattice_type.compare("cubic") == 0){
    //In cubic lattice conventional cell equals to unit cell
    _prim_miller = cut_planes / (_lattice_constant[0]);
  }

  else if (_lattice_type.compare("bcc") == 0){
    _prim_miller = (inv(prim_miller_basis) * cut_planes / (_lattice_constant[0]));
    scale_to_int(_prim_miller);
  }

  else if (_lattice_type.compare("fcc") == 0){
    _prim_miller=inv(prim_miller_basis) * cut_planes / (_lattice_constant[0]);
    scale_to_int(_prim_miller);
  }

  else if (_lattice_type.compare("hexagonal") == 0){
    _prim_miller(1,1) = cut_planes(1,1) / _lattice_constant[0]; _prim_miller(2,1) = cut_planes(2,1) / _lattice_constant[1]; _prim_miller(3,1) = cut_planes(3,1) / _lattice_constant[2];
    _prim_miller(1,2) = cut_planes(1,2) / _lattice_constant[0]; _prim_miller(2,2) = cut_planes(2,2) / _lattice_constant[1]; _prim_miller(3,2) = cut_planes(3,2) / _lattice_constant[2];
    _prim_miller(1,3) = cut_planes(1,3) / _lattice_constant[0]; _prim_miller(2,3) = cut_planes(2,3) / _lattice_constant[1]; _prim_miller(3,3)=cut_planes(3,3) / _lattice_constant[2];
    scale_to_int(_prim_miller);
  }

  else {
    _prim_miller(1,1) = cut_planes(1,1) / _lattice_constant[0]; _prim_miller(2,1) = cut_planes(2,1) / _lattice_constant[1]; _prim_miller(3,1) = cut_planes(3,1) / _lattice_constant[2];
    _prim_miller(1,2) = cut_planes(1,2) / _lattice_constant[0]; _prim_miller(2,2) = cut_planes(2,2) / _lattice_constant[1]; _prim_miller(3,2) = cut_planes(3,2) / _lattice_constant[2];
    _prim_miller(1,3) = cut_planes(1,3) / _lattice_constant[0]; _prim_miller(2,3) = cut_planes(2,3) / _lattice_constant[1]; _prim_miller(3,3) = cut_planes(3,3) / _lattice_constant[2];
    scale_to_int(_prim_miller);
  };


};




void AtomisticGenerator::parse_parameters(Material* mat)
{

  Atom tmp;
  Tensor1 T;
  int i, j, n;

  if ( !(mat->is_alloy()) )
    {
      //WORKS ONLY FOR BULK, EXTEND TO ALLOY

      //lattice constant are expressed in Amstrong
      //_lattice_constant[0] = data("a", 0.0) * 10.0;

      Database& db = mat->get_database();
      db.set_section("lattice");
      _lattice_constant[0] = db.get("a", 0.0) * 10.0;
      if (_lattice_constant[0] == 0.0) std::cerr << "At least lattice constant a must be defined !!!!" << std::endl;

      /*_lattice_constant[1] = data("b", 0.0) * 10.0;*/
      _lattice_constant[1] = db.get("b", 0.0) * 10.0;
      if (_lattice_constant[1] == 0.0) _lattice_constant[1] = _lattice_constant[0];
      ////////////////////////////////////////
      _lattice_constant[2] = db.get("c", 0.0) * 10.0;
      if (_lattice_constant[2] == 0.0) _lattice_constant[2] = _lattice_constant[0];
      db.set_section("atomistic_structure");
      set_lattice_type(db.get("lattice_type", "none"));

      unsigned int n_basis_specie = db.get("n_basis_specie", 0);

      for (i = 1; i <= n_basis_specie; i++)
        {
          std::string record, s, n_s;
          std::stringstream out;

          out << i;
          s = out.str();

          record = "n_" + s;

          n = db.get(record.c_str(), 0);

          for (j = 1; j <= n; j++)
            {
              record.clear(); n_s.clear();
              record = "T_" + s + "_";
              out.str(std::string());
              out.clear(std::stringstream::goodbit);
              out << j;
              n_s = out.str();
              n_s = record + n_s;

              //Putting specie (defined by an integer) temporary in flag data
              tmp.set_flag(i);

              record = n_s + "_x";
              T(1) = db.get(record, 0.0);
              record = n_s + "_y";
              T(2) = db.get(record, 0.0);
              record = n_s + "_z";
              T(3) = db.get(record, 0.0);

              tmp.set_position(T);

              //Insert tmp in basis
              _crystal_basis.push_back(tmp);
            }
        }
    }

  if (mat->is_alloy())
    {
      //Cannot act dynamic cast on mat itself because constant
      Alloy* mat_alloy = dynamic_cast<Alloy*>(mat);

      //Get database for alloy (db) and for parental materials
      //NOTE: IMPLEMENTATION IS GOOD ONLY FOR BINARY COMPOUNDS
      //Nota: usiamo solo un pointer perche' per tutti i materiali viene istanziato solo un oggetto
      //database, a cui di volta in volta (ogni volta che chiamiamo un get) viene associato un data file.
      //Quindi non possiamo inizializzare 3 oggetti database e portarceli appresso, perche' saranno tutti
      //collegati al datafile settato dall'ultima assegnazione. Questa cosa va cambiata nella classe Database
      // (TODO)
      Database* db = &(mat_alloy->get_database());

      //std::cout << "component A is  " << mat_alloy->get_component_A()->get_name() << std::endl;
      //std::cout << "component B is  " << mat_alloy->get_component_B()->get_name() << std::endl;
      //Database& db1 = mat_alloy->get_component_A()->get_database();
      //Database& db2 = mat_alloy->get_component_B()->get_database();

      if (db->get("alloy_type", 2) == 2)
        {
          double ax_1, ay_1, az_1, ax_2, ay_2, az_2;
          db = &(mat_alloy->get_component_A()->get_database());
          db->set_section("lattice");

          //We express lattice constant in Amstrong
          ax_1 = db->get("a", 0.0) * 10.0;
          if (ax_1 == 0.0) std::cerr << "At least lattice constant a must be defined !!!!" << std::endl;

          ay_1 = db->get("b", 0.0) * 10.0;
          if (ay_1 == 0.0) ay_1 = ax_1;

          az_1 = db->get("c", 0.0) * 10.0;
          if (az_1 == 0.0) az_1 = ax_1;

          db = &(mat_alloy->get_component_B()->get_database());
          db->set_section("lattice");

          ax_2 = db->get("a", 0.0) * 10.0;
          if (ax_2 == 0.0) std::cerr << "At least lattice constant a must be defined !!!!" << std::endl;

          ay_2 = db->get("b", 0.0) * 10.0;
          if (ay_2 == 0.0) ay_2 = ax_2;

          az_2 = db->get("c", 0.0) * 10.0;
          if (az_2 == 0.0) az_2 = ax_2;

          //Setting lattice parameters for the alloy
          double molar_fraction = mat->get_options().get_option("x", 1.0);
          _lattice_constant[0] = ax_1 * molar_fraction + ax_2 * (1.0 - molar_fraction);
          _lattice_constant[1] = ay_1 * molar_fraction + ay_2 * (1.0 - molar_fraction);
          _lattice_constant[2] = az_1 * molar_fraction + az_2 * (1.0 - molar_fraction);

        }

      db = &(mat_alloy->get_component_A()->get_database());
      db->set_section("");
      db->set_section("atomistic_structure");
      set_lattice_type(db->get("lattice_type", "none"));

      db = &(mat_alloy->get_database());
      db->set_section("atomistic_structure");

      unsigned int n_basis_specie = db->get("n_basis_specie", 0);

      db = &(mat_alloy->get_component_A()->get_database());
      db->set_section("atomistic_structure");

      for (i = 1; i <= n_basis_specie; i++)
        {
          std::string record("");
          std::string s("");
          std::stringstream out;
          out << i;
          s = out.str();
          record = "n_" + s;

          unsigned int n_x = (db->get(record, 0));
          for (j = 1; j <= n_x; j++)
            {
              std::string s2;
              record = "T_" + s + "_";
              out.str(std::string());
              out.clear(std::stringstream::goodbit);
              out << j;
              s2 = out.str();
              s2 = record + s2;

              //Putting specie (defined by an integer) temporary in flag data
              // ???????????? CHECK IT , WHY i IS SET AS FLAG???????????
              tmp.set_flag(i);

              record = s2 + "_x";
              T(1) = db->get(record, 0.0);
              record = s2 + "_y";
              T(2) = db->get(record, 0.0);
              record = s2 + "_z";
              T(3) = db->get(record, 0.0);
              tmp.set_position(T);

              //Insert tmp in basis
              _crystal_basis.push_back(tmp);
            }
        }


    }


};




void AtomisticGenerator::make_conv_cell()
{
  //Calculate conventional cell vectors in the directions given by cut planes (conventional growth cell)
  Tensor1 m1,m2,m3,select_vect(0);
  Tensor1 conv1, conv2, conv3;
  int i;
  Tensor2Gen prim_miller_basis;

  prim_miller_basis = reciprocal(_prim_vec);
  _conv_prim = inv(_prim_vec) * prim_miller_basis * _prim_miller;
  scale_to_int(_conv_prim);
  _conv_vect = _prim_vec * _conv_prim;

  for (i = 1; i <= 3; i++) {conv1(i) = _conv_vect(i, 1); conv2(i) = _conv_vect(i, 2); conv3(i) = _conv_vect(i, 3);};

  //At least x growth index must be specified (growth direction in 1D structure)
  assert (norm(conv1) > tol);

  // If no other indexes are specified, set y growth direction by default (orthogonal to x index)
  if ((norm(conv2) < tol)&&(norm(conv3) < tol)) {

    std::cerr << "Warning: only x growht direction is defined. Building other direction orthogonal " << std::endl;

    conv2(1) = - conv1(3); conv2(2) = 0.0; conv2(3) = conv1(1);
    conv3 =  vectorProduct(conv1, conv2);
    conv3 = conv3 / norm(conv3);
  }

  //If only one index is not specified, build it by default (orthogonal)
  else if (norm(conv2) < tol){

    std::cerr << "Warning: only x and z growht direction is defined. Building y direction orthogonal " << std::endl;

    conv2 =  vectorProduct(conv1, conv3);
    conv2 = conv2 / norm(conv2);
  }

  else if (norm(conv3) < tol){

    std::cerr << "Warning: only x and y growht direction is defined. Building z direction orthogonal " << std::endl;

    conv3 =  vectorProduct(conv1, conv2);
    conv3 = conv3 / norm(conv3);
  };


  //   //If not all miller indexes are specified, build others as a default
  //   if (norm(conv1) < tol) {conv1(1) = - conv3(3); conv1(2) = 0.0; conv1(3) = conv3(1);}

  //   //If only z and x miller indexes are specifed, get y miller index
  //   if (norm(conv2) < tol) {
  //     conv2 =  vectorProduct(conv1, conv3);
  //     conv2 = conv2 / norm(conv2);
  //   };


  for (i = 1; i <= 3; i++) {_conv_vect(i, 1) = conv1(i); _conv_vect(i, 2) = conv2(i); _conv_vect(i, 3) = conv3(i);};


  _conv_prim = inv(_prim_vec) * _conv_vect;
  scale_to_int(_conv_prim);
  _conv_vect = _prim_vec * _conv_prim;

  // Calculate distance between equivalent planes for every cut plane
  //prim_miller = inv(prim_miller_basis) * prim_vec * conv_prim;
  //for (i = 1; i <= 3; i++) {m1(i) = prim_miller(i, 1); m2(i) = prim_miller(i, 2); m3(i) = prim_miller(i, 3);};
  //planar_distance[0] = 1 / norm(prim_miller_basis * m1);
  //planar_distance[1] =1 / norm(prim_miller_basis * m2);
  //planar_distance[2] = 1 / norm(prim_miller_basis * m3);

};


void AtomisticGenerator::make_conv_basis()
{
  //Fill the conventional growth cell with atomic basis
  int lower_1, lower_2, lower_3, upper_1, upper_2, upper_3, i;
  Tensor1 prim_position, tmp_check;
  Tensor1 tmp_position;
  Tensor1 vec_x,vec_y,vec_z;
  std::vector<Atom>::iterator basis_iterator;
  Tensor2Gen rot_tmp;

  //Make a preliminar rotation only if conventional cell vectors are orthogonal

  vec_x(1) = _conv_vect(1,1); vec_x(2) = _conv_vect(2,1); vec_x(3) = _conv_vect(3,1);
  vec_y(1) = _conv_vect(1,2); vec_y(2) = _conv_vect(2,2); vec_y(3) = _conv_vect(3,2);
  vec_z(1) = _conv_vect(1,3); vec_z(2) = _conv_vect(2,3); vec_z(3) = _conv_vect(3,3);


  //If z and y vectors are not orthogonal, build rotation vector y as a vector orthogonal to x
  if ( (vec_z * vec_y ) > 1e-10) {

    if ( (vec_y * vec_x) < 1e-10 ) {vec_z = vectorProduct(vec_x, vec_y);}
    else if  ( (vec_z * vec_x) < 1e-10)  {vec_y = vectorProduct(vec_x, vec_z);}
    else {std::cout << "Warning: at least x or y growth direction should be orthogonal to z growth direction" << std::endl;}

  }


  //vec_y = vectorProduct(vec_x, vec_z);
  if (((vec_x * vec_y) < 1e-10) && ((vec_x * vec_z) < 1e-10) && ((vec_y * vec_z) < 1e-10)) {

    assert(norm(vec_x) > 1e-10);
    assert(norm(vec_y) > 1e-10);
    assert(norm(vec_z) > 1e-10);

    for ( int i = 1; i <=3; i++ )  rot_tmp(1,i) = vec_x(i)/norm(vec_x);
    for ( int i = 1; i <=3; i++ )  rot_tmp(2,i) = vec_y(i)/norm(vec_y);
    for ( int i = 1; i <=3; i++ )  rot_tmp(3,i) = vec_z(i)/norm(vec_z);
  }
  else {
    rot_tmp(1,1) = 1.0; rot_tmp(1,2) = 0.0; rot_tmp(1,3) = 0.0;
    rot_tmp(2,1) = 0.0; rot_tmp(2,2) = 1.0; rot_tmp(2,3) = 0.0;
    rot_tmp(3,1) = 0.0; rot_tmp(3,2) = 0.0; rot_tmp(3,3) = 1.0;
    std::cout << "Warning: no rotation has been done, orientation of atoms in space depends on primitive vectors definitions" << std::endl;
  }

  //Also a user-defined rotation is allowed, which put the axes in a direction different
  //from canonical basis
  for (i = 1; i <= 3; i++){vec_x(i) = _rotation(i,1); vec_y(i) = _rotation(i,2); vec_z(i) = _rotation(i,3);}
  vec_x =vec_x / norm(vec_x); vec_y = vec_y / norm(vec_y); vec_z = vec_z / norm(vec_z);
  for (i = 1; i <= 3; i++){_rotation(i,1) = vec_x(i); _rotation(i,2) = vec_y(i); _rotation(i,3) = vec_z(i);}
  _rotation = _rotation * rot_tmp;


  //Define a box including conventional cell
  lower_1 = int(std::min(0.0,std::min(_conv_prim(1,1),std::min(_conv_prim(1,2),_conv_prim(1,3)))));
  upper_1 = int(std::max(0.0,std::max(_conv_prim(1,1),std::max(_conv_prim(1,2),_conv_prim(1,3)))));
  lower_2 = int(std::min(0.0,std::min(_conv_prim(2,1),std::min(_conv_prim(2,2),_conv_prim(2,3)))));
  upper_2 = int(std::max(0.0,std::max(_conv_prim(2,1),std::max(_conv_prim(2,2),_conv_prim(2,3)))));
  lower_3 = int(std::min(0.0,std::min(_conv_prim(3,1),std::min(_conv_prim(3,2),_conv_prim(3,3)))));
  upper_3 = int(std::max(0.0,std::max(_conv_prim(3,1),std::max(_conv_prim(3,2),_conv_prim(3,3)))));

  for (int i = lower_1 - 1; i <= upper_1 + 1; i++){
    for (int j = lower_2 - 1; j <= upper_2 + 1; j++){
      for (int l = lower_3 - 1; l <= upper_3 + 1; l++){

        prim_position(1) = double(i); prim_position(2) = double(j); prim_position(3) = double(l);

        tmp_position = prim_position;

        bool check_boundary;

        tmp_check = inv(_conv_prim) * tmp_position;

        check_boundary= ((tmp_check(1) >= -tol) && (tmp_check(1) < (1.0 - tol)))&&
        ((tmp_check(2) >= -tol) && (tmp_check(2) < (1.0 - tol))) &&
        ((tmp_check(3) >= -tol) && (tmp_check(3) < (1.0 - tol)));

        if (check_boundary){
          tmp_position = _rotation * _prim_vec * prim_position;

          _conv_lattice_basis.push_back(tmp_position);

        }

      };
    };
  };
  _conv_vect = _rotation * _conv_vect;

};



//Bond map generation (cluster)
void  AtomisticGenerator::bond_map_gen(std::vector<Atom> &basis){

  //use internal member, if already used delete it
  if (_bondmapobject == NULL) _bondmapobject = new BondMap;
  else
    {
      delete _bondmapobject;
      _bondmapobject = new BondMap;
    }

  _bondmapobject->do_init(basis, _period);
  std::cout << "Solving bond map " << std::endl;
  _bondmapobject->do_solve(basis);
  //std::cout << "Getting and returning ";
  //return _bondmapobject->get_bond_map();

};


void AtomisticGenerator::passivate()
{
  unsigned int ** bond_map;
  Tensor1 position;
  double hydrogen_distance = 1.0;
  Atom* bonded_atom;

  std::cout << "Starting passivate " << std::endl;

  if (_bondmapobject == NULL)
    {
      bond_map_gen(_structure_basis);
    }

  bond_map = _bondmapobject->get_bond_map();


  //Warning: cycle end must be defined before as size will change dynamically during cycle
  //and we need acting only on already defined structure
  unsigned int size_before_passivating = _structure_basis.size();

  for (unsigned int i = 0; i < size_before_passivating; i++)
    {
      if (_structure_basis[i].belong_to_structure)
        {

          for (unsigned int j = 0; j < bond_map[i][8]; j++)
            {
              if (bond_map[i][8] != 4) std::cout << "Warning, atom has not 4 neighbours " << std::endl;

              bonded_atom = &(_structure_basis[bond_map[i][j]]);
              if (!((*bonded_atom).belong_to_structure))
                {
                  //TODO: using default copy constructor, with further modifications to
                  //Atom class it could not work anymore!
                  //Position must be modified in order to put Hydrogen atom near,
                  //and also as we cannot have hydrogen bonded to more than one atom,
                  //so in some cases we cannot keep crystal positions
                  Atom tmp(*bonded_atom);
                  tmp.set_specie("H");
                  tmp.belong_to_structure = true;

                  position = _structure_basis[i].get_position() +
                  ( ( bonded_atom->get_position() - _structure_basis[i].get_position()) /
                      (norm(bonded_atom->get_position() - _structure_basis[i].get_position() ) ) ) *
                      hydrogen_distance;
                  tmp.set_position(position);

                  _structure_basis.push_back(tmp);

                }

            }
        }

    }

  std::cout << "Passivate done " << std::endl;


}
//void AtomisticGenerator::passivate_cluster(std::vector<Atom> &basis){
//
//  Tensor1 u, u1, u2, r1, r2, r3, O;
//  double R1, R2;
//  const double sq3 = sqrt(3.0);
//  const double d = 0.64;
//  unsigned int i;
//  std::vector<Atom> hydrogens;
//  Atom tmp;
//  double sin109, cos109;
//  unsigned int ** bond_map;
//
//  sin109 = sin ( ( 180.0 / ( asin(1.0) * 2.0 ) ) * 109.471 );
//  cos109 = cos ( ( 180.0 / ( asin(1.0) * 2.0 ) ) * 109.471 );
//
//
//
//  if (_bondmapobject == NULL) {
//    bond_map_gen(basis);
//  }
//
//  bond_map = _bondmapobject->get_bond_map();
//
//  for (i = 0; i < basis.size(); i++){
//
//    if ( bond_map[i][8] == 3 ){
//
//      O = basis[i].get_position();
//      r1 = basis[ bond_map[i][0] ].get_position() - O;
//      r2 = basis[ bond_map[i][1] ].get_position() - O;
//      r3 = basis[ bond_map[i][2] ].get_position() - O;
//
//      u1 = r1 + r2 + r3;
//      R1 = norm(u1);
//      u1 = -d * (u1/R1);
//
//      tmp.set_specie("H");
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      tmp.set_position ( O + u1 );
//      hydrogens.push_back(tmp);
//    }
//
//    else if  ( bond_map[i][8] == 2 ){
//
//      O = basis[i].get_position();
//      r1 = basis[ bond_map[i][0] ].get_position() - O;
//      r2 = basis[ bond_map[i][1] ].get_position() - O;
//
//      u1 = r1 + r2;
//      u2 = vectorProduct(r1, r2);
//      R1 = norm(u1); R2 = norm(u2);
//
//      u = - (u1 / R1) - sq3 * (u2 / R2);
//      u = d * (u / 2.0);
//
//      tmp.set_specie("H");
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//      tmp.set_position ( O + u );
//      hydrogens.push_back(tmp);
//
//      u = - (u1 / R1) + sq3 * (u2 / R2);
//      u = d * (u / 2.0);
//
//      tmp.set_specie("H");
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      tmp.set_position ( O + u );
//      hydrogens.push_back(tmp);
//    }
//
//    else if (bond_map[i][8] == 1){
//
//      O = basis[i].get_position();
//      r1 = basis[ bond_map[i][0] ].get_position() - O;
//
//      Tensor2Gen vect_rot(0);
//      vect_rot(1,1) = cos109; vect_rot(1,2) = sin109; vect_rot(1,3) = 0.0;
//      vect_rot(2,1) = sin109; vect_rot(2,2) = cos109; vect_rot(2,3) = 0.0;
//      vect_rot(3,1) = 0.0; vect_rot(3,2) = 0.0; vect_rot(3,3) = 1.0;
//
//      u = vect_rot * r1;
//      R1 = norm(r1);
//      r2 = u * R1;
//      u = d * (u / norm(u));
//
//      tmp.set_specie("H");
//
//      // TEMPORARY SOLUTION FOR WURTZITE
//      //u = (-r1)/norm(r1); u(1) = u(1) + 0.1;
//      //////////////////////////////
//
//
//
//      tmp.set_position ( O + u );
//      tmp.set_flag ( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      hydrogens.push_back(tmp);
//
//      u1 = r1 + r2;
//      u2 = vectorProduct(r1, r2);
//      R1 = norm(u1); R2 = norm(u2);
//
//      u = - (u1 / R1) - sq3 * (u2 / R2);
//      u = d * (u / 2.0);
//
//      tmp.set_specie("H");
//
//      // TEMPORARY SOLUTION FOR WURTZITE
//           // u = (-r1)/norm(r1); u(2) = u(2) + 0.1;
//            //////////////////////////////
//
//      tmp.set_position ( O + u );
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      hydrogens.push_back(tmp);
//
//      u = - (u1 / R1) + sq3 * (u2 / R2);
//      u = d * (u / 2.0);
//
//      tmp.set_specie("H");
//      tmp.set_flag( basis[i].get_flag() );
//
//      // Up to now hydrogen is considered part of the same material of bonded atom
//      tmp.set_region_ID(basis[i].get_region_ID());
//
//      // TEMPORARY SOLUTION FOR WURTZITE
//                // u = (-r1)/norm(r1); u(2) = u(2) - 0.1; u(1)=u(1) - 0.1;
//                 //////////////////////////////
//
//      tmp.set_position ( O + u );
//      hydrogens.push_back(tmp);
//    }
//
//    else if (bond_map[i][8] != 4) {std::cout << "Warning! atom " << i
//      << " is bonded to " << bond_map[i][8] << " atoms" << std::endl;}
//
//  }
//  for (i = 0; i < hydrogens.size(); i++) {basis.push_back( hydrogens[i] );}
//
//
//
//};


//Some data manipulation function useful only in this class


Tensor2Gen
AtomisticGenerator::reciprocal(Tensor2Gen real_basis)
{

  //Build the reciprocal basis related to input 2-rank tensor
  Tensor1 a1,a2,a3;
  Tensor1 b1,b2,b3;
  Tensor1 select_vect(0);
  Tensor2Gen reciprocal;

  //Select vector a1
  select_vect(1) = 1.0; a1 = real_basis * select_vect;

  //Select vector a2
  select_vect(1) = 0.0; select_vect(2) = 1.0; a2 = real_basis * select_vect;

  //Select vector a3
  select_vect(2) = 0.0; select_vect(3) = 1.0; a3 = real_basis * select_vect;

  const double volume = a1 * vectorProduct(a2,a3);
  assert (volume != 0);
  b1 = vectorProduct(a2,a3) / volume;
  b2 = vectorProduct(a3,a1) / volume;
  b3 = vectorProduct(a1,a2) / volume;

  reciprocal(1,1)=b1(1); reciprocal(2,1)=b1(2);reciprocal(3,1)=b1(3);
  reciprocal(1,2)=b2(1);reciprocal(2,2)=b2(2);reciprocal(3,2)=b2(3);
  reciprocal(1,3)=b3(1);reciprocal(2,3)=b3(2);reciprocal(3,3)=b3(3);

  assert(det(reciprocal) != 0);

  return reciprocal;
};


int
AtomisticGenerator::compare_tol(double a, double b)
{
  //Comparison routine with a tolerance defined as internal constant.
  //If absolute value of difference between a and b is minor than tolerance,
  //a and b are considered equal
  if (std::fabs(a-b) < tol) return 1;
  else return 0;
};


int
AtomisticGenerator::double_to_int_cast_checked(double a)
{
  //Convert a double to the nearest integer, within a certain tolerance
  int n;
  if (std::abs(std::floor(a)-a) < std::abs(std::ceil(a) - a)) n = int(std::floor(a));
  else n = int(std::ceil(a));
  assert (std::abs(double(n) - a) < tol);
  return n;
};


double
AtomisticGenerator::double_to_int_value_checked(double a)
{
  //Gives the double number equal to the integer nearest to a, within a certain tolerance
  double b;
  b = double(double_to_int_cast_checked(a));
  return b;
};


void
AtomisticGenerator::double_to_int_value_checked(Tensor1& a)
{
  double tmp;;

  tmp = double_to_int_value_checked(a(1)); a(1) = tmp;
  tmp = double_to_int_value_checked(a(2)); a(2) = tmp;
  tmp = double_to_int_value_checked(a(3)); a(3) = tmp;
};


int
AtomisticGenerator::gcd(int a, int b)
{
  //Calculate greater common denominator between integers
  //(return 0 if gcd(a,0) or gcd(0,a)

  //Added: if illegal operation gcd(a,0) performed return abs(a) (useful for reduce_vector routine)
  if ((a == 0) || (b == 0)) return std::max( std::abs(a), std::abs(b) );

  //by Derek Chandler, MEng, MIEE
  int reminder;
  do{
    reminder = a % b;
    if (reminder != 0)
      {
        a = b;
        b = reminder;
      }
  } while (reminder);
  return b;
};


Tensor1
AtomisticGenerator::reduce_vector(Tensor1 v)
{
  //Reduce a vector of double containing integer values to its minimal form
  Tensor1 v_tmp;
  int gcd_tmp,gcd_value;

  if (norm(v) < tol) return v;

  //Find the maximimum common denominator
  gcd_tmp = gcd(double_to_int_cast_checked(v(1)),double_to_int_cast_checked(v(2)));
  gcd_value = gcd(gcd_tmp,double_to_int_cast_checked(v(3)));
  v_tmp = v / double(gcd_value);
  double_to_int_value_checked(v_tmp);
  return v_tmp;

};


void AtomisticGenerator::scale_to_int(Tensor1& a)
{
  //Expand a double vector to a vector having same direction but integer values
  //If vector is not in a reduced form (having integer values with gcd > 1) it's reduced
  //If input is a zero vector, a zero vector is returned

  //Check if a is already integer
  if ( (fabs(a(1) - round(a(1))) < tol) && (fabs(a(2) - round(a(2))) < tol) && (fabs(a(3) - round(a(3))) < tol) )
    {
      a(1) = double_to_int_value_checked(a(1));
      a(2) = double_to_int_value_checked(a(2));
      a(3) = double_to_int_value_checked(a(3));
    }

  else

    {

      int i = 0;
      Tensor1 a_tmp;

      do{
        i = i + 1;
        a_tmp = a * i;
      }while (  ( fabs (a_tmp(1) - round(a_tmp(1))) >= tol)  ||  (fabs(a_tmp(2) - round(a_tmp(2))) >= tol)  || (fabs(a_tmp(3) - round(a_tmp(3))) >= tol) );

      a(1) = a_tmp(1); a(2) = a_tmp(2); a(3) = a_tmp(3);

    }

  a = reduce_vector(a);

}



void
AtomisticGenerator::scale_to_int(Tensor2Gen& a)
{
  //Same of scale_to_int with Tensor1 argument, considering the columns of a 2-rank tensor
  Tensor1 tmp;
  Tensor1 select_vect(0);

  select_vect(1) = 1.0;  tmp = a * select_vect;
  scale_to_int(tmp);
  a(1,1) = tmp(1); a(2,1) = tmp(2); a(3,1) = tmp(3);

  select_vect(1) = 0.0; select_vect(2) = 1.0; tmp = a * select_vect;
  scale_to_int(tmp);
  a(1,2) = tmp(1); a(2,2) = tmp(2); a(3,2) = tmp(3);

  select_vect(2) = 0.0; select_vect(3) = 1.0; tmp = a * select_vect;
  scale_to_int(tmp);
  a(1,3) = tmp(1); a(2,3) = tmp(2); a(3,3) = tmp(3);
};


Tensor2Gen
AtomisticGenerator::reduce_vector(Tensor2Gen a)
{
  //Same of reduce_vector with Tensor1 argument, considering the columns of a 2-rank tensor
  Tensor1 tmp1,tmp2;
  Tensor1 select_vect(0);
  Tensor2Gen b(0);

  //Select vector a1
  select_vect(1) = 1.0;  tmp1 = a * select_vect;
  tmp2 = reduce_vector(tmp1);
  b(1,1) = tmp2(1); b(2,1) = tmp2(2); b(3,1) = tmp2(3);

  //Select vector a2
  select_vect(1) = 0.0; select_vect(2) = 1.0; tmp1 = a * select_vect;
  tmp2 = reduce_vector(tmp1);
  b(1,2) = tmp2(1); b(2,2) = tmp2(2); b(3,2) = tmp2(3);

  //Select vector a3
  select_vect(2) = 0.0; select_vect(3) = 1.0; tmp1 = a * select_vect;
  tmp2 = reduce_vector(tmp1);
  b(1,3) = tmp2(1); b(2,3) = tmp2(2); b(3,3) = tmp2(3);
  return b;

};

