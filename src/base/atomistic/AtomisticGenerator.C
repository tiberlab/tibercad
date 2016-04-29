// $Id$

#include "AtomisticGenerator.h"
#include "AtomisticStructure.h"
#include "AtomisticGenerator1D.h"
#include "AtomisticGenerator2D.h"
#include "AtomisticGenerator3D.h"
#include "BondMap.h"
#include "Messages.h"
#include "MeshUtils.h"
#include "Specie.h"
#include "Utils.h"
#include "RotatedCrystal.h"
#include "Atom.h"
#include "BulkCrystal.h"
#include <plane.h>

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>
#include <ctime>
#include <tr1/random>


using namespace std;


AtomisticGenerator::AtomisticGenerator(void)
:_bondmap(NULL),
_reference_material(NULL),
_conv_vect(0),
_conv_prim(0),
_local_origin(0),
_translation(0),
_cell_translation(0),
_period(0),
_bulk(NULL)
{

}

AtomisticGenerator::~AtomisticGenerator(void)
{
  if (_bondmap != NULL) delete _bondmap;
}


const double AtomisticGenerator::tol = 5e-2;

AtomisticGenerator*
AtomisticGenerator::create(AtomisticStructure* const as, unsigned int dimension)
{
  AtomisticGenerator* ag =  NULL;
  if (dimension == 1) ag = AtomisticGenerator1D::create(as);
  if (dimension == 2) ag = AtomisticGenerator2D::create(as);
  if (dimension == 3) ag = AtomisticGenerator3D::create(as);

  return ag;
}


void AtomisticGenerator::print_basis(std::vector<Atom> &basis, const std::string filename){

  std::ofstream output_file;

  std::vector<Atom>::iterator basis_iterator = basis.begin();

  output_file.open(filename.c_str());
  output_file << basis.size() << std::endl << std::endl;

  do{

    output_file << std::setw(2) << (*basis_iterator).get_specie()
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).get_position(0))
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).get_position(1))
    << std::setw(20) << std::setprecision(10)<< std::fixed  << double((*basis_iterator).get_position(2)) << "\n";

    ++basis_iterator;

  }while(basis_iterator != basis.end());

  output_file.close();

};



void
AtomisticGenerator::init_commons()
{

  std::ostringstream os;
  
  //Set dimensional scale
  scale = _as->get_scale();
  // Set material informations
  //-----------------------------------------------------------------------------------------
  std::string ref_region = _as->get_options().get_option("reference_region", "");
  std::set<ID> ids;
  _as->get_device()->get_active_region_ids(ref_region, ids);
  if (ids.size() != 0)
  {
    _reference_material = _as->get_device()->get_material(*ids.begin());
  }
  else
  {
    ModelOptions::const_submodel_iterator it =
        _as->get_options().submodels_begin("reference_material");
    if (it != _as->get_options().submodels_end("reference_material"))
    {
      const ModelOptions& refmat_opts = it->second;
      Material* refmat = Material::create(refmat_opts.get_name(), refmat_opts);
      refmat->init();
      _reference_material = refmat;
    }
  }

  if (_reference_material == NULL)
    throw InitFailedException("Reference region/material missing or badly "
        "defined for structure " + _as->get_name() );

  //Build the right BulkCrystal object
  //Additional options respect to the region should be specified here
  Messages::info("Create reference lattice: "+_reference_material->get_name());

  _bulk = BulkCrystal::create(_reference_material); 
  
  _bulk->do_init();

  //_bulk->print_xyb("Ref.xyb");
  //-----------------------------------------------------------------------
  
  //A translation vector can be specified to modify supercell alignment
  if ( _as->get_options().find_option("conventional_cell_shift") )
  {
    std::vector<double> translation(3, 0.0);
   _as->get_options().get_option("conventional_cell_shift", translation);
   _cell_translation(1) = translation[0]; 
   _cell_translation(2) = translation[1]; 
   _cell_translation(3) = translation[2];
   if (abs(_cell_translation(1)) > 1.0 || 
       abs(_cell_translation(2)) > 1.0 || 
       abs(_cell_translation(3)) > 1.0 )
      throw InitFailedException("cell_translation in relative coordinates must be < 1.0");
  }

  if ( _as->get_options().find_option("translation") )
  {
    std::vector<double> translation(3, 0.0);
   _as->get_options().get_option("translation", translation);
   _translation(1) = translation[0]; 
   _translation(2) = translation[1]; 
   _translation(3) = translation[2];
  }
  
  
}


void
AtomisticGenerator::do_init()
{
  std::ostringstream os;

  init_commons();

  // Set the vector of elements covered by structure regions,
  // useful for change specie and cut
  MeshBase::element_iterator el(_as->get_device()->get_mesh().elements_begin());
  const MeshBase::element_iterator el_end(_as->get_device()->get_mesh().elements_end());
  //number of elements in atomistic regions
  unsigned int num_elem = 0;
  for ( ; el != el_end; ++el)
  {
    Elem* elem = *el;
    if (_as->get_IDset().count(elem->subdomain_id()))
      ++num_elem;
  }
  _structure_elements.reserve(num_elem);
  for (el = _as->get_device()->get_mesh().elements_begin(); el != el_end; ++el)
  {
    Elem* elem = *el;
    if (_as->get_IDset().count(elem->subdomain_id()))
      _structure_elements.push_back(elem);
  }

  //Build up supercell structure with proper options
  Messages::info("Build a first oversized structure");
  build();
  
  //os<<"Initial structure with "<<_super_basis.size()<<" atoms"<<std::endl;
  //Messages::info(os.str());

  std::string preserve="none";
  //preserve = _as->get_options().get_option("preserve", "none");
 
  // associate the right element to each atom
  associate_elements(_as->get_IDset());
  
  // cut the structure (only flags atoms)
  cut(_as->get_IDset(), preserve);

  // iterates on _super_basis and assign species
  assign_species();

  //eventually enlarge along dummy supercell directions
  //(build bondmap if not present)
  check_periodic();

  //rebuild bondmap with new periodicity
  Messages::info("Rebuild bondmap");
  delete _bondmap; _bondmap = NULL;
  bond_map_gen(_super_basis);

  //Passivate structure (build bondmap if not present)
  if (_as->get_options().get_option("passivation", false))
    passivate();

  // remove unflagged atoms
  remove_atoms();

  // assign random-alloy species (only cations for the moment)
  if (_as->is_random_alloy())
  {
    if (_as->get_options().get_option("old_random_routine",false))
    {
      build_random_old();
    }
    else
      build_random_alloy();
  }

  // delete generator bondmap
  delete _bondmap;
  _bondmap = NULL;

}


//Eliminate not included atoms from structure
//Check structure to eliminate unincluded atoms (using swap in another vector)
void
AtomisticGenerator::remove_atoms(void)
{
  _structure_basis.clear();
  _structure_basis.reserve(_super_basis.size());

  for (unsigned int i = 0; i < _super_basis.size(); i++)
  {
    if ( _belong_to_structure[i] )
    {
      _structure_basis.push_back(_super_basis[i]);
    }
  }

  // We destroy the array that is now inconsistent with new structure. 
  // NOTE: if any use is needed later it must be properly defined here
  _belong_to_structure.clear();

}

void 
AtomisticGenerator::finalize(void)
{

  //Pass data to AtomisticStructure
  //--------------------------------------------------------------------------------------------------

  Atom tmp_atom;
  //TODO:not safe, better swap arrays, so then we can delete AtomisticGenerator instance
  _as->set_structure_atoms(_structure_basis);
  
  _as->set_ttype_lattice_vectors(_period);
 
  //_as->set_periodicity(_periodicity);

  _as->set_N_atoms( _structure_basis.size() );

  std::set<std::string> atom_types;

  for (unsigned int i = 0; i < _structure_basis.size(); i++)
  {
    atom_types.insert(_structure_basis[i].get_specie().get_string());
  }

  _as->set_N_types( atom_types.size() );

  _as->clear_atom_types();

  _as->set_atom_types(atom_types);


};


void 
AtomisticGenerator::associate_elements(const std::set<ID>& reg_ids)
{

  // the tensor grid to real mesh mapper for fast association atom->Elem
  // NOTE: we pass the relevant ID set, since Quantum contacts could be present
  //       which have to be included in the atomistic structure
  MeshUtils::GridMapper& mapper =
      MeshUtils::GridMapper::get_mapper(_as->get_device()->get_mesh(), reg_ids);


  Utils::Progress prog("Assign elements", _super_basis.size());
  unsigned int progress = 0;

  std::vector<Atom>::iterator atom = _super_basis.begin();

  for ( ; atom != _super_basis.end(); ++atom)
  {
    
    Point p((*atom).get_position());
    p *= 1.0 / scale;
    
    // set unneeded dim to 0, so atoms are associated to the correct elements
    switch (_dim)
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
      (*atom).set_elem(elem);
    }
    else
    {
      //cout<<endl<<"(debug) atom coords: "<<p(0)<<", "<<p(1)<<", "<<p(2)<<"  have no elem."<<endl;
    }

    progress++;
    prog.progress_message(progress);

  }
}


void
AtomisticGenerator::cut(const std::set<ID>& reg_ids, const std::string preserve)
{

  //Different strategies if preserving conventional cell or preserving basis are needed
  if (preserve.compare("none") == 0)
  {
    _belong_to_structure.reserve(2 * _super_basis.size());

    for (unsigned int i=0; i<  _super_basis.size(); i++)
    {
      Atom& atom = _super_basis[i];

      const Elem* elem = atom.get_elem();

      if (elem == NULL)
      {
        _belong_to_structure.push_back(false);
      }
      else
      {
        ID el_reg = elem->subdomain_id();
        _belong_to_structure.push_back( (reg_ids.count(el_reg) > 0) ?  true : false );
      }                            
    }
  }
  else
  {
    Messages::error("preserve conventional and primitive cell not implemented yet");
  }

}

// Assign species. 
//
void
AtomisticGenerator::assign_species(void)
{
  std::set<ID> reg_ids(_as->get_IDset());
  std::map<ID, std::map<Atom::label_t, Specie> > assign;

  // we will need the reference region for atoms falling outside of
  // the atomistic structure's regions
  std::set<ID>::iterator reg(reg_ids.begin());
  unsigned int ref_rgn_id = *reg;

  for ( ; reg != reg_ids.end(); ++reg)
  {

    const Material* mat = _as->get_device()->get_material( (*reg) );

    const Database& db = mat->get_database();

    db.set_section("atomistic_structure");
    
    //Build up conversion map from file
    for (unsigned int i = 1; i <= db.get("n_basis_specie", 0); i++)
    {
      std::string record;
      std::string s;
      std::stringstream out;

      out << i;
      s = out.str();
      record = "specie_" + s;
      //note: db gets the species of the first alloy component
      std::string db_record = db.get(record.c_str(),"none");
      assign[*reg][static_cast<Atom::label_t>(i)] = Specie(db_record);
    }

    //No more reading from section atomistic_structure in database are needed
    db.set_section("");

  }

  //
  // Cycle upon all atoms and change specie according to assign map
  //
  int n_super_basis = _super_basis.size();
 
  Utils::Progress prog("Assign Species", n_super_basis);

  unsigned int progress_counter = 0;
  
  for (unsigned int i=0; i < _super_basis.size(); i++)
  {
    Atom& atom = _super_basis[i];

    if (_belong_to_structure[i])
    {    
      ID el_reg = atom.get_elem()->subdomain_id();
      Specie tmp =  assign[el_reg][atom.get_label()];
      atom.set_specie(tmp);
    }
    else
    {
      //Species of atoms outside the structure are assigned for the bondmap
      Specie tmp = assign[ref_rgn_id][atom.get_label()];
      atom.set_specie(tmp);
    }
    
    progress_counter += 1;
    prog.progress_message(progress_counter);
  }


}


void
AtomisticGenerator::dorestrict(bool passivation)
{

  _super_basis = _as->get_structure_atoms();

  _period = _as->get_ttype_lattice_vectors();


  cut(_as->get_IDset());

  // We must re-attach passivation Hydrogens
  

  if (passivation) passivate();  

  remove_atoms();

}


double
AtomisticGenerator::substitution_probability(size_t id, const Specie& sp)
{
  const BondMap& bm = *_bondmap;

  int same_species = 0;
  int n_neigh = 0;

  std::set<size_t> visited;
  visited.insert(id);
  const std::vector<unsigned int>& neigh = bm[id];
  for (unsigned int i = 0; i < neigh.size(); ++i)
  {
    const std::vector<unsigned int>& nn = bm[neigh[i]];
    for (unsigned int j = 0; j < nn.size(); ++j)
    {
      if (!visited.count(nn[j]))
      {
        n_neigh++;
        visited.insert(nn[j]);
        if (_structure_basis[nn[j]].get_specie() == sp)
          same_species += 1;
      }

      const std::vector<unsigned int>& nn2 = bm[nn[j]];
      for (unsigned int ii = 0; ii < nn2.size(); ++ii)
      {
        const std::vector<unsigned int>& nn3 = bm[nn2[ii]];
        for (unsigned int jj = 0; jj < nn3.size(); ++jj)
        {
          if (!visited.count(nn3[jj]))
          {
            n_neigh++;
            visited.insert(nn3[jj]);
            if (_structure_basis[nn3[jj]].get_specie() == sp)
              same_species++;
          }
        }
      }
    }
  }


  // for zb, wz:
  // # nn = 12
  // # nn2 = 56
  if (n_neigh == 0) n_neigh = 1;
  double ratio = static_cast<double>(same_species) / n_neigh;
  ratio = (1 - cos(ratio*M_PI/2.0));
  return ratio;
}






void AtomisticGenerator::make_supercell(double l1, double l2, double l3)
{

  //Build a supercell, defined by the lenght of conventional growth cell vectors
  int i,j,l;
  int start_i = 0, start_j = 0, start_l = 0;
  double conv_l1, conv_l2, conv_l3;
  Atom basis_atom;
  Tensor1 lattice_point;
  Tensor2Gen supercell_vect,inv_supercell_vect;
  Tensor1 tmp_check, tmp_conv;
  std::ostringstream os;

  std::vector<Atom> basis = _bulk->get_rotated_basis();
  std::vector<Atom>::const_iterator basis_iterator;
  std::vector<Tensor1>::iterator conv_iterator;
  
  Messages::debug("Running make_supercell...");

  //Check values. l1,l2,l3 cannot be unwisely large (no more than (1um)^3)
  assert((l1*l2*l3) < 1e+12);

  //Find lenght of conventional cell sides
  conv_l1 = sqrt(_conv_vect(1,1) * _conv_vect(1,1) + _conv_vect(2,1) * _conv_vect(2,1) 
                                                   + _conv_vect(3,1) * _conv_vect(3,1));
  conv_l2 = sqrt(_conv_vect(1,2) * _conv_vect(1,2) + _conv_vect(2,2) * _conv_vect(2,2) 
                                                   + _conv_vect(3,2) * _conv_vect(3,2));
  conv_l3 = sqrt(_conv_vect(1,3) * _conv_vect(1,3) + _conv_vect(2,3) * _conv_vect(2,3) 
                                                   + _conv_vect(3,3) * _conv_vect(3,3));

  int n1 = std::max(1, static_cast<int>(floor(l1 / conv_l1)));
  int n2 = std::max(1, static_cast<int>(floor(l2 / conv_l2)));
  int n3 = std::max(1, static_cast<int>(floor(l3 / conv_l3)));

  // override with integer factors from input file
  n1 = _as->get_options().get_option("cells_along_x", n1);
  n2 = _as->get_options().get_option("cells_along_y", n2);
  n3 = _as->get_options().get_option("cells_along_z", n3);

  os << "Conventional cells along x, y, z: " << n1 << " " << n2  << " " << n3;
  Messages::info(os.str());

  //Set supercell periodical vectors
  Tensor2Gen lmat(0);
  lmat(1,1) = (n1);
  lmat(2,2) = (n2);
  lmat(3,3) = (n3);

  _period = _conv_vect * lmat;

  if (!_as->is_periodic(0))
  {
    start_i = -2;
    n1 = n1 + 2;
  }
  if (!_as->is_periodic(1))
  {
    start_j = -2;
    n2 = n2 + 2;
  }
  if (!_as->is_periodic(2))
  {
    start_l = -2;
    n3 = n3 + 2;
  }

  //Definition of number of conventional cells, useful for reserving arrays
  unsigned int max_number_of_cells = n1*n2*n3;
  _super_basis.reserve(max_number_of_cells * _conv_lattice.size() * basis.size());
  
  Tensor1 origin(_local_origin +  _translation);

   cout<<"(debug) basis size: "<<basis.size()<<endl;
   cout<<"(debug) conv_latt size: "<<_conv_lattice.size()<<endl;

  //Need to construct a redundant supercell (for passivation purposes)
  //Note that it must be redundant only in non periodic directions
  for (i = start_i; i < n1; i++){
    for (j = start_j; j < n2; j++){
      for (l = start_l; l < n3; l++){

        conv_iterator = _conv_lattice.begin();

        //Fill conventional edges basis (super_conv)
        tmp_conv(1) = (i * _conv_vect(1,1)) + (j * _conv_vect(1,2)) + (l * _conv_vect(1,3));
        tmp_conv(2) = (i * _conv_vect(2,1)) + (j * _conv_vect(2,2)) + (l * _conv_vect(2,3));
        tmp_conv(3) = (i * _conv_vect(3,1)) + (j * _conv_vect(3,2)) + (l * _conv_vect(3,3));
        
        basis_iterator = _conv_basis.begin();

        do
        {
           basis_atom = (*basis_iterator);
           basis_atom.set_position( origin + tmp_conv + 
                (*basis_iterator).get_ttype_position() );
            _super_basis.push_back(basis_atom);

          ++basis_iterator;
        }
        while(basis_iterator != _conv_basis.end());
       
      }
    }
  }
  cout<<"(debug) super_basis size: "<<_super_basis.size()<<endl;

}



void AtomisticGenerator::make_conv_cell()
{
  Messages::info("Make conventional cell");

  //Calculate conventional cell vectors in the directions given by cut planes (conventional growth cell)
  Tensor1 vec_x(0),vec_y(0),vec_z(0);
  Tensor2Gen rotated_prim_vec = _bulk->get_rotated_prim_vec();

  // In general:
  // vj = e1 * _conv_prim(1,j) + e2 * _conv_prim(2,j) + e3 * _conv_prim(3,j)   
  //
  // for e1, e2, e3 rotated primitive vectors:
  // (|  |  | )(n11 n12 n13) (a  0  0)
  // (e1 e2 e3)(n21 n22 n23)=(0  b  0) 
  // (|  |  | )(n31 n32 n33) (0  0  c)
  // 
  // the matrik nij is computed first by computing  U = e^-1  
  // then the columns of U are scaled to integers 
  //
  _conv_prim = inv(rotated_prim_vec);
  scale_to_int(_conv_prim);
 
  // now the matrix nij = _conv_prim
  // so the orthogonal cell vectors are 
  _conv_vect = rotated_prim_vec * _conv_prim;

  //Note: we don't know if the vectors conv_vect are positively or
  //negatively oriented along the standard basis x,y,z.
  //For the way we build the supercell (going from edge_min to edge_max 
  //in positive x,y,z direction) we need positive conv_vect. 
  //IF the supercell is built along standard basis, we can simply take the
  //absolute value of the component of conv_vect
  for (int i = 1; i <=3; i++)
  {
     if (_conv_vect(i,i)<0)
     { 
       _conv_vect(i,i) = -_conv_vect(i,i);
       for (int j = 1; j <=3; j++)
          _conv_prim(j,i) = -_conv_prim(j,i);
     }
  }

  vec_x(1) = _conv_vect(1,1); vec_x(2) = _conv_vect(2,1); vec_x(3) = _conv_vect(3,1);
  vec_y(1) = _conv_vect(1,2); vec_y(2) = _conv_vect(2,2); vec_y(3) = _conv_vect(3,2);
  vec_z(1) = _conv_vect(1,3); vec_z(2) = _conv_vect(2,3); vec_z(3) = _conv_vect(3,3);
  
  //Check orthogonality of the final conventional vectors
  assert(((vec_x * vec_y) < 1e-10) && 
      ((vec_x * vec_z) < 1e-10) && 
      ((vec_y * vec_z) < 1e-10) &&
    (norm(vec_x) > 1e-10) &&
    (norm(vec_y) > 1e-10) &&
    (norm(vec_z) > 1e-10)); 

  // make a right-hended basis (for later use in fold to cell) 
  if (det(_conv_vect)<0)
  {
     //swap vec_y <-> vec_z
     _conv_vect(1,2) = vec_z(1); _conv_vect(2,2) = vec_z(2); _conv_vect(3,2) = vec_z(3); 
     _conv_vect(1,3) = vec_y(1); _conv_vect(2,3) = vec_y(2); _conv_vect(3,3) = vec_y(3); 
  } 

};


// -------------------------------------------------------------------------------
// Minimal conv_cell for 1D heterostructures.
// (|  |  | )(n11 n12 n13) (a |  | )
// (e1 e2 e3)(n21 n22 n23)=(0 v1 v2) 
// (|  |  | )(n31 n32 n33) (0 |  | )
// 
//  v1 = Sum_k{ek*nk1};  
//  v2 = Sum_k{ek*nk2};  
//
//  We minimize the area product: A = |v1 x v2|
//
//  A = | sum_kl{(ek x el) nk2 nl3} |
//  with constrains:
//  v3*y = 0  v3*z = 0 => fix all numbers nk1
//  
//  v1*x = 0  v2*x = 0
//  x*ek = ck : v1*x=0 =>  n12 = -(c2 n22 + c3 n32)/c1  
//  x*ek = ck : v2*x=0 =>  n13 = -(c2 n23 + c3 n33)/c1  
//  if c1==0 ?? 
//
//
void AtomisticGenerator::minimal_conv_cell()
{
  Tensor1 v1, v2, vF, prim_pos1, prim_pos2;
  double c[4], F, Fmin, AngMax, Ang, F1, F2;
  int i,j, k,l, ilow,jlow,klow,llow, lower_2, lower_3, upper_2, upper_3;
  Tensor2Gen rotated_prim_vec = _bulk->get_rotated_prim_vec();

  for (i = 1; i<=3; i++)
    c[i] = rotated_prim_vec(1,i);

  //Define a box including conventional cell (this should be the maximal area cell)
  Point p(_conv_vect(1,1), _conv_vect(2,2), _conv_vect(3,3));
  Tensor1 vect(p);  
  Tensor1 indxs = inv(rotated_prim_vec) * vect;
  
  lower_2 = int(min(indxs(2), min(_conv_prim(2,1), min(_conv_prim(2,2),_conv_prim(2,3)))));
  upper_2 = int(max(indxs(2), max(_conv_prim(2,1), max(_conv_prim(2,2),_conv_prim(2,3)))));
  lower_3 = int(min(indxs(3), min(_conv_prim(3,1), min(_conv_prim(3,2),_conv_prim(3,3)))));
  upper_3 = int(max(indxs(3), max(_conv_prim(3,1), max(_conv_prim(3,2),_conv_prim(3,3)))));
 
  lower_2 = int(min(0, lower_2));
  upper_2 = int(max(0, upper_2));
  lower_3 = int(min(0, lower_3));
  upper_3 = int(max(0, upper_3));
  // ----------------------------------------------------------------------------
  
  Fmin = 1e6;
  AngMax = 0.0; 

  for (int i = lower_2 ; i <= upper_2  ; i++)
  {
    for (int j = lower_3 ; j <= upper_3  ; j++)
    {
      prim_pos1(1) = double(-c[2]*i-c[3]*j)/c[1]; 
      prim_pos1(2) = double(i); 
      prim_pos1(3) = double(j);
      v1 = rotated_prim_vec * prim_pos1;
      F1 = norm(v1);

      for (int k = lower_2 ; k <= upper_2  ; k++)
      {
        for (int l = lower_3 ; l <= upper_3  ; l++)
        {
          prim_pos2(1) = double(-c[2]*k-c[3]*l)/c[1]; 
          prim_pos2(2) = double(k); 
          prim_pos2(3) = double(l); 
          
          v2 = rotated_prim_vec * prim_pos2;
          F2 = norm(v2);
          vF = v1^v2;
          F = norm(vF);
          Ang = asin(F/(F1*F2));

          if (vF(1)>0)
          {
            if (F<=Fmin) Fmin = F;
            if (F<=Fmin && Ang>AngMax)
            {
              ilow=i; jlow=j;klow=k;llow=l;
              AngMax = Ang;
              Fmin = F;
            }
          }
          else
          {
            continue;
          }
        }
      }    
      
    }
  }

  prim_pos1(1) = double(-c[2]*ilow-c[3]*jlow)/c[1];
  prim_pos1(2) = double(ilow); 
  prim_pos1(3) = double(jlow); 
  prim_pos2(1) = double(-c[2]*klow-c[3]*llow)/c[1];
  prim_pos2(2) = double(klow); 
  prim_pos2(3) = double(llow); 
  
  v1 = rotated_prim_vec * prim_pos1;
  v2 = rotated_prim_vec * prim_pos2;

  _conv_vect(1,2) = v1(1); _conv_vect(2,2) = v1(2); _conv_vect(3,2) = v1(3); 
  _conv_vect(1,3) = v2(1); _conv_vect(2,3) = v2(2); _conv_vect(3,3) = v2(3); 

  cout<<"conv_vect"<<endl;
  cout<<_conv_vect<<endl;
}

void AtomisticGenerator::make_conv_lattice()
{
  //Fill the conventional growth cell with atomic basis
  int lower[4],upper[4],i;
  Tensor1 prim_position, tmp_position;
  Tensor2Gen rotated_prim_vec = _bulk->get_rotated_prim_vec();

  /* IMPORTANT NOTE:
   * The lattice points covering the conventional cells are obtained by looking at the indeces
   * _conv_vect = rotated_prim_vec * _conv_prim
   * Quindi, poiche' _conv_vect= (ax,ay,az) dove ax = (a,0,0); ay=(0,b,0); az=(0,0,c)
   * ax = e1 * _conv_prim(1,1) + e2 * conv_prim(2,1) + e3 * conv_prim(3,1)   
   * ay = e1 * _conv_prim(1,2) + e2 * conv_prim(2,2) + e3 * conv_prim(3,2)   
   * az = e1 * _conv_prim(1,3) + e2 * conv_prim(2,3) + e3 * conv_prim(3,3)   
   * 
   * lo span della conv_cell si ha per i = min(0,_conv_prim(1,*)) : max(0,_conv_prim(1,*))
   *                                   j = min(0,_conv_prim(2,*)) : max(0,_conv_prim(2,*))
   *                                   k = min(0,_conv_prim(3,*)) : max(0,_conv_prim(3,*))
   *
   * Pero' in questo modo non e' garantito che anche il vertice (a,b,c) sia incluso!
   * Quindi vengono calcolati gli indici di questo da 
   * indxes = inv(rotated_prim_vec) * (a,b,c)
   * 
   * A questo punto, per costruzione anche (a,b,0) (0,b,c) e (a,0,c) 
   * sono inclusi. 
   */

  Tensor1 vect;
  vect(1) = _conv_vect(1,1)+_conv_vect(1,2)+_conv_vect(1,3);
  vect(2) = _conv_vect(2,1)+_conv_vect(2,2)+_conv_vect(2,3);
  vect(3) = _conv_vect(3,1)+_conv_vect(3,2)+_conv_vect(3,3);
  Tensor1 indxs = inv(rotated_prim_vec) * vect;

  Point a1,a2,a3,orig(-tol,-tol,-tol);
  a1(0) = _conv_vect(1,1); a1(1) = _conv_vect(2,1); a1(2) = _conv_vect(3,1);
  a2(0) = _conv_vect(1,2); a2(1) = _conv_vect(2,2); a2(2) = _conv_vect(3,2);
  a3(0) = _conv_vect(1,3); a3(1) = _conv_vect(2,3); a3(2) = _conv_vect(3,3);

  cout<<a1<<endl;
  cout<<a2<<endl;
  cout<<a3<<endl;
  //Define a box including conventional cell
 
  for (i=1; i<=3; i++)
  {
    double vals[]={0.0, indxs(i), _conv_prim(i,1), _conv_prim(i,2), _conv_prim(i,3)};
    lower[i] = int( *min_element(vals, vals+5) );
    upper[i] = int( *max_element(vals, vals+5) );
  }
  
  for (int i = lower[1]; i <= upper[1]; i++){
    for (int j = lower[2]; j <= upper[2]; j++){
      for (int l = lower[3]; l <= upper[3]; l++){

        prim_position(1) = double(i); 
        prim_position(2) = double(j); 
        prim_position(3) = double(l);
        tmp_position = rotated_prim_vec * prim_position;
        Point p(tmp_position(1), tmp_position(2), tmp_position(3));

        if (fold_in_cell(p,orig,a1,a2,a3,false)){
           _conv_lattice.push_back(tmp_position);
        }

      }
    }
  }

}

void AtomisticGenerator::move_origin()
{
  // Now we place the basis and translate the origin of the conv_cell in order to have all atoms
  // with positive coordinates. This is done to prevent loosing folded copies in periodic structures
  std::vector<Atom> basis = _bulk->get_rotated_basis();
  std::vector<Atom>::const_iterator basis_iterator;
  std::vector<Tensor1>::iterator conv_iterator = _conv_lattice.begin();      
  Tensor1 min_pos(1.0e10);
 
  do
  {
    basis_iterator = basis.begin();
    do
    {
      Tensor1 pos = (*basis_iterator).get_ttype_position() +
                   *conv_iterator;
      if (pos(1) < min_pos(1)) min_pos(1) = pos(1); 
      if (pos(2) < min_pos(2)) min_pos(2) = pos(2); 
      if (pos(3) < min_pos(3)) min_pos(3) = pos(3); 
      ++basis_iterator;
    }
    while(basis_iterator != basis.end());
    ++conv_iterator;
  }
  while(conv_iterator != _conv_lattice.end());
  
  _local_origin = _conv_vect * _cell_translation - min_pos;
 
  ostringstream os;
  os<< "New origin: "<<_local_origin;
  Messages msg;
  msg.info(os.str()); 

}

void AtomisticGenerator::make_conv_basis()
{
  Tensor1 atm_coords;
  std::vector<Atom> basis = _bulk->get_rotated_basis();
  std::vector<Tensor1>::iterator conv_iterator(_conv_lattice.begin());
  std::vector<Atom>::iterator basis_iterator;
  _conv_basis.reserve(_conv_lattice.size() * basis.size() );
  Atom basis_atom;
  Point a1, a2, a3, origP(0.0);
  a1(0) = _conv_vect(1,1); a1(1) = _conv_vect(2,1); a1(2) = _conv_vect(3,1);
  a2(0) = _conv_vect(1,2); a2(1) = _conv_vect(2,2); a2(2) = _conv_vect(3,2);
  a3(0) = _conv_vect(1,3); a3(1) = _conv_vect(2,3); a3(2) = _conv_vect(3,3);

  do
  {

    // add the basis to lattice point
    basis_iterator = basis.begin();
    do
    {
      basis_atom = (*basis_iterator);
      atm_coords = _local_origin + *conv_iterator +
                 basis_atom.get_ttype_position();
      basis_atom.set_position(atm_coords);

      // fold atom within the conventional cell
      fold_in_cell(basis_atom, origP, a1, a2, a3);
      //if (atm_coords(1)<0.0)             atm_coords(1) += _conv_vect(1,1);
      //if (atm_coords(1)>_conv_vect(1,1)) atm_coords(1) -= _conv_vect(1,1);

      //if (atm_coords(2)<0.0)             atm_coords(2) += _conv_vect(2,2);
      //if (atm_coords(2)>_conv_vect(2,2)) atm_coords(2) -= _conv_vect(2,2);

      //if (atm_coords(3)<0.0)             atm_coords(3) += _conv_vect(3,3);
      //if (atm_coords(3)>_conv_vect(3,3)) atm_coords(3) -= _conv_vect(3,3);
      // --------------------------------------
      
      _conv_basis.push_back(basis_atom);
      ++basis_iterator;

    }
    while(basis_iterator != basis.end());

    ++conv_iterator;
  }
  while(conv_iterator != _conv_lattice.end());

}


bool AtomisticGenerator::fold_in_cell(Point& position,
    const Point& orig, const Point& a1, const Point& a2, const Point& a3, bool fold)
{
  bool is_inside=true;
  bool pfold;
  // The planes vectors are oriented to point outside the cell
  // NOTE: below_surface() includes on the surface
  Point origin(orig);

  Plane p(origin, origin + a2, origin + a1);
  pfold = p.below_surface(position);
  if (!pfold && fold) position += a3;
  is_inside &= pfold;

  p.create_from_three_points(origin, origin + a3, origin + a2);
  pfold = p.below_surface(position);
  if (!pfold && fold) position += a1;
  is_inside &= pfold;

  p.create_from_three_points(origin, origin + a1, origin + a3);
  pfold = p.below_surface(position);
  if (!pfold && fold) position += a2;
  is_inside &= pfold;

  Point corner(origin + a1 + a2 + a3);

  p.create_from_three_points(corner, corner + a1, corner + a2);
  pfold = p.below_surface(position);
  if (!pfold && fold) position -= a3;
  is_inside &= pfold;

  p.create_from_three_points(corner, corner + a2, corner + a3);
  pfold = p.below_surface(position);
  if (!pfold && fold) position -= a1;
  is_inside &= pfold;

  p.create_from_three_points(corner, corner + a3, corner + a1);
  pfold = p.below_surface(position);
  if (!pfold && fold) position -= a2;
  is_inside &= pfold;

  return is_inside;
}

bool AtomisticGenerator::fold_in_cell(Atom& atom,
    const Point& orig, const Point& a1, const Point& a2, const Point& a3, bool fold)
{
  Point p = atom.get_position();
  bool is_in = fold_in_cell(p, orig, a1, a2, a3, fold);
  if (!is_in && fold) atom.set_position(p);
  return is_in;
}


// Increase periodicity in non-periodic directions 'periodic' flag
void AtomisticGenerator::check_periodic(void)
{
  Messages msg;
  msg.info("Check periodicity ... ");
  msg.indent();

  if (_as->is_periodic())
  {
    ostringstream os;
    os << "structure is periodic along";
    if (_as->is_periodic(0))
      os << " x";
    if (_as->is_periodic(1))
      os << " y";
    if (_as->is_periodic(2))
      os << " z";
    msg.info(os.str());
  }

  Tensor2Gen periods(0);
  periods(1,1) = 1;
  periods(2,2) = 1;
  periods(3,3) = 1;

  switch (_dim)
  {
    case 3:
      if (!_as->is_periodic(2))
        periods(3,3) *= 2;

    case 2:
      if (!_as->is_periodic(1))
        periods(2,2) *= 2;

    case 1:
      if (!_as->is_periodic(0))
        periods(1,1) *= 2;

    default:
      break;
  }

  _period = _period * periods;
  _as->set_ttype_lattice_vectors(_period);
  if (_bondmap == NULL) bond_map_gen(_super_basis);

  BondMap& bond_map = *_bondmap;

  Messages::info("Adding needed atoms for correct passivation");

  // maybe need map pair<i,j>
  std::vector<std::set<unsigned int>> added(3);

  //Warning: cycle end must be defined before as size will change dynamically during cycle
  //and we need acting only on already defined structure
  unsigned int size_before_passivating = _super_basis.size();

  for (unsigned int i = 0; i < size_before_passivating; i++)
  {
    if (_belong_to_structure[i])
    {

      for (unsigned int j = 0; j < bond_map[i].size(); j++)
      {

        if (( _belong_to_structure[bond_map[i][j]]) )
        {
          Atom* bonded_atom = &(_super_basis[bond_map[i][j]]);
          // is it a periodic one?
          Point trnsl(_bondmap->get_translation()[i][j]);
          Tensor1 shift;
          shift(1) = trnsl(0); shift(2) = trnsl(1); shift(3) = trnsl(2);
          double norm = shift(1)*shift(1) + shift(2)*shift(2) + shift(3)*shift(3);
          if (norm > 1e-9)
          {
            // check periodicity
            for (int d = 0; d < 3; ++d)
            {
              if (!_as->is_periodic(d))
              {
                if ((fabs(_period(d+1) * shift) > 1e-6) && !added[d].count(bond_map[i][j]))
                {
                  Atom tmp(*bonded_atom);
                  _belong_to_structure[bond_map[i][j]] = false;

                  Tensor1 position(bonded_atom->get_ttype_position() + shift);
                  tmp.set_position(position);
                  //std::cerr << "adding atom at " << position << "\n";
                  _super_basis.push_back(tmp);
                  added[d].insert(bond_map[i][j]);

                }

                // adjust bond map
              }
            }
          }
        }
      }
    }

  }


  msg.unindent();
  Messages::info("done");


}

//Bond map generation (cluster)
void  AtomisticGenerator::bond_map_gen(const std::vector<Atom>& basis){

  std::ostringstream os;

  //use internal member, if already used delete it
  if (_bondmap == NULL) _bondmap = new BondMap(basis.size());
  else
  {
    delete _bondmap;
    _bondmap = new BondMap(basis.size());
  }

  //--------------------------------------------------------------------------
  //os << "calling bond map with period "
  //<< _period(1,1)<<" "<<_period(2,2)<<" "<< _period(3,3) << std::endl;
  //Messages::debug(os.str());
  //os.str(std::string());
  //---------------------------------------------------------------------------

  _bondmap->do_solve(basis, _period);

};


void AtomisticGenerator::passivate()
{
  double hydrogen_distance = _as->get_options().get_option("hydrogen_distance", 1.2);
  Atom* bonded_atom;

  Messages::info("Starting passivation...");

  if (_bondmap == NULL) bond_map_gen(_super_basis);

  const BondMap& bond_map = *_bondmap;


  //Warning: cycle end must be defined before as size will change dynamically during cycle
  //and we need acting only on already defined structure
  unsigned int size_before_passivating = _super_basis.size();

  for (unsigned int i = 0; i < size_before_passivating; i++)
  {
    if ( _belong_to_structure[i] )
    {

      for (unsigned int j = 0; j < bond_map[i].size(); j++)
      {
     
        if ( ! _belong_to_structure[bond_map[i][j]] )
        {
          //TODO: using default copy constructor, with further modifications to
          //Atom class it could not work anymore!
          //Position must be modified in order to put Hydrogen atom near,
          //and also as we cannot have hydrogen bonded to more than one atom,
          //so in some cases we cannot keep crystal positions
          bonded_atom = &(_super_basis[bond_map[i][j]]);
          Atom tmp(*bonded_atom);
          tmp.set_specie("H");
          // Assign the element of the atom that belong to structure 
          tmp.set_elem(_super_basis[i].get_elem());  
          //NB: set label==0 as mark of passivation atom
          tmp.set_label(0);

          Point bonded_rel_position = bonded_atom->get_position() +
              _bondmap->get_translation()[i][j] - _super_basis[i].get_position();

          Point position = _super_basis[i].get_position() +
                         ( ( bonded_rel_position) /
                             bonded_rel_position.size() *
                                         hydrogen_distance); 

          tmp.set_position(position);

          // Passivation atoms are added to structure (not substituted).
          _super_basis.push_back(tmp);

          _belong_to_structure.push_back(true); 

        }

      }
    }

  }

  Messages::info("Cut & Passivation done");


}


//Some data manipulation function useful only in this class

Tensor2Gen
AtomisticGenerator::reciprocal(Tensor2Gen real_basis)
{

  //Build the reciprocal basis related to input 2-rank tensor
  Tensor1 a1(0),a2(0),a3(0);
  Tensor1 b1(0),b2(0),b3(0);
  Tensor1 select_vect(0);
  Tensor2Gen reciprocal(0);

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

  reciprocal(1,1)=b1(1);reciprocal(2,1)=b1(2);reciprocal(3,1)=b1(3);
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
  v_tmp = v / double(fabs(gcd_value));
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

  
void
AtomisticGenerator::build_random_alloy()
{
  using namespace std;
  // map between species and nominal concentrations 
  typedef map<Specie, double> specie_fraction;
  typedef map<Specie, unsigned int> specie_number;

  Messages::newline();
  Messages::info("Building random alloy ...");

  if (_as->is_random_alloy() == false)
    Messages::error("build_random_alloy is called but AtomisticStructure is not "
        "a random alloy");

  if (_bondmap == NULL)
  {
    std::cout<< "yet another bondmap..."<<std::endl;
    bond_map_gen(_structure_basis);
  }

  Messages m;
  m.indent();

  vector<string> cluster;
  _as->get_options().get_option("clustering", cluster);
  bool clustering = (cluster.size() > 0);

  vector<double> cluster_seeds(cluster.size(), 0.01);
  if (clustering)
    _as->get_options().get_option("cluster_seeds", cluster_seeds);

  bool fix_mean_alloy_concentration =
      _as->get_options().get_option("fix_mean_alloy_concentration", true);

  // build a map associating species to seeds
  map<Specie, double> rand_percentage; 
  for (int i=0; i < cluster_seeds.size(); i++)
  {
    Specie tmp(cluster[i]);
    rand_percentage[tmp] = cluster_seeds[i];
  }

  // A random starting seed is needed to actually have different sequences
  // we try to use something that is different also if launching simulations
  // at the same time
  int seed = _as->get_options().get_option("random_generator_seed",
      static_cast<int>(time(NULL) * std::tr1::random_device()()));
  {
    std::ostringstream os;
    os << "Initializing  MT19937 random generator with seed " << seed;
        //std::ios::hex << seed;
    Messages::info(os.str());
  }
  std::tr1::mt19937 generator(seed);


  // First build up a map with info per region and label between species and fractions
  vector<vector<specie_fraction>> frac;
  // define two vectors   num_to_substitute[region][label][Specie]
  vector<vector<specie_number>> num_to_substitute;
  vector<vector<specie_number>> num_substituted;
  // define vectors atm_to_substitute[region][label]
  vector<vector<unsigned int>> atm_to_substitute;
  vector<vector<bool>> done;

  // the number of regions is counted in excess as max between region ids + 1 
  int numregions = *(max_element(_as->get_IDset().begin(), _as->get_IDset().end())) + 1;
  frac.resize(numregions);
  num_to_substitute.resize(numregions);
  num_substituted.resize(numregions);
  atm_to_substitute.resize(numregions);
  done.resize(numregions);

  for (set<ID>::iterator reg = _as->get_IDset().begin(); reg != _as->get_IDset().end(); ++reg)
  {
    const Material* mat = _as->get_device()->get_material( (*reg) );
    unsigned int num_labels = mat->count_labels();
    num_to_substitute[*reg].resize(num_labels+1); //+1 to start from 1
    num_substituted[*reg].resize(num_labels+1);
    atm_to_substitute[*reg].resize(num_labels+1);
    frac[*reg].resize(num_labels+1);
    done[*reg].resize(num_labels+1);

    if (mat->is_alloy())
    {
      const Alloy* alloy = dynamic_cast<const Alloy*>(mat);
      for (unsigned int lb=1; lb<=num_labels; lb++)
      {
        const specie_fraction& speciemap = alloy->get_species_map(lb);
        frac[*reg][lb] = speciemap;
        done[*reg][lb] = false; 
      }
    }
    else
    {
      for (unsigned int lb=1; lb<=num_labels; lb++)
      {
        done[*reg][lb] = false; 
        Material::crystal_species_iterator it = mat->species_begin(lb);
        Material::crystal_species_iterator itend = mat->species_end(lb);
        for ( ; it != itend; ++it)
          frac[*reg][lb][*it] = 1.0;
      }
    }
  }

  /* Debug messages 
  for (set<ID>::iterator reg = _as->get_IDset().begin(); reg != _as->get_IDset().end(); ++reg)
  {
    std::cout<<"reg id: "<<*reg<<std::endl;
    for (unsigned int lb=1; lb<frac[*reg].size(); lb++)
    {
      specie_fraction::iterator it    = frac[*reg][lb].begin();
      specie_fraction::iterator itend = frac[*reg][lb].end();
      for( ; it != itend; ++it)
      {
        std::cout<<"lb: "<<lb<<" "<<(it->first).get_string()
                 <<" frac: "<<frac[*reg][lb][it->first]<<std::endl;
      }
    }
  }
  */
  
  // now we count how many atms with a given label 
  // may be substituted in each region
  // temporary use atm flag virtual_type to flag mutable atoms 
  for (unsigned int i = 0; i < _structure_basis.size(); i++)
  {
    Atom& atm = _structure_basis[i];
    unsigned int lb = static_cast<unsigned int>(atm.get_label());

    atm.set_type(0); 
    unsigned int regid = atm.get_region_ID();
    done[regid][lb] = true;
    // we skip passivation atoms 
    if (lb == 0) continue;

    const Material* mat = _as->get_material(atm);
    if (mat->is_alloy())
    {
      const Alloy* alloy = dynamic_cast<const Alloy*>(mat);
      if (alloy->is_mutable(lb))
      {
        atm_to_substitute[regid][lb] += 1;
        atm.set_type(1);
        done[regid][lb] = false;
        
        if (clustering)
        {
          specie_fraction::iterator it = frac[regid][lb].begin();
          for(; rand_percentage.count(atm.get_specie()) > 0; )
          {     
            atm.set_specie(it->first);
            if ( (++it) == frac[regid][lb].end() ) break;
          }
          //std::cout<<"Clustering, initial ion: "<<atm.get_specie().get_string()<<std::endl;
        }
      }
    }
  }
  
  // Loop to actually compute the number of atoms that need to be substituted.
  // this set is used to check if we have to do something in a region
  std::set<ID> not_finished;
  int total_subst = 0;

  for (set<ID>::iterator reg = _as->get_IDset().begin(); reg != _as->get_IDset().end(); ++reg)
  {
    for (unsigned int lb=1; lb < atm_to_substitute[*reg].size(); ++lb)
    {
      if (atm_to_substitute[*reg][lb] > 0)
      {
        // Compute number of atoms in each specie, ok if (fix_mean_alloy_concentration)
        int count = 0;
        specie_fraction::iterator it = frac[*reg][lb].begin(); 
        specie_fraction::iterator itend = frac[*reg][lb].end(); 
        for(; it != itend; ++it)
        { 
          // floor is used. Later the last ion is adjusted to ensure sum.
          num_to_substitute[*reg][lb][it->first] = 
                                     std::floor(it->second * atm_to_substitute[*reg][lb]);
          count += num_to_substitute[*reg][lb][it->first];
        }

        // here the number of ions can fluctuate in a multinomial way. 
        // Approximate by a normal distribution with 
        // mean=x*N and sigma=sqrt(N*x*(1-x))
        if (!fix_mean_alloy_concentration)
        {
          count = 0;
          specie_fraction::iterator it = frac[*reg][lb].begin(); 
          for(; it != itend; ++it)
          { 
             const double eps = std::numeric_limits<double>::min();
             double u1;
             do{  
               u1 = static_cast<double>(generator()) / generator.max();
             }while (u1 <= eps);
             double u2 = static_cast<double>(generator()) / generator.max();
             double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
             unsigned int N = num_to_substitute[*reg][lb][it->first];
             double x = frac[*reg][lb][it->first];
             double sigma = sqrt(N*x*(1-x));
             num_to_substitute[*reg][lb][it->first] = z * sigma + N;
             count += num_to_substitute[*reg][lb][it->first];
          }
        }
        
        // correct last specie to ensure sum rule
        specie_fraction::reverse_iterator rit = frac[*reg][lb].rbegin();  
        num_to_substitute[*reg][lb][rit->first] += (atm_to_substitute[*reg][lb] - count);

        // Write some output 
        double n_tot = atm_to_substitute[*reg][lb]; 
        total_subst += n_tot;
        it = frac[*reg][lb].begin(); 
        for(; it != itend; ++it)
        { 
          Specie sp(it->first);
          std::ostringstream os;
          os << "Region " << *reg << " Label "<<lb<<" Specie "
             << sp.get_string()<<": x = " << it->second << " -> " 
             << num_to_substitute[*reg][lb][sp] << " atoms out of " << n_tot 
             << " to be substituted (x_eff = " << std::setprecision(3) 
             << static_cast<float>(num_to_substitute[*reg][lb][sp]) / n_tot << ")" << std::endl;
          Messages::info(os.str());
        }
        not_finished.insert(*reg);
      }
    }    
  }

  
  Messages::info("Substituting atoms randomly");

  //
  // Now we extract random numbers between 0 and _structure_basis.size() - 1
  // to randomly pick an atom. If it is flagged as 1, and is in a region where
  // atoms need to be substituted, and is not already substituted it will be changed.
  // This is repeated until in all regions we have substituted the required number
  // of atoms.
  //
  std::tr1::uniform_int<size_t> random(0, _structure_basis.size() - 1);

  // We keep drowing numbers until a region is done (not_finished)
  size_t ctr = 0;
  size_t id = 0;
  for (; !not_finished.empty(); ++ctr)
  {
    id = random(generator);
    Atom& atm = _structure_basis[id];
    unsigned int regid = atm.get_region_ID();
    unsigned int lb = atm.get_label();

    //std::cout << "atom: "<<id<<" reg: "<<regid<<" lb: "<<lb<<"; fg: "
    //          << static_cast<unsigned int>(atm.get_type())
    //          << " is done? "<<done[regid][lb]<<std::endl;

    // NOTE: random numbers may repeat, so we have to check if this atom
    // has already been substituted: we use atom.type to flag substituted atoms
    if ( atm.get_type() && !done[regid][lb] )
    {

      map<Specie, double> new_frac(frac[regid][lb]);
      specie_fraction::iterator it = new_frac.begin(); 
      unsigned int nsp = new_frac.size();
      Specie clust_sp(it->first); // clustering specie 

      //for ( ; it != new_frac.end(); ++it)
      //{
      //  if (rand_percentage.count(it->first))
      //      std::cout<< (it->first).get_string()<<": "<<rand_percentage[it->first]<<" ";
      //}
      //std::cout<< std::endl;

      // CLUSTERING
      // rr is used in case of clustering to increase the prob
      // of the clustering specie w.r.t. to the others:
      //
      // x' = rr
      // y' = (1-rr)/(1-x) y
      // z' = (1-rr)/(1-x) z
      //
      // x'+y'+z' = 1
      //
      // NOTE 1: Clustering occurs only after some atoms have been placed at rnd. 
      // NOTE 2: First we need to place all clustering ions 
      if (clustering)
      {
        double rr=0.0;
        double ss=0.0;
        //changed=true when frac are actually redefined for clustering
        bool changed = false; 
        it = new_frac.begin();
        Specie sp(it->first);
        // first sets rr and ss for the clustering specie
        for ( ; it != new_frac.end(); ++it)
        {
          sp = it->first;
          // note: rand_percentage is used to see if a specie is a clustering specie
          if (rand_percentage.count(sp))  clust_sp = sp;

          if (rand_percentage.count(sp) && num_substituted[regid][lb][sp] > 
                      rand_percentage[sp] * num_to_substitute[regid][lb][sp])
          {
             rr = substitution_probability(id, sp);
             ss = it->second; 
             //std::cout<<"rr: "<<rr<<std::endl;
             new_frac[sp] = rr;
             changed = true;
             //std::cout<<"clustering: "<<sp.get_string()<<" "<<num_substituted[regid][lb][sp]
             //      <<"/"<<rand_percentage[sp] * num_to_substitute[regid][lb][sp]
             //      <<" : "<<num_to_substitute[regid][lb][sp] <<std::endl;
          }
        }

        it = new_frac.begin();
        for ( ; it != new_frac.end(); ++it)
        {
          sp = it->first;
          if (!rand_percentage.count(sp) &&  changed)
             new_frac[sp] = (1.0 - rr)/(1.0 - ss) * new_frac[sp]; 
        }
      }

      // Assign a specie according to fraction
      // we set one of the possible species at random
      // e.g., for quaternaries set sp according to fraction:
      //    x        y        1-x-y
      // |------|-------*-|-----------|
      // 0           rnd^             1
      // 
      it = new_frac.begin(); 
      Specie sp(it->first);
      double x = it->second;

      double rnd = static_cast<double>(generator()) / generator.max();
      //std::cout<<"nsp "<<nsp<<" x "<<x<<" "<<rnd<<std::endl;
      for (unsigned int i=0; i < nsp; i++)
      {
        if (rnd <= x)
        {
          sp = it->first;
          break;
        }
        ++it;
        x += it->second; 
      }

      // if (clustering)  => a specie can be substituted ONLY IF 
      // all clustering species have been already substituted.
      // This is done to avoid that non-clustering ions take away all slots,
      // restraining the clustering ions into forced positions 
      if ( !clustering || sp == clust_sp || num_substituted[regid][lb][clust_sp] 
                             >= num_to_substitute[regid][lb][clust_sp] )
      {
        if (num_substituted[regid][lb][sp] < num_to_substitute[regid][lb][sp])
        {
          atm.set_specie(sp);
          //std::cout<<sp.get_string()<<" subs "<<num_substituted[regid][lb][sp]
          //              <<"/"<<num_to_substitute[regid][lb][sp] << std::endl;
          num_substituted[regid][lb][sp] += 1;
          atm.set_type(0); // flag atom as substituted
        }
      }

      // checks whether a region/label is done
      if (num_substituted[regid][lb][sp] >= num_to_substitute[regid][lb][sp])
      {
         specie_number::iterator it = num_substituted[regid][lb].begin(); 
         specie_number::iterator itend = num_substituted[regid][lb].end(); 
         done[regid][lb] = true;
         for ( ; it != itend; ++it)
         {
            if (num_substituted[regid][lb][it->first] < 
                                  num_to_substitute[regid][lb][it->first])
            {
                done[regid][lb] = false;
                break;
            }
         }
         // remove a region when all of a lb are substituted
         bool alldone = true; 
         for (unsigned int i=1; i < done[regid].size(); i++)
         {
           //std::cout<<"reg: "<<regid<<" "<<done[regid][i]<<";"; 
           if (done[regid][i] == false){ alldone=false; break;}
         }
         //std::cout<<std::endl;
         if (alldone) not_finished.erase(regid);
      }  

    }
  }

  std::ostringstream os;
  os << "Needed " << ctr << " random number extractions to substitute " << total_subst << " atoms";
  Messages::info(os.str());
  Messages::newline();
          
}


void
AtomisticGenerator::build_random_old()
{
  Messages::newline();
  Messages::info("Building random alloy ...");

  if (_bondmap == NULL) bond_map_gen(_structure_basis);

  Messages m;
  m.indent();

  std::map<ID, std::map<unsigned char, Specie> > assignA;
  std::map<ID, std::map<unsigned char, Specie> > assignB;
  std::map<ID, double> a_to_b_prob;

  if (_as->is_random_alloy() == false)
    Messages::error("build_random_alloy is called but AtomisticStructure is not "
        "a random alloy");

  bool clustering = _as->get_options().get_option("clustering", false);
  double rand_percentage;
  if (clustering)
    rand_percentage = _as->get_options().get_option("cluster_seeds", 0.02);
  bool fix_mean_alloy_concentration =
      _as->get_options().get_option("fix_mean_alloy_concentration", true);

  // By default in VCA the specie assigned is the one of parent A
  // so we first swap all atoms and then change back to species A
  // I do this because then the rest of the code becomes
  // more intuitive.

  for (std::set<ID>::iterator reg = _as->get_IDset().begin(); reg != _as->get_IDset().end(); ++reg)
  {
    const Material* mat = _as->get_device()->get_material( (*reg) );

    if (mat->is_alloy())
    {
      const Alloy* alloy = dynamic_cast<const Alloy*>(mat);
      const Material* matA = alloy->get_component_A();
      const Material* matB = alloy->get_component_B();

      Database dbA = matA->get_database();
      Database dbB = matB->get_database();

      dbA.set_section("atomistic_structure");
      dbB.set_section("atomistic_structure");

      // Build up conversion map from file
      for (unsigned int i = 1; i <= dbA.get("n_basis_specie", 0); i++)
      {
        std::string record;
        std::string s;
        std::stringstream out;

        out << i;
        s = out.str();

        record = "specie_" + s;
        std::string db_record = dbA.get(record.c_str(),"none");
        assignA[*reg][static_cast<unsigned char>(i)] = Specie(db_record);
        // Note: probability to switch specie is 1-x
        a_to_b_prob[*reg] = mat->get_options().get_option("x", 1.0);
      }

      for (int i = 1; i <= dbB.get("n_basis_specie", 0); i++)
      {
        std::string record;
        std::string s;
        std::stringstream out;

        out << i;
        s = out.str();

        record = "specie_" + s;
        std::string db_record = dbB.get(record.c_str(),"none");
        assignB[*reg][static_cast<unsigned char>(i)] = Specie(db_record);
      }

    }
  }

  // A random starting seed is needed to actually have different sequences
  // we try to use something that is different also if launching simulations
  // at the same time
  int seed = _as->get_options().get_option("random_generator_seed",
      static_cast<int>(time(NULL) * std::tr1::random_device()()));
  {
    std::ostringstream os;
    os << "Initializing  MT19937 random generator with seed " << seed;
        //std::ios::hex << seed;
    Messages::info(os.str());
  }
  std::tr1::mt19937 generator(seed);


  //
  // First we count for each region the number of atoms that have to be substituted
  //

  // these vectors have space for all regions, so we can access them in a fast way
  std::vector<int> num_to_substitute(
      *(std::max_element(_as->get_IDset().begin(), _as->get_IDset().end())) + 1, 0);
  std::vector<int> num_substituted(num_to_substitute.size(), 0);


  for (unsigned int i = 0; i < _structure_basis.size(); i++)
  {
    Atom& atm = _structure_basis[i];

    // we skip passivation atoms
    if (atm.get_specie() == Specie::H)
      continue;

    // if it is flagged as first atom in the basis it may be substituted
    if (atm.get_label() == 1)
    {
      ID regid = atm.get_region_ID();
      const Material* mat = _as->get_material(atm);
      if (mat->is_alloy())
      {
        num_to_substitute[regid] += 1;
        // swap species
        atm.set_specie(assignB[regid][atm.get_label()]);
      }
    }
  }

  // this set is used to check if we have to do something in a region
  std::set<ID> not_finished;

  for (int i = 0; i < num_to_substitute.size(); ++i)
  {
    if (num_to_substitute[i] > 0)
    {
      // if we want to fix the number of atmos to be substituted,
      // we calculate this number now
      if (fix_mean_alloy_concentration)
      {
        double x = a_to_b_prob[i];
        // NOTE: we use floor() here to not have any fluctuation due to numerical
        // roundoff errors
        int n_tot = num_to_substitute[i];
        num_to_substitute[i] = std::floor(x * num_to_substitute[i]);
        std::ostringstream os;
        os << "Region " << i << ": x = " << x << " -> " << num_to_substitute[i] <<
            " atoms out of " << n_tot << " to be substituted (x_eff = " <<
            std::setprecision(3) <<
            static_cast<float>(num_to_substitute[i]) / n_tot << ")" << std::endl;
        Messages::info(os.str());
      }

      not_finished.insert(i);
    }
  }



  //
  // Now we extract random numbers between 0 and _structure_basis.size() - 1
  // to randomly pick an atom. If it is flagged as 1, and is in a region where
  // atoms need to be substituted, and is not already substituted it will be changed.
  // This is repeated until in all regions we have substituted the required number
  // of atoms.
  //
  std::tr1::uniform_int<size_t> random(0, _structure_basis.size() - 1);


  size_t ctr = 0;
  size_t id = 0;
  for (; !not_finished.empty(); ++ctr)
  {
    if (fix_mean_alloy_concentration)
      id = random(generator);
    else
    {
      // in this case we visit every atom
      id = ctr;
      if (id >= _structure_basis.size())
        break;
    }

    Atom& atm = _structure_basis[id];

    if ((atm.get_label() == 1) &&
        (atm.get_specie() != Specie::H))
    {
      ID regid = atm.get_region_ID();

      if (num_substituted[regid] < num_to_substitute[regid])
      {
        // NOTE: random numbers may repeat, so we have to check if this atom
        // has already been substituted!!
        Specie sp(assignA[regid][atm.get_label()]);

        if (atm.get_specie() != sp)
        {
          double prob = 1.0;
          double rnd = 0.0;
          // the first X% will be distributed randomly
          if (clustering && (num_substituted[regid] > rand_percentage * num_to_substitute[regid]))
          {
            prob = substitution_probability(id, sp);
            rnd = static_cast<double>(generator()) / generator.max();
          }
          if (!fix_mean_alloy_concentration)
          {
            prob = a_to_b_prob[regid];
            rnd = static_cast<double>(generator()) / generator.max();
          }

          if (rnd <= prob)
          {
            atm.set_specie(sp);
            ++num_substituted[regid];
          }
        }
      }
      else
      {  
        not_finished.erase(regid);
      }  

    }

  }

  if (fix_mean_alloy_concentration)
  {
    size_t subst = 0;
    for (int i = 0; i < num_substituted.size(); ++i)
      subst += num_substituted[i];

    std::ostringstream os;
    os << "Needed " << ctr << " random number extractions to substitute " << subst << " atoms";
    Messages::info(os.str());
    Messages::newline();
  }
  else
  {
    for (int i = 0; i < num_substituted.size(); ++i)
    {
      if (num_to_substitute[i] > 0)
      {
        std::ostringstream os;
        os << _as->get_device()->get_region_name(i) << " : substituted "
            << num_substituted[i] << " out of " << num_to_substitute[i]
            << " ( " << static_cast<double>(100 * num_substituted[i])
                             / num_to_substitute[i]
            << "%, nominally " << 100 * a_to_b_prob[i] << "%)";
        Messages::info(os.str());
      }
    }
    Messages::newline();
  }

}

