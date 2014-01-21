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



AtomisticGenerator::AtomisticGenerator(void)
:_bondmap(NULL),
_reference_material(NULL),
_conv_vect(0),
_conv_prim(0),
_local_origin(0),
_period(0),
_bulk(NULL)
{

}

AtomisticGenerator::~AtomisticGenerator(void)
{
  if (_bondmap != NULL) delete _bondmap;
}


const double AtomisticGenerator::tol = 1e-2;

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

    basis_iterator++;

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
  //if (!(_as->get_options().find_option("reference_region")))
  //{
  //  os << "No material could be set: reference_region is mandatory "
  //     <<  "in Atomistic section when no structure path is specified " << std::endl;
  //  Messages::warning(os.str());
  //}
  std::set<ID> ids;
  std::string ref_region;
  ref_region = _as->get_options().get_option("reference_region", "");
  _as->get_device()->get_active_region_ids(ref_region, ids);
  if (ids.size() != 0)
  {
    _reference_region_id = *ids.begin();
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
    _reference_region_id = *(_as->get_IDset().begin());
  }

  if (_reference_material == NULL)
    throw InitFailedException("Reference region/material badly "
        "defined for structure " + _as->get_name() );

  //Build the right BulkCrystal object
  //Additional options respect to the region should be specified here
  Messages::debug("(AG) Create bulk "+_reference_material->get_name());

  _bulk = BulkCrystal::create(_reference_material); 
  
  _bulk->do_init();

  //_bulk->print_xyb("Ref.xyb");
  //-----------------------------------------------------------------------
  
  //A translation vector can be specified to modify supercell alignment
  std::vector<double> translation (3,0.0);
  if ( _as->get_options().find_option("translation") )
  {
   _as->get_options().get_option("translation", translation);
   _local_origin(1) += translation[0]; 
   _local_origin(2) += translation[1]; 
   _local_origin(3) += translation[2];
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
  for ( ; el != el_end; el++)
  {
    Elem* elem = *el;
    if (_as->get_IDset().count(elem->subdomain_id()))
      ++num_elem;
  }
  _structure_elements.reserve(num_elem);
  for (el = _as->get_device()->get_mesh().elements_begin(); el != el_end; el++)
  {
    Elem* elem = *el;
    if (_as->get_IDset().count(elem->subdomain_id()))
      _structure_elements.push_back(elem);
  }

  //Build up supercell structure with proper options
  Messages::debug("(AG) Build supercell");
  build();
  
  //os<<"Initial structure with "<<_super_basis.size()<<" atoms"<<std::endl;
  //Messages::info(os.str());

  std::string preserve;
  preserve = _as->get_options().get_option("preserve", "none");
 
  // assign the right element to each atom
  assign_elements(_as->get_IDset());
  
  // cut the structure (only flags atoms)
  cut(_as->get_IDset(), preserve);

  // iterates on _super_basis and assign species
  assign_species();

  //eventually enlarge along dummy supercell directions
  check_periodic();

  if (_as->get_options().get_option("passivation", false))
    passivate();

  // remove unflagged atoms
  remove_atoms();

  // assign random-alloy species (only cations for the moment)
  if (_as->is_random_alloy())
    build_random_alloy();

  // delete generator bondmap
  delete _bondmap;
  _bondmap = NULL;

};


//Eliminate not included atoms from structure
//Check structure to eliminate unincluded atoms (using swap in another vector)
void
AtomisticGenerator::remove_atoms(void)
{
  _structure_basis.reserve(_super_basis.size());

  for (unsigned int i = 0; i < _super_basis.size(); i++)
  {
    if ((_super_basis[i].belong_to_structure))
    {
      _structure_basis.push_back(_super_basis[i]);
    }
  }

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
AtomisticGenerator::assign_elements(const std::set<ID>& reg_ids)
{

  // the tensor grid to real mesh mapper for fast association atom->Elem
  // NOTE: we pass the relevant ID set, since Quantum contacts could be present
  //       which have to be included in the atomistic structure
  MeshUtils::GridMapper& mapper =
      MeshUtils::GridMapper::get_mapper(_as->get_device()->get_mesh(), reg_ids);


  Utils::Progress prog("Assign elements",_super_basis.size());
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
    
    if (elem != NULL) (*atom).set_elem(elem);

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
    std::vector<Atom>::iterator atom = _super_basis.begin();
    for ( ; atom != _super_basis.end(); ++atom)
    {
      
      const Elem* elem = (*atom).get_elem();

      if (elem == NULL)
      {
        (*atom).belong_to_structure = false;
      }
      else
      {
        ID el_reg = elem->subdomain_id();
        
        (*atom).belong_to_structure = ( (reg_ids.count(el_reg) > 0) ?  true : false ); 
      }                            
    }
  }
  else
  {
    Messages::error("preserve conventional and primitive cell not implemented yet");
  }

}

void
AtomisticGenerator::assign_species(void)
{
  std::set<ID> reg_ids(_as->get_IDset());
  std::map<ID, std::map<unsigned char, Specie> > assign;

  // we will need the reference region for atoms falling outside of
  // the atomistic structure's regions
  reg_ids.insert(_reference_region_id);


  std::set<ID>::iterator reg(reg_ids.begin());
  for ( ; reg != reg_ids.end(); ++reg)
  {

    const Material* mat = _as->get_device()->get_material( (*reg) );

    Database db = mat->get_database();

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
      std::string db_record = db.get(record.c_str(),"none");
      assign[*reg][static_cast<unsigned char>(i)] = Specie(db_record);

    }

    //No more reading from section atomistic_structure in database are needed
    db.set_section("");

    //If some doping is present, a strategy must be studied and implemented here,
    // as we can choose arbitrarily species name (e.g. calling one doping Silicon Si1 and another one Si2)

  }


  //
  // Cycle upon all atoms and change specie according to assign map
  //
  int n_super_basis = _super_basis.size();
  //_structure_basis.clear();
  //_structure_basis.reserve(n_super_basis);
 
  Utils::Progress prog("Assign Species", n_super_basis);

  unsigned int progress_counter = 0;
  
  std::vector<Atom>::iterator atom = _super_basis.begin();
  for ( ; atom != _super_basis.end(); ++atom)
  {

    if ((*atom).belong_to_structure)
    {    
      ID el_reg = (*atom).get_elem()->subdomain_id();
      Specie tmp =  assign[el_reg][(*atom).get_label()];
      (*atom).set_specie(tmp);
      //_structure_basis.push_back(*atom);
    }
    else
    {
      //Species of atoms outside the structure are assigned for the bondmap
      Specie tmp = assign[_reference_region_id][(*atom).get_label()];
      (*atom).set_specie(tmp);
      //_structure_basis.push_back(*atom);
    }
    
    progress_counter += 1;
    prog.progress_message(progress_counter);
  }


}


void
AtomisticGenerator::restrict(bool passivation)
{
      

  _super_basis = _as->get_structure_atoms();

  _period = _as->get_ttype_lattice_vectors();


  cut(_as->get_IDset());

  // We must re-attach passivation Hydrogens
  

  if (passivation) passivate();  

  remove_atoms();

}





void
AtomisticGenerator::build_random_alloy()
{
  Messages::newline();
  Messages::info("Building random alloy ...");

  if (_bondmap == NULL) bond_map_gen(_structure_basis);

  Messages m;
  m.indent();

  std::map<ID, std::map<unsigned char, Specie> > assignA;
  std::map<ID, std::map<unsigned char, Specie> > assignB;
  std::map<ID, double> a_to_b_prob;
  bool done;

  if (_as->is_random_alloy() == false)
    Messages::error("build_random_alloy is called but AtomisticStructure is not "
        "a random alloy");

  bool clustering = _as->build_clusters();
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
    if (atm.belong_to_structure && (atm.get_label() == 1))
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

  // this we need for the substitution probability
  std::tr1::uniform_real<double> random2;

  size_t ctr = 0;
  size_t id = 0;
  for (; !not_finished.empty(); ++ctr)
  {
    if (fix_mean_alloy_concentration)
      id = random(generator);
    else
      id = ctr;

    Atom& atm = _structure_basis[id];

    if (atm.belong_to_structure &&
        (atm.get_label() == 1) &&
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


double
AtomisticGenerator::substitution_probability(size_t id, const Specie& sp)
{
  const BondMap& bm = *_bondmap;

  int same_species = 1;
  int n_neigh = 1;

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

        const std::vector<unsigned int>& nn2 = bm[neigh[j]];
        for (unsigned int ii = 0; ii < nn2.size(); ++ii)
        {
          const std::vector<unsigned int>& nn3 = bm[neigh[ii]];
          for (unsigned int jj = 0; jj < nn2.size(); ++jj)
          {
            if (!visited.count(nn[jj]))
            {
              n_neigh++;
              visited.insert(nn[jj]);
              if (_structure_basis[nn[jj]].get_specie() == sp)
                same_species++;
            }
          }
        }
      }
    }
  }

  n_neigh = 13;
  double ratio = static_cast<double>(same_species) / n_neigh;
  ratio = (ratio < 1) ? ratio : 1;
  //ratio = 1 - (1 - ratio)*(1 - ratio)*(1 - ratio);
  return ratio;
}


void AtomisticGenerator::make_supercell(double l1, double l2, double l3){

  //Build a supercell, defined by the lenght of conventional growth cell vectors
  std::vector<Tensor1>::iterator conv_iterator;
  int i,j,l;
  int n1,n2,n3,start_i,start_j,start_l;
  double conv_l1, conv_l2, conv_l3;
  Atom basis_atom;
  Tensor1 lattice_point;
  Tensor2Gen supercell_vect,inv_supercell_vect;
  Tensor1 tmp_check, tmp_conv;
  std::ostringstream os;

  std::vector<Atom> basis = _bulk->get_rotated_basis();
  std::vector<Atom>::const_iterator basis_iterator = basis.begin();


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

  n1 = int(floor(l1 / conv_l1)); n2 = int(floor(l2 / conv_l2)); n3 = int(floor(l3 / conv_l3));

  _conv_cells_supercell_lenght[0] = n1 + 1;
  _conv_cells_supercell_lenght[1] = n2 + 1;
  _conv_cells_supercell_lenght[2] = n3 + 1;

  l1 = (n1 + 1) * conv_l1;
  l2 = (n2 + 1) * conv_l2;
  l3 = (n3 + 1) * conv_l3;

  //Set supercell periodical vectors
  Tensor2Gen lmat(0);

  lmat(1,1) = (n1 + 1);
  lmat(2,2) = (n2 + 1);
  lmat(3,3) = (n3 + 1);

  _period = _conv_vect * lmat;

  //----------------------------------------------------
  //os << "in make_conv_cell period is "
  //<< _period(1,1)<<" "<< _period(2,2)<<" "<< _period(3,3) << std::endl;
  //Messages::info(os.str());
  //os.str(std::string());
  //----------------------------------------------------

  //Define vectors with same direction of conventional cell vectors, but with size specifed by l1,l2,l3
  supercell_vect(1,1) = _conv_vect(1,1) * (l1 / conv_l1);
  supercell_vect(2,1) = _conv_vect(2,1) * (l1 / conv_l1);
  supercell_vect(3,1) = _conv_vect(3,1) * (l1 / conv_l1);
  supercell_vect(1,2) = _conv_vect(1,2) * (l2 / conv_l2);
  supercell_vect(2,2) = _conv_vect(2,2) * (l2 / conv_l2);
  supercell_vect(3,2) = _conv_vect(3,2) * (l2 / conv_l2);
  supercell_vect(1,3) = _conv_vect(1,3) * (l3 / conv_l3);
  supercell_vect(2,3) = _conv_vect(2,3) * (l3 / conv_l3);
  supercell_vect(3,3) = _conv_vect(3,3) * (l3 / conv_l3);
  inv_supercell_vect = inv(supercell_vect);

  if (_dim == 1) {start_i = -2; start_j = 0; start_l = 0; n1 = n1 + 2;}
  if (_dim == 2) {start_i = -2; start_j =-1; start_l = 0; n1 = n1 + 2; n2 = n2 + 2;}
  if (_dim == 3) {start_i = -2; start_j =-2; start_l = -2; n1 = n1 + 2; n2 = n2 + 2; n3 = n3 + 2;}

  //Definition of number of conventional cells, useful for reserving arrays
  unsigned int max_number_of_cells = n1 + n2 + n3 + 6;
  _super_conv.reserve(max_number_of_cells);
  _super_lattice.reserve(max_number_of_cells * _conv_lattice_basis.size());
  _super_basis.reserve(max_number_of_cells * _conv_lattice_basis.size() * basis.size());

  //Need to construct a redundant supercell (for passivation purposes)
  //Note that it must be redundant only in non periodic directions
  for (i = start_i; i <= n1; i++){
    for (j = start_j; j <= n2; j++){
      for (l = start_l; l <= n3; l++){

        conv_iterator = _conv_lattice_basis.begin();

        //Fill conventional edges basis (super_conv)
        //Note: we don't know if the vectors conv_vect are positively or
        //negatively oriented along the standard basis x,y,z.
        //For the way we build the supercell (going from edge_min to edge_max 
        //in positive x,y,z direction) we need positive conv_vect. 
        //IF the supercell is built along standard basis, we can simply take the
        //absolute value of the component of conv_vect
        tmp_conv(1) = (i * _conv_vect(1,1)) + (j * _conv_vect(1,2)) + (l * _conv_vect(1,3));
        tmp_conv(2) = (i * _conv_vect(2,1)) + (j * _conv_vect(2,2)) + (l * _conv_vect(2,3));
        tmp_conv(3) = (i * _conv_vect(3,1)) + (j * _conv_vect(3,2)) + (l * _conv_vect(3,3));
        _super_conv.push_back(tmp_conv + _local_origin);

        //        }

        do
        {
          //Assign lattice point position
          lattice_point(1) = (*conv_iterator)(1) + (i * _conv_vect(1,1)) 
                                                 + (j * _conv_vect(1,2)) 
                                                 + (l * _conv_vect(1,3));
          lattice_point(2) = (*conv_iterator)(2) + (i * _conv_vect(2,1)) 
                                                 + (j * _conv_vect(2,2)) 
                                                 + (l * _conv_vect(2,3));
          lattice_point(3) = (*conv_iterator)(3) + (i * _conv_vect(3,1)) 
                                                 + (j * _conv_vect(3,2)) 
                                                 + (l * _conv_vect(3,3));

          //Put lattice point into supercell lattice points array
          _super_lattice.push_back(lattice_point + _local_origin);
          basis_iterator = basis.begin();

          do
          {
            basis_atom = (*basis_iterator);
            basis_atom.set_position ( _local_origin + lattice_point +
                (*basis_iterator).get_ttype_position() );
            _super_basis.push_back(basis_atom);
            ++basis_iterator;

          }
          while(basis_iterator != basis.end());


          conv_iterator++;

        }
        while(conv_iterator != _conv_lattice_basis.end());
        
      };
    };
  };

};



void AtomisticGenerator::make_conv_cell()
{
  //Calculate conventional cell vectors in the directions given by cut planes (conventional growth cell)
  Tensor1 m1,m2,m3,select_vect(0);
  Tensor1 conv1, conv2, conv3;
  Tensor2Gen rotated_prim_vec = _bulk->get_rotated_prim_vec();
  
  //std::cout<<"rotated_prim_vec:"<<std::endl;
  //std::cout<<rotated_prim_vec<<std::endl;

  _conv_prim = inv(rotated_prim_vec);

  //std::cout<<"_conv_prim:"<<std::endl;
  //std::cout<<_conv_prim<<std::endl;

  scale_to_int(_conv_prim);
  _conv_vect = rotated_prim_vec * _conv_prim;

  //Note: we don't know if the vectors conv_vect are positively or
  //negatively oriented along the standard basis x,y,z.
  //For the way we build the supercell (going from edge_min to edge_max 
  //in positive x,y,z direction) we need positive conv_vect. 
  //IF the supercell is built along standard basis, we can simply take the
  //absolute value of the component of conv_vect. Conv_vect defines
  //a lattice vectors, therefore if magnitude and angles are preserved, the
  //definition is still legit
  for (int i = 1; i <=3; i++)
  {
     for (int j = 1; j <=3; j++)
     {
       _conv_vect(i,j) = fabs(_conv_vect(i,j));
     }
  }

};


void AtomisticGenerator::make_conv_basis()
{
  //Fill the conventional growth cell with atomic basis
  int lower_1, lower_2, lower_3, upper_1, upper_2, upper_3, i;
  Tensor1 prim_position, tmp_check;
  Tensor1 tmp_position;
  Tensor1 vec_x(0),vec_y(0),vec_z(0);
  std::vector<Atom>::const_iterator basis_iterator;
  Tensor2Gen rotated_prim_vec = _bulk->get_rotated_prim_vec();
  
  //Make a preliminar rotation only if conventional cell vectors are orthogonal

  vec_x(1) = _conv_vect(1,1); vec_x(2) = _conv_vect(2,1); vec_x(3) = _conv_vect(3,1);
  vec_y(1) = _conv_vect(1,2); vec_y(2) = _conv_vect(2,2); vec_y(3) = _conv_vect(3,2);
  vec_z(1) = _conv_vect(1,3); vec_z(2) = _conv_vect(2,3); vec_z(3) = _conv_vect(3,3);

  //Check orthogonality
  assert(((vec_x * vec_y) < 1e-10) && 
      ((vec_x * vec_z) < 1e-10) && 
      ((vec_y * vec_z) < 1e-10) &&
    (norm(vec_x) > 1e-10) &&
    (norm(vec_y) > 1e-10) &&
    (norm(vec_z) > 1e-10)); 

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

        prim_position(1) = double(i); 
        prim_position(2) = double(j); 
        prim_position(3) = double(l);
        tmp_position = prim_position;
        bool check_boundary;
        tmp_check = inv(_conv_prim) * tmp_position;
        check_boundary= ((tmp_check(1) >= -tol) && (tmp_check(1) < (1.0 - tol)))&&
            ((tmp_check(2) >= -tol) && (tmp_check(2) < (1.0 - tol))) &&
            ((tmp_check(3) >= -tol) && (tmp_check(3) < (1.0 - tol)));

        if (check_boundary){
          tmp_position = rotated_prim_vec * prim_position;
          _conv_lattice_basis.push_back(tmp_position);

        }

      };
    };
  };

};

// Increase periodicity in non-periodic directions 'periodic' flag
void AtomisticGenerator::check_periodic(void)
{

  Tensor2Gen periods(0);

  periods(1,1) = 1;
  periods(2,2) = 1;
  periods(3,3) = 1;

  bool periodic = _as->get_options().get_option("periodic", false);
  if (_as->get_options().get_option("passivation", false)) periodic = false;

  if ( !periodic )
  {
    switch (_dim)
    {
      case 3:
        periods(3,3) *= 2;

      case 2:
        periods(2,2) *= 2;

      case 1:
        periods(1,1) *= 2;

      default:
        break;
    }
  }

  _period = _period * periods;

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

  Messages::info("Building Bond Map...");
  _bondmap->do_solve(basis, _period);
  Messages::info("Bond Map completed");

};


void AtomisticGenerator::passivate()
{
  Tensor1 position;
  double hydrogen_distance = 1.2;
  Atom* bonded_atom;

  Messages::info("Starting passivation...");

  if (_bondmap == NULL) bond_map_gen(_super_basis);

  const BondMap& bond_map = *_bondmap;


  //Warning: cycle end must be defined before as size will change dynamically during cycle
  //and we need acting only on already defined structure
  unsigned int size_before_passivating = _super_basis.size();

  for (unsigned int i = 0; i < size_before_passivating; i++)
  {
    if (_super_basis[i].belong_to_structure)
    {

      for (unsigned int j = 0; j < bond_map[i].size(); j++)
      {
     
        bonded_atom = &(_super_basis[bond_map[i][j]]);
        if (!((*bonded_atom).belong_to_structure))
        {
          //TODO: using default copy constructor, with further modifications to
          //Atom class it could not work anymore!
          //Position must be modified in order to put Hydrogen atom near,
          //and also as we cannot have hydrogen bonded to more than one atom,
          //so in some cases we cannot keep crystal positions
          Atom tmp(*bonded_atom);
          tmp.set_specie("H");
          tmp.set_label((*bonded_atom).get_label());
          tmp.belong_to_structure = true;
          Tensor1 bonded_rel_position = bonded_atom->get_ttype_position() +
              _bondmap->get_translation()[i][j] - _super_basis[i].get_ttype_position();

          position = _super_basis[i].get_ttype_position() +
              ( ( bonded_rel_position) /
                  norm(bonded_rel_position ) *
                  hydrogen_distance);
          tmp.set_position(position);

          _super_basis.push_back(tmp);

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




/*
void
AtomisticGenerator::cut_and_change_specie(std::string preserve)
{


  Messages::debug("Running cut_and_change_specie");
  
  std::map<ID, std::map<unsigned char, Specie> > assign;
  assign.clear();

  std::set<ID> reg_ids(_as->get_IDset());

  // we will need the reference region for atoms falling outside of
  // the atomistic structure's regions
  reg_ids.insert(_reference_region_id);

  std::set<ID>::iterator reg(reg_ids.begin());
  for ( ; reg != reg_ids.end(); ++reg)
  {

    const Material* mat = _as->get_device()->get_material( (*reg) );

    Database db = mat->get_database();

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
      std::string db_record = db.get(record.c_str(),"none");
      assign[*reg][static_cast<unsigned char>(i)] = Specie(db_record);

    }


    //No more reading from section atomistic_structure in database are needed
    db.set_section("");

    //If some doping is present, a strategy must be studied and implemented here,
    // as we can choose arbitrarily species name (e.g. calling one doping Silicon Si1 and another one Si2)

  }

  //
  // Cycle upon all atoms and change specie according to assign map
  //
  int n_super_basis = _super_basis.size();
  _structure_basis.clear();
  _structure_basis.reserve(n_super_basis);


  // the tensor grid to real mesh mapper for fast association atom->Elem
  // NOTE: we pass the relevant ID set, since Quantum contacts could be present
  //       which have to be included in the atomistic structure
  MeshUtils::GridMapper& mapper =
      MeshUtils::GridMapper::get_mapper(_as->get_device()->get_mesh(),
          _as->get_IDset());

  Utils::Progress progress("Atomistic Generator", n_super_basis);
  //unsigned int progress_step = (n_super_basis > 100) ? n_super_basis / 100 : 1;
  unsigned int progress_counter = 0;
 
  //Different strategies if preserving conventional cell or preserving basis are needed
  if (preserve.compare("none") == 0)
  {
    for (std::vector<Atom>::iterator atom = _super_basis.begin();
        atom != _super_basis.end(); ++atom)
    {

      Point p((*atom).get_position());
      p *= 1.0 / scale;

      // set unneeded dimensions to 0, so the atoms are associated to the correct
      // elements
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
      bool done = false;

      if (elem != NULL)
      {
        ID el_reg = elem->subdomain_id();

        if (reg_ids.count(el_reg))
        {

          Specie tmp =  assign[el_reg][(*atom).get_label()];
          (*atom).set_specie(tmp);
          (*atom).belong_to_structure = true;
          (*atom).set_elem(elem);
          _structure_basis.push_back(*atom);
          done = true;
        }
      }

      if (!done)
      {
        Specie tmp =  assign[_reference_region_id][(*atom).get_label()];
        (*atom).set_specie(tmp);
        (*atom).belong_to_structure = false;
        (*atom).set_elem(elem);
        _structure_basis.push_back(*atom);
      }

      progress_counter += 1;
      progress.progress_message(progress_counter);
    }
  }


  if (preserve.compare("lattice") == 0)
  {
    std::vector<Atom> basis = _bulk->get_rotated_basis();
    Atom tmp_atom;

    std::vector<Tensor1>::iterator lattice = _super_lattice.begin();
    for ( ; lattice != _super_lattice.end(); lattice++)
    {

      Point p((*lattice)(1), (*lattice)(2), (*lattice)(3));
      p *= 1.0 / scale;

      // set unneeded dimensions to 0, so the atoms are associated to the correct
      // elements
      switch (_dim)
      {
        case 0:
          p(0) = 0;
        case 1:
          p(1) = 0;
        case 2:
          p(2) = 0;
        default:
          break;
      }



      const Elem* elem = mapper.get_element(p);
      bool done = false;

      if (elem != NULL)
      {
        ID el_reg = elem->subdomain_id();

        if (reg_ids.count(el_reg))
        {

          for ( std::vector<Atom>::const_iterator atom = basis.begin();
              atom != basis.end(); atom++)
          {

            tmp_atom.set_position((*lattice) + (*atom).get_ttype_position());

            Specie tmp =  assign[el_reg][(*atom).get_label()];
            tmp_atom.set_specie(tmp);
            tmp_atom.belong_to_structure = true;
            // TODO what if it falls out of the current element?
            tmp_atom.set_elem(elem);

            _structure_basis.push_back(tmp_atom);
          }

          done = true;
        }
      }

      if (!done)
      {
        for ( std::vector<Atom>::const_iterator atom = basis.begin();
            atom != basis.end(); atom++)
        {

          tmp_atom.set_position((*lattice) + (*atom).get_ttype_position());

          Specie tmp = assign[_reference_region_id][(*atom).get_label()];
          tmp_atom.set_specie(tmp);
          tmp_atom.belong_to_structure = false;
          _structure_basis.push_back(tmp_atom);
        }
      }
    }
  }

  if (preserve.compare("conventional") == 0)
  {
    std::vector<Atom> basis = _bulk->get_rotated_basis();
    Atom tmp_atom;
    for (std::vector<Tensor1>::iterator conv = _super_conv.begin();
        conv != _super_conv.end(); conv++)
    {

      Point p((*conv)(1), (*conv)(2), (*conv)(3));
      p *= 1.0 / scale;

      // set unneeded dimensions to 0, so the atoms are associated to the correct
      // elements
      switch (_dim)
      {
        case 0:
          p(0) = 0;
        case 1:
          p(1) = 0;
        case 2:
          p(2) = 0;
        default:
          break;
      }



      const Elem* elem = mapper.get_element(p);
      bool done = false;

      if (elem != NULL)
      {
        ID el_reg = elem->subdomain_id();

        if (reg_ids.count(el_reg))
        {
          done = true;

          for ( std::vector<Tensor1>::iterator conv_lattice_basis_it = _conv_lattice_basis.begin();
              conv_lattice_basis_it != _conv_lattice_basis.end(); conv_lattice_basis_it++)
          {
            for ( std::vector<Atom>::const_iterator atom = basis.begin();
                atom != basis.end(); atom++)
            {
              tmp_atom.set_position( (*conv) + (*conv_lattice_basis_it) +
                  (*atom).get_ttype_position());
              Specie tmp =  assign[el_reg][(*atom).get_label()];
              tmp_atom.set_specie(tmp);
              tmp_atom.belong_to_structure = true;
              // TODO what if it falls out of the current element?
              tmp_atom.set_elem(elem);


              _structure_basis.push_back(tmp_atom);
            }
          }
        }
      }

      if (!done)
      {
        for ( std::vector<Tensor1>::iterator conv_lattice_basis_it = _conv_lattice_basis.begin();
            conv_lattice_basis_it != _conv_lattice_basis.end(); conv_lattice_basis_it++)
        {
          for ( std::vector<Atom>::const_iterator atom = basis.begin();
              atom != basis.end(); atom++)
          {
            tmp_atom.set_position( (*conv) + (*conv_lattice_basis_it) +
                (*atom).get_ttype_position());

            Specie tmp =  assign[_reference_region_id][(*atom).get_label()];
            tmp_atom.set_specie(tmp);
            tmp_atom.belong_to_structure = false;
            _structure_basis.push_back(tmp_atom);
          }
        }
      }
    }
  }

  Messages::newline();
  Messages::debug("Finished cut_and_change_specie");

}
*/
