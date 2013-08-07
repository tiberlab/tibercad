#include "StrainLattice.h"
#include "TypeDefs.h"
#include "Messages.h"
#include "Material.h"
#include "Alloy.h"
#include "AtomisticStructure.h"
#include "BulkCrystal.h"
#include "RuntimeException.h"
#include "AtomisticBasis.h"
#include "Specie.h"
#include "Atom.h"
#include "mesh.h"

StrainLattice::StrainLattice()
:_as(NULL)
{
}


void
StrainLattice::init(AtomisticStructure* as)
{ 
  Messages::info("Initializing atomistic strain projection");
  _as = as;
  fill_materials_set();
  //! Build reference for every material
  std::set<const Material*>::iterator it;
  for (it = _materials.begin(); it != _materials.end(); ++it)
  {
  BulkCrystal* bulk = BulkCrystal::create(*it);
  bulk->do_init();
  const AtomisticBasis* structure =  bulk; 
  assert(bulk->get_structure_atoms().size() > 0);

  //Center the reference on the cation
  unsigned int ref_center = 0;
  Atom ref_atm = bulk->get_structure_atoms()[0]; 
  //Note that you dereference the iterator to get a pointer
  if ((*it)->is_cation(ref_atm.get_specie()))
    _reference[*it] = build_tetraedron(structure, 0);
  else 
  {
    ref_atm = bulk->get_structure_atoms()[1]; 
    if ((*it)->is_cation(ref_atm.get_specie()))
        _reference[*it] = build_tetraedron(structure, 1);
    else
      Messages::error("Could not find cation in StrainLattice");
  }

  delete bulk; 
  }
}


void 
StrainLattice::fill_materials_set(void)
{
  std::set<ID>::iterator begin = _as->get_IDset().begin(); 
  std::set<ID>::iterator end = _as->get_IDset().end();
  const Material* mat = NULL;
  const Alloy* alloy = NULL;

  for (std::set<ID>::iterator it = begin; it != end; ++it)
  {
    mat = _as->get_device()->get_material(*it);
    _materials.insert(mat);
    //
    //Up to now I don't use the parent materials. I leave these line
    //commented for reference, I can delete them in a while
    //
    //if (mat->is_alloy())
    //{  
    //  alloy = dynamic_cast<const Alloy*>(mat);
    //  _materials.insert(alloy->get_component_A());
    //  _materials.insert(alloy->get_component_B());
    //}
  }
}


void
StrainLattice::do_solve(void)
{
  Tensor2Gen strain;
  Tensor2Gen identity(1);
  Tetra tet, tmp;
  const Bondmap& bondmap = _as->get_bond_map();
  const std::vector<Atom>& atoms = _as->get_structure_atoms();
  _solution.reserve(_as->get_N_without_H() / 2 + 1);
  TensorField sol;

  //!Do the job on any atom which is compatible with reference tetraedra
  for (unsigned int ind = 0; ind < _as->get_N_without_H(); ind++)
  {
    const Atom& atm = atoms[ind];
    const Tetra& ref = _reference[_as->get_material(atm)];
    //If bonded to hydrogen or not 4-coordinate, discharge
    if (bondmap[ind].size() != 4) continue;
    if ((atoms[bondmap[ind][0]].get_specie() == Specie::H) ||
       (atoms[bondmap[ind][1]].get_specie() == Specie::H) ||
       (atoms[bondmap[ind][2]].get_specie() == Specie::H) ||
       (atoms[bondmap[ind][3]].get_specie() == Specie::H)) continue;

    //Check if it's the right specie, i.e. not reference vertex
    if (atm.get_specie() == ref.vertex_sp) continue;

    //If previuos conditions were not fullfilled, we should have a valid one, we calculate the strain.
    tmp = build_tetraedron(_as, ind);
    tet = rearrange(ref, tmp);
    strain = tet.edges * inv(ref.edges) - identity;
    //Dirty hack to deal with WZ double basis (you have 2 tetraedra there!!)
    //ONLY TEMPORARY, and only good in the direction I'm using for the calculations!!!
    if (norm(strain) > 0.5)
    {
    if (_as->get_device()->get_mesh().mesh_dimension() == 3)
    {
      tmp.bonds[0](1) = -1.0 * tmp.bonds[0](1);
      tmp.bonds[1](1) = -1.0 * tmp.bonds[1](1);
      tmp.bonds[2](1) = -1.0 * tmp.bonds[2](1);
      tmp.bonds[3](1) = -1.0 * tmp.bonds[3](1);
    }
    if (_as->get_device()->get_mesh().mesh_dimension() == 1)
    {
      tmp.bonds[0](2) = -1.0 * tmp.bonds[0](2);
      tmp.bonds[1](2) = -1.0 * tmp.bonds[1](2);
      tmp.bonds[2](2) = -1.0 * tmp.bonds[2](2);
      tmp.bonds[3](2) = -1.0 * tmp.bonds[3](2);
    }

      tet = rearrange(ref, tmp);
      strain = tet.edges * inv(ref.edges) - identity;
    }  
    //END HACK


    sol.tensor = strain;
    sol.atom_p = &atm;
    _solution.push_back(sol);
  }

}


StrainLattice::Tetra
StrainLattice::build_tetraedron(const AtomisticBasis* as, unsigned int atm)
{
  StrainLattice::Tetra tet;
  Tensor1 bond(0);
  Tensor1 edge;
  const std::vector<Atom> atoms = as->get_structure_atoms();

  //!Store first specie as reference (it should be the cation 
  //! in III-V alloys)
  tet.sp = atoms[atm].get_specie();

  //Bulid bonds and tetraedron
  const Bondmap& bondmap = as->get_bond_map();
  const std::vector<std::vector<Tensor1> > translation =
      as->get_neighbor_translation();
  unsigned int n_bonds = bondmap[atm].size();
  if (n_bonds != 4)
    Messages::error("Bulk badly defined in StrainLattice. An atom has not 4 bonds.");
  //All the vertices are supposed to have same specie, we pick the first neighbor
  tet.vertex_sp = atoms[bondmap[atm][0]].get_specie();
  //We check that all neghbors are really of the same specie
  for (unsigned int i = 0; i < n_bonds; i++)
  {
    if (atoms[bondmap[atm][i]].get_specie() != tet.vertex_sp)
      Messages::error("Error in StrainLattice. Try to build tetraedron on ill conditioned structure");
  } 

  for (unsigned int i = 0; i < n_bonds; i++)
  {
   //calculate bonds and put in tet.bonds R01 = R1 - R0
   bond = atoms[bondmap[atm][i]].get_ttype_position() + translation[atm][i] 
    - atoms[atm].get_ttype_position();
   tet.bonds[i] = bond; 
  }  
   //do similar job with edges R12, R23, R34 and store them by column
   //Note that R12 = R02 - R01, we can just use the bonds
   for (unsigned int i =1; i < n_bonds; i++)
   {
   edge = tet.bonds[i] - tet.bonds[i-1];
   tet.edges(1, i) = edge(1); 
   tet.edges(2, i) = edge(2);
   tet.edges(3, i) = edge(3);
   }

   return tet;
}

StrainLattice::Tetra::Tetra(void)
:bonds(4, Tensor1(0))
{
}

StrainLattice::Tetra
StrainLattice::rearrange(const Tetra& ref, const Tetra& tet)
{
  Tetra out = tet;
  double max, prod;
  std::vector<unsigned int> swap(4, 0);

  //Define a swap array by maximizing dot products
  for (unsigned int i = 0; i < 4; i++)
  {
    max = 0.0;
    for (unsigned int j = 0; j < 4; j++)
    {
      if (tet.bonds[j] * ref.bonds[i] > max)
        swap[j] = i;
    }
  }  
  for (unsigned int i = 0; i < 4; i++)
  {
    out.bonds[swap[i]] = tet.bonds[i];
  }
 //Tetraedron is built from bonds, which are now consistent
 for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = 0; j < 3; j++)
    {
      out.edges(i+1, j+1) = out.bonds[j + 1](i+1) - out.bonds[j](i + 1);
    }
  }
  return out;
}
