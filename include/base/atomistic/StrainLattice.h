#ifndef _STRAINLATTICE_H_
#define _STRAINLATTICE_H_

#include "TiberModelObject.h"
#include "tensor.h"
#include "Specie.h"
#include "point.h"
#include "Atom.h"

//Forward declarations
class AtomisticStructure;
class Material;
class AtomisticBasis;

//! This class implement a strain projector based on 
//! deformation of single tetraedra of fcc and wz structures
//! The algorithm is taken from J. Appl.Phys. 83, 2548 (1998)
class StrainLattice
{

  public:

  StrainLattice(void);

  void init(AtomisticStructure* as);

  //!Calculate the solution vector
  void do_solve(void);

  //! A solution point is identified by a coordinate and a tensor
  class TensorField
  {
    public:
    //TODO: coord is ging to be deprecated, all we need is in Atom
    Point coord;
    //------------------------------------------------------------
    const Atom* atom_p;
    Tensor2Gen tensor;
  };

  //! Get solution vector
  const std::vector<TensorField>& get_solution(void) const;

  protected:

  private:    

    class Tetra
    {
      public:
      Tetra(void);
      //!Edges R12, R23, R34 stored by columns
      Tensor2Gen edges;
      //! Bonds R01, R02, R03, R04 
      //Note: bonds are redundant, geometry is defined by edges, 
      //but it's handy to have them separated if the overload 
      //is not too high
      std::vector<Tensor1> bonds;
      //!Specie of central atom
      Specie sp;
      //!Specie of vertices (useful to discriminate central atom
      //in alloys)
      Specie vertex_sp;
    };

    //! Build tetraedron for a given atom of a given structure
    Tetra build_tetraedron(const AtomisticBasis* as,
        unsigned int atm);

    //! Reference tetraedron for every material
    std::map<const Material*, Tetra> _reference;

    //! Structure to be processed
    AtomisticStructure* _as;

    //! Fill the set of all possible needed materials, 
    //to build corresponding reference bulk and 
    //tetraedra
    void fill_materials_set(void);

    //! Contains all needed material (including parent)
    std::set<const Material*> _materials;

    //! Return a tetraedron with vertices arranged as reference.
    //! Decision is made by maximizing dot product of bond vectors
    Tetra rearrange(const Tetra& ref, const Tetra& tet);

    //! Solution vector
    std::vector<TensorField> _solution;

};


////////////////////////
// Inline functions  //
//////////////////////

inline
const std::vector<StrainLattice::TensorField>& 
StrainLattice::get_solution(void) const
{
  return _solution;
}

#endif // _STRAINLATTICE_H_
