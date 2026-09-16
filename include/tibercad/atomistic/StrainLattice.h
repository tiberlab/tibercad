/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file StrainLattice.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef TC_STRAINLATTICE_H
#define TC_STRAINLATTICE_H

#include "tibercad/module/TiberModelObject.h"
#include "tibercad/atomistic/Specie.h"
#include "tibercad/math/Tensor1.h"
#include "tibercad/math/Tensor2.h"
#include "tibercad/atomistic/Atom.h"

#include "libmesh/point.h"

//Forward declarations
class AtomisticStructure;
class Material;
class AtomisticBasis;

/*!
 * \brief This class implement a strain projector based on
 * deformation of single tetraedra of fcc and wz structures
 *
 * The algorithm is taken from J. Appl.Phys. 83, 2548 (1998)
 */
class StrainLattice
{

  public:

  StrainLattice(void);

  void init(AtomisticStructure* as);

  //! Calculate the solution vector
  void do_solve(void);

  //! A solution point is identified by a coordinate and a tensor
  class TensorField
  {
    public:
    //TODO: coord is ging to be deprecated, all we need is in Atom
    //Point coord;
    //------------------------------------------------------------
    const Atom* atom_p;
    Tensor2 tensor;
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
      Tensor2 edges;

      //! Bonds R01, R02, R03, R04 
      /*! Note: bonds are redundant, geometry is defined by edges,
       * but it's handy to have them separated if the overhead
       * is not too high.
       */
      std::vector<Tensor1> bonds;

      //! Specie of central atom
      //Specie sp;

      //! label of central atom
      unsigned char central_atom_label;

      //! Specie of vertices
      /*! useful to discriminate central atom in alloys */
      Specie vertex_label;
    };

    //! Build tetraedron for a given atom of a given structure
    //Tetra build_tetraedron(const AtomisticBasis* as,
    //    unsigned int atm);
    void build_tetraedron(const AtomisticBasis* as,
        unsigned int atm, Tetra& tet);

    //! Reference tetraedron for every material
    std::map<const Material*, Tetra> _reference;

    //! Additional reference tetraedron for wurtztite
    std::map<const Material*, Tetra> _reference2;
    
    //! Structure to be processed
    AtomisticStructure* _as;

    //! Fill the set of all possible needed materials, 
    //to build corresponding reference bulk and 
    //tetraedra
    void fill_materials_set(void);

    //! Contains all needed material (including parent)
    std::set<const Material*> _materials;

    //! Return a tetraedron with vertices arranged as reference.
    /*!
     * Decision is made by maximizing dot product of bond vectors
     */
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

#endif // TC_STRAINLATTICE_H
