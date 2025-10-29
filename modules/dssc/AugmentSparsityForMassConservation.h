// $Id$

#ifndef AUGMENTSPARSITYFORMASSCONSERVATION_H_
#define AUGMENTSPARSITYFORMASSCONSERVATION_H_


#include "libmesh/dof_map.h"

namespace libMesh {
  class ImplicitSystem;
}

using libMesh::DofMap;
using libMesh::ImplicitSystem;
using libMesh::dof_id_type;

/*!
 * \brief Augment the sparsity pattern for mass conservation
 */
class AugmentSparsityForMassConservation : public DofMap::AugmentSparsityPattern
{

  public:

    /*!
     * constructor
     */
    AugmentSparsityForMassConservation(ImplicitSystem& sys,
        const std::map<dof_id_type, std::set<unsigned int>>& coupled_vars);


    /*!
     * Our function to do the actual augment
     */
    virtual void augment_sparsity_pattern(libMesh::SparsityPattern::Graph & ,
                                          std::vector<dof_id_type> & n_nz,
                                          std::vector<dof_id_type> & n_oz);


  private:

    /*!
     * The equation systems object
     */
    ImplicitSystem& _sys;


    /*!
     * The other variables the DOF should couple to
     */
    std::map<unsigned int, std::set<unsigned int>> _coupled_vars;
};

#endif /* AUGMENTSPARSITYFORMASSCONSERVATION_H_ */
