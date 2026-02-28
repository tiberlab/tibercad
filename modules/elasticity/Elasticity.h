/*  
 * This file is part of the tiberCAD module elasticity.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file Elasticity.h
 * \brief tiberCAD elasticity module header.
 *
 * \note This file is part of module elasticity.
 */


#ifndef TC_ELASTICITY_H
#define TC_ELASTICITY_H

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/solver/TiberLinearSystem.h"
#include "tibercad/base/tiber_dll.h"
#include "tibercad/base/Device.h"

#include "libmesh/tensor_value.h"

class Tensor4DSym;

/*!
 * 
 * \brief This implements the linear elasticity model.
 *
 *
 */
class TC_DLLOCAL Elasticity : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Elasticity(void);

    //! The assembly function
    void assemble(void);


  protected:

    //! The constructor
    Elasticity(const ModelOptions &options);

    //! The initialization
    virtual void do_init(void);

    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);


    //! Solve the MyPoisson equation
    virtual void do_solve(void);


    //! Print some useful information
    virtual void do_print_info(void);


    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
        const Material* mat) const;

    //! We need to create boundary condition model
    virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
        const MaterialBoundary* boundary) const;


    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);



  private:

    // A local helper class to be used to access assembly routine
    class MyAssembly : public TiberLinearSystem::Assembly
    {
      public:
      MyAssembly(Elasticity* obj) : _obj(obj) {};

      void assemble() override
      {
        _obj->assemble();
      }

      private:
      Elasticity* _obj;

    };

    MyAssembly _my_assembly;


    //! We need this to store the accumulated elemental strain
    class SymTensor
    {
      public:
      SymTensor(void)
      {
        _val[0] = _val[1] = _val[2] = _val[3] = _val[4] = _val[5] = 0.0;
      }

      SymTensor(const libMesh::RealTensor& t)
      {
        _val[0] = t(0,0);
        _val[1] = t(1,1);
        _val[2] = t(2,2);
        _val[3] = t(0,1);
        _val[4] = t(1,2);
        _val[5] = t(0,2);
      }

      void get_tensor(libMesh::RealTensor& t)
      {
        t(0,0) = _val[0];
        t(1,1) = _val[1];
        t(2,2) = _val[2];
        t(0,1) = t(1,0) = _val[3];
        t(1,2) = t(2,1) = _val[4];
        t(0,2) = t(2,0) = _val[5];
      }

      SymTensor& operator+=(libMesh::RealTensor& t)
      {
        _val[0] += t(0,0);
        _val[1] += t(1,1);
        _val[2] += t(2,2);
        _val[3] += t(0,1);
        _val[4] += t(1,2);
        _val[5] += t(0,2);
        return(*this);
      }

      private:
      double _val[6];
    };

    //! The accumulated elemental strain, used for shape deformation
    HashMap<const Elem *, SymTensor>::Type _accumulated_strain;

    //! A pointer to device
    Device *_device;

    //! Compute the elastic energy
    Real compute_elastic_energy(void);

    typedef std::map<const Elem *, libMesh::RealGradient> force_vector;

   
    struct options
    {
      ID magnification; //!< Magnification
      ID edge;
      bool non_linear_strain;
      bool deformation;
      double shape_error;                   //!< Max tollerance for shape deformation
      unsigned int shape_iterations;        //! Max number of shape itarations
      std::string structure_to_be_strained; //! Atomistic strucure to be strained
    };

    options myopt;

    //! These are the known solution variables
    /*!
     * This is an enum, but we use the string representation of
     * the enum values to refer to solutions for plotting or
     * for data exchange with other modules.
     *
     * \note Do \em not use (\c INVALID_ID - 1) or the strings \c RegionIDs
     * or \c materials as they are used to plot the materials/region IDs.
     *
     * \note The name "all" is used to plot all solutions
     */
    enum Solutions
    {
      Strain,            /*!< the strain */
      StrainCell,        /*!< the strain */
      StrainCrystal,     /*!< the strain in the crystal system*/
      Stress,            /*!< the total stress */
      RelativeStrain,    /*!< the strain relative to the structure (= reference material) */
      StressCrystal,     /*!< the stress in the crystal system*/
      Displacement,      /*!< the displacement */
      ForceSource,       /*!< force source */
      StrainSource,      /*!< strain source */
      StressSource,      /*!< stress source */
      HydrostaticStrain, /*!< the hydrostatic strain */
      EnergyDensity,     /*!< elastic energy density */
      Energy             /*!< Force */
    };

    //! The total displacement
    std::unique_ptr<libMesh::NumericVector<Number>> sol;

    force_vector internal_force;

    std::vector<unsigned int> uvar;

    //! Apply the deformation
    void apply_shape_deformation();

    //! Apply correction for internal strain to atoms
    void internal_strain_correction(AtomisticStructure *as);

    //! Update the aqccumulated strain
    void accumulate_strain(void);

    //! Restore the shape
    void restore_shape();

    libMesh::RealTensor get_subtensor(const Tensor4DSym &C_calc, unsigned int i, unsigned int k);
};




#endif // TC_ELASTICITY_H
