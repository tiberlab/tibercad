/*  
 * This file is part of the tiberCAD module quantumespresso.
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
 * \file QEInterface.h
 * \brief tiberCAD quantumespresso module header.
 *
 * \note This file is part of module quantumespresso.
 */

#ifndef TC_QEINTERFACE_H
#define TC_QEINTERFACE_H

#include "tibercad/module/SimulationInterface.h"

/*!
 * \brief Module to interface to QuantumEspresso
 *
 * This module contains functionality to interface to QuantumEspresso
 * calculations. For the moment, it can produce input files for QE
 * based on control parameters given from the input file, and using 
 * the atomistic structure provided.
 */
class QEInterface : public SimulationInterface
{

  public:

    //! Destructor
    ~QEInterface(void);

    //! Create an instance of QEInterface
    static QEInterface* create(const ModelOptions& options);


  protected:

    /*! 
     * \brief Constructor
     * 
     * There is no sense in deriving from this module, presumably.
     */
    QEInterface(const ModelOptions& options);



    /*!
     * \brief Initialize the module
     *
     * There is no sense in deriving from this module, presumably.
     */
    virtual void do_init(void) final;


    /*!
     * \brief Call the solver
     *
     * For now, this just prepares the files
     */
    virtual void do_solve(void) final;


  private:

    //! The directory with the pseudopotentials
    std::string _qe_pseudo_dir;

    //! The wavefunction cutoff energy
    double _qe_ecutwfc;

    //! Convergence threshold for QE
    double _qe_conv_thr;

    //! Mixing parameters for electrons
    double _qe_mixing_beta;

    //! Number of bands for nscf
    int _qe_nbnd;

    //! The k_points as strings, for the scf file
    std::vector<std::string> _qe_k_points_scf;

    //! The k_points as strings, for the nscf file
    std::vector<std::string> _qe_k_points_nscf;

    //! Pseudopotentials to be used
    std::map<std::string, std::string> _qe_pseudos;

    //! Output directory for QE
    std::string _qe_outdir;

    //! degauss value for pdos
    double _qe_degauss; //

    //! DeltaE value for pdos
    double _qe_DeltaE; //

    //! Emin value for pdos
    double _qe_Emin; //

    //! Emax value for pdos
    double _qe_Emax; //

    //! filplot option for pp.el1
    std::string _qe_filplot; //

    //! fileout for pp.el1
    std::string _qe_fileout; //

    //! plot_num for pp.el1
    int _qe_plot_num; //

    //! kpoint for pp.el1
    int _qe_kpoint; //

    //! nfile option for pp.el1
    int _qe_nfile; //

    //! iflag for pp.el1
    int _qe_iflag; //

    //! nx for pp.el1
    int _qe_nx; //

    //! ny for pp.el1
    int _qe_ny; //

    //! nz for pp.el1
    int _qe_nz; //

    //! output format for pp.el1
    int _qe_pp_outputformat; //

    //! weight for pp.el1
    double _qe_weight; //

};


#endif // TC_QEINTERFACE_H
