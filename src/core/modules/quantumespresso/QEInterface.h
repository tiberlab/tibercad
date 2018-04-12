// $Id$

#include "SimulationInterface.h"

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

    /*! 
     * \brief Constructor
     * 
     * There is no sense in deriving from this module, presumably.
     */
    QEInterface(const ModelOptions& options);


    //! The directory with the pseudopotentials
    std::string _qe_pseudo_dir;

    //! The wavefunction cutoff energy
    double _qe_ecutwfc;

    //! Convergence threshold for QE
    double _qe_conv_thr;

    //! Number of bands for nscf
    int _qe_nbnd;

    //! The k_points as strings
    std::vector<std::string> _qe_k_points;

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


