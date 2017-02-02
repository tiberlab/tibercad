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
    double _eq_conv_thr;

    //! The k_points as strings
    std::vector<std::string> _qe_k_points;

    //! Pseudopotentials to be used
    std::map<std::string, std::string> _qe_pseudos;

};


