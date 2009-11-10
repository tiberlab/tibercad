#ifndef _DFTB_H_
#define _DFTB_H_

//---------------------------------------------------------------------
#include "TightBinding.h"

#include <map>


//! A class for Tight Binding simulations with DFTB+ code
/*!
 *This class provides methods for invoking DFTB+ (Density Functional
 *Tight Binding) library, in order to make calculations at equilibrium
 */
class Dftb : public TightBinding{

public:

  //! A class for Dftb options
  class DftbOptions
  {
  public:

    DftbOptions(void);
    ~DftbOptions(void);

    char* skNames;
    char* specieNames;
    std::vector <std::string> specieNameStrings ;
    int nAtom;
    int nType;
    double* coords;
    int* species;
    int iPeriodic;
    double* latVecs;
    double eTemp;
    int orbResolved;
    int skInterp;
    int nkPoints;
    double *kPoints;
    double* kWeights;
    int* mAngs;
    bool supersampling;
    double* samplingcoeffs;
    double* samplingshift;
  };


  class DftbSolverOptions{
  public:
    std::string solver;
  };


  class DftbpCouplingOptions
  {
  public:
    //! Simulation which provides external potential
    std::string potential_sim;
  };

  //! Constructor
  Dftb(void);

  //! Destructor
  ~Dftb(void);

  //! Create object
  static Dftb* create();

  virtual PhysicalModel* create_physical_model(const ModelOptions &options,
      const Material* mat) const throw (ModelErrorException);

private:

  //! A map containing which kind of shell parametrization is used for any specie (maximum angular momentum + 1)
  std::map <std::string, int> _shell;

  //! Maximum angular momentum in parametrization of all species
  int _max_shell;

  //! Get options suited for DFTB+ tight binding builder and solver
  void get_dftbp_options();

  //! A function for building a DFTB compatible char of SK and species names
  void build_names(void);

  //! Function for building options from atomistic structure
  void build_structure_options(void);

  //! Function for building options from input
  void build_input_options(void);

  //! Structure containing options for DFTB+ tight binding builder
  DftbOptions _dftb_options;

  //! Structure containing options for DFTB+ tight binding solver
  DftbSolverOptions _dftb_solver_options;

  //! Print all _dftb_options for debugging purposes
  void print_dftb_options(void);

  //! Dftb instance associated to the simulation
  DftbpWrapper* inst;

  //! Structure containing options for coupling DFTB with other calculations
  DftbpCouplingOptions _dftbp_coupling_options;

  //! Set external potential provided by other simulations using TightBinding method
  void add_shifts();

protected:

  virtual void do_init(void);

  virtual void do_solve (void);

  virtual void  parse_options(void);

  //! Gives Hubbard parameters and put them in TightBinding member
  /*!
   * do_init() and parse_parameter() must have been done, note that if updatecoords runs
   * this function must be refreshed
   */
  void obtain_hubbard_parameters(void);

  void read_kpoints(void);

};



inline
Dftb* Dftb::create()
{
  return new Dftb;
}





#endif
