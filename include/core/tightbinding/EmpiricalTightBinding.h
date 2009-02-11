#ifndef _EMPIRICALTIGHTBINDING_H_
#define _EMPIRICALTIGHTBINIDNG_H_


#include "TightBinding.h"


class ETB : public TightBinding
{

public:
  //! A class for Dftb options
  class UptOptions
  {
  public:

    UptOptions(void);
    ~UptOptions(void);

    int verbose;
    int max_TB_order;
    bool harrison_flag;
    bool relat_flag;
    bool potential_flag;
    bool opt_flag;
    int poldir;
    char *database_path;
    char *work_path;
    char *upt_filename;
    char *gen_outfile;
    double *c_axis;
      
  };

  class UptSolverOptions
  {
  public:
    UptSolverOptions(void);
    ~UptSolverOptions(void);

    std::string solver;
    int n_vb;
    int n_cb;
    int min_iter;
    int long_iter;
    int max_iter;
    double guess_vb;
    double guess_cb;
    double fast_tol;
    double long_tol;
    double ort_tol;

  };

  //! Constructor
  ETB(void);

  //! Destructor
  ~ETB(void);

  //! Create object
  static ETB* create();

  virtual PhysicalModel* create_physical_model(const ModelOptions &options,
      const Material* mat) const throw (ModelErrorException);



private:

  //! Get options suited for DFTB+ tight binding builder and solver
  void get_upt_options();

  //! Function for building options from input
  void build_input_options(void);

  //! Print all _dftb_options for debugging purposes
  void print_upt_options(void);

  //! Structure containing options for DFTB+ tight binding builder
  UptOptions _upt_options;

  //! Structure containing options for DFTB+ tight binding solver
  UptSolverOptions _upt_solver_options;
  
  //! Uptight instance associated to the simulation
  UptWrapper* inst;

protected:

  virtual void do_init(void);

  virtual void do_solve (void);

  virtual void parse_options(void);

  void read_kpoints(void);

};

// Inline members definition
//----------------------------

inline
ETB* ETB::create()
{
  return new ETB;
}










#endif
