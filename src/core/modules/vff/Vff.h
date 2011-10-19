
#ifndef _VFF_H_
#define _VFF_H_

#include "SimulationInterface.h"

class TBDLLOCAL Vff : public SimulationInterface

{

public:

  //! Destructor
  /*!
   * We do not declare it virtual here, as we will not allow
   * to derive from this class anyway.
   */
  ~Vff(void);

  //! We need a public static creator function
  static Vff* create(const ModelOptions& options);

protected:

  //! The initialization
  virtual void do_init(void){};


  //! Parse the options from the input file
  virtual void parse_options(void){};


  //! Setup the available variables
  virtual void do_setup_solution_variables(void){};


  //! Solve the MyPoisson equation
  virtual void do_solve(void){};


  //! Print some useful information
  virtual void do_print_info(void){};


  //! We need to create a physical model
  virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
      const Material* mat) const{};

  //! We need to create boundary condition model
  virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
      const Material* material_A, const Material* material_B) const{};


  //! We have to provide somehow our solution variables
  virtual void get_solution_secure(const Elem* elem,
      std::map<ID, std::vector<double> >& values,
      const std::vector<Point>& p){};




private:

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
    Displacement     /*!< the atom displacement */
  };

  //! The constructor
  /*!
   * Being private disables further inheritance.
   */
  Vff(const ModelOptions& options);

  //! A static pointer to this
  static Vff* _this;

  class Options
  {
  public:
    Options(void);

    std::string boundary_conditions;
    std::string substrate_plane;
  };

  Options _options;

};





#endif
