// $Id$

#ifndef _EIGENVALUEPROBLEM_H_
#define _EIGENVALUEPROBLEM_H_

#include "SimulationInterface.h"


//! Abstract class to solve complex valued eigenvalue problem
class EigenvalueProblem : public SimulationInterface
{

  public:

    //! Constructor
    EigenvalueProblem(void) { };

    //! Destructor
    ~EigenvalueProblem(void) { };

    struct eigen_problem_solution
    {
      //! particle type ("electron", "hole", "photon", "exciton", ...) 
      std::string particle;
      //! eigenvalue
      double eigen_energy; 
      //! eigenvector
      std::vector< std::complex<double>  > eigen_vector; 
      //! statistic type ("Bose" or "Fermi")
      std::string statistics;
      //! electro-chemical potential [eV] \f$ \langle \psi |\mu({\bf r} | \psi \rangle \f$
      double electro_chem_pot;
      //! Level temperature
      double Temperature; 
    };

    const std::vector<eigen_problem_solution>& get_eigen_solution() const
                                                           { return _solution;};

    void get_eigenvalues(std::vector<double>& values) const;

    void get_populations(std::vector<double>& values) const;

    virtual void assemble(const ModelOptions& options){};
    
    virtual std::complex<double> calculate_matrix_element(const std::string& i_particle,
							  unsigned int i, 
							  const std::string& j_particle,
							  unsigned int j){}; 
   
  protected:

    std::vector<eigen_problem_solution> _solution;
   

  private:

    double Fermi(double Energy, double Fermi_energy, double Temperature) const;

    double Bose(double Energy, double elec_chem, double Temperature) const;
};


#endif // _EIGENVALUEPROBLEM_H_
