// $Id$

#ifndef _EIGENVALUEPROBLEM_H_
#define _EIGENVALUEPROBLEM_H_

#include "SimulationInterface.h"


//! Abstract class to solve complex valued eigenvalue problem
class EigenvalueProblem : public SimulationInterface
{

  public:

    //! Constructor
    EigenvalueProblem(const ModelOptions& options);

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
      double temperature; 
    };

    //! 
    void set_k_vector(std::vector<double> k_point);

    const std::vector<double>& get_k_vector(void) const;
  
    bool has_new_k(void) const;

    void k_is_old(void);

    //! get the eigensolution vector
    const std::vector<eigen_problem_solution>& get_eigen_solution() const
      {return _solution;};

    void get_eigenvalues(const std::string& particle, std::vector<double>& values) const;

    void get_populations(const std::string& particle, std::vector<double>& values) const;

    double get_population(int i) const;

    virtual void assemble(const ModelOptions& options){};
    
    //! computes matrix elements between state i of particle_i and state j of particle_j
    virtual std::complex<double> calculate_matrix_element(const std::string& i_particle,
							  unsigned int i, 
							  const std::string& j_particle,
							  unsigned int j){}; 
    //! get number of states of a given particle type
    unsigned int get_num_states(const std::string& particle) const;

    /*! Note: for the moment calculate_matrix_element relays on the fact that the first
     *  n_vb states are for valence, then there are all the electron states.
    */

    void write_states(void) const;

    void write_states(const std::string& filename) const;

  protected:

    std::vector<eigen_problem_solution> _solution;

    double Fermi(double Energy, double Fermi_energy, double Temperature) const;

    double Bose(double Energy, double elec_chem, double Temperature) const;

  private:

  std::vector<double> _k_vector;
  
  bool _new_k;

};


inline
EigenvalueProblem::EigenvalueProblem(const ModelOptions& options)
 : SimulationInterface(options)
{
  _k_vector.reserve(3);
  _k_vector[0]=0.0;   _k_vector[1]=0.0;   _k_vector[2]=0.0; 
}

//inline 
//const std::vector<eigen_problem_solution>& EigenvalueProblem::get_eigen_solution() const
//{ 
//   return _solution;
//}

inline
const std::vector<double>& EigenvalueProblem::get_k_vector() const
{
  return _k_vector;
}

inline 
void EigenvalueProblem::set_k_vector(std::vector<double> k_point)
{
  _k_vector = k_point;
  _new_k = true;
}

inline 
bool EigenvalueProblem::has_new_k(void) const
{
  return _new_k;
}

inline
void EigenvalueProblem::k_is_old(void)
{
  _new_k=false;
}


#endif // _EIGENVALUEPROBLEM_H_
