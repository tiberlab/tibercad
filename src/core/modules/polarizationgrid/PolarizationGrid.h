#ifndef _POLAGRID_H_
#define _POLAGRID_H_

#include "SimulationInterface.h"
#include "TensorGrid.h"
#include <tr1/random>
#include "mesh.h"
#include "elem.h"

using namespace std;

class TBDLLOCAL PolarizationGrid : public SimulationInterface
{
 public:

  enum Solution
  {
    Polarization //
  };

  static PolarizationGrid* create(const ModelOptions& options); 
   
  virtual ~PolarizationGrid(void){};

  //! set number of monte Carlo steps
  void set_steps(unsigned int);

  //! set dipole moment in units defined by str
  void set_dipole(double p, string str); 

  void set_kT(double);
 
  void set_cutoff(void);
  
  //! single dipole energy contribution
  double energy1(unsigned int);

  //! dipole-dipole contribution
  double energy2(unsigned int);

  void setup_grid(const Point& p0, const Point& p1, int, int, int);

  void set_efield(const vector<double> &Efield);  

  void rnd_orientation(Point &p);

  void set_rnd_dipoles(void);

  void write_dipoles(void);

  void read_dipoles(void);

  void kmc(void); 
 
  void write_antiferro(void);

 protected:
  
  PolarizationGrid(const ModelOptions& opt);

  void do_init(void); 

  void do_solve(void);

  void parse_options(void);

  void get_solution_secure(const Elem* elem,
                           map<ID, std::vector<double> >& values,
                           const vector<Point>& p);
   
  void setup_mesh(void);

 private:  
  
  unsigned int Nstep;
  unsigned int IOstep;

  const double angle_res = 2;

  double kbT;

  double eps_r;
  
  double cutoff;

  double dipole;
 
  int seed;

  TensorGrid grid;
  
  vector<Point> pp;

  vector<Point> efield;

  mt19937 generator;

  unsigned int nx, ny, nz;
  int dk, dl, dm;

};

inline
void
PolarizationGrid::set_steps(unsigned int N)
{
  Nstep = N;
}


inline
void
PolarizationGrid::set_kT(double kT)
{
  kbT = kT;
}


#endif
