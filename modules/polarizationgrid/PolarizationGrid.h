/*  
 * This file is part of the tiberCAD module polarizationgrid.
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
 * \file PolarizationGrid.h
 * \brief tiberCAD polarizationgrid module header.
 *
 * \note This file is part of module polarizationgrid.
 */

#ifndef TC_POLAGRID_H
#define TC_POLAGRID_H

#include "tibercad/module/SimulationInterface.h"
#include "tibercad/geom/TensorGrid.h"
#include <random>
#include "mesh.h"
#include "elem.h"


class TC_DLLOCAL PolarizationGrid : public SimulationInterface
{
 public:

  enum Solution
  {
    Polarization //
  };

  virtual ~PolarizationGrid(void){};

  //! set number of monte Carlo steps
  void set_steps(unsigned int);

  //! set dipole moment in units defined by str
  void set_dipole(double p, std::string str); 

  void set_kT(double);
 
  void set_cutoff(void);
  
  //! single dipole energy contribution
  double energy1(unsigned int);

  //! dipole-dipole contribution
  double energy2(unsigned int);

  void setup_grid(const Point& p0, const Point& p1, int, int, int);

  void set_efield(const std::vector<double> &Efield);  

  void rnd_orientation(Point &p);
  
  void rnd_orientation_discrete(Point &p);

  void set_rnd_dipoles(void);

  void write_dipoles(void);

  void read_dipoles(void);

  void kmc(void); 
 
  void write_dipoles_man(void);

  void write_antiferro(void);

  void average_dipole(void);


 protected:
  
  explicit PolarizationGrid(const ModelOptions& opt);

  virtual void do_setup_solution_variables(void) override;

  virtual void do_init(void) override; 
  
  virtual void do_reinit(void) override; 

  virtual void do_solve(void) override;

  void parse_options(void);

  virtual void get_solution_secure(const Elem* elem,
                           std::map<ID, std::vector<double> >& values,
                           const std::vector<Point>& p) override;
  
  //NumericVector<double>& do_get_solution_vector(void);
    
  virtual void setup_mesh(void) override;
 
  void dot_product(unsigned int i, std::vector<double>& bins, double dr);

  void autocorrelation(int dir, unsigned int l1, unsigned int l2, std::vector<double>& c);

  virtual void plot_globaldata(void) override;

 private:  
  
  unsigned int Nstep;
  unsigned int IOstep;

  const double angle_res = 2;

  double kbT;

  double eps_r;
  
  double cutoff;

  double dipole;

  double scaling;
  
  double scaling_bkp;

  int seed;

  TensorGrid grid;
  
  std::vector<Point> pp;

  std::vector<Point> efield;

  std::mt19937 generator;

  unsigned int _nx[3];

  int _dx[3];

  unsigned int _ncell[3];
  
  bool _write_file;

  int _periodic[3];

  bool _discrete;

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
