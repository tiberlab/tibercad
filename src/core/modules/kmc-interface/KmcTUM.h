#ifndef _KMCTUM_H_
#define _KMCTUM_H_

#include "KmcInterface.h"
#include "TensorGrid.h"
#include "SimulationInterface.h"

//class Mesh;

using namespace std;

class TBDLLOCAL KmcTUM : public KmcInterface
{
 public:

  static KmcTUM* create(const ModelOptions& options); 
   
  // Class Destructor
  virtual ~KmcTUM() {};
 
 protected:
 
  // All reimplemented from SimulationInterface: 
  void do_setup_solution_variables(void);

  void do_init(void); 
  
  void do_reinit(void); 

  void do_solve(void);

  void parse_options(void);

  void get_solution_secure(const Elem* elem,
                           map<ID, std::vector<double> >& values,
                           const vector<Point>& p);
  
  void setup_mesh(void);

  void plot_globaldata(void);

 private:  

  //! private constructor. can be created only from 'create'
  KmcTUM(const ModelOptions& opt);
  
  //! write input file
  void write_input(void);

  void write_boundary_potential(void);
  
  void write_potential(void);

  bool read_density(void);
  
  bool read_current(void);

  double linear_interpolation(const Point& p);

  TensorGrid grid;

  double _ecurr_bottom;
  double _ecurr_top;
  double _hcurr_bottom;
  double _hcurr_top;

  Point _p0;

  string _exe_name;
  string _kmc_scratch;
 
  double kbT;
  double eps_r;

  bool _create_blend; 
  unsigned int _generation_steps;
  double _HOMO_acceptor;
  double _LUMO_acceptor;
  double _HOMO_donor;
  double _LUMO_donor;
  double _mue_cathode; 
  double _muh_cathode; 
  double _mue_anode; 
  double _muh_anode; 
  double _sigma;
  double _x_v0;
  double _x_sep;
  double _x_rec;
  double _el_v0;
  double _hl_v0;
  double _el_inj;
  double _hl_inj;
  double _alpha;
  double _el_hl_R;
  double _el_coll;
  double _hl_coll;
  double _bias;
  double _volume_ratio;     
  unsigned int _coll_el_cathode;
  unsigned int _coll_el_anode;
  unsigned int _coll_hl_cathode;
  unsigned int _coll_hl_anode;
  double _kmc_total_time;
  double _kmc_stat_time;
  double _const_gen_rate;
  double _buffer_top;
  double _buffer_bottom;

   
  map<Boundary*, double>  _mue;
  map<Boundary*, double>  _muh;

};

#endif
