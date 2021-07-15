#include "PolarizationGrid.h"
#include "InitFailedException.h"
#include "MeshUtils.h"
#include "Constants.h"
#include "mesh_generation.h" 
#include "SimulationEnvironment.h"
#include "Messages.h"
#include <sstream> 
#include <fstream> 
#include "TiberModule.h"
#include "gmsh_io.h"

using std::string;
using std::vector;
using std::map;
using std::endl;
using std::cout;
using std::uniform_int_distribution;
using std::uniform_real_distribution;

PolarizationGrid::PolarizationGrid(const ModelOptions& options) :
 SimulationInterface(options),
  eps_r(1.0),
  cutoff(100.0),
  dipole(0.0),
  Nstep(100000),
  kbT(0.026),
  _write_file(false)
{
   has_solution_vector(false);
}

PolarizationGrid*
PolarizationGrid::create(const ModelOptions& opt)
{
  return new PolarizationGrid(opt);
}


void
PolarizationGrid::do_setup_solution_variables(void)
{
  declare_solution(Polarization, VECTOR, CELL, "C/m^2");
}



void
PolarizationGrid::do_init()
{ 
  parse_options();  
 
  generator.seed(seed);

  if (get_options().find_option("read_file"))
  {
    read_dipoles();
  }
  else
  {
    set_rnd_dipoles();
  }

  set_cutoff();

  if (get_options().find_option("write_file"))
    _write_file = true;
}
 
void
PolarizationGrid::setup_mesh(void)
{
  // the grid must be created first (before any other option)
  _nx[0] = get_option("Nx",10); 
  _nx[1] = get_option("Ny",10); 
  _nx[2] = get_option("Nz",1); 

  double gs = get_option("grid_step",0.65); // in nanometers
  double aa = get_option("grid_step_a",gs); // in nanometers
  double bb = get_option("grid_step_b",aa); // in nanometers
  double cc = get_option("grid_step_c",gs); // in nanometers

  Point p1(0,0,0);
  Point p2(_nx[0]*aa,_nx[1]*bb,_nx[2]*cc);
  setup_grid(p1,p2,_nx[0],_nx[1],_nx[2]);
  //create internal rectangular mesh 
  //
  //Mesh* mymesh = new Mesh(3);
  //ElemType type(HEX8);

  //MeshTools::Generation::build_cube(*mymesh, 
	//			       _nx[0], _nx[1], _nx[2],
  //             0.0, _nx[0] * aa, 
  //             0.0, _nx[1] * bb,
  //             0.0, _nx[2] * cc,
  //             type);

  // set in environment
  //get_environment().set_mesh(mymesh);

  // set in my object
  //set_mesh(mymesh);
 
  SimulationInterface::setup_mesh();

  //GmshIO(*mymesh).write("mymesh.msh");

}

void
PolarizationGrid::do_reinit()
{
}

void
PolarizationGrid::do_solve()
{

  //for (unsigned int i = 0; i < grid.num_elements(); ++i)
  //{
  //   pp[i] = pp[i] * scaling;
  //}


  // If the electric field is taken from an external module we
  // have to reload the fields
  string poisson = get_option("poisson_model", "");
  if (!poisson.empty())
  {
    SolutionProvider poisson_sol = this->find_solution_provider(poisson);
    // check

    SimulationInterface* poisson_mod = poisson_sol.first;
    ID sol_id = poisson_sol.second;

    MeshUtils::GridMapper& mapper =
        MeshUtils::GridMapper::get_mapper(poisson_mod->get_mesh());

    vector<double> values(3, 0.0);

    for (unsigned int i = 0; i < grid.num_elements(); ++i)
    {
      Point p(grid.get_centroid(i));
      const Elem* elem = mapper.get_element(p);
      if ((elem != NULL) && poisson_mod->get_solution(elem, sol_id, values, p))
      {
        efield[i](0) = values[0] * 1e-7; // V/nm
        efield[i](1) = values[1] * 1e-7;
        efield[i](2) = values[2] * 1e-7;
      }
    }
  }

  this->kmc();

  if (_write_file)
  {
    write_dipoles();
    //write_dipoles_man();
    write_antiferro();
  }


  //for (unsigned int i = 0; i < grid.num_elements(); ++i)
  //{
  //   pp[i] = pp[i] / scaling;
  //}

}

void
PolarizationGrid::parse_options(void)
{
  Nstep = get_option("Nstep",10);
  IOstep = get_option("IOstep",1000);
  eps_r = get_option("eps_r",1.0);
  scaling = 1.0;
  scaling_bkp = 1.0;
  get_parameter("scaling_factor", scaling);
  string units = get_option("dipole_units","Debye");
  double dip = get_option("dipole", 2.29);
  set_dipole(dip,units);

  std::ostringstream ostr;
  ostr<<"dipole= "<<dipole<<" e*nm "<<endl;
  Messages::info(ostr.str());  

  vector<double> E(3,0.0);
  get_option("external_field",E);
  set_efield(E);

  kbT = get_option("kbT", 0.0);
  
  seed = get_option("seed",
      static_cast<int>(time(NULL) * std::random_device()()));

  // faster (I think) to use array than std::vector;
  vector<bool> periodic;
  get_option("periodic",periodic);

  for (int i=0; i<3; i++)
  {
    _periodic[i] = 0;
    if (periodic[i]) _periodic[i] = 1;
  }

  _discrete=get_option("discrete_orientations",false);

}



void
PolarizationGrid::set_cutoff(void)
{
  // e^2/(4 pi eps0) = 1.44 eV nm

  double tmp = 4 * Constants::EE * dipole * dipole / eps_r;
 
  // cutoff in nm to get 0.1 meV energy difference
  cutoff = pow(tmp/0.0001,1.0/3.0);
  
  std::ostringstream os;
  os<<"Computed Cutoff= "<<cutoff<<endl;
  Messages::info(os.str());

  cutoff = get_option("cutoff",cutoff);
  
  for(unsigned int i=0; i<3; i++) 
  {
    _dx[i] =static_cast<int>( ceil(cutoff/grid.grid_step(i)) );

    if(_periodic[i]) 
    {
      _ncell[i] = cutoff / _nx[i] + 1;
    }
    else
    {
      _ncell[i]=0;
    }
  }

  os.str("");
  os<<"Cutoff= "<<cutoff<<endl;
  Messages::info(os.str());
}

void
PolarizationGrid::set_dipole(double p, string str)
{
  // The dipole must be transformed in e*nm
  if (str == "Debye")
  {
    // 1 Debye = 3.336 10^-30 C m = 0.0208 e * nm
    dipole = p * 0.0208;
  }
  else if (str == "C*m")
  {
    // 1 C m =  e/(1.60e-19) *  nm/1e-9 
    dipole = p / 1.60e-28;
  }
  else if (str == "C/m^2")
  {
    // P*vol = C/m^2 * nm*nm^2 = C/m^2 * nm*1e-18 m^2 = 1e-18/1.60e-19 e*nm;
    double vol = grid.grid_step(0)*grid.grid_step(1)*grid.grid_step(2);
    //cout<<"vol= "<<vol<<"nm^3"<<endl;
    dipole = p * vol / 0.1602;
    //cout<<"dipole= "<<dipole<<" e*nm"<<endl;
  }
  else if (str == "e*nm")
  {
    dipole = p;
  }
  else
  {
    throw InitFailedException("set_dipole: unrecongnized units");
  }
}

void
PolarizationGrid::setup_grid(const Point& p0, const Point& p1, int nx, int ny, int nz)
{ 
  grid.setup(p0,p1,nx,ny,nz);
  pp.resize(grid.num_elements());
}

void 
PolarizationGrid::set_efield(const vector<double>& E)
{
  // field is in V/cm
  efield.resize(grid.num_elements());

  for(int i=0; i < efield.size(); i++)
  {
    efield[i](0) = E[0] * 1e-7; // V/nm
    efield[i](1) = E[1] * 1e-7; 
    efield[i](2) = E[2] * 1e-7; 
  }
}

void
PolarizationGrid::rnd_orientation(Point &p)
{

   uniform_real_distribution<float> random1(0.0, 1.0);

   double teta;
   if (grid.num_elements(2) > 1)
   {
     teta = acos(2.0 * random1(generator) - 1.0); 
   }
   else
   {
     teta = M_PI/2.0; 
   }
  
   double phi = random1(generator) * 2.0 * M_PI; 
   
   p(0) = sin(teta) * cos(phi);
   p(1) = sin(teta) * sin(phi);    
   p(2) = cos(teta);

   p *= dipole;

}


void
PolarizationGrid::rnd_orientation_discrete(Point &p)
{
  uniform_int_distribution<int> random1(0, 5);

  int rnd = random1(generator); 
    
  p(0) = 0.0; 
  p(1) = 0.0; 
  p(2) = 0.0; 

  switch(rnd)
  {
  case 0:
    p(0) = dipole; 
    break;
  case 1:
    p(0) = -dipole; 
    break;
  case 2:
    p(1) = dipole; 
    break;
  case 3:
    p(1) = -dipole; 
    break;
  case 4:
    p(2) = dipole; 
    break;
  case 5:
    p(2) = -dipole; 
    break;
  }
}


void
PolarizationGrid::write_dipoles()
{
   string out_path = get_output_directory();
   string file_name = out_path+"/"+get_option("filename","polarization.dat");
   std::ofstream file; 

   file.open(file_name.c_str());
   for (unsigned int i=0; i<pp.size(); i++)
   {
     file<<pp[i](0)<<" "<<pp[i](1)<<" "<<pp[i](2)<<" ";
   }
   file.close();
}

void
PolarizationGrid::read_dipoles()
{
   string file_name = get_option("read_file","polarization.dat");
  
   std::ifstream file; 

   file.open(file_name.c_str());
   for (unsigned int i=0; i<pp.size(); i++)
   {
     file>>pp[i](0);
     file>>pp[i](1);
     file>>pp[i](2);
   }
   file.close();
}

void
PolarizationGrid::write_antiferro()
{
   string out_path = get_output_directory();
   string file_name = out_path+"/"+get_option("filename","antiferro.dat");
   std::ofstream file; 

   unsigned int nx = grid.num_elements(0);
   unsigned int ny = grid.num_elements(1);

   file.open(file_name.c_str());

   for (unsigned int j=0; j<ny; j++)
     for (unsigned int i=0; i<nx; i++)
     {
       file<<dipole*pow(-1.0,j)<<" "<<0.0<<" "<<0.0<<" ";
     }

   file.close();
}
  
void
PolarizationGrid::write_dipoles_man()
{
   int nx = 150;
   int ny = 50;
   int nzs = get_option("nz_start",20);
   int nze = get_option("nz_end",24);
   int nz = nze - nzs + 1;
   double aa = get_option("grid_step",0.65); // in nanometers

   Point p1(0,0,0);
   Point p2(nx*aa,ny*aa,nz*aa);
   TensorGrid grid2;
   grid2.setup(p1,p2,nx,ny,nz);
   vector<Point> pp2;
   pp2.resize(grid2.num_elements());

   // copy pp on pp2
   for (int k1 = 0; k1 <= 49; k1++)
   {
     for (int l1 = 0; l1 <= 49; l1++)
     { 
       for (int m1 = nzs; m1 <= nze; m1++) 
       {
         int j = grid.index_to_element(k1,l1,m1);
      
         int j1 = grid2.index_to_element(k1,l1,m1-nzs);
         pp2[j1] = pp[j];

         j1 = grid2.index_to_element(k1 + 50,l1,m1-nzs);
         pp2[j1] = pp[j];
         

         j1 = grid2.index_to_element(k1 + 100,l1,m1-nzs);
         pp2[j1] = pp[j];
        
       }
     }
   }


   string out_path = get_output_directory();
   string file_name = out_path+"/"+get_option("filename","polarization_man.dat");
   std::ofstream file; 

   file.open(file_name.c_str());
   for (unsigned int i=0; i<pp2.size(); i++)
   {
     file<<pp2[i](0)<<" "<<pp2[i](1)<<" "<<pp2[i](2)<<" ";
   }
   file.close();
}

void
PolarizationGrid::set_rnd_dipoles(void)
{

  if (_discrete)
  { 
    for (unsigned int i=0; i < pp.size(); i++) 
    { 
      this->rnd_orientation_discrete(pp[i]);
    }
  }
  else
  {
    for (unsigned int i=0; i < pp.size(); i++) 
    { 
      this->rnd_orientation(pp[i]);
    }
  }  
}


double
PolarizationGrid::energy1(unsigned int i)
{
  // e*nm  * V/nm  = eV
  return -pp[i]*efield[i];
}

double 
PolarizationGrid::energy2(unsigned int i)
{

  double en = 0.0;
  int ii[3];
  int start[3], end[3];

  grid.element_to_index(i,ii[0],ii[1],ii[2]);

  // way of treating periodicity: 
  // define start-end grid points and fold to central cell.
  for (unsigned int l = 0; l < 3; l++)
  {
    if (_periodic[l])
    {
      start[l] = ii[l] - _dx[l];
      end[l] = ii[l] + _dx[l];
    }
    else
    {
      // start = min(0, ii[l] - _dx[l]); 
      start[l] = ii[l] - _dx[l] > 0 ? ii[l] - _dx[l] : 0;
      // start = max(ii[l] + _dx[l], _nx[l]-1); 
      end[l] = ii[l] + _dx[l] < _nx[l]-1 ? ii[l] + _dx[l] :  _nx[l]-1;
    }
  }

  for (int k1 = start[0]; k1<=end[0]; k1++)
  {
    for (int l1 = start[1]; l1<=end[1]; l1++)
    { 
      for (int m1 = start[2]; m1<=end[2]; m1++) 
      {
        if (k1==ii[0] && l1==ii[1] && m1==ii[2]) continue;

        // note: distance also works for folded cells
        Point rr = grid.distance(ii[0],ii[1],ii[2], k1,l1,m1);
        double dd = rr.size();
        rr /= dd; 

        // fold to central cell (nx folds to 0 and so on)
        unsigned int kk1 = (k1 + _ncell[0] * _nx[0])%_nx[0];
        unsigned int ll1 = (l1 + _ncell[1] * _nx[1])%_nx[1];
        unsigned int mm1 = (m1 + _ncell[2] * _nx[2])%_nx[2];
        int j = grid.index_to_element(kk1,ll1,mm1);
        
        en += (pp[i]*pp[j] - 3.0*(pp[i]*rr) * (pp[j]*rr))/(dd*dd*dd); 
      }
    }
  } 
  en *= Constants::EE/eps_r;

  return en;
}

void
PolarizationGrid::dot_product(unsigned int i, std::vector<double>& bins, double dr)
{
  int ii[3];
  int start[3], end[3];

  grid.element_to_index(i,ii[0],ii[1],ii[2]);

  // way of treating periodicity: 
  // define start-end grid points and fold to central cell.
  for (unsigned int l = 0; l < 3; l++)
  {
    if (_periodic[l])
    {
      start[l] = ii[l] - _dx[l];
      end[l] = ii[l] + _dx[l];
    }
    else
    {
      // start = min(0, ii[l] - _dx[l]); 
      start[l] = ii[l] - _dx[l] > 0 ? ii[l] - _dx[l] : 0;
      // start = max(ii[l] + _dx[l], _nx[l]-1); 
      end[l] = ii[l] + _dx[l] < _nx[l]-1 ? ii[l] + _dx[l] :  _nx[l]-1;
    }
  }

  for (int k1 = start[0]; k1<=end[0]; k1++)
  {
    for (int l1 = start[1]; l1<=end[1]; l1++)
    { 
      for (int m1 = start[2]; m1<=end[2]; m1++) 
      {
        // note: distance also works for folded cells
        Point rr = grid.distance(ii[0],ii[1],ii[2], k1,l1,m1);
        double dd = rr.size();

        if (dd < 0.2){ continue; }
        // fold to central cell (nx folds to 0 and so on)
        unsigned int kk1 = (k1 + _ncell[0] * _nx[0])%_nx[0];
        unsigned int ll1 = (l1 + _ncell[1] * _nx[1])%_nx[1];
        unsigned int mm1 = (m1 + _ncell[2] * _nx[2])%_nx[2];
        int j = grid.index_to_element(kk1,ll1,mm1);
        
        // define the proper bin 
        int kk = (int) (dd/dr);

        bins[kk] += (pp[i]*pp[j]) /(pp[i].size()*pp[j].size());
        
        //cout<<dd<<"  "<<kk<<"  "<<bins[kk]<<endl;

      }
    }
  } 
 
}



void
PolarizationGrid::autocorrelation(int dir, unsigned int l1, unsigned int l2, vector<double>& c)
{
  int ii[3];
  int start[3], end[3];
  int i,j;

  for (unsigned int l=0; l < _nx[dir]; l++)
  {

    switch(dir)
    {
     case 0:
       ii[0]=l; ii[1]=l1; ii[2]=l2;
       i = grid.index_to_element(l,l1,l2);
       break;
     case 1:
       ii[0]=l1; ii[1]=l; ii[2]=l2;
       i = grid.index_to_element(l1,l,l2);
       break;
     case 2:
       ii[0]=l1; ii[1]=l2; ii[2]=l;
       i = grid.index_to_element(l1,l2,l);
    }

    // way of treating periodicity: 
    // define start-end grid points and fold to central cell.
    if (_periodic[dir])
    {
       start[dir] = ii[dir] - _dx[dir];
       end[dir] = ii[dir] + _dx[dir];
    }
    else
    {
      start[dir] = ii[dir] - _dx[dir] > 0 ? ii[dir] - _dx[dir] : 0;
      end[dir] = ii[dir] + _dx[dir] < _nx[dir]-1 ? ii[dir] + _dx[dir] :  _nx[dir]-1;
    }
    
    for (int k = start[dir]; k <= end[dir]; k++)
    { 
       // array index (c[kk])
       unsigned int kk = (k + _ncell[dir] * _nx[dir])%_nx[dir];
       
       // cell index for pp
       unsigned int k1 = (k + l + _ncell[dir] * _nx[dir])%_nx[dir];

       if (k1 > _nx[dir]-1){ continue; }

       switch(dir)
       {
        case 0:
          j = grid.index_to_element(k1,ii[1],ii[2]);
          c[kk] += (pp[i](0)*pp[j](0)) / (pp[i].size()*pp[j].size());
          break;
        case 1:
          j = grid.index_to_element(ii[0],k1,ii[2]);
          c[kk] += (pp[i](1)*pp[j](1)) / (pp[i].size()*pp[j].size());
          break;
        case 2:
          j = grid.index_to_element(ii[0],ii[1],k1);
          c[kk] += (pp[i]*pp[j]) / (pp[i].size()*pp[j].size());
       }

          
    }

  }

}


void
PolarizationGrid::kmc(void)
{
  uniform_int_distribution<int> random1(0, grid.num_elements()-1);
  uniform_int_distribution<int> random2(0, 1);
  uniform_real_distribution<double> rrand(0.0, 1.0);    
  std::ostringstream ostr;
  std::ofstream file;
  string filename = get_output_directory()+"/energy.dat";
  file.open(filename.c_str());
  file<<"#iteration  Energy"<<endl;
  file.precision(8);
  ostr.precision(8);

  ostr<<"Computing total energy"<<endl;
  Messages::info(ostr.str());
  ostr.str("");
  double En_before=0.0, En_after=0.0, En_total=0.0;
  if (Nstep>=0) 
  {
    for (unsigned int l=0; l<grid.num_elements(); l++)
      En_total += energy1(l) + 0.5 * energy2(l); 
  }

  for(unsigned int l=0; l<Nstep; l++)
  {

    unsigned int i = random1(generator);
    Point pp_tmp(pp[i]);
    // NOTE: in this energy difference the factor 1/2 is NOT needed !
    En_before = energy1(i) + energy2(i);
    if (_discrete) 
    {
      rnd_orientation_discrete(pp[i]);
    }
    else
    {
      rnd_orientation(pp[i]);
    }

    En_after = energy1(i) + energy2(i);
    
    // Accept or Reject with a probability
    if (En_after < En_before)
    {
      En_total += En_after - En_before; 
    }
    else if (En_after + kbT * log(rrand(generator)) < En_before )
    {
      En_total += En_after - En_before; 
    }
    else
    {
      pp[i] = pp_tmp;
    }


    if (l%IOstep == 0)
    {
      ostr<<"Iteration "<<l<<": Energy = "<<En_total<<" dE = "<<En_after - En_before<<endl;
      Messages::info(ostr.str());
      ostr.str("");
      file<<l<<"   "<<En_total<<endl; 
    } 
  }
 
  ostr<<"Iteration "<<Nstep<<": Energy = "<<En_total<<endl;
  Messages::info(ostr.str());  
  file<<Nstep<<"   "<<En_total<<endl; 
  file.close();
}

void 
PolarizationGrid::get_solution_secure(const Elem* elem,
                           map<ID, vector<double> >& values,
                           const vector<Point>& p)
{

   if (values.count(Polarization))
   {  
     int i = grid.find_element(elem->centroid());
     double vol = grid.grid_step(0)*grid.grid_step(1)*grid.grid_step(2);
     
     for (int k=0; k<p.size(); k=k+3)
     { 
       values[Polarization][k] = pp[i](0) * Constants::e/vol*1e18;  // C/m^2
       values[Polarization][k+1] = pp[i](1) * Constants::e/vol*1e18;
       values[Polarization][k+2] = pp[i](2) * Constants::e/vol*1e18;
     }
   }

}

void 
PolarizationGrid::average_dipole()
{
   Point p;
   p(0) = 0.0; p(1) =0.0; p(2)=0.0;

   for (unsigned int i = 0; i < grid.num_elements(); ++i)
   {
     p = p + pp[i];
   }

   p = p / grid.num_elements();

   std::ostringstream os;
   os<<"Average Dipole P= "<<p(0)<<" "<<p(1)<<" "<<p(2)<<endl;
   Messages::info(os.str());

}

void
PolarizationGrid::plot_globaldata(void)
{
   double gs = get_option("grid_step",0.65); // in nanometers
   double aa = get_option("grid_step_a",gs); // in nanometers
   double bb = get_option("grid_step_b",aa); // in nanometers
   double cc = get_option("grid_step_c",gs); // in nanometers

   if (plot_solution("radial_correlation"))  
   {
     double rmax = get_option("rmax",2.0*cutoff);
     double dr = get_option("dr",0.2);
     int nint = (int) rmax/dr;
     cout<<"computing radial "<<dr<<"  "<<nint<<endl;

     std::vector<double> bins(nint,0.0);

     for (unsigned int l=0; l < grid.num_elements(); l++)
     {
        dot_product(l,bins,dr);    
     }

     double rho2d = aa*bb;
     double rho3d = aa*bb*cc;

     for (unsigned int k=1; k < nint; k++)
     {
        double r = k*dr;
        if (_nx[3]>1)
        {
           bins[k] /= grid.num_elements() * 2.0 * M_PI * r * dr * rho2d;
        }
        else
        {
           bins[k] /= grid.num_elements() * 2.0 * M_PI * r * r * dr * rho3d;
        }    
     }
      
     std::ofstream of(TiberCad::get_output_dir() + "/" + get_name() +
          "_radial.dat");

     of << "# r  radial_distribution"<<endl; 
     
     for (unsigned int k=0; k < nint; k++)
     {
        of << dr*k<< "  " << bins[k] << endl;
     }

   }

   if (plot_solution("autocorrelation"))  
   {
     vector<double> cx(_nx[0], 0.0);
     vector<double> cy(_nx[1], 0.0);
     vector<double> cz(_nx[2], 0.0);

     cout<<"computing autocorrelation x "<<_nx[0]<<endl;
     for (unsigned int l1 = 0; l1<_nx[1]; l1++)
     {
       for (unsigned int l2 = 0; l2<_nx[2]; l2++)
       {
         autocorrelation(0, l1, l2, cx);
       }
     }
     std::ofstream ofx(TiberCad::get_output_dir() + "/" + get_name() +
          "_autocorr_x.dat");
     ofx << "# x autocorrelation_x"<<endl; 
     for (unsigned int k=0; k < _nx[0]; k++)
     {
        ofx << aa*k<< "  " << cx[k]/(_nx[0]*_nx[1]*_nx[2]) << endl;
     }
     
     
     cout<<"computing autocorrelation y "<<_nx[1]<<endl;
     for (unsigned int l1 = 0; l1<_nx[0]; l1++)
     {
       for (unsigned int l2 = 0; l2<_nx[2]; l2++)
       {
         autocorrelation(1, l1, l2, cy);
       }
     }
     std::ofstream ofy(TiberCad::get_output_dir() + "/" + get_name() +
          "_autocorr_y.dat");
     ofy << "# y autocorrelation_y"<<endl; 
     for (unsigned int k=0; k < _nx[1]; k++)
     {
        ofy << bb*k<< "  " << cy[k]/(_nx[0]*_nx[1]*_nx[2])  << endl;
     }

     cout<<"computing autocorrelation z "<<_nx[2]<<endl;
     for (unsigned int l1 = 0; l1<_nx[0]; l1++)
     {
       for (unsigned int l2 = 0; l2<_nx[1]; l2++)
       {
         autocorrelation(2, l1, l2, cz);
       }
     }
     std::ofstream ofz(TiberCad::get_output_dir() + "/" + get_name() +
          "_autocorr_z.dat");
     ofz << "# z autocorrelation_z"<<endl; 
     for (unsigned int k=0; k < _nx[2]; k++)
     {
        ofz << cc*k<< "  " << cz[k]/(_nx[0]*_nx[1]*_nx[2])  << endl;
     }

   }

}

