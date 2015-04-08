#include "PolarizationGrid.h"
#include "InitFailedException.h"
#include "Constants.h"
#include "mesh_generation.h" 
#include "SimulationEnvironment.h"
#include "Messages.h"
#include <sstream> 
#include <fstream> 
#include "TiberModule.h"
#include "gmsh_io.h"

PolarizationGrid::PolarizationGrid(const ModelOptions& options) :
 SimulationInterface(options),
  eps_r(1.0),
  cutoff(100.0),
  dipole(0.0),
  Nstep(100000),
  kbT(0.026)
{
}

PolarizationGrid*
PolarizationGrid::create(const ModelOptions& opt)
{
  return new PolarizationGrid(opt);
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

  declare_solution(Polarization, VECTOR, CELL, "C/m^2");

}
 
void
PolarizationGrid::setup_mesh(void)
{
  // the grid must be created first (before any other option)
  nx = get_option("Nx",10); 
  ny = get_option("Ny",10); 
  nz = get_option("Nz",1); 

  double aa = get_option("grid_step",0.65); // in nanometers

  Point p1(0,0,0);
  Point p2(nx*aa,ny*aa,nz*aa);
  setup_grid(p1,p2,nx,ny,nz);
  //create internal rectangular mesh 
  //
  Mesh* mymesh = new Mesh(3);
  ElemType type(HEX8);

  MeshTools::Generation::build_cube(*mymesh, 
				       nx, ny, nz,
               0.0, nx * aa, 
               0.0, ny * aa,
               0.0, nz * aa,
               type);

  // set in environment
  get_environment().set_mesh(mymesh);

  // set in my object
  set_mesh(mymesh);

  GmshIO(*mymesh).write("mymesh.msh");

}


void
PolarizationGrid::do_solve()
{
  this->kmc();

  if (get_options().find_option("write_file"))
        write_dipoles();

  write_antiferro();
}

void
PolarizationGrid::parse_options(void)
{
  Nstep = get_option("Nstep",10);
  IOstep = get_option("IOstep",1000);
  eps_r = get_option("eps_r",1.0);
  string units = get_option("dipole_units","Debye");
  double dip = get_option("dipole", 2.29);
  set_dipole(dip,units);

  vector<double> E(3,0.0);
  get_option("external_field",E);
  set_efield(E);

  kbT = get_option("kbT", 0.0);
  
  seed = get_option("seed",
      static_cast<int>(time(NULL) * random_device()()));

}



void
PolarizationGrid::set_cutoff(void)
{
  // e^2/(4 pi eps0) = 1.44 eV nm

  double tmp = 4 * Constants::EE * dipole * dipole / eps_r;
 
  // cutoff in nm to get 1 meV energy difference
  cutoff = pow(tmp/0.001,1.0/3.0);

  cutoff = get_option("cutoff",cutoff);
  
  dk =static_cast<int>( ceil(cutoff/grid.grid_step(0)) );
  dl =static_cast<int>( ceil(cutoff/grid.grid_step(1)) );
  dm =static_cast<int>( ceil(cutoff/grid.grid_step(2)) );
}

void
PolarizationGrid::set_dipole(double p, string str)
{
  // The dipole must be transformed in e*nm
  if (str == "Debye")
  {
    // 1 Debye = 3.336 e-30 C m = 0.0208 e * nm
    dipole = p * 0.0208;
  }
  else if (str == "Cm")
  {
    // 1 C m =  e/(1.60e-19) *  nm/1e9 
    dipole = p / 1.60e-10;
  }
  else if (str == "C/m^2")
  {
    // P*vol = C/m^2 * nm*nm^2 = C/m^2 * nm*1e-18 m^2 = 1e-18/1.60e-19 e*nm;
    double vol = grid.grid_step(0)*grid.grid_step(1)*grid.grid_step(2);
    dipole = p * vol / 0.1602;
  }
  else if (str == "eA")
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
   uniform_int_distribution<int> random1(0, 179);
   uniform_int_distribution<int> random2(0, 89);

   double teta;

   if (grid.num_elements(2) > 1)
   {
     teta = static_cast<double>(random2(generator)) * 2.0; 
   }
   else
   {
     teta = 90.0;
   }

   double phi = static_cast<double>(random1(generator)) * 2.0; 

   p(0) = cos(phi * M_PI/180.0) * sin(teta * M_PI/180.0);
   p(1) = sin(phi * M_PI/180.0) * sin(teta * M_PI/180.0);    
   p(2) = cos(teta * M_PI/180.0);

   p *= dipole; 
}

void
PolarizationGrid::write_dipoles()
{
   string out_path = get_output_directory();
   string file_name = out_path+"/"+get_option("filename","polarization.dat");
   ofstream file; 

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
   ifstream file; 

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
   ofstream file; 

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
PolarizationGrid::set_rnd_dipoles(void)
{

  for (unsigned int i=0; i < pp.size(); i++) 
  {
    this->rnd_orientation(pp[i]);
  }
    
}


double
PolarizationGrid::energy1(unsigned int i)
{
  // e*nm  * V/nm  = eV
  return pp[i]*efield[i];
}

double 
PolarizationGrid::energy2(unsigned int i)
{

  double en = 0.0;
  int k,l,m;

  grid.element_to_index(i,k,l,m);

  unsigned int kstart = (k-dk >= 0 ? k-dk : 0);
  unsigned int kend = (k+dk <= nx-1 ? k+dk : nx-1);
  unsigned int lstart = (l-dl >= 0 ? l-dl : 0);
  unsigned int lend = (l+dl <= ny-1 ? l+dl : ny-1);  
  unsigned int mstart = (m-dm >= 0 ? m-dm : 0);
  unsigned int mend = (m+dm <= nz-1 ? m+dm : nz-1);

  for (unsigned int k1 = kstart; k1<=kend; k1++)
  {
    for (unsigned int l1 = lstart; l1<=lend; l1++)
    { 
      for (unsigned int m1 = mstart; m1<=mend; m1++) 
      {
        int j = grid.index_to_element(k1,l1,m1);
        if (j==i) continue;
        Point rr = grid.distance(k,l,m,k1,l1,m1);
        double dd = rr.size();
        rr /= dd; 
        en += (pp[i]*pp[j] + 3.0*(pp[i]*rr) * (pp[j]*rr))/(dd*dd*dd); 
      }
    }
  } 
  en *= Constants::EE/eps_r;

  return en;
}


void
PolarizationGrid::kmc(void)
{
  uniform_int_distribution<int> random1(0, grid.num_elements()-1);
  uniform_int_distribution<int> random2(0, 1);
  uniform_real_distribution<double> rrand(0.0, 1.0);    
  ostringstream ostr;
  ofstream file;
  string filename = get_output_directory()+"/energy.dat";
  file.open(filename.c_str());
  file<<"#iteration  Energy"<<endl;
  file.precision(8);
  ostr.precision(8);

  double En_before=0.0, En_after=0.0, En_total=0.0;
  for (unsigned int l=0; l<grid.num_elements(); l++)
    En_total += energy1(l) + 0.5 * energy2(l); 

  for(unsigned int l=0; l<Nstep; l++)
  {

    unsigned int i = random1(generator);
    Point pp_tmp(pp[i]);
    // NOTE: in this energy difference the factor 1/2 is NOT needed !
    En_before = energy1(i) + energy2(i);
    rnd_orientation(pp[i]);
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
      ostr<<"Iteration "<<l<<": Energy = "<<En_total<<endl;
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
       values[Polarization][k] = pp[i](0);
       values[Polarization][k+1] = pp[i](1);
       values[Polarization][k+2] = pp[i](2);
     }
   }

}


