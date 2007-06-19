#include "ParticleThermalConductivity.h"
#include "getpot.h"
#include "Material.h"
#include "Database.h"
#include "Constants.h"


//-------------------------------------------------------------------------//

void  ParticleThermalConductivity::copy_from(const PhysicalModelInterface *rhs)
{
  const ParticleThermalConductivity* mod = dynamic_cast<const ParticleThermalConductivity*> (rhs);

  _kappa_e = mod -> _kappa_e;

  _kappa_h = mod -> _kappa_h;


}

//-------------------------------------------------------------------------//


void   ParticleThermalConductivity::calculate_VCA (const PhysicalModelInterface *comp_A, 
                                                const PhysicalModelInterface *comp_B, double xa) 
{ 
  const  ParticleThermalConductivity* modA = dynamic_cast<const  ParticleThermalConductivity*>(comp_A);

  const ParticleThermalConductivity* modB = dynamic_cast<const  ParticleThermalConductivity*>(comp_B);


   alloy(_kappa_e,modA->_kappa_e, modB->_kappa_e, xa);  

   alloy(_kappa_h,modA->_kappa_h, modB->_kappa_h, xa); 
  
}

//-------------------------------------------------------------------------//

//--------------------------------------------------------//
void   ParticleThermalConductivity::read_database(void)
{

  const ModelOptions& options = get_options();

  _particle_kappa_model = options.get_option("model","constant");

  if  ( _particle_kappa_model.compare("constant") == 0 )
  {
    const Material* mat = get_material();
    GetPot data((mat->get_database()).get_data_file());
    
    _kappa_e_x = data("kappa_e", 0.0);
    _kappa_h_x = data("kappa_h", 0.0); 
    
  }	
  
}

//---------------------------------------------------------//

void  ParticleThermalConductivity::do_init(void)
{

     
  if (  _particle_kappa_model.compare("constant") == 0)
  {
    
    ModelOptions& options = get_options ();
    
    _kappa_e_x = options.get_option("kappa_e", _kappa_e_x );
    
    _kappa_h_x = options.get_option("kappa_h", _kappa_h_x );
    
  }
  else if  ( _particle_kappa_model.compare("Wiedemann_Franz") == 0 ) 
  {
    //Here we assume a scalar electronic conducibility

    _kappa_e_x =  _Tloc * _sigma_e * Constants::Lorenz_Number;

    _kappa_h_x =  _Tloc * _sigma_h * Constants::Lorenz_Number; 
    
    
    _kappa_e(1,1) =  _kappa_e_x;
    _kappa_e(2,2) =  _kappa_e_x;
    _kappa_e(3,3) =  _kappa_e_x;
    
    _kappa_h(1,1) =  _kappa_h_x;
    _kappa_h(2,2) =  _kappa_h_x;
    _kappa_h(3,3) =  _kappa_h_x;  


  }
  else
  {


  throw InitFailedException("Unrecongnition the particle thermal conductivity model:  " +  _particle_kappa_model);
  

  } 

 
  
}

void  ParticleThermalConductivity::re_init(void)
{
 
  this->do_init();
 

  
}







//-------------------------------------------------------------------------//
