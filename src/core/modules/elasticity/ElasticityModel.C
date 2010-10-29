// $Id$

#include "ElasticityModel.h"



using namespace std;




void
ElasticityModel::do_init(void)
{

 

   PhysicalModelInterface::SubmodelIterator  it;
   it = submodels_begin("stiffness");
   StiffnessModel* _sm = dynamic_cast<StiffnessModel*> ((*it).second);
   _stiffness =  _sm->get_stiffness();


   //BodyForceModel
   it = submodels_begin("body_force");
   const PhysicalModelInterface::SubmodelIterator  it_end(submodels_end("body_force"));
   for ( ; it != it_end ; ++it)
     _bfm.push_back(dynamic_cast<BodyForceModel*> ((*it).second));
  
}



void
ElasticityModel::calculate(const Elem* elem, const Point& point)
{
 
  //Update Body force model
  _force = 0;
  _strain = 0;
  _stress = 0;
  for (ID n = 0 ; n <_bfm.size() ; n++)
  {
    _bfm[n]->calculate(elem,point);
    _force  +=  _bfm[n]->get_force_source();
    _strain +=  _bfm[n]->get_strain_source();
    _stress +=  _bfm[n]->get_stress_source();
  }

}

void
ElasticityModel::create_submodels(void)
{

  //Heat Transport Default
  if (!get_options().has_submodel("stiffness"))
  {
   
    ModelOptions opts;
    opts.set_option("type","anisotropic");
    get_options().add_submodel("stiffness",opts);  

  }


}
