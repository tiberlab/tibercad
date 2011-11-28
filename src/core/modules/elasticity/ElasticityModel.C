// $Id$

#include "ElasticityModel.h"



using namespace std;


ElasticityModel*
ElasticityModel::create(const Material* mat, const ModelOptions& options)
{

  return PhysicalModelInterface::create<ElasticityModel>(_create, _destroy, mat, options);

}


void
ElasticityModel::do_init(void)
{
  PhysicalModelInterface::SubmodelIterator  it;
  it = submodels_begin("stiffness");
  StiffnessModel* _sm = dynamic_cast<StiffnessModel*> ((*it).second);
  _stiffness =  _sm->get_stiffness();


}

//! Print some useful information
void
ElasticityModel::do_print_info(void)
{
Messages::info("Stiffness Constant:");
ostringstream os;
os <<"  C11: "<<_stiffness(1,1,1,1)<<" GPa\n";
os <<"  C12: "<<_stiffness(2,2,1,1)<<" GPa\n";
os <<"  C13: "<<_stiffness(3,3,1,1)<<" GPa\n";
os <<"  C44: "<<_stiffness(3,2,3,2)<<" GPa\n";
Messages::info(os.str());
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
ElasticityModel::prepare_submodels(void)
{
  ModelOptions opts;
  opts.set_option("type", "anisotropic");
  PhysicalModelInterface* pm;
  create_submodel(pm, "stiffness", opts);

  create_submodels(_bfm, "body_force");
}
