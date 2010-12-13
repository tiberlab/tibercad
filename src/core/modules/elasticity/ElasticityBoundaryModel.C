// $Id$

#include "ElasticityBoundaryModel.h"
#include "BoundaryModel.h"


using namespace std;




void
ElasticityBoundaryModel::do_init(void)
{

 
  PhysicalModelInterface::SubmodelIterator  it;
  it = submodels_begin("boundary"); 
  if ((*it).second != NULL)
  {
    _bm = dynamic_cast<BoundaryModel*> ((*it).second);
    _bm->get_coefficients(_H, _R);
  }
  

}

//! Calculate properties 
void 
ElasticityBoundaryModel::calculate(const Elem* elem, unsigned int side,
		   const Point& point)
{

  if (_bm != NULL)
  {
    _bm->set_normal(_normal);
    _bm->calculate(elem,side,point);
    _bm->get_coefficients(_H, _R);
  }


}

// void
// ElasticityBoundaryModel::create_submodels(void)
// {

//   //Heat Transport Default
//   if (!get_options().has_submodel("boundary"))
//   {
   
//     ModelOptions opts;
//     opts.set_option("type","clamp");
//     get_options().add_submodel("boundary",opts);  

//   }


// }
