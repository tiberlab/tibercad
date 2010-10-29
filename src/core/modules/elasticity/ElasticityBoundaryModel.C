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
    BoundaryModel* _bm = dynamic_cast<BoundaryModel*> ((*it).second);
    _bm->get_coefficients(_H, _A, _R);
  }
  cout<<_R<<endl;

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
