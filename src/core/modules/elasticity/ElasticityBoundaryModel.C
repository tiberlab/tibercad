// $Id$

#include "ElasticityBoundaryModel.h"
#include "BoundaryModel.h"


using namespace std;




// void
// ElasticityBoundaryModel::do_init(void)
// {

 
//   PhysicalModelInterface::SubmodelIterator  it;
//   it = submodels_begin("boundary"); 
//   if (it != submodels_end("boundary"))
//   {
//     _bm = dynamic_cast<BoundaryModel*> ((*it).second);
//     _bm->get_coefficients(_H_tens, _R_vec);
//   }
  

// }

//! Calculate properties 

//void 
//ElasticityBoundaryModel::calculate(const Elem* elem, unsigned int side,
//		   const Point& point)
//{

//  if (_bm != NULL)
//  {
//    _bm->set_normal(_normal);
//    _bm->calculate(elem,side,point);
//    _bm->get_coefficients(_H_tens, _R_vec);
//  }
//}


ElasticityBoundaryModel*
ElasticityBoundaryModel::create(const ModelOptions& options)
{
 

   std::string type = options.get_option("type", "clamp");

   std::cout<<"ebnd_" + type<<std::endl;
  ElasticityBoundaryModel* mod = dynamic_cast<ElasticityBoundaryModel*>(
      PhysicalModelInterface::create("ebnd_" + type, options));

  if (mod == NULL)
  {
    ostringstream os;
    os << "Elasticity boundary model \'" << type << "\' cannot be found.";
    throw InitFailedException(os.str());
  }

  return mod;

}


