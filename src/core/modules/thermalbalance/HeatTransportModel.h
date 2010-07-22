// $Id$

#ifndef _HEATTRANSPORTMODEL_H_
#define _HEATTRANSPORTMODEL_H_

#include "PhysicalModelInterface.h"


#define TIBER_MODULE_PREFIX thermal_model

class Elem;
class Point;

using namespace std;

//! The base class for Poisson boundary conditions
//class TBDLEXPORT HeatTransportModel : public PhysicalModel
class HeatTransportModel : public PhysicalModelInterface
{

  public:

    //! Destructor
    ~HeatTransportModel(void) {};

     //! Creator function
    static HeatTransportModel* create(const ModelOptions& options);

  //HeatTransportType
  enum Type
    {
      Fourier = 0,
      Gray = 1
    };

 Type get_type(void) const;

  protected:

    //! Constructor
    HeatTransportModel(const ModelOptions& options);

  void set_type(Type type);
 
  private:

  
 Type type;

};

inline 
HeatTransportModel::Type
HeatTransportModel::get_type(void) const
{
  return type;
}

inline 
void 
HeatTransportModel::set_type(Type type_in)
{
  type = type_in;
}



inline
HeatTransportModel::HeatTransportModel(const ModelOptions& options) :
  PhysicalModelInterface(options)
 {
 }



// inline
// HeatTransportModel* 
// HeatTransportModel::create(const ModelOptions& options)
//  {

   
//    std::string name = options.get_option("type", "fourier");

//    HeatTransportModel* mod = dynamic_cast<HeatTransportModel*>(
//        PhysicalModelInterface::create("heat_transport_" + name, options));

//    if (mod == NULL)
//    {
//      ostringstream os;
//      os << "Heat transport model \'" << name << "\' cannot be found.";
//      throw InitFailedException(os.str());
//    }

//    return mod;
//  }






#endif // _HEATTRANSPORTMODEL_H_
