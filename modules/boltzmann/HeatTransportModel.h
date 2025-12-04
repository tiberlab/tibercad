/*  
 * This file is part of the tiberCAD module boltzmann.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file HeatTransportModel.h
 * \brief tiberCAD boltzmann module header.
 *
 * \note This file is part of module boltzmann.
 */


#ifndef _HEATTRANSPORTMODEL_H_
#define _HEATTRANSPORTMODEL_H_

#include "tibercad/physics/PhysicalModel.h"



using namespace std;

//! The base class for Poisson boundary conditions
//class TBDLEXPORT HeatTransportModel : public PhysicalModel
class HeatTransportModel : public PhysicalModel
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
  PhysicalModel(options)
 {
 }



// inline
// HeatTransportModel* 
// HeatTransportModel::create(const ModelOptions& options)
//  {

   
//    std::string name = options.get_option("type", "fourier");

//    HeatTransportModel* mod = dynamic_cast<HeatTransportModel*>(
//        PhysicalModel::create("heat_transport_" + name, options));

//    if (mod == NULL)
//    {
//      ostringstream os;
//      os << "Heat transport model \'" << name << "\' cannot be found.";
//      throw InitFailedException(os.str());
//    }

//    return mod;
//  }






#endif // _HEATTRANSPORTMODEL_H_
