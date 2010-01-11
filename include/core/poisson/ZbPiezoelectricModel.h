#ifndef _ZBPIEZOELECTRICMODEL_H_
#define _ZBPIEZOELECTRICMODEL_H_

#include "PiezoelectricModel.h"

class ZbPiezoelectricModel : public PiezoelectricModel
{

 public:

  //!Empty constructor
  ZbPiezoelectricModel(const ModelOptions& options) : PiezoelectricModel(options) {};
  ~ZbPiezoelectricModel(){};  


 inline static ZbPiezoelectricModel* create(const ModelOptions& options);

 void calculate_piezopolarization(const Elem* elem,const Point& point);

 protected:

  //reads database
  virtual void read_database ( );
  
  virtual  PhysicalModelInterface* create_new(void) const;

  virtual void do_init (void);


private:
  
  double e14;

  
  enum strain_variables
  {
     EXX = 0, 
     EXY,  
     EXZ,
     EYY,
     EYZ,
     EZZ
   };

   //!Strain variables 
   std::set<ID> ID_set;

   //!Variable map
   std::map<ID,ID> var_map;
};


inline  PhysicalModelInterface* ZbPiezoelectricModel::create_new(void) const
{
  return ( new ZbPiezoelectricModel(get_options()) ) ;
}


inline ZbPiezoelectricModel* ZbPiezoelectricModel::create(const ModelOptions& options)
{
   return new ZbPiezoelectricModel(options) ;
}

#endif
