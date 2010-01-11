#ifndef _WZPIEZOELECTRICMODEL_H_
#define _WZPIEZOELECTRICMODEL_H_

#include "PiezoelectricModel.h"

class WzPiezoelectricModel : public PiezoelectricModel
{

 public:

  //!Empty constructor
  WzPiezoelectricModel(const ModelOptions& options) : PiezoelectricModel(options) {};
  ~WzPiezoelectricModel(){};  


 inline static WzPiezoelectricModel* create(const ModelOptions& options);

 void calculate_piezopolarization(const Elem* elem, const Point& p);

 protected:

  //reads database
  virtual void read_database ( );
  
  virtual  PhysicalModelInterface* create_new(void) const;

  virtual void do_init (void);


private:
  
  double e33;
  
  double e31;
  
  double e15;

  
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


inline  PhysicalModelInterface* WzPiezoelectricModel::create_new(void) const
{
  return ( new WzPiezoelectricModel(get_options()) ) ;
}


inline WzPiezoelectricModel* WzPiezoelectricModel::create(const ModelOptions& options)
{
   return new WzPiezoelectricModel(options) ;
}

#endif
