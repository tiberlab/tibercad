#ifndef _ALLOYMODEL_H_
#define _ALLOYMODEL_H_


#include <string>

class AlloyModel
{
 public:
  AlloyModel() {};
  ~AlloyModel() {};
  
  
  void set_model(const std::string& model);  

  void set_x_min(const double  x_min);

  void set_x_max(const double  x_max);

  void set_x_constant(const double  x_constant);

  const std::string&  get_model();

  const double  get_x_min();

  const double  get_x_max();


  
 private:

  std::string _model_type;
  double  _x_min, _x_max;


};

//
// inline  methods 
//

inline void
AlloyModel::set_model(const std::string& model) 
{
  
  _model_type = model;
  
}

inline void
AlloyModel::set_x_min(const double  x_min ) 
{
  
  _x_min = x_min;
  
}


inline void
AlloyModel::set_x_max(const double  x_max)  
{
  
  _x_max = x_max;
  
}

inline void
AlloyModel::set_x_constant(const double  x_constant)  
{
  
  _x_max = x_constant;
  _x_min = x_constant;
  
}


inline const std::string& 
AlloyModel::get_model()
{
  
  return _model_type;
  
}

inline const double  
AlloyModel::get_x_min() 
{
  
  return _x_min;
  
}


inline const double
AlloyModel::get_x_max()
{
  
  return _x_max;
  
}



#endif /*_ALLOYMODEL_H_*/








