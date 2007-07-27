#ifndef _WZOPTDIELECTRICCONSTANT_H_  
#define _WZOPTDIELECTRICCONSTANT_H_ 


#include  "OptDielectricConstant.h"

class  WzOptDielectricConstant: public OptDielectricConstant

{
 public:
  //!constructor
WzOptDielectricConstant() {};
 
 //!destructor
 ~WzOptDielectricConstant() {};


 inline  static WzOptDielectricConstant* create();

  virtual void update_tensor();


 private:

  std::string _eps_model;

  double _eps_a_x; 
  double _eps_b_x; 
  double _eps_c_x;
  double _eps_x;

  double _eps_a_z; 
  double _eps_b_z; 
  double _eps_c_z;
  double _eps_z;

 protected:
 
  virtual void read_database(void);


  virtual void do_init(void);

 inline  virtual PhysicalModelInterface*  create_new (void) const;

};




inline
WzOptDielectricConstant* WzOptDielectricConstant::create()

{
  return (new WzOptDielectricConstant());
}

inline
PhysicalModelInterface*  WzOptDielectricConstant::create_new (void) const
{
  return (new WzOptDielectricConstant() ); 
}




#endif
