// $Id$

#ifndef _WZOPTDIELECTRICCONSTANT_H_  
#define _WZOPTDIELECTRICCONSTANT_H_ 


#include  "OptDielectricConstant.h"

class  WzOptDielectricConstant: public OptDielectricConstant

{
 public:
  //!constructor
  WzOptDielectricConstant() {};
 
 //!destructor
 virtual ~WzOptDielectricConstant() {};


 inline  static WzOptDielectricConstant* create();

 


 private:

 
 //!epsilon along [10-10] direction
 double _eps_a; 
 
 //!epsilon along [0001] direction
 double _eps_c; 
  
 

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
