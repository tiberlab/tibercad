// $Id$

#ifndef _ZBOPTDIELECTRICCONSTANT_H_  
#define _ZBOPTDIELECTRICCONSTANT_H_ 


#include  "OptDielectricConstant.h"

class  ZbOptDielectricConstant: public OptDielectricConstant

{
 public:
  //!constructor
  ZbOptDielectricConstant() {};
 
  //!destructor
  virtual ~ZbOptDielectricConstant() {};


 inline  static ZbOptDielectricConstant* create();

 


 private:

 
 //!epsilon along [100] direction
 double _eps; 
 

  
 

 protected:
 
  virtual void read_database(void);


  virtual void do_init(void);

  
  inline  virtual PhysicalModelInterface*  create_new (void) const;

};




inline
ZbOptDielectricConstant* ZbOptDielectricConstant::create()

{
  return (new ZbOptDielectricConstant());
}

inline
PhysicalModelInterface*  ZbOptDielectricConstant::create_new (void) const
{
  return (new ZbOptDielectricConstant() ); 
}




#endif
