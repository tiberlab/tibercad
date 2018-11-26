// $Id$

#ifndef _ZBOPTDIELECTRICCONSTANT_H_  
#define _ZBOPTDIELECTRICCONSTANT_H_ 


#include  "OptDielectricConstant.h"

class  ZbOptDielectricConstant: public OptDielectricConstant

{
 public:
 
  //!destructor
  virtual ~ZbOptDielectricConstant() {};


  static ZbOptDielectricConstant* create(const ModelOptions& options);

 


 private:

 
 //!epsilon along [100] direction
 double _eps; 
 

  
 

 protected:
 
  //!constructor
  ZbOptDielectricConstant(const ModelOptions& options);

  virtual void read_database(void);


  virtual void do_init(void);

  
  inline  virtual PhysicalModelInterface*  create_new (void) const;

};



inline
ZbOptDielectricConstant::ZbOptDielectricConstant(const ModelOptions& options) :
  OptDielectricConstant(options)
{
}

inline
ZbOptDielectricConstant* ZbOptDielectricConstant::create(const ModelOptions& options)

{
  return (new ZbOptDielectricConstant(options));
}

inline
PhysicalModelInterface*  ZbOptDielectricConstant::create_new (void) const
{
  return (new ZbOptDielectricConstant(get_options()) );
}




#endif
