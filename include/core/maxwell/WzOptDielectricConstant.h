// $Id$

#ifndef _WZOPTDIELECTRICCONSTANT_H_  
#define _WZOPTDIELECTRICCONSTANT_H_ 


#include  "OptDielectricConstant.h"

class  WzOptDielectricConstant: public OptDielectricConstant

{
 public:
 
 //!destructor
 virtual ~WzOptDielectricConstant() {};


 inline  static WzOptDielectricConstant* create(const ModelOptions& options);

 


 private:

 
 //!epsilon along [10-10] direction
 double _eps_a; 
 
 //!epsilon along [0001] direction
 double _eps_c; 
  
 

 protected:
 
  //!constructor
  WzOptDielectricConstant(const ModelOptions& options);

  virtual void read_database(void);


  virtual void do_init(void);

  
  inline  virtual PhysicalModelInterface*  create_new (void) const;

};


inline
WzOptDielectricConstant::WzOptDielectricConstant(const ModelOptions& options) :
  OptDielectricConstant(options)
{
}


inline
WzOptDielectricConstant* WzOptDielectricConstant::create(const ModelOptions& options)

{
  return new WzOptDielectricConstant(options);
}

inline
PhysicalModelInterface*  WzOptDielectricConstant::create_new (void) const
{
  return new WzOptDielectricConstant(get_options());
}




#endif
