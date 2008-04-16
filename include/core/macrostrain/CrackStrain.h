#ifndef _CRACKSTRAIN_H_
#define _CRAINSTRAIN_H_
#include "StrainSimulation.h"

class CrackStrain : public StrainSimulation
{
 public:

 
  // virtual void get_solution_secure(const Elem* elem,
  //				  const std::vector<Point>& p, const std::set<ID>& ids,
  //				   std::vector<std::map<ID, double> >& values) {} ;


  static CrackStrain* create(void);

 protected:

  virtual void parse_options(void);
  
  virtual void do_init(void) ;
 
  virtual void do_solve(void);


 private:

  double _Ki;

  double _x0;

  double _y0;

};


inline
CrackStrain* CrackStrain::create()
{
  return new CrackStrain();
}

#endif
