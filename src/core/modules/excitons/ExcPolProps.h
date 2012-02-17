#ifndef _EXCPOLPROPS_H_
#define _EXCPOLPROPS_H_

#include <string>
#include <map>

class ExcPolProps {
  public:
    static const double a0 = 1e3;
    static const double b0 = 1e-9;
    static const double pol_tau0 = 1e-11;
    
    double a;
    double b;
    double pol_tau;
    double density_renormalization;

    ExcPolProps() {
      a = a0;
      b = b0;
      pol_tau = pol_tau0;

      density_renormalization = 1;
    }
};





#endif /* _EXCPOLPROPS_H_*/
