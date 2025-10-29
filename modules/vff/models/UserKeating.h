#ifndef _USERKEATING_H_
#define _USERKEATING_H_

#include "tiber_dll.h"
#include "Keating.h"

//! User defined Keating model parameters
class TBDLLOCAL UserKeating : public Keating
{
public:

  //! Destructor
  ~UserKeating(void) {};

  //! Creator function
  static UserKeating* create(const ModelOptions& options);

  double get_alpha();

  //! Assign value to parameters
  void do_init(void);

protected:

private:

  UserKeating(const ModelOptions& options);

  void assign_alpha(void);

  void assign_beta(void);

  //! Parse parent material parameters. To be used if the models belongs to an alloy
  void assign_alpha_parents(void);

  void assign_beta_parents(void);


};

inline
UserKeating*
UserKeating::create(const ModelOptions& options)
{
  return new UserKeating(options);
}



#endif
