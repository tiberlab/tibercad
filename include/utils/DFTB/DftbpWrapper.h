#ifndef _DFTBPWRAPPER_H_
#define _DFTBPWRAPPER_H_

#include "dftbp.h"
#include <iostream>

//-----------------------------------------------------------------------

class DftbpWrapper
{

public:

//!Wrapper class for callings to DFTB+ library (libdftbp.so)



//!Constructor
//!Assign an handler to DFTB+ instance, transparent to programmer


  DftbpWrapper();
  
  ~DftbpWrapper();

  static DftbpWrapper* create(); 

  void fill_param (int nAtoms, int nType, double eTemp, int iPeriodic,
                   char *speciesNames, int *species);

  void addskdata (char *skNames, int *mAngs, int orbResolved, int skInterp,
                  int nType);

  void initdftb ();

  void up_coords (int nAtom, double *newCoords);

  void get_energy (double &energy);

  void addlattice (double *latVecs);

  void addkpoints (int nKPoint, double *kPoints, double *kWeights);

private:
  int _handler[DFTBP_HSIZE];

};




#endif
