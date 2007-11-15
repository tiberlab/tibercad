#include "DftbpWrapper.h"

//---------------------------------------------------------------------


DftbpWrapper::DftbpWrapper(){
  std::cout << "Constructing DFTB instance... ";
    f77_dftbp_init(_handler);
    std::cout << "done." << std::endl;
    std::cout << "Received handler: ";
    for  (int ii = 0; ii < DFTBP_HSIZE; ++ii) {
      std::cout << _handler[ii] << " ";
    }
    std::cout << std::endl;
}


DftbpWrapper::~DftbpWrapper(){
  std::cout << "Destructing DFTB instance... ";
    f77_dftbp_destruct(_handler);
    std::cout << "done." << std::endl;
}


DftbpWrapper* DftbpWrapper::create()
{
  return new DftbpWrapper();
}



//!Assign simulation parameters to DFTB+ instance
void DftbpWrapper::fill_param(int nAtoms, int nType, double eTemp, int iPeriodic,
			char *speciesNames, int *species) {

 f77_dftbp_fillbasicparameters(_handler, nAtoms, nType, speciesNames,
                                species, iPeriodic, eTemp);

} 



//!Set Slater Koster parameters (read from file)
void DftbpWrapper::addskdata (char *skNames, int *mAngs, int orbResolved, int skInterp,
                  int nType) {
    f77_dftbp_addskdatafromfile(_handler, nType, skNames, mAngs, orbResolved,
                              skInterp);
  }


//!Initialize DFTB+ instance (allocations)
void DftbpWrapper::initdftb () {
    f77_dftbp_initdftb(_handler);
  }


//!Set atoms coordinates
void DftbpWrapper::up_coords (int nAtom, double *newCoords) {
    f77_dftbp_updatecoords(_handler, nAtom, newCoords);
  }


//!Solve total system energy 
void DftbpWrapper::get_energy (double &energy) {
    f77_dftbp_gettotalenergy(_handler, energy);
  }


void DftbpWrapper::addlattice (double *latVecs) {
    f77_dftbp_addlattice(_handler, latVecs);
  }


//!Set k points for calculations
void DftbpWrapper::addkpoints (int nKPoint, double *kPoints, double *kWeights) {
    f77_dftbp_addkpoints(_handler, nKPoint, kPoints, kWeights);
  }


void DftbpWrapper::getchargesperatom(int nAtom, double* charges) {
	f77_dftbp_getchargesperatom (_handler, nAtom, charges);
}
