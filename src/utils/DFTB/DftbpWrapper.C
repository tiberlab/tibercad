#include "DftbpWrapper.h"

//---------------------------------------------------------------------


DftbpWrapper::DftbpWrapper()
{
  std::cout << "Constructing DFTB instance... ";
  f77_dftbp_init(_handler);

  for  (int ii = 0; ii < DFTBP_HSIZE; ++ii) {
    std::cout << _handler[ii] << " ";
  }
  std::cout << std::endl;
}


DftbpWrapper::~DftbpWrapper()
{
  std::cout << "Destructing DFTB instance... ";
  f77_dftbp_destruct(_handler);
  std::cout << "done." << std::endl;
}


DftbpWrapper* DftbpWrapper::create()
{
  return new DftbpWrapper();
}



//!Assign simulation parameters to DFTB+ instance
void
DftbpWrapper::fill_param(int nAtoms, int nType, double eTemp, int iPeriodic,
    char *speciesNames, int *species)
{
  f77_dftbp_fillbasicparameters(_handler, nAtoms, nType, speciesNames,
      species, iPeriodic, eTemp);
}



//!Set Slater Koster parameters (read from file)
void
DftbpWrapper::addskdata (char *skNames, int *mAngs, int orbResolved, int skInterp,
    int nType)
{
  f77_dftbp_addskdatafromfile(_handler, nType, skNames, mAngs, orbResolved,
      skInterp);
}


//!Initialize DFTB+ instance (allocations)
void
DftbpWrapper::initdftb ()
{
  f77_dftbp_initdftb(_handler);

}


//!Set atoms coordinates
void
DftbpWrapper::up_coords (int nAtom, double *newCoords)
{
  f77_dftbp_updatecoords(_handler, nAtom, newCoords);
}


//!Solve total system energy
void
DftbpWrapper::get_energy (double &energy)
{
  f77_dftbp_gettotalenergy(_handler, energy);
}


void
DftbpWrapper::addlattice (double *latVecs)
{
  f77_dftbp_addlattice(_handler, latVecs);
}


//!Set k points for calculations
void
DftbpWrapper::addkpoints (int nKPoint, double *kPoints, double *kWeights)
{
  f77_dftbp_addkpoints(_handler, nKPoint, kPoints, kWeights);
}


//! Set automatic kpoints, based on supercell folding
//! See DFTB+ 1.0 manual for further informations.
void
DftbpWrapper::addsupersampling(double *coeffs, double *shifts, int noinv)
{
  f77_dftbp_addsupersampling(_handler, coeffs, shifts, noinv);
}


void
DftbpWrapper::getchargesperatom(int nAtom, double* charges)
{
  f77_dftbp_getchargesperatom (_handler, nAtom, charges);
}


void
DftbpWrapper::getnetchargesperatom(int nAtom, double* charges)
{
  f77_dftbp_getnetchargesperatom (_handler, nAtom, charges);
}


void DftbpWrapper::getmatrix(int &nrow, int &ncol, int &nzval, int &isreal, int* &colind, int* &rowpnt, double* &val, std::string matrix, double *kPoint)
{
  f77_dftbp_recreatecsrbuffer(_handler, nrow, ncol, nzval, isreal);
  int nrow1 = nrow + 1;

  if (kPoint == NULL) {

    std::cerr << "Pointer in Wrapper " << colind << std::endl;
    colind = new int[nzval];
    rowpnt = new int[nrow1];
    val = new double[nzval];
    std::cerr << "Pointer in Wrapper " << colind << std::endl;

    if (isreal == 0) std::cerr << "ERROR: requested a real H or S while complex calculation is computed" << std::endl;

    if (matrix.compare("H") == 0) {
      std::cerr << "Getting real Hamiltonian " << std::endl;
      f77_dftbp_getrcsrhamiltonian (_handler, nrow1, nzval, colind, rowpnt, val);
      std::cout << "First value of H " << colind[0] << " " << rowpnt[0] << " " << val[0] << " " << std::endl;
      std::cout << "Second value of H " << colind[1] << " " << rowpnt[0] << " " << val[1] << " " << std::endl;
      std::cout << "Number of rows is " << nrow << std::endl;

    }

    else if (matrix.compare("S") == 0) {

      f77_dftbp_getrcsroverlap (_handler, nrow1, nzval, colind, rowpnt, val);

    }

    else std::cerr << "Invalid matrix selection in DftbpWrapper::getmatrix(...). Choose between H and S" << std::endl;

  }

  else {

    colind = new int[nzval];
    rowpnt = new int[nrow1];

    //Note: to val is associated a fortran complez vector, that's why allocation is double-sized.
    //val should be read as a [re1, im1, re2, im2.....] array
    val = new double[nzval*2];

    if (isreal == 1) std::cerr << "ERROR: requested a complex H or S while complex calculation is computed" << std::endl;
    assert(isreal==0);

    if (matrix.compare("H") == 0) {

      f77_dftbp_getzcsrhamiltonian (_handler, kPoint, nrow1, nzval, colind, rowpnt, val);

    }

    else if (matrix.compare("S") == 0) {

      f77_dftbp_getzcsroverlap (_handler, kPoint, nrow1, nzval, colind, rowpnt, val);

    }

    else std::cerr << "Invalid matrix selection in DftbpWrapper::getmatrix(...). Choose between H and S" << std::endl;

  }


  std::cout << "First value of H " << colind[0] << " " << rowpnt[0] << " " << val[0] << " " << std::endl;
  std::cout << "Second value of H " << colind[1] << " " << rowpnt[0] << " " << val[1] << " " << std::endl;
  std::cout << "Last value of H " << colind[nzval - 1] << " " << rowpnt[nrow1] << " " << val[nzval - 1] << " " << std::endl;
  std::cout << "Number of rows is " << nrow << std::endl;
}


void
DftbpWrapper::setexternalshift(int nAtom, double *pot)
{
  f77_dftbp_setexternalshift (_handler, nAtom, pot);
}


void DftbpWrapper::gethubbards(int nAtom, int max_shell, double* u_hubbard)
{
  f77_dftbp_gethubbards (_handler, nAtom, max_shell, u_hubbard);
}
