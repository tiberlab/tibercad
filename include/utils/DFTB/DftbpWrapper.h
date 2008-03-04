#ifndef _DFTBPWRAPPER_H_
#define _DFTBPWRAPPER_H_

#include "dftbp.h"
#include <iostream>
#include <assert.h>

//-----------------------------------------------------------------------

class DftbpWrapper
{

public:

  //!Wrapper class for callings to DFTB+ library (libdftbp.so)



  //!Constructor
  /*!Assign an handler to DFTB+ instance, transparent to programmer
   *
   */
  DftbpWrapper();


  //! Destructor  
  ~DftbpWrapper();


  //!Static method to create a Dftb wrapper instance
  static DftbpWrapper* create(); 


  //!Function to fill Dftb parameters:
  /*!
   * \param nAtoms (in) total number of atoms
   * \param nType (in) total number of atomic species
   * \param eTemp (in) electronic temperature
   * \param iPeriodic (in) 1 if structure is periodic, 0 if not
   * \param speciesNames[nType] (in) char array containing atomic species names (sized nType * DFTBP_MC * sizeof(char)),
   * where DFTBP_MC is defined in dftbp.h
   * \param species[nAtom] (in) integer array containing the specie index for any atom
   */
  void fill_param (int nAtoms, int nType, double eTemp, int iPeriodic,
                   char *speciesNames, int *species);


  //!Pass the path of sk parameters files to Dftb
  /*!
   * \param skNames[ ] (in) sk files path (sized nType * DFTBP_LC * sizeof(char)),
   * where DFTBP_MC is defined in dftbp.h
   * \param mAngs[nType] (in) maximum angular momentum for each specie
   * \param orbResolved (in) tells if Hubbard U for all shells should be read
   * \param skInterp (in) interpolation method for SK tabel points (1 or 2, see DFTB+ documentation for further info)
   * \param nType (in) total number of atomic species
   */
  void addskdata (char *skNames, int *mAngs, int orbResolved, int skInterp,
                  int nType);


  //!Initialize the Dftb instance
  void initdftb ();


  //!Update atom positions
  /*!
   * \param nAtom (in) total number of atoms
   * \param newCoords[nAtom * 3] (in) coordinates (x1, y1, z1, x2, y2, z2, ...)
   */
  void up_coords (int nAtom, double *newCoords);


  //!Get system energy
  /*!
   * \param energy (out) value of energy
   */
  void get_energy (double &energy);


  //!Add periodicity vectors 
  /*!
   * \param latVecs[9] (in) periodicity vectors (x1, y1, z1, x2, y2, x3, y3)
   */
  void addlattice (double *latVecs);


  //!Add explicitly k points for k space integration
  /*!
   * \param nKpoint (in) total number of k points
   * \param kPoints[nKpoint * 3] (in) array containing k points (k1x, k1y, k1z, k2x, k2y,...)
   * \param kWeights[nKpoint] (in) array containing k points weights 
   */
  void addkpoints (int nKPoint, double *kPoints, double *kWeights);


  //!Specify a supercell sampling for k points calculations (please refer to Dftb+ documentation)
  void addsupersampling (double *coeffs, double *shifts, int noinv);

  
  //!Get electronic charge for each atom
  /*!
   * \param nAtom (in) total number of atoms
   * \param charges[nAtom] (out) charges per atom
   */
  void getchargesperatom (int nAtom, double* charges);


  //!Get Real/Complex Hamiltonian/Overlap in CSR sparse format
  /*!
   * \param nrow (out) number of rows
   * \param ncol (out) number of columns
   * \param nzval (out) number of non zero values
   * \param isreal (out) 0 if matrix is complex, 1 if matrix is real
   * \param colind[nzval] (out) column index (starting from 1)
   * \param rowpnt[nrow] (out) pointer to rows
   * \param val[nzval if real, 2*nzval if complex] (out) non zero values (if complex, (re1, im1, re2, im2,...))
   * \param matrix (in) choose which matrix is needed (H or S)
   * \param kPoint[3] (in) if not specified, a real matrix in Gamma point is expected. If specified, H or S for that K point is given  
   * WARNING: colind, rowpnt and val are internally allocated
   */
  void getmatrix(int nrow, int ncol, int nzval, int isreal, int *colind, int *rowpnt, double *val, std::string matrix, double *kPoint = NULL);


private:
  int _handler[DFTBP_HSIZE];

};




#endif
