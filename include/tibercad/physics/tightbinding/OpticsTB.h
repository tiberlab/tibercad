/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file OpticsTB.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef _OPTICSTB_H_
#define _OPTICSTB_H_

#include "tibercad/physics/optics/Optics.h"
#include "tibercad/physics/schroedinger/EigenvalueProblem.h"

//! A base class of optics calculation
/*!
 * The task of this class is to calculate the spectrum at a certain given
 * k-point
 */
class OpticsTB : public Optics
{

  public:

    //! The constructor
    OpticsTB(const ModelOptions& options);

    //! The destructor
    virtual ~OpticsTB(void);

 
    static OpticsTB* create(const ModelOptions& options);


  protected:

    virtual void do_init(void);

    virtual void do_assemble(const ModelOptions& options);

    //! Assemble the P-matrix and compute its matrix elements.
    virtual void do_compute_matrix_elements(void);

    
    virtual void calculate_matrix_bulk(void){};


  private:


    //! checks that states for optics are really there
    void check_states(void);

};


inline OpticsTB* OpticsTB::create(const ModelOptions& options)
{
  return (new OpticsTB(options));
}


#endif // _OPTICSTB_H_
