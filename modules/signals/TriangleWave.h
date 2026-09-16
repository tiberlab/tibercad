/*  
 * This file is part of the tiberCAD module signals.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file TriangleWave.h
 * \brief tiberCAD signals module header.
 *
 * \note This file is part of module signals.
 */


#include "tibercad/math/SignalGenerator.h"

/*!
 * \brief Triangle wave
 *
 * Represent a triangle wave based on its Fourier series,
 * up to a given harmonic.
 */
class TriangleWave : public SignalGenerator
{

  public:

    //! Destructor
    virtual ~TriangleWave(void);


  protected:

    //! Constructor
    explicit TriangleWave(const ModelOptions& options);

    //! Initialize
    void do_init(void) override;

    //! Recalculate dependent variables
    void do_update_dependent_variables(void) override;


  private:

    //! The output value
    double _output;

    //! The variable name of the output variable
    std::string _output_str;

    //! The periodicity
    double _period;

    //! The mean value
    double _mean;

    //! The amplitude
    double _amplitude;

    //! The highest harmonic
    unsigned int _max_harmonic;

};
