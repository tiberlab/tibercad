// $Id$

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

    //! Create an instance
    static TriangleWave* create(const ModelOptions& options);


  protected:

    //! Constructor
    TriangleWave(const ModelOptions& options);

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
