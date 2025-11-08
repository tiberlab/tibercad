// $Id: DriftDiffusionProperties.h 4145 2015-10-02 11:53:20Z maufder $

#ifndef _MASTEREQUATIONSPROPERTIES_H_
#define _MASTEREQUATIONSPROPERTIES_H_



#include "PhysicalModel.h"

#include "TemperatureInterface.h"
#include "MEBandProperties.h"
#include "SimulationOptions.h"
#include "MasterEquationsDefs.h"
#include "TiberCad.h"
#include "Constants.h"
#include "TypeDefs.h"

#include "tensor.h"

#include "tensor_value.h"
#include "vector_value.h"
#include "point.h"


#include <vector>
#include <set>
#include <map>


// forward declarations
class Elem;
class Point;
class SimulationInterface;
class MasterEquations;


class MasterEquationsProperties : public PhysicalModel //TBDLEXPORT
{

  public:


    //! A default (empty) destructor
    virtual ~MasterEquationsProperties(void);


    //! Create a named master-equations model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static MasterEquationsProperties* create(const std::string& name,
        const Material* mat, const ModelOptions& options = ModelOptions());


    //! Set the coupling type
    void set_coupling_type(MasterEquationsDefs::Coupling coupling)
      { _coupling = (int) coupling; };


    //! Set the coupling type
    void set_coupling_type(int coupling)
      { _coupling = coupling; };

    //! Get the coupling type
    MasterEquationsDefs::Coupling get_coupling_type(void) const
      { return (MasterEquationsDefs::Coupling) _coupling; };



    //! Set the statistics to be used
    /*!
     * \param statistics the statistics
     */
    void set_statistics(TiberCad::Statistics statistics);

    //! Get the statistics to be used
    /*!
     * \return the statistics
     */
    TiberCad::Statistics get_statistics(void) const;


    //! (Re-)Initialize for the given element
    /*!
     * \c reinit() calls \c prepare_element_data() which needs to be
     * implemented in derived classes
     *
     * \param elem the current element
     * \param dd_prop a pointer to the semiconductor model
     */
    void reinit(const Elem* elem);


    //! Set the coordinates
    void set_coordinates(const Point& pt);


    //! Set the carrier temperatures
    /*!
     * The lattice, electron and hole temperatures have to be set before
     * calling any of the \c calculate_xxx methods
     *
     * \param T_lat the lattice temperature
     * \param T_e the electron temperature
     * \param T_h the hole temperature
     */
    void set_carrier_temperatures(double T_e, double T_h);



    //! Set the carrier densities
    /*!
     * This method can be used if the densities are the independent variables.
     * The densities have to be set before calling any of the
     * \c calculate_xxx methods.
     *
     * \param n the electron density
     * \param p the hole density
     */
    void set_densities(double n, double p);



    //! Get the conduction band properties as pointer
    //void set_conduction_band(BandProperties* cb);


    //! Get the valence band properties as pointer
    //void set_valence_band(BandProperties* vb);


    //! Get the element we are currently working on
    const Elem* get_element(void) const;

    //! Get the coordinates of the point we are currently working on
    const Point& get_coordinates(void) const;


    //! Set the lattice temperature (in K)
    void set_lattice_temperature(double T);


    //! Get the lattice temperature (in units of eV)
    double get_lattice_temperature(void) const;



    //! Calculate electron and hole densities and derivatives
    /*!
     * This method calculates electron and hole densities and their
     * derivatives with respect to the electric and the
     * electro-chemical potentials.
     */
    void calculate_densities(void);



    //! Get the electron density
    /*!
     * Get the electron density as calculated by \c calculate_all(...)
     *
     * \return the electron density
     */
    double get_electron_density(void) const
      { return electron_density; };


    //! Get the hole density
    /*!
     * Get the hole density as calculated by \c calculate_all(...)
     *
     * \return the hole density
     */
    double get_hole_density(void) const
      { return hole_density; };


    //! Get the lowest conduction band edge
    double get_conduction_band_edge(void) const
      { return _conduction_band_edge->get_band_edge(); };


    //! Get the highest valence band edge
    double get_valence_band_edge(void) const
      { return _valence_band_edge->get_band_edge(); };


    //! Get the conduction band properties
    MEBandProperties& get_conduction_band(void) const //??? cost or not const
      { return *_conduction_band; };


    //! Get the valence band properties
    MEBandProperties& get_valence_band(void) const //??? cost or not const
      { return *_valence_band; };


    //! Get the band gap
    double get_band_gap(void) const
      { return _conduction_band->get_band_edge() - _valence_band->get_band_edge(); };


    //! Get the square of the intrinsic density
    double get_intrinsic_density_squared(void) const
      { return _intrinsic_density * _intrinsic_density; };


    //! Get the intrinsic density
    double get_intrinsic_density(void) const
      { return _intrinsic_density; };


    //! Get the equilibrium electron density
    double get_equilibrium_electron_density(void) const
      { return _equilibrium_n; };


    //! Get the equilibrium hole density
    double get_equilibrium_hole_density(void) const
      { return _equilibrium_p; };


    //! Get equilibrium fermi level
    double get_equilibrium_fermi_level(void) const
      { return _equilibrium_fermi_level; };



  protected:



    //! The empty constructor.
    MasterEquationsProperties(const ModelOptions& options);


    //! Initialize this model
    /*!
     * This reads the database and calls init for all submodels
     * A derived class which reimplements this method has to call
     * explicitly the one of this class!
     */
    virtual void do_init(void);


    //! Create some of the submodels
    //virtual void prepare_submodels(void);


    //! Set the element we are currently working on
    void set_element(const Elem* elem);


    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     * This method can be overriden by derived classes.
     */
    virtual void prepare_element_data(void) {};



    //! Get the temperature interface
    TemperatureInterface& get_temperature_interface(void);


    //! Get the conduction band properties as pointer
    void set_conduction_band(MEBandProperties* cb);


    //! Get the valence band properties as pointer
    void set_valence_band(MEBandProperties* vb);


    //! Set the intrinsic density
    void set_intrinsic_density(double ni)
      { _intrinsic_density = ni; };


    //! Set the equilibrium electron density
    void set_equilibrium_n(double n)
      { _equilibrium_n = n; };


    //! Set the equilibrium hole density
    void set_equilibrium_p(double p)
      { _equilibrium_p = p; };


    //! Set the equilibrium Fermi level
    void set_equilibrium_fermi_level(double Ef)
      { _equilibrium_fermi_level = Ef; }


    //! The lattice temperature in eV (\f$= k_B T_{lat} / e\f$)
    double lattice_vt;


    //! The electron density
    double electron_density;


    //! The hole density
    double hole_density;


    //! The conduction band edge
    double conduction_edge;


    //! The valence band edge
    double valence_edge;


    //! The electron temperature in eV (\f$= k_B T_e / e\f$)
    double electron_vt;

    //! The hole temperature in eV (\f$= k_B T_h / e\f$)
    double hole_vt;



  private:

    //! The copy constructor is disabled
    MasterEquationsProperties(const MasterEquationsProperties& rhs);

    //! The assignment operator is disabled
    MasterEquationsProperties& operator=(const MasterEquationsProperties& rhs);

    //! The interface to the lattice temperature simulation
    TemperatureInterface _lattice_temp;

    //! The element we are currently working on
    const Elem* _elem;

    //! The coordinates of the point we are working on
    const Point* _coord;

    //! The conduction band properties
    /*!
     * Band properties are assumed to be elemental data, \em not nodal data
     */
    MEBandProperties* _conduction_band;

    //! The conduction band properties
    /*!
     * Band properties are assumed to be elemental data, \em not nodal data
     */
    MEBandProperties* _valence_band;


    MEBandProperties* _valence_band_edge;


    MEBandProperties* _conduction_band_edge;

    //! The lattice temperature in eV (\f$= k_B T_{lat} / e\f$)
    double _lattice_vt;

    //! The statistics to be used
    TiberCad::Statistics _statistics;


    int _coupling;


    double _equilibrium_fermi_level;

    //! The intrinsic density
    double _intrinsic_density;

    //! The equilibrium electron density
    double _equilibrium_n;

    //! The equilibrium electron density
    double _equilibrium_p;


};


//
// inline members
//



inline
void
MasterEquationsProperties::set_statistics(TiberCad::Statistics statistics)
{
  _statistics = statistics;
}


inline
TiberCad::Statistics
MasterEquationsProperties::get_statistics(void) const
{
  return _statistics;
}



inline
void
MasterEquationsProperties::set_coordinates(const Point& pt)
{
  _coord = &pt;
}



inline
void
MasterEquationsProperties::set_densities(double n, double p)
{
  electron_density = n;
  hole_density = p;
}



inline
void
MasterEquationsProperties::set_element(const Elem* elem)
{
  _elem = elem;
}



inline
const Elem*
MasterEquationsProperties::get_element(void) const
{
  return _elem;
}



inline
const Point&
MasterEquationsProperties::get_coordinates(void) const
{
  return *_coord;
}


inline
void
MasterEquationsProperties::set_lattice_temperature(double T)
{
  _lattice_vt = T * Constants::k_B;
}


inline
double
MasterEquationsProperties::get_lattice_temperature(void) const
{
  return _lattice_vt;
}




inline
TemperatureInterface&
MasterEquationsProperties::get_temperature_interface(void)
{
  return _lattice_temp;
}


inline
void
MasterEquationsProperties::set_carrier_temperatures(double T_e, double T_h)
{
  electron_vt = T_e;
  hole_vt = T_h;
}




#endif /* _MASTEREQUATIONSPROPERTIES_H_ */
