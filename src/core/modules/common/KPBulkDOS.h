// $Id: KPBulkDOS.h 2117 2010-10-20 20:56:13Z maufder $

#ifndef _KPBULKDOS_H_
#define _KPBULKDOS_H_


#include "DensityOfStates.h"
#include "StrainInterface.h"
#include "HashMap.h"

class DDsemiconductor;

/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT KPBulkDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~KPBulkDOS(void) {};


    //! Creator function
    static KPBulkDOS* create(const ModelOptions& options);



  protected:

    //! Constructor
    KPBulkDOS(const ModelOptions& options);

    //! Prepare the submodel
    virtual void prepare_submodels(void);

    //! Initialize
    virtual void do_init(void);

    //! Reinitialize for elem
    virtual void do_reinit(const Elem* elem);

    //! Print some info
    virtual void do_print_info(void);

    //! Get occupied states and derivative w.r.t. Ef
    virtual std::pair<double, double>
    calculate_density_and_derivative(double Ef, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;

    //overloading for Trap.C
    virtual std::pair<double, double>
    calculate_density_and_derivative(double Ef, double Epot, double kT, double kTlattice) const;

  private:

    //! A class for the cache
    class Cache
    {
      public:

      struct Data
      {
        int strain_seq_num = -1;
        double ref_energy;
        double dos_mass;
        double eff_dos;
      };

      private:

      HashMap<const Elem*, Data>::Type _cache;

      public:

      bool get_data(const Elem* elem, int strain_seq, Data& data)
      {
        bool valid = false;
        HashMap<const Elem*, Data>::Type::iterator it(_cache.find(elem));
        if ((it != _cache.end()) && ((it->second).strain_seq_num >= strain_seq))
        {
          data = it->second;
          valid = true;
        }

        return valid;
      }

      void set_data(const Elem* elem, const Data& data)
      {
        _cache[elem] = data;
      }
    };

    //! The underlying model
    DDsemiconductor* _bulk_model;

    //! The interface to the strain simulation
    StrainInterface _strain_if;

    //! We allow for several subbands
    std::vector<double> _ref_energies;

    //! The effective masses
    std::vector<double> _dos_mass;

    //! The degeneracies
    std::vector<int> _degeneracy;

    //! If recompute is true, for each element we need to reinit all data
    bool _recompute;

    //! Data cache for element data
    Cache _cache;

    //! The DOS factor
    /*!
     * In 3D the effective DOS is
     * \f[ N_{3D} = \f]
     */
    double _dos_factor;


    std::vector<double>& band_edge(void) { return _ref_energies; }
    std::vector<double>& dos_mass(void) { return _dos_mass; }
    std::vector<int>& degeneracy(void) { return _degeneracy; }

    //! Do the kp calculation
    void _solve_kp(void);


};

//
// inline methods
//

inline
KPBulkDOS*
KPBulkDOS::create(const ModelOptions& options)
{
  return new KPBulkDOS(options);
}


#endif // _BULKDOS_H_
