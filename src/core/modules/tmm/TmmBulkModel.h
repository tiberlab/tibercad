// $Id: TmmBulkModel.h 4391 2017-04-07 11:16:58Z pecchia $

#ifndef _TMMBULKMODEL_H_
#define _TMMBULKMODEL_H_

#include "PhysicalModel.h"
#include "TmmDipoleSource.h"
// #include "Tmm.h"

// class TmmDipoleSource

//! This is the base class for the TMM bulk physical model
class TmmBulkModel : public PhysicalModel
{

  public:


    //! Destructor
    virtual ~TmmBulkModel(void) {};

    //! Creator function
    static TmmBulkModel* create(const Material* mat, const ModelOptions& options);
	
    void  calculate(const libMesh::Elem* elem, const libMesh::Point& point, double lambda);


    //! Get the relative permittivity at a given wavelength
    libMesh::Complex get_permittivity(double lambda) const;


    //! Get the complex refractive index
    libMesh::Complex get_refractive_index(double lambda) const;

    double get_coherent_index(void) const;
	
	
    double get_emission_power(void) const;






  protected:

    //! Constructor
    TmmBulkModel(const ModelOptions& options);

    /*!
     * \brief Reads complex refractive index from database
     *
     * Optical refractive index is provided in an extra ASCII file, and in the
     * \c permittivity block of the database this file is indicated by
     * \c optical_data, which has to be a relative path (relative to the directory where
     * the main database file resides)
     */
    virtual void read_database(void) override;

    /*!
     * \brief Setup the optical constants
     *
     * The data file is read, and the internal data structure containing complex
     * refractive index (n, k) is set up.
     */
    virtual void do_init(void) override;

    virtual void prepare_submodels(void) override;


    virtual void do_print_info(void) override;

    /*!
     * \brief Interpolate (n, k) to a given wavelength
     *
     * If \c wavelength is outside the data range, the nearest value will be taken
     */
    std::pair<double, double> interpolate(double wavelength) const;

  private:

    //! The constructor method
    static TiberModelObject* _create(const ModelOptions& options, const void*);

    //! The destructor method
    static void _destroy(TiberModelObject* p);

    //! The file where we read optical data from
    std::string _datafile;

    //! The wavelengths
    std::vector<double> _wavelengths;

    //! n
    std::vector<double> _n_data;

    //! k
    std::vector<double> _k_data;

    //! An addressing array for faster access during interpolation
    std::vector<int> _addressing;

    double _incoherent_index {0};
	
	std::vector<TmmDipoleSource* > _DS;
	
	double _emission_power;



};




#endif // _TMMBULKMODEL_H_
