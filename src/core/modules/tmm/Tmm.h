// $Id$

#ifndef _TMM_H_
#define _TMM_H_

#include "SimulationInterface.h"

using namespace std;



/*!
 * \brief Implementation of Transfer Matrix Method for electromagnetic fields
 *
 * This class implements standard Transfer Matrix Method (TMM) to calculate the
 * electromagnetic field distribution in a 1D layered structure.
 *
 * Author:
 * Contributors:
 */

class TBDLLOCAL Tmm : public SimulationInterface
{

  public:

    //! Destructor
    /*!
     * We do not declare it virtual here, as we will not allow
     * to derive from this class anyway.
     */
    ~Tmm(void);

    //! We need a public static creator function
    static Tmm* create(const ModelOptions& options);





  protected:

    //***********************function to print a matrix**********************************
    void show_matrix(vector<vector<complex<double>>> matrix);

    //***This function print the matrix
    //***"matrix" is input matrix

    //********************function to calculate angle in all layers*************************************
    vector<double> theta_cal(vector<double> n_real , double incident_angle);

    //***This function calculate propagating angle of the light using snell's law in all layers
    //***"n_real" is a Refractive index vector
    //***"incident_angle" is angle of light in the first layer


    //**************************function to calculate M matrix********************************************
    vector<vector<complex<double>>> get_M(double n_real,double n_imag,double length,double lambda, double theta);

    //***This function calculate the M matrix (propagation matrix)
    //***"n_real" and "n_imag" are refractive index of the layer
    //***"length" is the length of the layer
    //***"lambda" is the wavelength of light
    //***"theta" is the light's angle


    //******************************function to calculate D matrix******************************************
    vector<vector<complex<double>>> get_D(double n1_real,double n1_imag,double n2_real,double n2_imag,double theta_layer1, double theta_layer2);

    //***This function calculate the D or I matrix (transition matrix)
    //***"n1_real" and "n1_imag" are refractive index of the first layer
    //***"n2_real" and "n2_imag" are refractive index of the second layer
    //***"theta_layer1" is the light's angle in first layer
    //***"theta_layer2" is the light's angle in second layer


    //*******************************function to calculate matrixs product*******************************
    vector<vector<complex<double>>> matrix_product(vector<vector<complex<double>>> A,vector<vector<complex<double>>> B);

    //***This function calculate multiplication of two inputs
    //***"A" is first matrix
    //***"A" is second matrix

    //! The initialization
    virtual void do_init(void);


    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);


    //! Solve the Poisson equation
    virtual void do_solve(void);


    //! Print some useful information
    virtual void do_print_info(void);


    //! We need to create a physical model
    virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
        const Material* mat) const;

    //! We need to create boundary condition model
    virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
        const MaterialBoundary* boundary) const;


    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);



  private:


    /*!
     * \brief The known solution variables
     */
    enum Solutions
    {
      EField,
      HField,
      Intensity,
      GenerationRate
    };
    
    /*!
     * \brief Constructor
     *
     * Being private disables further inheritance.
     */
    Tmm(const ModelOptions& options);

    /*!
     * \brief The wavelengths
     */
    std::vector<double> _wavelengths;

    /*!
     * \brief The incident angle
     */
    std::vector<double> _incident_angle;





};


#endif // _TMM_H_
