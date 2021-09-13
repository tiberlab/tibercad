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
    


    //! function to calculate light's propagation angle in all layers
    /*!
     * This function calculate propagating angle of the light using snell's law in all layers
     * "n_real" is a Refractive index vector
     * "incident_angle" is angle of light in the first layer
     * normal incident is equal to incident_angle = 0
     */
    vector<double> theta_cal(vector<double> n_real , double incident_angle);



    //! defining a class to work with 2*2 matrics 
    class sqr {
    	public:
    	
        sqr();
        sqr(double a00, double a01, double a10, double a11);
        
        //! defining operator  for two matrics product
        sqr operator*(sqr);
        
        //! definig a function to print matrix elements
        void print();
        
        //! defining a function to return M matrix
        /*!
         * "n_real" is real part of refractive index
         * "n_imag" is imaginary part of refractive index
         * "length" is the lenth of the layer
         * "lambda" is the light's wavelength
         * "theta" is the light's traveling angle(normal incident is equal to '0')
         */           
        void get_M(double n_real,double n_imag,double length,double lambda, double theta);
        
        
        //! defining a function to return D matrix
        /*!
         * "n1_real" is first layer real part of the refractive index
         * "n1_imag" is first layer imaginary part of the refractive index
         * "n2_real" is second layer real part of the refractive index
         * "n2_imag" is second layer imaginary part of the refractive index
         * "theta_layer1" is the first layer light's traveling angle(normal incident is equal to '0')
         * "theta_layer2" is the second layer light's traveling angle(normal incident is equal to '0')
         */  
        void get_D(double n1_real,double n1_imag,double n2_real,double n2_imag,double theta_layer1, double theta_layer2);
       
        //! variables corresponding to elements of the matrix
        std::complex<double> m00;
        std::complex<double> m01;
        std::complex<double> m10;
        std::complex<double> m11;
        
    };


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
