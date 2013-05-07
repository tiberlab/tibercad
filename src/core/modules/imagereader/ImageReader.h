// $Id$

#ifndef _IMAGEREADER_H_
#define _IMAGEREADER_H_

#include "SimulationInterface.h"

/*!
 * 
 * \brief Extract data from image files
 *
 *
 * The scope of this module is to import an image file,
 * like a TEM image or similar, to extract the data and to
 * make it accessible to other modules
 *
 */
class TBDLLOCAL ImageReader : public SimulationInterface
{

  public:

    //! Destructor
    ~ImageReader(void);

    //! We need a public static creator function
    static ImageReader* create(const ModelOptions& options);



  protected:

    //! The initialization
    virtual void do_init(void);


    //! Parse the options from the input file
    virtual void parse_options(void);


    //! Setup the available variables
    virtual void do_setup_solution_variables(void);


    //! Solve the MyPoisson equation
    virtual void do_solve(void);


    //! Print some useful information
    //virtual void do_print_info(void);


    //! We need to create a physical model
    //virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
    //    const Material* mat) const;

    //! We need to create boundary condition model
    //PhysicalModel* create_boundary_model(const ModelOptions& options,
    //    const MaterialBoundary* boundary) const;

    //! We have to provide somehow our solution variables
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);



  private:

    //! These are the known solution variables
    enum Solutions
    {
      Data        /*!< the imported data */
    };

    
    ImageReader(const ModelOptions& options);


    //! Do the actual read
    void _import_picture(void);

    //! The image file
    std::string _file;

    //! The origin (lower left corner) of the image
    Point _origin;

    //! The pixel size
    double _pix_size;

    //! The data vector
    std::vector<double> _data;

};





#endif // _IMAGEREADER_H_
