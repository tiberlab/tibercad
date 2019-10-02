// $Id$

#ifndef _DATAIMPORTER_H_
#define _DATAIMPORTER_H_

#include <string>

#include "SimulationInterface.h"

class TensorGrid;

/*!
 *
 * \brief Module to read 1d, 2d or 3d data from files
 * 
 */
class TBDLLOCAL DataImporter : public SimulationInterface
{
  public:
    //! Destructor
    virtual ~DataImporter(void);

    //! Creator function as in exmaple Poisson module
    static DataImporter* create(const ModelOptions& options);

    //! Constructor
    DataImporter(const ModelOptions& options);


  protected:

    //! Initialisation
    virtual void do_init(void) override;

    //! Option parser
    void parse_options(void);

    //! Setup available variables
    virtual void do_setup_solution_variables(void) override;

    //! Solver function
    virtual void do_solve(void) override;

    //! Print information
    virtual void do_print_info(void) override;

    //! Provide generation profile to outside world
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p) override;

  private:

    //! Known solution variables
    enum Solutions
    {   
      Data //! Generation rate from file
    };

    //! Assembly function
    //static void assemble(EquationSystems& es, const std::string& system_name);

    //! Real assembly function
    //void do_assemble(EquationSystems& es, const std::string& system_name);

    //! Static pointer to this
    static DataImporter* _this;

    //! Wrapper for reading the right file type
    void _read_file(void);

    //! HDF Reader
    void _read_hdf5(void);
    
    //! Image Reader
    void _read_image(void);
    
    //! n-dimensional CSV Reader
    void _read_csv(void);
    
    //! VTK Reader
    void _read_vtk(void);

    //! COMSOL Reader
    void _read_comsol(void);

    //! Create an unstructured grid from a point set
    void _create_mesh_from_points(const std::vector<double>* x = nullptr,
                                  const std::vector<double>* y = nullptr,
                                  const std::vector<double>* z = nullptr);

    static const std::string valid_filetypes[];
    static const int num_valid_filetypes;

    unsigned int _dims;
    size_t _size_x, _size_y, _size_z;

    //! The file to read from
    std::string _filename;

    //! The file type
    std::string _filetype, _variable_name,
    _unit, _variable_alias, _dataset_name, _num_dimensions, _sizes, _print_data;

    //! Delimiters
    std::string _delimiter;

    //! Characters interpreted as comments
    std::set<char> _comment_chars;

    double* _data;
    TensorGrid* _tensorgrid;
    Point _origin,_bound;

};



#endif
