// $Id$

#ifndef _DATAIMPORTER_H_
#define _DATAIMPORTER_H_

#include<string>

#include "SimulationInterface.h"
#include "TensorGrid.h"



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
    virtual void do_init(void);

    //! Option parser
    virtual void parse_options(void);

    //! Setup available variables
    virtual void do_setup_solution_variables(void);

    //! Solver function
    virtual void do_solve(void);

    //! Print information
    virtual void do_print_info(void);

    //! Provide generation profile to outside world
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p);

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

    struct Options
    {
      std::string filename;
      std::string filetype;
      std::string variable_name;
      std::string unit;
      std::string variable_alias;
      std::string dataset_name;
      std::string num_dimensions;
      std::string sizes;
      std::string delimiter;
      std::string print_data;
      int int_num_dimensions;
      std::vector<int> int_sizes;
    };

    //! Sanity check of supplied options
    void _check_options(void);

    //! Information printout
    void _print_module_info(void);

    //! Wrapper for reading the right file type
    void _read_file(void);

    //! HDF Reader
    void _read_hdf5(void);
    
    //! Image Reader
    void _read_image(void);
    
    //! 1-dimensional CSV Reader
    void _read_csv1d(void);
    
    //! 2-dimensional CSV Reader
    void _read_csv2d(void);
    
    //! 3-dimensional CSV Reader
    //! Each line has the form:
    //! x;y;z;value
    void _read_csv3d(void);

    //! VTK Reader
    void _read_vtk(void);

    //! COMSOL Reader
    void _read_comsol(void);

    //! Index conversion from 2d to 1d
    int _at(int pos_x, int pos_y);

    //! Index conversion from 3d o 1d
    int _at(int pos_x, int pos_y, int pos_z);


    static const std::string valid_filetypes[];
    static const int num_valid_filetypes;

    int _dims, _size_x, _size_y, _size_z;

    std::string _filename, _filetype, _variable_name,
    _unit, _variable_alias, _dataset_name, _num_dimensions, _sizes, _print_data;

    //! Delimiters
    std::string _delimiter;

    //! Characters interpreted as comments
    std::set<std::string> _comment_chars;

    double* _data;
    TensorGrid* _tensorgrid;
    Point _origin,_bound;

    Options myopts;
};



#endif
