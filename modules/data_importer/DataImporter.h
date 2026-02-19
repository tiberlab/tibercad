/*  
 * This file is part of the tiberCAD module data_import.
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
 * \file DataImporter.h
 * \brief tiberCAD data_import module header.
 *
 * \note This file is part of module data_import.
 */



/*!
 * \file DataImporter.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */


#ifndef TC_DATAIMPORTER_H
#define TC_DATAIMPORTER_H

#include <string>

#include "tibercad/module/SimulationInterface.h"

class TensorGrid;

/*!
 *
 * \brief Module to read 1d, 2d or 3d data from files
 * 
 */
class TC_DLLOCAL DataImporter : public SimulationInterface
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

    //! Setup the mesh, read from file
    virtual void setup_mesh(void) override;

    //! Solver function
    virtual void do_solve(void) override;

    //! Print information
    virtual void do_print_info(void) override;

    //! Provide generation profile to outside world
    virtual void get_solution_secure(const Elem* elem,
        std::map<ID, std::vector<double> >& values,
        const std::vector<Point>& p) override;

  private:


    //! Known file types
    enum Filetype
    {
      csv,      //! comma or space separated text file
      comsol,   //! comsol generated data file
      unknown
    };

    //! Assembly function
    //static void assemble(EquationSystems& es, const std::string& system_name);

    //! Real assembly function
    //void do_assemble(EquationSystems& es, const std::string& system_name);

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


    //! Spatial dimensions
    unsigned int _dims;

    //! The file to read from
    std::string _filename;

    //! The filetype
    //Filetype _filetype;

    //! The file type
    std::string _filetype;

    //! Delimiters
    std::string _delimiter;

    //! Characters interpreted as comments
    std::set<char> _comment_chars;

    TensorGrid* _tensorgrid;

    //! coordinate translation
    Point _translate;

};



#endif
