// $Id$

#ifndef _PVMODULEMODEL_H_
#define _PVMODULEMODEL_H_

#include "PhysicalModel.h"



class Elem;
class ElementaryCell;


//! This is the basic physical modle of PVModule
class PVModuleModel : public PhysicalModel
{

  public:

    //! The type of region
    enum RegionType
    {
      ACTIVE, /*!< active region */
      P1,     /*!< P1, bottom contact layer interrupted */
      P2,     /*!< P2, top-to-bottom contact */
      P3      /*!< P3, isolation from P2 */
    };

    //! Destructor
    ~PVModuleModel(void) {};

    //! Creator function
    static PVModuleModel* create(const Material* mat, const ModelOptions& options);

    //! Get region type
    RegionType get_region_type(void) const;

    //! Get the top and bottom sheet resistances for the given point
    /*!
     * The top sheet resistivity is given in the first argument, the bottom one
     * in the second.
     */ 
    std::pair<double, double> get_sheet_resistances(const Elem* elem,
                                                    const Point& point) const;

    //! Get the vertical connection resistance
    double get_connection_resistance(const Elem* elem,
                                     const Point& point) const;
    
    //! Get the elementary cell model
    /*!
     * If the area is zero, an empty representation is returned,
     * otherwise the elementary cell scaled to the given area.
     * 
     * TODO need to define what to return, and what else is needed.
     * For example, we might have to pass top and bottom node
     * indices, and next free node index.
     */
    ElementaryCell* get_elementary_cell(void) const;


  protected:

    //! Constructor
    PVModuleModel(const ModelOptions& options);

    virtual void do_init(void) final;

    virtual void prepare_submodels(void) final;

    virtual void do_print_info(void) final;


  private:

    //! The region type
    RegionType _region_type {ACTIVE};

    //! The top sheet resistivity
    /*!
     * Units are Ohm/square
     */
    double _top_rsheet = 0.001;

    //! The bottom sheet resistivity
    /*!
     * Units are Ohm/square
     */
    double _bottom_rsheet = 0.001;

    //! The top-to-bottom connection resistivity
    /*!
     * Units are Ohm*cm^2, corresponding to the conductivity of
     * the connection metal times the thickness of the layer stack.
     */
    double _connection_res = 0.001;

    //! The elementary cell model
    ElementaryCell* _elementary_cell = nullptr;

    //! The constructor method
    static TiberModelObject* _create(const ModelOptions& options, const void*);

    //! The destructor method
    static void _destroy(TiberModelObject* p);

};



inline
PVModuleModel::RegionType
PVModuleModel::get_region_type(void) const
{
  return(_region_type);
}


inline
ElementaryCell*
PVModuleModel::get_elementary_cell(void) const
{
  return _elementary_cell;
}

#endif // _PVMODULEMODEL_H_
