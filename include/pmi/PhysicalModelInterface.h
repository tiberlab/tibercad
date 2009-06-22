// $Id$

#ifndef _PHYSICALMODELINTERFACE_H_
#define _PHYSICALMODELINTERFACE_H_

#include "tiber_config.h"
#include "TiberModelObject.h"
#include "TypeDefs.h"
#include "Database.h"
#include "ModelOptions.h"
#include "InitFailedException.h"

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#include "tensor.h"
#include "xtensor.h"

#include "TiberModule.h"


class Material;


//! Base class for the different physical models
/*!
 * This is the base class for all implementations of any kind of models
 * used in simulations.
 */
class PhysicalModelInterface : public TiberModelObject
{

  public:

    //! Destructor
    virtual ~PhysicalModelInterface(void);


    //! Get the unique ID of this model
    ID get_id(void) const;


    //! Get the ID of the simulator this model is used with
    /*!
     * This ID is needed by the simulator to get the right models
     */
    ID get_simulator_id(void) const;


    //! Set the ID of the simulator this model is used with
    /*!
     * This ID is needed by the simulator to get the right models
     */
    void set_simulator_id(ID id);


    //! Get the ID of the model with name \c model_name
    /*!
     * If the model was not already registered, it gets inserted in the
     * model list and assigned a new ID.
     *
     * \param model_name the name of the model
     * \return the ID of the model or 0 if it does not exist
     *
     * The template parameter identifies the model family which the
     * model should belong to.
     */
    template <typename T>
    static ID get_id_from_name(const std::string& model_name);

    //! Creates a new named model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param options the options as given in the input file
     * \return a pointer to the newly created object
     */
    static PhysicalModelInterface* create(const std::string& name,
        const ModelOptions& options = ModelOptions());


    //! Creates a new model from a given creator function
    /*!
     * Use this method only in special cases, e.g. if you don't have single
     * libraries for the different models.
     */
    static PhysicalModelInterface* create(create_t create_fnc,
        destroy_t destroy_fnc, const ModelOptions& options = ModelOptions());


    //! Create a new model as a copy of this
    /*!
     * \return a pointer to the newly created model, \c NULL if
     * creation failed
     *
     * \note copy will usually \em not copy any members that can be set up
     * during do_init()
     */
    PhysicalModelInterface* copy(void) const;


    //! Create a new submodel model as a copy of another one
    /*!
     *\return a pointer to the newly created model, \c NULL if
     * creation failed
     * \param other the model to be copied
     *
     * The new model will be associated to the same material as this
     */
    template <typename T>
    T* create_submodel_copy(const T* other) const;


    //! Create a new submodel as an alloy model
    template <typename T>
    T* create_submodel_alloy(const T* comp_A, const T* comp_B, double xa);


    //! Get a reference to the database
    /*!
     * the database will already be setup for the material this model is
     * associated to
     */
    Database& get_database(void);

    //! Set a reference to the material this model belongs to
    void set_material(Material* material);


    //! Get a reference to the material this model belongs to
    const Material* get_material(void) const;


    //! Get a writeable reference to the material this model belongs to
    Material* get_material(void);


    //! Get the user defined name of this model
    const std::string& get_name(void) const;


    //! Get the type of this model
    /*!
     * The type is the identifying string which defines at creation time
     * which simulation to create. It's the same string one writes in the
     * input file.
     */
    const std::string& get_type(void) const;


    //! Get the default name for this model
    std::string get_default_name(void) const;


    //! Initialize this model
    /*!
     * It calls read_database() and then do_init()
     */
    void init(void);


    //! Initialize this model as an alloy with two components
    /*!
     * It calls read_database_alloy() and then do_init_alloy()
     */
    void init_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);


    //! Print some info
    void print_info(void);



  protected:

    //! Empty constructor
    PhysicalModelInterface(void);


    //! Set the name of a model
    /*!
     * Use with caution as it could break standard behaviour!
     */
    void set_name(const std::string& name);


    //! Print some info
    /*!
     * The implementation should add 4 spaces at the beginning of
     * each line.
     */
    virtual void do_print_info(void){};


    //! Initialize the model
    /*!
     * This method should set all model options and call
     * \c init() of any associated model.
     *
     * It will be called \em only for non-alloy models!
     *
     * This method should be implemented in derived classes.
     */
    virtual void do_init(void) {};


    //! Create a new model of the same type
    virtual PhysicalModelInterface* create_new(void) const = 0;


    //! Copy data from another model to this one
    /*!
     * This method should copy \em class members and data structures from
     * rhs which cannot be setup at initialization from scratch,
     * such as members that were set up at creation time.
     *
     * If you reimplement this in a derived class, call the method
     * of the base class, too. If not, you may not copy important
     * things.
     */
    virtual void copy_from(const PhysicalModelInterface* rhs) {};


    //! Read the properties from the database
    /*!
     * Reads all needed physical properties from the database for \em non alloy
     * models.
     * (The default behaviour is to do nothing at all.)
     * If you reimplement this in a derived class, call the method
     * of the base class, too.
     */
    virtual void read_database(void) {};


    //! Read alloy parameters from the database
    /*!
     * Read all parameters for an alloy model from the database.
     * (The default behaviour is to do nothing at all.)
     * If you reimplement this in a derived class, call the method
     * of the base class, too.
     */
    virtual void read_database_alloy(void) {};


    //! Initialize an alloy
    /*!
     * Calculates all parameters of an alloy \f$A_xB_{x-1}C\f$.
     *
     * You have to do everything that in a normal material would be done in
     * do_init()
     */
    virtual void do_init_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);


    /*! \copydoc TiberModelObject::override_parameter_string()
     *
     * This method looks first in the ModelOptions object, and then in the
     * ModelOptions of the material. It will try in this order:
     * \li \c name
     * \li \c modelname.name
     * \li \c simulationname.modelname.name
     */
    virtual void override_parameter_string(const std::string& name,
        std::string& s) const;


    //! calculate an alloy parameter in VCA approximation
    /*!
     * In a ternary compound semiconductor
     * \f$Q = A_xB_{1-x}C\f$ the value of a
     * material parameter can (in the virtual crystal approximation) be
     * calculated as
     * \f[\alpha_Q = x\alpha_{AC} + (1-x)\alpha_{BC} - bx(1-x)\f]
     * where \em b is called bowing parameter and describes deviation
     * from the nonlinear behaviour.
     *
     * \param val_a the value for material \f$AC\f$
     * \param val_b the value for material \f$BC\f$
     * \param xa the molar fraction of \f$AC\f$
     * \param bowing the bowing parameter
     */
    double alloy(double val_a, double val_b, double xa, double bowing = 0.0);


    /*!
     * \copydoc alloy(double, double, double, double)
     *
     * \param result the value for material \f$  A_xB_{1-x}C\f$
     */
    void alloy(Tensor4DSym& result, const Tensor4DSym& val_a,
        const Tensor4DSym& val_b, double xa,
        const Tensor4DSym& bowing = Tensor4DSym(0));


    /*! \copydoc alloy(double, double, double, double) */
    void alloy(Tensor2Sym& result, const Tensor2Sym& val_a,
        const Tensor2Sym& val_b, double xa,
        const Tensor2Sym& bowing = Tensor2Sym(0));



  private:

    //! An iterator for the models
    typedef std::map<const std::string, ID>::iterator model_id_iterator;


    //! Disable copy constructor
    PhysicalModelInterface(const PhysicalModelInterface&);


    //! Disable assignment operator
    PhysicalModelInterface& operator=(const PhysicalModelInterface&);


    //! The unique ID of this model
    ID _id;


    //! The ID of the simulator this model is used for
    /*!
     * This ID is needed by the simulator to get the right models
     */
    ID _simulator_id;


    //! A user defined name for this model
    /*!
     * The name is assigned from the ModelOptions.
     */
    std::string _name;


    //! The identifying string for the type of this model
    std::string _type;


    //! The material this properties belong to
    /*!
     * This pointer can be used by this or associated models
     */
    Material* _material;


    //! A map with ID/model name pairs
    /*!
     * Models are counted starting from 1. 0 means undefined model.
     */
    static std::map<const std::string, ID> _model_ids;


    //! Register a new model
    /*!
     * This method registers every new model that gets created and assigns
     * it a unique model ID.
     */
    static void register_model(PhysicalModelInterface* model);


    //! Set the model type (= identifier)
    /*!
     * The identifier is used at creation time to know which type of
     * model to create.
     */
    void set_type(const std::string& type);


};


//
// inline methods
//

inline
PhysicalModelInterface::PhysicalModelInterface(void)
  : _id(INVALID_ID),
    _simulator_id(INVALID_ID),
    _name(""),
    _material(NULL)
{
}


inline
void
PhysicalModelInterface::set_type(const std::string& type)
{
  _type = type;
}



inline
const std::string&
PhysicalModelInterface::get_type(void) const
{
  return _type;
}



inline
double
PhysicalModelInterface::alloy(double val_a, double val_b,
    double xa, double bowing)
{
  return val_b + (val_a - val_b) * xa - bowing * xa * (1 - xa);
}



inline
void
PhysicalModelInterface::alloy(Tensor4DSym& result, const Tensor4DSym& val_a,
    const Tensor4DSym& val_b, double xa, const Tensor4DSym& bowing)
{
 result = (1 - xa) * val_b + xa * val_a - xa * (1 - xa) * bowing ;
}


inline
void
PhysicalModelInterface::alloy(Tensor2Sym& result, const Tensor2Sym& val_a,
    const Tensor2Sym& val_b, double xa, const Tensor2Sym& bowing)
{
  result = (1 - xa) * val_b + xa * val_a - xa * (1 - xa) * bowing ;
}



inline
ID
PhysicalModelInterface::get_id(void) const
{
  return _id;
}



inline
ID
PhysicalModelInterface::get_simulator_id(void) const
{
  return _simulator_id;
}



inline
void
PhysicalModelInterface::set_simulator_id(ID id)
{
  _simulator_id = id;
}



inline
const std::string&
PhysicalModelInterface::get_name(void) const
{
  return _name;
}



inline
void
PhysicalModelInterface::set_name(const std::string& name)
{
  _name = name;
}



inline
void
PhysicalModelInterface::set_material(Material* material)
{
  _material = material;
}




inline
Material*
PhysicalModelInterface::get_material(void)
{
  return _material;
}



inline
const Material*
PhysicalModelInterface::get_material(void) const
{
  return _material;
}




inline
void
PhysicalModelInterface::init(void)
{
  read_database();
  do_init();
}



template <typename T>
T*
PhysicalModelInterface::create_submodel_copy(const T* other) const
{
  assert(other != NULL);
  T* newmod = static_cast<T*>(other->copy());
  assert(newmod != NULL);
  newmod->set_material(_material);
  return newmod;
}



template <typename T>
T*
PhysicalModelInterface::create_submodel_alloy(const T* comp_A,
    const T* comp_B, double xa)
{
  assert((comp_A != NULL) && (comp_B != NULL));
  T* newmod = create_submodel_copy(comp_A);
  newmod->init_alloy(comp_A, comp_B, xa);
  return newmod;
}




inline
void
PhysicalModelInterface::init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  assert(typeid(*comp_A) == typeid(*comp_B));
  read_database_alloy();
  do_init_alloy(comp_A, comp_B, xa);
}




inline
void
PhysicalModelInterface::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  ignore_unused_variable(comp_A);
  ignore_unused_variable(comp_B);
  ignore_unused_variable(xa);
}


inline
void
PhysicalModelInterface::print_info(void)
{
  do_print_info();
}


template <typename T>
ID
PhysicalModelInterface::get_id_from_name(const std::string& name)
{
  ID id = 0;

  PhysicalModelInterface* rec = T::create(name);

  if (rec != NULL)
    id = rec->get_id();

  // rec is either valid or NULL
  destroy(rec);


  return id;
}


#endif // _PHYSICALMODELINTERFACE_H_
