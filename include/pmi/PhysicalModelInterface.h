// $Id$

#ifndef _PHYSICALMODELINTERFACE_H_
#define _PHYSICALMODELINTERFACE_H_

#include "tiber_config.h"
#include "TypeDefs.h"
#include "ModelOptions.h"
#include "InitFailedException.h"

// For debugging
#include "reference_counted_object.h"

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#include "tensor.h"
#include "xtensor.h"


#ifndef TIBER_MODULE_NAME
# define TIBER_MODULE_NAME
#endif

#ifdef BUILD_TIBER_MODULES
/*!
 * \def TIBER_MODULE(classname, libname)
 *
 * \brief Creates methods to create and destroy a simulation object
 * 
 * In each implementation derived from SimulationInterface, put
 * this macro somewhere in the source file to be able to compile
 * it as TiberCad module.
 *
 * \param name the name of the class that should be 'creatable'
 * \param libname the name for this module
 *
 * \c libname will be used to create the library name, and the model
 * will have to be referred to in the input file by \c libname
 */
# ifndef TIBER_MODULE
#  define TIBER_MODULE(classname, libname) \
  extern "C" { \
    void destroy(PhysicalModelInterface* p) { \
      delete p; \
    } \
    classname* create(void) { \
      return new classname(); \
    } \
    const char* _tiber_module_ ## libname = #libname; \
    const char* library_name(void) { \
      return _tiber_module_ ## libname; \
    } \
  }
# endif
#else
# ifndef TIBER_MODULE
#  define TIBER_MODULE(classname, libname)
# endif
#endif



class Material;


//! Base class for the different physical models
/*!
 * This is the base class for all implementations of any kind of models
 * used in simulations.
 */
class PhysicalModelInterface
  : public ReferenceCountedObject<PhysicalModelInterface>
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

    
    //! Deletes a model
    /*!
     * \param p the pointer to the model to be destroyed
     */
    static void destroy(PhysicalModelInterface* p);

    
    //! Create a new model as an exact copy of this
    /*!
     * \return a pointer to the newly created model, \c NULL if
     * creation failed
     */
    PhysicalModelInterface* copy(void) const;
    

    //! Set a reference to the material this model belongs to
    void set_material(Material* material);


    //! Get a reference to the material this model belongs to
    const Material* get_material(void) const;

    
    //! Get a writeable reference to the material this model belongs to
    Material* get_material(void);


    //! Get the user defined name of this model
    const std::string& get_name(void) const;


    //! Get the default name for this model
    std::string get_default_name(void) const;

    
    //! Initialize this model
    /*!
     * It calls read_database() and then do_init()
     */
    void init(void) throw (InitFailedException);

    
    //! Build parameters for an alloy
    void build_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);


    //! Print some info
    void print_info(void);
    


  protected:

    //! Empty constructor
    PhysicalModelInterface(void);
 
    
    //! Set options for this model
    /*!
     * The options are stored internally and are accessible through
     * special methods.
     * Options have to be specified at creation time.
     */
    void set_options(const ModelOptions& options);

    
    //! Get a reference to the model options
    ModelOptions& get_options(void);
   

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


    //! Get the value of a parameter from the input file
    /*!
     *
     * This method looks first in the ModelOptions object, and then in the
     * ModelOptions of the material. It will try in this order:
     * \li \c name
     * \li \c modelname.name
     * \li \c simulationname.modelname.name
     *
     * \param name the name of the option
     * \param default_value the default value, which also defines
     * the type of the option
     * \return the value
     */
    template <typename T>
    T get_parameter(const std::string& name, T default_value) const;


    //! Get a parameter which is a vector of values (of the same type)
    /*!
     *
     * This method looks first in the ModelOptions object, and then in the
     * ModelOptions of the material. It will try in this order:
     * \li \c name
     * \li \c modelname.name
     * \li \c simulationname.modelname.name
     *
     * \param name the name of the option
     * \param vec the vector, where the values will be stored. \c vec can
     * contain default values, but it's size will be changed according to
     * the vector found in the options.
     */
    template <typename T>
    void get_parameter(const std::string& name, std::vector<T>& vec) const;
    
    
    //! Initialize the model
    /*!
     * This method should set all model options and call
     * \c init() of any associated model
     *
     * This method should be implemented in derived classes.
     */
    virtual void do_init(void) {};

    
    //! Create a new model of the same type
    virtual PhysicalModelInterface* create_new(void) const = 0;

    
    //! Copy all data from another model to this one
    /*!
     * If you reimplement this in a derived class, call the method
     * of the base class, too. If not, you may not copy important
     * things.
     */
    virtual void copy_from(const PhysicalModelInterface* rhs) = 0;

    
    //! Read the properties from the database
    /*!
     * Reads all needed physical properties from the database.
     * (The default behaviour is to do nothing at all.)
     * If you reimplement this in a derived class, call the method
     * of the base class, too.
     */
    virtual void read_database(void) {};

    
    //! Read the bowing parameters from the database
    /*!
     * Read the bowing parameters for an alloy material from the database.
     * (The default behaviour is to do nothing at all.)
     * If you reimplement this in a derived class, call the method
     * of the base class, too.
     */
    virtual void read_bowing_parameters(void) {};

    
    //! Calculate parameters for an alloy
    /*!
     * Calculates all parameters of an alloy \f$A_xB_{x-1}C\f$ in
     * virtual crystal approximation
     */
    virtual void calculate_VCA(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

    
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

    
    //! The creation method signature
    typedef PhysicalModelInterface* (*create_t)(void);

    
    //! The destruction method signature
    typedef void (*destroy_t)(PhysicalModelInterface*);


    //! The type for library handles
    typedef void* libhandle_t;

    
    //! Disable copy constructor
    PhysicalModelInterface(const PhysicalModelInterface&);

    
    //! Disable assignement operator
    PhysicalModelInterface& operator=(const PhysicalModelInterface&);


    //! The library handle for this model
    libhandle_t _libhandle;

    
    //! The creation method for this model
    create_t _create;

    
    //! The destruction method for this model
    destroy_t _destroy;

    
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

    
    //! The material this properties belong to
    /*!
     * This pointer can be used by this or associated models
     */
    Material* _material;

    
    //! The options for this model as read from the input file
    ModelOptions _options;

    
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
};


//
// inline methods
// 

inline
PhysicalModelInterface::PhysicalModelInterface(void)
  : _libhandle(NULL),
    _create(NULL),
    _destroy(NULL),
    _material(NULL)
{
}


inline
PhysicalModelInterface::~PhysicalModelInterface(void)
{
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
PhysicalModelInterface::set_options(const ModelOptions& options)
{
  _options = options;
}

inline
ModelOptions&
PhysicalModelInterface::get_options(void)
{
  return _options;
}


inline
void
PhysicalModelInterface::init(void) throw (InitFailedException)
{
  read_database();
  do_init();
}


inline
void
PhysicalModelInterface::build_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{
  read_bowing_parameters();
  calculate_VCA(comp_A, comp_B, xa);
}


inline
void
PhysicalModelInterface::calculate_VCA(const PhysicalModelInterface* comp_A,
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


#endif // _PHYSICALMODELINTERFACE_H_
