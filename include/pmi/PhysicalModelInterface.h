// $Id$

#ifndef _PHYSICALMODELINTERFACE_H_
#define _PHYSICALMODELINTERFACE_H_

#include "tiber_config.h"
#include "TiberModelObject.h"
#include "TypeDefs.h"
#include "ModelOptions.h"
#include "InitFailedException.h"
#include "ModelErrorException.h"
#include "Material.h"

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#include "tensor.h"
#include "xtensor.h"


#ifndef MODULE_NAME
#define MODULE_NAME
#endif

class PhysicalObject;
class MaterialBoundary;
class EdgeObject;
class NodeObject;
class Database;
class Material;


//! Base class for the different physical models
/*!
 * This is the base class for all implementations of any kind of models
 * used in simulations.
 *
 * Any class derived from this one has to implement not only the needed
 * virtual methods, but also
 *  static PhysicalModelInterface* create(const ModelOptions&)
 * This method could use the options for instantiation of specialized
 * classes.
 *
 *
 * Every physical model is owned by a PhysicalObject object. In the case of
 * interface, edge and node models, they get the according bulk material assigned
 * during calculations which can be obtained with get_material().
 */
class PhysicalModelInterface : public TiberModelObject
{

  public:

    //! Destructor
    virtual ~PhysicalModelInterface(void);


    //! Get the unique ID of this model type
    ID get_id(void) const;

    /*//! Get the unique ID of this model type
    ID get_unique_id (void) const;

    ID get_unique_id (std::string& name, std::string& material) const;
    */

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
    //template <typename T>
    //static ID get_id_from_name(const std::string& model_name);

    //! Creates a new named model
    /*!
     * The model is created according to the given model name.
     * If it is not known, the NULL pointer is returned.
     *
     * \param name the model name
     * \param owner the PhysicalObject this model is associated with (can be \c NULL)
     * \param options the options as given in the input file
     * \param module the module this model belongs to (in most cases
     *   found automatically)
     * \return a pointer to the newly created object
     */
    template <typename T = PhysicalModelInterface>
    static T* create(const std::string& name,
        const PhysicalObject* owner,
        const ModelOptions& options = ModelOptions(),
        const std::string& module = xstr(MODULE_NAME));


    //! Creates a new model from a given creator function
    /*!
     * Use this method only in special cases, e.g. if you don't have single
     * libraries for the different models.
     */
    template <typename T = PhysicalModelInterface>
    static T* create(create_t create_fnc, destroy_t destroy_fnc,
        const PhysicalObject* owner,
        const ModelOptions& options = ModelOptions(),
        const std::string& module = xstr(MODULE_NAME));



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
    const Database& get_database(void);


    //! Set the bulk material
    /*!
     * A model, even defined on a lower dimensional region, should
     * always have a reference to a bulk material. This reference
     * can be set using this method explicitly, but it will be set
     * automatically when obtaining the model using one of the
     * methods provided in SimulationInterface.
     */
    void set_material(const Material* mat);


    //! Get the bulk material
    /*!
     * A model, even defined on a lower dimensional region, should
     * always have a reference to a bulk material. This reference
     * can be set using this method explicitly using set_bulk_material(),
     * but it will be set automatically when obtaining the model
     * using one of the methods provided in SimulationInterface.
     */
    const Material* get_material(void) const;


    //! Set a reference to the physical object this model belongs to
    void set_owner(const PhysicalObject* owner);


    //! Get a reference to the physical object this model belongs to
    //PhysicalObject* get_owner(void);


    //! Get a reference to the physical object this model belongs to
     const PhysicalObject* get_owner(void) const;



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
     * It calls read_database(), do_init() and create_submodels()
     *
     * Note that do_init() will also be called for alloys to allow for
     * correct parameter override and to assure complete initialization.
     */
    void init(void);


    //! Reinitialize the model
    /*!
     * This method will be called before any solve of the associated
     * simulation module.
     */
    void reinit(void);


    //! Call this to prepare model data
    void reinit(const Elem* elem);


    //! Initialize this model as an alloy with two components
    /*!
     * It calls read_database_alloy(), do_init_alloy() and do_init()
     */
    void init_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);


    //! Initialize this model as interface model
    /*!
     * It calls read_interface_database() and do_init_interface()
     */
    void init_interface(const Material* comp_A,
        const Material* comp_B);


    //! Print some info
    void print_info(void);


    //! Get the name of the module this object belong to
    const std::string& get_module_name(void) const;



  protected:


    //! The type of the submodel list
    typedef std::multimap<std::string, PhysicalModelInterface*> SubmodelMap;

    //! An iterator for the submodels
    typedef SubmodelMap::iterator SubmodelIterator;

    //! An const iterator for the submodels
    typedef SubmodelMap::const_iterator ConstSubmodelIterator;


    //! The constructor
    PhysicalModelInterface(const ModelOptions& options);


    //! Print some info
    virtual void do_print_info(void) {};


    //! Reinitialize before solving the module
    /*!
     * May be reimplemented if necesary
     */
    virtual void do_reinit(void) {};


    //! Reinitialize for a certain element to prepare model data
    /*!
     * May be reimplemented if necesary
     */
    virtual void do_reinit(const Elem* elem);


    //! Initialize the model
    /*!
     * This method should be implemented in derived classes.
     */
    virtual void do_init(void) {};


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
    virtual void copy_from(const PhysicalModelInterface* rhs);


    //! Create a new instance of the same model
    /*!
     * This is needed for creating alloy models. Ususally, the method does
     * not have to be implemented.
     */
    virtual PhysicalModelInterface* create_new(void) const;


    //! Read the properties from the database
    /*!
     * Reads all needed physical properties from the database.
     * In the case of alloys it will do also the mixing according
     * to the chosen rule (VCA for now)
     * If you reimplement this in a derived class, call the method
     * of the base class, too.
     */
    virtual void read_database(void) {};


    //! Read alloy parameters from the database
    /*!
     * Read parameters for an alloy model from the database.
     * This method will read exclusively from the alloy database and not
     * attempt to do any mixing.
     * If you reimplement this in a derived class, call the method
     * of the base class, too.
     */
    virtual void read_database_alloy(void) {};


    //! Read interface parameters from the database
    /*!
     * Read parameters for an interface model from the database.
     * If you reimplement this in a derived class, call the method
     * of the base class, too.
     */
    virtual void read_interface_database(void) {};


    //! Initialize an alloy
    /*!
     * Calculates all parameters of an alloy \f$A_xB_{x-1}C\f$.
     *
     * \note do_init() will be called after do_init_alloy()
     */
    virtual void do_init_alloy(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);


    //! Initialize an interface model
    /*!
     * Calculates all parameters of an interface.
     *
     * \note do_init() will \emph not be called automatically
     *      after do_init_interface()
     *
     * \note The default implementation calls do_init()
     */
    virtual void do_init_interface(const Material* comp_A,
        const Material* comp_B);


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


    //! Create a submodel
    /*!
     * If no submodel of type \c type is provided in the input file,
     * \c NULL will be returned.
     * An exception is thrown, if more than one specification is found
     * in the input file.
     */
    template <typename T>
    void create_submodel(T*& model, const std::string& type);


    //! Create a submodel
    /*!
     * If no submodel of type \c type is provided in the input file,
     * the default options are used to create an instance of the model.
     */
    template <typename T>
    void create_submodel(T*& model, const std::string& type,
        const ModelOptions& default_opts);


    //! Create multiple instances of submodel \c type
    template <typename T>
    void create_submodels(std::vector<T*>& models,
        const std::string& type);


    //! Create multiple instances of submodel \c type
    /*!
     * If no submodel of type \c type is provided in the input file,
     * the default options are used to create an instance of the model.
     */
    template <typename T>
    void create_submodels(std::vector<T*>& models,
        const std::string& type, const ModelOptions& default_opts);


    //! Create submodels
    /*!
     * This method is to be used to create submodels, if they are needed for all bulk,
     * interface, edges and nodes objects. Otherwise, submodels can be created inside
     * one of the \c init_xxx methods.
     *
     * Submodels can be created be using one of the \c create_submodel() or
     * \c create_submodels() methods or by calling directly PhysicalModelInterface::create()
     * and using \c add_submodel() subsequently.
     *
     * Subsequent operations on the submodels assume that they are ordered exactly
     * the same way in models associated to alloy or interface components.
     * This is assured as long as all submodels are created in create_submodels() in
     * a way independent of the type of "owner" (PhysicalObject)
     */
    virtual void prepare_submodels(void) {};


    //! Add an externally created submodel
    /*!
     * \param key the name of the model
     * \param pm the pointer to the model
     */
    void add_submodel(const std::string& key, PhysicalModelInterface* pm);


    //! Destroy a submodel
    void delete_submodel(const std::string& key);


    //! Get the iterator to the first submodel
    SubmodelIterator submodels_begin(void);


    //! Get the iterator to the first submodel
    ConstSubmodelIterator submodels_begin(void) const;


    //! Get the submodels past-the-end iterator
    SubmodelIterator submodels_end(void);


    //! Get the submodels past-the-end iterator
    ConstSubmodelIterator submodels_end(void) const;


    //! Get the iterator for a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return the const iterator for the first appearance of the model
     */
    SubmodelIterator submodels_begin(const std::string& name);


    //! Get the iterator for a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return the const iterator for the first appearance of the model
     */
    ConstSubmodelIterator submodels_begin(const std::string& name) const;


    //! Get the past-the-end iterator for a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return the past-the-end iterator for the model
     */
    SubmodelIterator submodels_end(const std::string& name);


    //! Get the past-the-end iterator for a certain submodel
    /*!
     * \param name the name of the model to look for
     * \return the past-the-end iterator for the model
     */
    ConstSubmodelIterator submodels_end(const std::string& name) const;

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

    /*//! used if a model needs a unique id to set _has_unique_id = true
    void has_unique_id(bool flag = false)
    { _has_unique_id = true; };*/


  private:

    //! An iterator for the models
    typedef std::map<const std::string, ID>::iterator model_id_iterator;

    /*//! An iterator for models with unique ID
    typedef std::map< std::pair<const std::string,
                                const std::string>,
                      ID >::iterator unique_model_id_iterator;*/


    //! Disable copy constructor
    PhysicalModelInterface(const PhysicalModelInterface&);


    //! Disable assignment operator
    PhysicalModelInterface& operator=(const PhysicalModelInterface&);


    //! The unique ID of this model type
    ID _id;

    /*//! The unique ID of this model
    ID _unique_id;

    //! False by default. If the model needs a unique ID it has to be set to true in derived classes
    bool _has_unique_id = false;*/

    //! The ID of the simulator this model is used for
    /*!
     * This ID is needed by the simulator to get the right models
     */
    ID _simulator_id;


    //! The identifying string for the type of this model
    std::string _type;


    //! The physical object this properties belong to
    /*!
     * This pointer can be used by this or associated models.
     */
    const PhysicalObject* _owner;


    //! Lower dimensional objects can have a bulk material assigned
    const Material* _bulk_material;


    //! The name of the module this object is part of
    std::string _module;


    //! A list of submodels
    SubmodelMap _submodels;


    //! A map with ID/model name pairs
    /*!
     * Models are counted starting from 1. 0 means undefined model.
     */
    static std::map<const std::string, ID> _model_ids;

    /*//! A map with unique ID/model name pairs
    static std::map< std::pair<const std::string,
                               const std::string>,
                     ID > _unique_model_ids;*/


    //! Register a new model
    /*!
     * This method registers every new model that gets created and assigns
     * it a unique model ID.
     */
    static void _register_model(PhysicalModelInterface* model) TBDLLOCAL;


    //! Set the model type (= identifier)
    /*!
     * The identifier is used at creation time to know which type of
     * model to create.
     */
    void _set_type(const std::string& type);


    //! Set the module name
    void _set_module_name(const std::string& module);


    //! Create automatically all submodels
    /*!
     * Calls prepare_submodels() which can be overridden by module developers.
     */
    void _create_submodels(void);

    //! Create a submodel
    void _create_submodel(PhysicalModelInterface*& model, const std::string& type);

    //! Create a submodel
    void _create_submodel(PhysicalModelInterface*& model, const std::string& type,
        const ModelOptions& default_opts);

    //! Create submodels
    void _create_submodels(std::vector<PhysicalModelInterface*>& models,
        const std::string& type);

    //! Create submodels
    void _create_submodels(std::vector<PhysicalModelInterface*>& models,
        const std::string& type, const ModelOptions& default_opts);

    //! The internal implementation of the create method
    static PhysicalModelInterface* _create(
        const std::string& name,
        const PhysicalObject* owner,
        const ModelOptions& options,
        const std::string& module);

    //! The internal implementation of the create method
    static PhysicalModelInterface* _create(
        create_t create_fnc, destroy_t destroy_fnc,
        const PhysicalObject* owner,
        const ModelOptions& options,
        const std::string& module);


};


//
// inline methods
//

inline
PhysicalModelInterface::PhysicalModelInterface(const ModelOptions& options)
  : TiberModelObject(options),
    _id(INVALID_ID),
    _simulator_id(INVALID_ID),
    _owner(NULL),
    _bulk_material(NULL),
    _module("")
{
}


inline
void
PhysicalModelInterface::_set_type(const std::string& type)
{
  _type = type;
}



inline
void
PhysicalModelInterface::_set_module_name(const std::string& module)
{
  _module = module;
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

/*inline
ID
PhysicalModelInterface::get_unique_id(void) const
{
  return _unique_id;
}

inline
ID
PhysicalModelInterface::get_unique_id(std::string& name, std::string& material) const
{
  unique_model_id_iterator it = _unique_model_ids.find(std::make_pair(name, material));

  if (it != _unique_model_ids.end())
    return _unique_model_ids[std::make_pair(name, material)];

  return INVALID_ID;
}
*/

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
const Material*
PhysicalModelInterface::get_material(void) const
{
  return _bulk_material;
}



/*
inline
PhysicalObject*
PhysicalModelInterface::get_owner(void)
{
  return _owner;
}
*/

inline
const PhysicalObject*
PhysicalModelInterface::get_owner(void) const
{
  return _owner;
}



inline
const std::string&
PhysicalModelInterface::get_module_name(void) const
{
  return _module;
}



template <typename T>
T*
PhysicalModelInterface::create_submodel_copy(const T* other) const
{
  assert(other != NULL);
  T* newmod = static_cast<T*>(other->copy());
  assert(newmod != NULL);
  newmod->set_owner(_owner);
  newmod->set_material(_bulk_material);
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




template <typename T>
T*
PhysicalModelInterface::create(const std::string& name, const PhysicalObject* owner,
    const ModelOptions& options, const std::string& module)
{
  PhysicalModelInterface* mod =
      PhysicalModelInterface::_create(name, owner, options, module);

  if (mod != dynamic_cast<T*>(mod))
  {
    throw ModelErrorException("Given model type \'" + name + "\' does not correspond to "
        "type of created model.");
  }
  return static_cast<T*>(mod);
}


template <typename T>
T*
PhysicalModelInterface::create(create_t create_fnc, destroy_t destroy_fnc,
        const PhysicalObject* owner,
        const ModelOptions& options,
        const std::string& module)
{
  PhysicalModelInterface* mod =
      PhysicalModelInterface::_create(create_fnc, destroy_fnc,
          owner, options, module);

  // No check, this would be a programmer error
  return static_cast<T*>(mod);
}



template <typename T>
void
PhysicalModelInterface::create_submodel(T*& model, const std::string& type)
{
  PhysicalModelInterface* mod = NULL;
  _create_submodel(mod, type);
  model = static_cast<T*>(mod);
}


template <typename T>
void
PhysicalModelInterface::create_submodel(T*& model, const std::string& type,
    const ModelOptions& default_opts)
{
  PhysicalModelInterface* mod = NULL;
  _create_submodel(mod, type, default_opts);
  model = static_cast<T*>(mod);
}



template <typename T>
void
PhysicalModelInterface::create_submodels(std::vector<T*>& models,
    const std::string& type)
{
  std::vector<PhysicalModelInterface*> mod;
  _create_submodels(mod, type);
  models.resize(mod.size());
  for (size_t i = 0; i < mod.size(); i++)
    models[i] = static_cast<T*>(mod[i]);
}


//! Create multiple instances of submodel \c type
template <typename T>
void
PhysicalModelInterface::create_submodels(std::vector<T*>& models,
    const std::string& type, const ModelOptions& default_opts)
{
  std::vector<PhysicalModelInterface*> mod;
  _create_submodels(mod, type, default_opts);
  models.resize(mod.size());
  for (size_t i = 0; i < mod.size(); i++)
    models[i] = static_cast<T*>(mod[i]);
}



inline
void
PhysicalModelInterface::do_init_alloy(const PhysicalModelInterface*,
    const PhysicalModelInterface*, double)
{
}


inline
void
PhysicalModelInterface::do_init_interface(const Material*,
    const Material*)
{
  do_init();
}


inline
void
PhysicalModelInterface::do_reinit(const Elem*)
{
}


inline
void
PhysicalModelInterface::print_info(void)
{
  do_print_info();
}

/*
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
*/

inline
void
PhysicalModelInterface::copy_from(const PhysicalModelInterface* rhs)
{
  ignore_unused_variable(rhs);
}



inline
PhysicalModelInterface::SubmodelIterator
PhysicalModelInterface::submodels_begin(void)
{
  return _submodels.begin();
}


inline
PhysicalModelInterface::SubmodelIterator
PhysicalModelInterface::submodels_end(void)
{
  return _submodels.end();
}


inline
PhysicalModelInterface::SubmodelIterator
PhysicalModelInterface::submodels_begin(const std::string& name)
{
  return _submodels.lower_bound(name);
}



inline
PhysicalModelInterface::SubmodelIterator
PhysicalModelInterface::submodels_end(const std::string& name)
{
  return _submodels.upper_bound(name);
}


inline
PhysicalModelInterface::ConstSubmodelIterator
PhysicalModelInterface::submodels_begin(void) const
{
  return _submodels.begin();
}


inline
PhysicalModelInterface::ConstSubmodelIterator
PhysicalModelInterface::submodels_end(void) const
{
  return _submodels.end();
}


inline
PhysicalModelInterface::ConstSubmodelIterator
PhysicalModelInterface::submodels_begin(const std::string& name) const
{
  return _submodels.lower_bound(name);
}



inline
PhysicalModelInterface::ConstSubmodelIterator
PhysicalModelInterface::submodels_end(const std::string& name) const
{
  return _submodels.upper_bound(name);
}



#endif // _PHYSICALMODELINTERFACE_H_
