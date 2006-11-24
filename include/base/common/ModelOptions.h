// $Id$

#ifndef _MODELOPTIONS_H_
#define _MODELOPTIONS_H_

#include <map>
#include <string>

//! A class to store model options
class ModelOptions
{

  public:

    //! Get the value of an option as double
    /*!
     * \param name the name of the option
     * \param default_value the default value
     */
    double get_option(const std::string& name, double default_value) const;

    //! Get the value of an option as int
    /*!
     * \param name the name of the option
     * \param default_value the default value
     */
    int get_option(const std::string& name, int default_value) const;

    //! Get the value of an option as bool
    /*!
     * \param name the name of the option
     * \param default_value the default value
     */
    bool get_option(const std::string& name, bool default_value) const;
    
    //! Get the value of an option as string
    /*!
     * \param name the name of the option
     * \param default_value the default value
     */
    const std::string& get_option(const std::string& name,
        const std::string& default_value) const;
    
    //! Get the value of an option as const char*
    /*!
     * \param name the name of the option
     * \param default_value the default value
     */
    const char* get_option(const std::string& name,
        const char* default_value) const;

    //! Check if an option is present
    bool find_option(const std::string& name) const;

    //! Set an option
    template <typename T>
    void set_option(const std::string& name, const T value);

    //! Set or get an option in string representation
    std::string& operator[](const std::string& name);


    //! Delete an option
    void delete_option(const std::string& name);

    //! Clear all options
    void clear(void);

    //! operator to add options
    ModelOptions& operator+=(const ModelOptions& rhs);


  private:

    //! typedef for convenience
    typedef std::map<const std::string, std::string> OptionsMap;
    
    //! The map holding all options
    OptionsMap _options;


};


inline
void
ModelOptions::delete_option(const std::string& name)
{
  _options.erase(name);
}

inline
std::string&
ModelOptions::operator[](const std::string& name)
{
  return _options[name];
}

inline
bool
ModelOptions::find_option(const std::string& name) const
{
  bool res = true;

  OptionsMap::const_iterator it(_options.find(name));
  
  if (it == _options.end())
    res = false;

  return res;
}

inline
void
ModelOptions::clear(void)
{
  _options.clear();
}


/*
template <typename T>
T
ModelOptions::get_option(const std::string& name, T default_value)
{
  return _options.get_option(name, default_value);
  
  T val(default_value);

  ModelOptions::const_iterator it = _options.find(name);

  if (it != _options.end())
  {
    std::istringstream s((it->second).c_str());
    s >> val;
  }

  return val;
}
*/


#endif // _MODELOPTIONS_H_
