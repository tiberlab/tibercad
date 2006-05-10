// $Id$

#ifndef _BOUNDARYDESCRIPTOR_H_
#define _BOUNDARYDESCRIPTOR_H_

// C++ includes
#include <vector>
#include <string>
#include <map>
#include <cmath>

// forward declarations

class BoundaryDescriptor
{

  public:

    enum BoundaryConditionType {DIRICHLET, NEUMANN, MIXED};

    typedef std::map<const std::string,
            std::vector<double> >::iterator iterator;

    typedef std::map<const std::string,
            std::vector<double> >::const_iterator const_iterator;

    explicit BoundaryDescriptor(const std::string identifier = "");

    const std::string& get_id(void) const;

    int get_size(void) const;

    const_iterator end(void) const;

    const_iterator begin(void) const;

    const_iterator find(const std::string& variable) const;

    iterator begin(void);

    iterator find(const std::string& variable);

    const std::vector<double>*
      get_coefficients(const std::string& variable) const;

    const std::vector<double>& operator[](const std::string& variable);

    void set_coefficients(const std::string& variable,
        const std::vector<double>& coefficients);

    BoundaryConditionType get_type(const std::string& variable) const;

    /**
     * returns the coefficients for the boundary condition formulated as
     *  du/dn = - a' * u + c'
     * where a' = a / b
     *       c' = c / b
     * 
     * For |b| < eps, a penalty function approach is used
     */
    const std::vector<double>
      get_scaled_to_normal_derivative(const std::string& variable) const;

    void print_info(void) const;
      
    double pinning;

  private:

    std::string _id;


    const static double _penalty_value;
    
    /**
     * the coefficients for the mixed boundary condition
     * a * u + b * du/dn = c
     * in the order [a b c] for each variable
     */
    std::map<const std::string, std::vector<double> > _data;
  
};

//
// inline member functions
// 

inline
BoundaryDescriptor::BoundaryDescriptor(const std::string identifier)
  : _id(identifier),
    pinning(-256.0)
{
}

inline
void
BoundaryDescriptor::set_coefficients(const std::string& variable,
        const std::vector<double>& coefficients)
{
  _data[variable] = coefficients;
}

inline
const std::string&
BoundaryDescriptor::get_id(void) const
{
  return _id;
}

inline
const std::vector<double>*
BoundaryDescriptor::get_coefficients(const std::string& variable) const
{
  const_iterator it = find(variable);
  if (it != end())
    return &(it->second);
  return NULL;
}

inline
const std::vector<double>&
BoundaryDescriptor::operator[](const std::string& variable)
{
  return _data[variable];
}

inline
BoundaryDescriptor::const_iterator
BoundaryDescriptor::find(const std::string& variable) const
{
  return _data.find(variable);
}

inline
BoundaryDescriptor::iterator
BoundaryDescriptor::find(const std::string& variable)
{
  return _data.find(variable);
}

inline
int
BoundaryDescriptor::get_size(void) const
{
  return _data.size();
}

inline
BoundaryDescriptor::const_iterator
BoundaryDescriptor::end(void) const
{
  return _data.end();
}

inline
BoundaryDescriptor::const_iterator
BoundaryDescriptor::begin(void) const
{
  return _data.begin();
}

inline
BoundaryDescriptor::iterator
BoundaryDescriptor::begin(void)
{
  return _data.begin();
}

inline
BoundaryDescriptor::BoundaryConditionType
BoundaryDescriptor::get_type(const std::string& variable) const
{
  double eps = 1e-12;

  BoundaryConditionType type = NEUMANN;
  const_iterator it = find(variable);
  if (it != end())
  {
    if (std::fabs(it->second[0]) > eps)
      if (std::fabs(it->second[1]) > eps)
        type = MIXED;
      else
        type = DIRICHLET;
  }
  return type;
}


#endif //_BOUNDARYDESCRIPTOR_H_
