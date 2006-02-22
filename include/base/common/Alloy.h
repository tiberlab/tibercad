#ifndef ALLOY_H_
#define ALLOY_H_


#include <vector>
#include "Material.h"
#include "PhysicalProperties.h"
#include <point.h>

using  namespace  std;



//! Description of  an  alloy  material.
/*!
 *  This class is derived  from  \c Material. It  extends Material with a vector of 
 *  pointers to the  component materials of  the  alloy. \c get_properties is  reimplemented to combine  phisical properties
 * of  component  materials.
 */
class Alloy : public Material
{
  public:
   //! Construct an  alloy material with a given structure 
  /*!
   * At construction one has to specify at least the material name and
   * optionally also the structure. The default for the structure is 
   * zincblende.
   *
   * \param name the name of the Material
   * \param structure the crystal structure
   */
    Alloy(const std::string& name, const std::string& structure = "zb") : Material(name, structure){};
#include <PhysicalProperties.h>
    //! Destructor
   /*!
    * Deletes all \c PhysicalProperties objects
   */
    virtual ~Alloy();
    
    //! Initialize the alloy material
    /*!
     * Read all needed bowing parameters  data from the database
     * It calls the \c read_database_bowing_parameters() of each \c PhysicalProperties
     * object
     *
     * \param database the database to read from
   */	
    void init(const Dummy& database);
    
	
     //! Used to define the component materials  of the  \c Alloy object
   /*!
     *  \param matpoint pointer to a \c Material object
    */
    void set_components(Material* matpoint);
	
	
    //  const  PhysicalProperties* Alloy::get_properties(const std::string& id,  const vector<double>& coord)  ; 
    //! Get physical properties of a given type
  /*!
     * Get the physical properties of type \c id.
     * Reimplements base class method to calculate and set \c PhysicalProperties 
     * of the \c Alloy  material. 
     * It will return the \c NULL pointer if the requested properties
     * are not in the list.
     *
     * \param id the identifier for the set of properties
     * \param coord  coordinates of the  real space point (for position-dependent properties)
     * \return a const pointer to the property object
   */
    const  PhysicalProperties*  get_properties(const std::string& id,  const Point& coord)  ; 
    
	
	
  private:
    
    //! Component material of  the   \c Alloy
    /*!
   * A  vector of  pointers  to  \c Material objects.
   */
    vector<Material*>  components; 
    
     //! Molar fraction  of  the   \c Alloy
    /*!
     *  Composition of  the  \c Alloy.  It is calculated by \c calculate_molar_fraction  
     *  in function of  the given point coordinates .
     */
    double  molar_fraction;
    //	double calculate_molar_fraction(const vector<double>& coord);
    
     //! Calculates Molar fraction  of  the   \c Alloy
    /*!
     *  Calculates Molar fraction  of  the   \c Alloy  according to  a  given  dependence (possibly constant)
     *  on the  current  point coordinates.
     */
    double calculate_molar_fraction(const Point& coord);
			
};

//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------
inline 
double 
Alloy::calculate_molar_fraction(const Point& coord) 
{
	return 0.0;
}




#endif /*ALLOY_H_*/
