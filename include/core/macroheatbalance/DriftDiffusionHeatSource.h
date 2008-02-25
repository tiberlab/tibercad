#ifndef _DRIFTDIFFUSIONHEATSOURCE_H_
#define _DRIFTDIFFUSIONHEATSOURCE_H_

#include "PhysicalModelInterface.h"
#include "HeatSourceInterface.h"       




//! This class handles the heat source


class DriftDiffusionHeatSource : public HeatSourceInterface
{

public:
  
  //!Constructor 
  DriftDiffusionHeatSource(){}

   //!Destructor
  ~DriftDiffusionHeatSource(){}

  //! set the electron conducibility;
  void get_total_heat_source();

  //! Costructor
  static  DriftDiffusionHeatSource* create();
 

  virtual void get_heat_source_output(std::vector<Point> h_point,const Elem* elem,
			       std::vector< std::vector< double > >& heat_source);


  virtual void get_heat_source(std::vector<Point> h_point,const Elem* elem,
				      std::vector< double >& heat_source);

  virtual void get_flux_heat_source(std::vector<Point> h_point,const Elem* elem,
			       std::vector<RealGradient>& heat_source);
     
  virtual std::vector<std::string>  get_legend();

private:

  //  std::vector<std::string> _legend;

  enum heat_variables
    {
      JNGRADPHIE,
      JPGRADPHIH,
      EJOULE,
      HJOULE,  
      PHIE,
      PHIH,  
      PN,
      PP,
      TOLD,
      SRHREC,
      PELTIERNX,
      PELTIERNY,
      PELTIERNZ,
      PELTIERPX,
      PELTIERPY,
      PELTIERPZ,
      PELTHOME,
      PELTHOMH
      
    };


   
  std::vector<std::string> _legend;

  SimulationInterface* _simul;


  //!Source variables
  std::set<ID> ID_set;

  std::vector<ID> ID_vector;



protected:

    

  virtual void do_init (void);

  virtual void copy_from(const PhysicalModelInterface *rhs);

  virtual void read_database(void){};

  virtual void read_bowing_parameters(void){};

  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 

  virtual PhysicalModelInterface* create_new (void) const;
 
 
 
};




inline
void    DriftDiffusionHeatSource::get_total_heat_source()
{
  
 
}



inline
PhysicalModelInterface*  DriftDiffusionHeatSource::create_new () const
{
  return (new    DriftDiffusionHeatSource() ); 
}

inline
DriftDiffusionHeatSource*
DriftDiffusionHeatSource::create()
{
  return new DriftDiffusionHeatSource(); 
}



inline
std::vector<std::string>
DriftDiffusionHeatSource::get_legend(void)
{

 return  _legend;

}

#endif
