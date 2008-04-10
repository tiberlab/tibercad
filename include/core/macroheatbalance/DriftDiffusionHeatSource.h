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
 
 
  virtual void get_heat_sources(std::vector<Point> h_point,
			       std::vector< std::vector< double > >& heat_source);


  virtual void get_power_fluxes(std::vector<Point> h_point,
				std::vector<std::vector<RealGradient> >& power_fluxes);
               
					
  //!Set the current element
   virtual void set_heat_model(HeatModel* heat_model);
    
  //!Init the heat source model
  // virtual void re_init(void){};
 
  virtual std::vector<std::string>  get_source_legend(const std::set<std::string>& variables);

  virtual std::vector<std::string>  get_flux_legend();

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
      TEMP,
      PELTHE,
      PELTHH,
      SRHREC,
      WNX,
      WNY,
      WNZ,
      WPX,
      WPY,
      WPZ
      
    };


  unsigned int n_s;
     
  unsigned int n_f;

  std::vector<std::string> _flux_legend;

  std::vector<std::string> _source_legend;

  SimulationInterface* _simul;

  HeatModel* _heat_model;

  //!Source variables
  std::set<ID> ID_set;

  std::vector<ID> ID_vector;

  //!Set the current element
  //void set_element(const Elem* elem);
  
  //!Set the current elemement side index
  //void  set_side(unsigned int side);

 
protected:

    

  virtual void do_init (void);

  virtual void copy_from(const PhysicalModelInterface *rhs);

  virtual void read_database(void){};

  virtual void read_bowing_parameters(void){};

  virtual void calculate_VCA (const PhysicalModelInterface *comp_A, const PhysicalModelInterface *comp_B, double xa); 

  virtual PhysicalModelInterface* create_new (void) const;
 
 
 
};






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
DriftDiffusionHeatSource::get_flux_legend(void)
{

 return  _flux_legend;

}

inline
void 
DriftDiffusionHeatSource::set_heat_model(HeatModel* heat_model)
{

 _heat_model = heat_model;

}

// inline
// void 
// DriftDiffusionHeatSource::set_element(const Elem* elem)
// {

//   _elem = elem;

// }



// inline
// void 
// DriftDiffusionHeatSource::set_side(unsigned int side)
// {

//   _side = side;

// }

#endif
