

#include "HostGuestRecombination.h"
#include "DriftDiffusionProperties.h"
#include "Database.h"
#include "SimulationInterface.h"
#include "SimulationEnvironment.h"
#include "Messages.h"

#include "DDBulkModel.h"

#include "mesh_base.h"
#include "quadrature.h"

#include "TiberModule.h"



using namespace std;

bool HostGuestRecombination::_coupled;

void
HostGuestRecombination::read_database(void)
{
  const Database& db = get_database();
  db.set_section("permittivity");

  _er = db.get("permittivity", _er);
}

void
HostGuestRecombination::do_init(void)
{
  get_parameter("gamma", _gamma);
  get_parameter("alpha", _alpha);
  _langevin = get_option("langevin", true);

  string csim = get_option("simulation", "");
  string old_sim = get_option("old_solution_simulation", "");

  _coupled_sim = SimulationInterface::find_simulation(csim);
  _saved_old_sim = SimulationInterface::find_simulation(old_sim);

  if (_coupled_sim == NULL)
  {
    std::string msg("Host/guest simulation " + std::string(csim) + " not found");
    throw InitFailedException(msg);
  }

  if (old_sim != "" ) {
    if (_saved_old_sim == NULL)
    {
      std::string msg("Saved solution simulation " + std::string(old_sim) + " not found");
    throw InitFailedException(msg);
    }
    _adiabatic = true;

    _old_erec = _saved_old_sim->get_solution_id("e" + get_default_name() );
    _old_hrec = _saved_old_sim->get_solution_id("h" + get_default_name() );

    //cout<<"_old_erec "<<_old_erec<<" _old_hrec "<<_old_hrec<<endl;

    if ( ( _old_erec == INVALID_ID ) || (_old_hrec == INVALID_ID ) )
    {
      std::string msg("HostGuestRecombination: simulation " + old_sim + " does not provide all necessary variables");
      throw InitFailedException(msg);
    }
  }

  _eDensity = _coupled_sim->get_solution_id("eDensity");
  _hDensity = _coupled_sim->get_solution_id("hDensity");
  _eDensity0 = _coupled_sim->get_solution_id("eDensity0");
  _hDensity0 = _coupled_sim->get_solution_id("hDensity0");
  _eMobility = _coupled_sim->get_solution_id("eMobility");
  _hMobility = _coupled_sim->get_solution_id("hMobility");

  _this_sim = SimulationInterface::get_simulation(get_simulator_id());

  HostGuestRecombination::_coupled = true;

}

void
HostGuestRecombination::get_net_recombination_rates(double& recomb_e, double& recomb_h)
{
  bool ok;
  if ( _coupled_sim->is_solved() && _coupled ) // && !(_this_sim->equilibrium() && _coupled_sim->equilibrium() ) )
  {

    const DriftDiffusionProperties& this_m = get_driftdiffusionproperties(); //host properties
    const Elem* el = this_m.get_element();

    DDBulkModel* coupled_m = _coupled_sim->get_bulk_model<DDBulkModel>(el);  //guest properties
    //coupled_m->set_coordinates(this_m.get_coordinates());

    //const Point& ct = this_m.get_coordinates();
    //const Point& cc = coupled_m->get_coordinates();

    if ( coupled_m == NULL )
    {
      recomb_e = recomb_h = 0.0;
      return;
    }

    double old_e_rec, old_h_rec = 0.0;

    double nt   = this_m.get_electron_density();
    double pt   = this_m.get_hole_density();
    double n0t  = this_m.get_quasi_equilibrium_n();
    double p0t  = this_m.get_quasi_equilibrium_p();
    double munt = this_m.get_electron_mobility();
    double mupt = this_m.get_hole_mobility();
    double Ect  = this_m.get_conduction_band().get_band_edge();
    double Evt  = this_m.get_valence_band().get_band_edge();
    double Nct  = this_m.get_conduction_band().get_effective_DOS();
    double Nvt  = this_m.get_valence_band().get_effective_DOS();
    double T    = this_m.get_lattice_temperature();

    double nc, pc, n0c, p0c, munc, mupc = 0.0;
    _coupled = false;

    _coupled_sim->get_solution(el, _eDensity, nc, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _hDensity, pc, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _eDensity0, n0c, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _hDensity0, p0c, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _eMobility, munc, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _hMobility, mupc, this_m.get_coordinates());

    _coupled = true;

    double Ecc  = coupled_m->get_conduction_band().get_band_edge();
    double Evc  = coupled_m->get_valence_band().get_band_edge();
    double Ncc  = coupled_m->get_conduction_band().get_effective_DOS();
    double Nvc  = coupled_m->get_valence_band().get_effective_DOS();

    //cout<<"nt = "<<nt<<" pt = "<<pt<<" munt = "<<munt<<" mupt = "<<mupt<<" Ect = "<<Ect<<" Evt = "<<Evt<<" Nct = "<<Nct<<" Nvt = "<<Nvt<<" x = "<<ct(0)<<" y = "<<ct(1)<<" z = "<<ct(2)<<endl;
    //cout<<"nc = "<<nc<<" pc = "<<pc<<" munc = "<<munc<<" mupc = "<<mupc<<" Ecc = "<<Ecc<<" Evc = "<<Evc<<" Ncc = "<<Ncc<<" Nvc = "<<Nvc<<endl;

    //cout<<"Sim name: "<<_sim->get_name()<<endl;

    if (Ect > Ecc) {
      recomb_e = nt * pc + nt * (Ncc - nc) - nc * (Nct - nt) * exp(-(Ect-Ecc)/T)  ;  //nt * pc +  nt * (Ncc - nc) - nc * (Nct - nt) * exp(-(Ect-Ecc)/T) 
      //recomb_e -= n0t * p0c + n0t * (Ncc - n0c) - n0c * (Nct - n0t) * exp(-(Ect-Ecc)/T);   //equilibrium term
      recomb_e *= _gamma;
      if (_langevin)
        recomb_e *= (munt + munc) * Constants::e * 100 / (_er * Constants::e0);

      //cout<<"nt = "<<nt<<" n0t = "<<n0t<<" nc = "<<nc<<" n0c = "<<n0c<<" recomb_e = "<<recomb_e<<endl;
      //cout<<"nt * (Ncc - nc) = "<<nt * (Ncc - nc)<<" n0t * (Ncc - n0c) = "<<n0t * (Ncc - n0c)<<" diff = "<<nt * (Ncc - nc) - n0t * (Ncc - n0c)<<endl;
      //cout<<"nc * (Nct - nt) = "<<nc * (Nct - nt)<<" n0c * (Nct - n0t) = "<<n0c * (Nct - n0t)<<" diff = "<<nc * (Nct - nt) - n0c * (Nct - n0t)<<endl;
      //cout<<"exp(-(Ect-Ecc)/T) = "<<exp(-(Ect-Ecc)/T)<<endl;
      //cout<<"Ect > Ecc recomb_e = "<<recomb_e<<" x = "<<this_m.get_coordinates()(0)<<" This sim name: "<<_this_sim->get_name()<<endl;
    }
    else {
     recomb_e = nt * pc + nt * (Ncc - nc) * exp(-(Ecc-Ect)/T) - nc * (Nct - nt)  ;  //nt * pc + nt * (Ncc - nc) * exp(-(Ecc-Ect)/T) - nc * (Nct - nt) 
     //recomb_e -= n0t * p0c + n0t * (Ncc - n0c) * exp(-(Ecc-Ect)/T) - n0c * (Nct - n0t)  ;   //equilibrium term
     recomb_e *= _gamma;
     if (_langevin)
       recomb_e *= (munt + mupc) * Constants::e * 100 / (_er * Constants::e0);

       //cout<<"Ect < Ecc recomb_e = "<<recomb_e<<" x = "<<this_m.get_coordinates()(0)<<" This sim name: "<<_this_sim->get_name()<<endl;
    }
    if (Evt < Evc) {
      recomb_h =  pt * nc + pt * (Nvc - pc) - pc * (Nvt - pt) * exp(-(Evc-Evt)/T)  ;  //pt * nc + pt * (Nvc - pc) - pc * (Nvt - pt) * exp(-(Evc-Evt)/T) 
      //recomb_h -= p0t * n0c + p0t * (Nvc - p0c) - p0c * (Nvt - p0t) * exp(-(Evc-Evt)/T)  ;   //equilibrium term
      recomb_h *= _gamma;
      if (_langevin)
        recomb_h *= (mupt + munc) * Constants::e * 100 / (_er * Constants::e0);
    }
    else {
      recomb_h = pt * nc + pt * (Nvc - pc) * exp(-(Evt-Evc)/T) - pc * (Nvt - pt)  ;  //pt * nc + pt * (Nvc - pc) * exp(-(Evt-Evc)/T) - pc * (Nvt - pt) 
      recomb_h -= p0t * n0c + p0t * (Nvc - p0c) * exp(-(Evt-Evc)/T) - p0c * (Nvt - p0t)  ;     //equilibrium term
      recomb_h *= _gamma; 
      if (_langevin)
        recomb_h *= (mupt + munc) * Constants::e * 100 / (_er * Constants::e0);
    }

    if (_adiabatic)
    {
      _saved_old_sim->get_solution(el, _old_erec, old_e_rec, this_m.get_coordinates());
      _saved_old_sim->get_solution(el, _old_hrec, old_h_rec, this_m.get_coordinates());

      //cout<<"old_e_rec "<<old_e_rec<<" old_h_rec "<<old_h_rec<<endl;
      cout<<"alpha = "<<_alpha<<endl;

      recomb_e *= _alpha;
      recomb_e += (1.0 - _alpha) * old_e_rec;
      recomb_h *= _alpha;
      recomb_h += (1.0 - _alpha) * old_h_rec;
    }
 
    //if ( recomb_e != 0.0 )
      //cout<<"recomb_e = "<<recomb_e<<endl;
    //if ( recomb_h != 0.0 )
      //cout<<"recomb_h = "<<recomb_h<<endl;
    //if ( old_e_rec != 0.0 )
      //cout<<"old_recomb_e = "<<old_e_rec<<endl;
    //if ( old_h_rec != 0.0 )
      //cout<<"old_recomb_h = "<<old_h_rec<<endl;
    
  }
  else
  {
    recomb_e = recomb_h = 0.0;
  }

}


void
HostGuestRecombination::get_net_recombination_rate_derivatives(std::vector<double>& recomb_e, std::vector<double>& recomb_h)
{

  if ( _coupled_sim->is_solved() && _coupled ) // && !(_this_sim->equilibrium() && _coupled_sim->equilibrium() ) )
  {
    const DriftDiffusionProperties& this_m = get_driftdiffusionproperties(); //host properties
    const Elem* el = this_m.get_element();

    DDBulkModel* coupled_m = _coupled_sim->get_bulk_model<DDBulkModel>(el);  //guest properties
    //coupled_m->set_coordinates(this_m.get_coordinates());

    if ( coupled_m == NULL )
    {
      recomb_e[0] = recomb_e[1] = recomb_h[0] = recomb_h[1] = 0.0;
      return;
    }

    double nt   = this_m.get_electron_density();
    double pt   = this_m.get_hole_density();
    double n0t  = this_m.get_quasi_equilibrium_n();
    double p0t  = this_m.get_quasi_equilibrium_p();
    double munt = this_m.get_electron_mobility();
    double mupt = this_m.get_hole_mobility();
    double Ect  = this_m.get_conduction_band().get_band_edge();
    double Evt  = this_m.get_valence_band().get_band_edge();
    double Nct  = this_m.get_conduction_band().get_effective_DOS();
    double Nvt  = this_m.get_valence_band().get_effective_DOS();
    double T    = this_m.get_lattice_temperature();
    double dnt_dphi = this_m.get_electron_density_derivative();
    double dpt_dphi = this_m.get_hole_density_derivative();
    double dmunt_dphi = this_m.get_electron_mobility_derivative_potential();
    double dmupt_dphi = this_m.get_hole_mobility_derivative_potential();
    
    double nc, pc, n0c, p0c, munc, mupc = 0.0;
    _coupled = false;

    _coupled_sim->get_solution(el, _eDensity, nc, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _hDensity, pc, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _eDensity0, n0c, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _hDensity0, p0c, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _eMobility, munc, this_m.get_coordinates());
    _coupled_sim->get_solution(el, _hMobility, mupc, this_m.get_coordinates());

    _coupled = true;

    double Ecc  = coupled_m->get_conduction_band().get_band_edge();
    double Evc  = coupled_m->get_valence_band().get_band_edge();
    double Ncc  = coupled_m->get_conduction_band().get_effective_DOS();
    double Nvc  = coupled_m->get_valence_band().get_effective_DOS();

    if (Ect > Ecc) {
      recomb_e[0] = ( pc + Ncc - nc + nc * exp(-(Ect-Ecc)/T) );   //( pc + Ncc - nc + nc * exp(-(Ect-Ecc)/T) )
      if (_langevin)
      {
        recomb_e[0] *= (munt + mupc);
        recomb_e[0] += dmunt_dphi / dnt_dphi * ( nt * pc +  nt * (Ncc - nc) - nc * (Nct - nt) * exp(-(Ect-Ecc)/T));  //(nt * pc +  nt * (Ncc - nc) - nc * (Nct - nt) * exp(-(Ect-Ecc)/T))
        //recomb_e[0] -= dmunt_dphi / dnt_dphi * ( n0t * p0c + n0t * (Ncc - n0c) - n0c * (Nct - n0t) * exp(-(Ect-Ecc)/T));   //equilibrium term
        recomb_e[0] *= Constants::e * 100 / (_er * Constants::e0);
      }
      recomb_e[0] *= _gamma;

      recomb_e[1] = 0.0;
      if (_langevin)
      {
        recomb_e[1] =  nt * pc + nt * (Ncc - nc) - nc * (Nct - nt) * exp(-(Ect-Ecc)/T);  //nt * pc +  nt * (Ncc - nc) - nc * (Nct - nt) * exp(-(Ect-Ecc)/T)
        //recomb_e[1] -= n0t * p0c + n0t * (Ncc - n0c) - n0c * (Nct - n0t) * exp(-(Ect-Ecc)/T);     //equilibrium term
        recomb_e[1] *= dmunt_dphi / dpt_dphi ;
        recomb_e[1] *= _gamma * Constants::e * 100 / (_er * Constants::e0);
      }
    }
    else {
      recomb_e[0] =  ( pc + (Ncc - nc) * exp(-(Ecc-Ect)/T) + nc )  ; //( pc + (Ncc - nc) * exp(-(Ecc-Ect)/T) + nc ) 
      if (_langevin)
      {
        recomb_e[0] *= (munt + mupc);
        recomb_e[0] += dmunt_dphi / dnt_dphi * ( nt * pc + nt * (Ncc - nc) * exp(-(Ecc-Ect)/T) - nc * (Nct - nt) ) ;  // ( nt * pc + nt * (Ncc - nc) * exp(-(Ecc-Ect)/T) - nc * (Nct - nt) ) 
        //recomb_e[0] -= dmunt_dphi / dnt_dphi * ( n0t * p0c + n0t * (Ncc - n0c) * exp(-(Ecc-Ect)/T) - n0c * (Nct - n0t) ) ;   //equilibrium term
        recomb_e[0] *= Constants::e * 100 / (_er * Constants::e0);
      }
      recomb_e[0] *= _gamma;

      recomb_e[1] = 0.0;
      if (_langevin)
      {
        recomb_e[1] =  nt * pc + nt * (Ncc - nc) * exp(-(Ecc-Ect)/T) - nc * (Nct - nt)  ; //nt * pc + nt * (Ncc - nc) * exp(-(Ecc-Ect)/T) - nc * (Nct - nt) 
        //recomb_e[1] -= n0t * p0c + n0t * (Ncc - n0c) * exp(-(Ecc-Ect)/T) - n0c * (Nct - n0t)  ;   //equilibrium term
        recomb_e[1] *= dmunt_dphi / dpt_dphi ;
        recomb_e[1] *= _gamma * Constants::e * 100 / (_er * Constants::e0);
      }
    }

    if (Evt < Evc) {
      recomb_h[0] = 0.0;
      if (_langevin)
      {
        recomb_h[0] =  pt * nc + pt * (Nvc - pc) - pc * (Nvt - pt) * exp(-(Evc-Evt)/T)  ;  //pt * nc + pt * (Nvc - pc) - pc * (Nvt - pt) * exp(-(Evc-Evt)/T)
        //recomb_h[0] -= p0t * n0c + p0t * (Nvc - p0c) - p0c * (Nvt - p0t) * exp(-(Evc-Evt)/T)  ;   //equilibrium term
        recomb_h[0] *= dmupt_dphi / dnt_dphi ;
        recomb_h[0] *= _gamma * Constants::e * 100 / (_er * Constants::e0);
       }

      recomb_h[1] = ( nc + Nvc - pc + pc * exp(-(Evc-Evt)/T) ) ; //( nc + Nvc - pc + pc * exp(-(Evc-Evt)/T) ) 
      if (_langevin) 
      {
        recomb_h[1] *= (mupt + munc) ;
        recomb_h[1] += dmupt_dphi / dpt_dphi * ( pt * nc + pt * (Nvc - pc) - pc * (Nvt - pt) * exp(-(Evc-Evt)/T) ) ; //( pt * nc + pt * (Nvc - pc) - pc * (Nvt - pt) * exp(-(Evc-Evt)/T) )
        //recomb_h[1] -= dmupt_dphi / dpt_dphi * ( p0t * n0c + p0t * (Nvc - p0c) - p0c * (Nvt - p0t) * exp(-(Evc-Evt)/T) ) ;   //equilibrium term
        recomb_h[1] *= Constants::e * 100 / (_er * Constants::e0);
      }
      recomb_h[1] *= _gamma;
    }
    else {
      recomb_h[0] = 0.0;
      if (_langevin)
      {
        recomb_h[0] =  pt * nc + pt * (Nvc - pc) * exp(-(Evt-Evc)/T) - pc * (Nvt - pt)  ; //pt * nc + pt * (Nvc - pc) * exp(-(Evt-Evc)/T) - pc * (Nvt - pt) 
        //recomb_h[0] -= p0t * n0c + p0t * (Nvc - p0c) * exp(-(Evt-Evc)/T) - p0c * (Nvt - p0t)  ;   //equilibrium term
        recomb_h[0] *= dmupt_dphi / dnt_dphi ;
        recomb_h[0] *= _gamma * Constants::e * 100 / (_er * Constants::e0);
      }

      recomb_h[1] = ( nc + (Nvc - pc) * exp(-(Evt-Evc)/T) + pc ) ;   //( nc + (Nvc - pc) * exp(-(Evt-Evc)/T) + pc )
      if (_langevin)
      {
        recomb_h[1] *= (mupt + munc) ;
        recomb_h[1] += dmupt_dphi / dpt_dphi * ( pt * nc + pt * (Nvc - pc) * exp(-(Evt-Evc)/T) - pc * (Nvt - pt) ) ;  //( pt * nc + pt * (Nvc - pc) * exp(-(Evt-Evc)/T) - pc * (Nvt - pt) ) 
        //recomb_h[1] -= dmupt_dphi / dpt_dphi * ( p0t * n0c + p0t * (Nvc - p0c) * exp(-(Evt-Evc)/T) - p0c * (Nvt - p0t) ) ;    //equilibrium term
        recomb_h[1] *= Constants::e * 100 / (_er * Constants::e0);
      }
      recomb_h[1] *= _gamma ;
    }

    if (_adiabatic)
    {
      recomb_e[0] *= _alpha;
      recomb_e[1] *= _alpha;
      recomb_h[0] *= _alpha;
      recomb_h[1] *= _alpha;
    }

  }
  else
  {
    recomb_e[0] = recomb_e[1] = recomb_h[0] = recomb_h[1] = 0.0;
  }
}

void
HostGuestRecombination::do_reinit(void)
{

}
