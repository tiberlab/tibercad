// $Id$

#include "DopingDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "getpot.h"

void
DopingDependentMobility::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  _formula = data("mobility_formula", _formula);
  
  if (_formula == 1)
  {
    // Model of Masetti et al.

    std::string s("mumin1_");
    s += get_carrier_type();
    _mumin = data(s.c_str(), _mumin);

    s = "mumin2_";
    s += get_carrier_type();
    _am = data(s.c_str(), _am);

    s = "mu1_";
    s += get_carrier_type();
    _mud = data(s.c_str(), _mud);

    s = "Cr_";
    s += get_carrier_type();
    _ad = data(s.c_str(), _ad);

    s = "Cs_";
    s += get_carrier_type();
    _N0 = data(s.c_str(), _N0);

    s = "alpha_";
    s += get_carrier_type();
    _an = data(s.c_str(), _an);

    s = "beta_";
    s += get_carrier_type();
    _a = data(s.c_str(), _a);

    s = "Pc_";
    s += get_carrier_type();
    _aa = data(s.c_str(), _aa);
  }
  else
  {
    // Model of Arora et al.

    std::string s("mumin_");
    s += get_carrier_type();
    _mumin = data(s.c_str(), _mumin);

    s = "am_";
    s += get_carrier_type();
    _am = data(s.c_str(), _am);

    s = "mud_";
    s += get_carrier_type();
    _mud = data(s.c_str(), _mud);

    s = "ad_";
    s += get_carrier_type();
    _ad = data(s.c_str(), _ad);

    s = "N0_";
    s += get_carrier_type();
    _N0 = data(s.c_str(), _N0);

    s = "an_";
    s += get_carrier_type();
    _an = data(s.c_str(), _an);

    s = "a_";
    s += get_carrier_type();
    _a = data(s.c_str(), _a);

    s = "aa_";
    s += get_carrier_type();
    _aa = data(s.c_str(), _aa);
  }

}



void
DopingDependentMobility::do_init(void)
{
  if (_formula == 1)
  {
    _const_mob = MobilityModelInterface::create("constant");
    if (_const_mob == NULL)
    {
      std::string msg("DopingDependentMobility: Could not ");
      msg += "create constant mobility model needed for formula of Masetti.";
      throw InitFailedException(msg);
    }

    _const_mob->set_driftdiffusionproperties(&get_driftdiffusionproperties());
    _const_mob->set_carrier_type(get_carrier_type());
    _const_mob->set_material(get_material());
    _const_mob->init();
  }

}



double
DopingDependentMobility::get_mobility(void)
{
  double mu;
  double T = get_driftdiffusionproperties().get_lattice_temperature() / T0;
  double N = get_material()->get_total_doping_density();

  if (_formula == 1)
  {
    assert(_const_mob != NULL);
    double mu_const = _const_mob->get_mobility();
    double mumin1 = _mumin;
    double Pc = _aa;
    double mumin2 = _am;
    double Cr = _ad;
    double alpha = _an;
    double mu1 = _mud;
    double Cs = _N0;
    double beta = _a;

    N = (N > 1.0) ? N : 1.0;
    // Model of Masetti et al.
    mu = mumin1 * std::exp(-Pc / N);
    mu += (mu_const - mumin2) / (1.0 + std::pow(N / Cr, alpha));
    mu -= mu1 / (1.0 + std::pow(Cs / N, beta));
  }
  else
  {
    // Model of Arora et al.
    double muminA = _mumin * std::pow(T, _am);
    double mudA = _mud * std::pow(T, _ad);
    double N00 = _N0 * std::pow(T, _an);
    double aa = _a * std::pow(T, _aa);
    mu = muminA + mudA / (1.0 + std::pow(N / N00, aa));
  }
  return mu;
}



void
DopingDependentMobility::get_mobility_derivatives(std::vector<double>& dm)
{
  dm[0] = dm[1] = dm[2] = 0.0;
}


void
DopingDependentMobility::copy_from(const PhysicalModelInterface* rhs)
{
  MobilityModelInterface::copy_from(rhs);

  const DopingDependentMobility* mod =
    dynamic_cast<const DopingDependentMobility*>(rhs);
  _mumin = mod->_mumin;
  _am = mod->_am;
  _mud = mod->_mud;
  _ad = mod->_ad;
  _N0 = mod->_N0;
  _an = mod->_an;
  _a = mod->_a;
  _aa = mod->_aa;
}



void
DopingDependentMobility::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const DopingDependentMobility* scA =
    dynamic_cast<const DopingDependentMobility*>(comp_A);
  const DopingDependentMobility* scB =
    dynamic_cast<const DopingDependentMobility*>(comp_B);

  _mumin = alloy(scA->_mumin, scB->_mumin , xa);
  _am = alloy(scA->_am, scB->_am , xa);
  _mud = alloy(scA->_mud, scB->_mud , xa);
  _ad = alloy(scA->_ad, scB->_ad , xa);
  _N0 = alloy(scA->_N0, scB->_N0 , xa);
  _an = alloy(scA->_an, scB->_an , xa);
  _a = alloy(scA->_a, scB->_a , xa);
  _aa = alloy(scA->_aa, scB->_aa , xa);
}

