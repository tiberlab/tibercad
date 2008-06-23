// $Id$

#include "DopingDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "getpot.h"


TIBER_MODULE(DopingDependentMobility, doping_dependent)




void
DopingDependentMobility::read_database(void)
{
  const Material* mat = get_material();
  GetPot data((mat->get_database()).get_data_file());

  formula_ = data("mobility_formula", formula_);
  
  if (formula_ == 1)
  {
    // Model of Masetti et al.

    std::string s("mumin1_");
    s += get_carrier_type();
    mumin_ = data(s.c_str(), mumin_);

    s = "mumin2_";
    s += get_carrier_type();
    am_ = data(s.c_str(), am_);

    s = "mu1_";
    s += get_carrier_type();
    mud_ = data(s.c_str(), mud_);

    s = "Cr_";
    s += get_carrier_type();
    ad_ = data(s.c_str(), ad_);

    s = "Cs_";
    s += get_carrier_type();
    N0_ = data(s.c_str(), N0_);

    s = "alpha_";
    s += get_carrier_type();
    an_ = data(s.c_str(), an_);

    s = "beta_";
    s += get_carrier_type();
    a_ = data(s.c_str(), a_);

    s = "Pc_";
    s += get_carrier_type();
    aa_ = data(s.c_str(), aa_);
  }
  else
  {
    // Model of Arora et al.

    std::string s("mumin_");
    s += get_carrier_type();
    mumin_ = data(s.c_str(), mumin_);

    s = "am_";
    s += get_carrier_type();
    am_ = data(s.c_str(), am_);

    s = "mud_";
    s += get_carrier_type();
    mud_ = data(s.c_str(), mud_);

    s = "ad_";
    s += get_carrier_type();
    ad_ = data(s.c_str(), ad_);

    s = "N0_";
    s += get_carrier_type();
    N0_ = data(s.c_str(), N0_);

    s = "an_";
    s += get_carrier_type();
    an_ = data(s.c_str(), an_);

    s = "a_";
    s += get_carrier_type();
    a_ = data(s.c_str(), a_);

    s = "aa_";
    s += get_carrier_type();
    aa_ = data(s.c_str(), aa_);
  }

}



void
DopingDependentMobility::do_init(void)
{
  if (formula_ == 1)
  {
    const_mob_ = MobilityModelInterface::create("constant");
    if (const_mob_ == NULL)
    {
      std::string msg("DopingDependentMobility: Could not ");
      msg += "create constant mobility model needed for formula of Masetti.";
      throw InitFailedException(msg);
    }

    const_mob_->set_driftdiffusionproperties(&get_driftdiffusionproperties());
    const_mob_->set_carrier_type(get_carrier_type());
    const_mob_->set_material(get_material());
    const_mob_->init();
  }

}



double
DopingDependentMobility::get_mobility(void)
{
  double mu;
  double T = get_driftdiffusionproperties().get_lattice_temperature() / T0;
  double N = get_material()->get_total_doping_density();

  if (formula_ == 1)
  {
    assert(const_mob_ != NULL);
    double mu_const = const_mob_->get_mobility();
    double mumin1 = mumin_;
    double Pc = aa_;
    double mumin2 = am_;
    double Cr = ad_;
    double alpha = an_;
    double mu1 = mud_;
    double Cs = N0_;
    double beta = a_;

    N = (N > 1.0) ? N : 1.0;
    // Model of Masetti et al.
    mu = mumin1 * std::exp(-Pc / N);
    mu += (mu_const - mumin2) / (1.0 + std::pow(N / Cr, alpha));
    mu -= mu1 / (1.0 + std::pow(Cs / N, beta));
  }
  else
  {
    // Model of Arora et al.
    double muminA = mumin_ * std::pow(T, am_);
    double mudA = mud_ * std::pow(T, ad_);
    double N00 = N0_ * std::pow(T, an_);
    double aa = a_ * std::pow(T, aa_);
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
  mumin_ = mod->mumin_;
  am_ = mod->am_;
  mud_ = mod->mud_;
  ad_ = mod->ad_;
  N0_ = mod->N0_;
  an_ = mod->an_;
  a_ = mod->a_;
  aa_ = mod->aa_;
}



void
DopingDependentMobility::calculate_VCA(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const DopingDependentMobility* scA =
    dynamic_cast<const DopingDependentMobility*>(comp_A);
  const DopingDependentMobility* scB =
    dynamic_cast<const DopingDependentMobility*>(comp_B);

  mumin_ = alloy(scA->mumin_, scB->mumin_ , xa);
  am_ = alloy(scA->am_, scB->am_ , xa);
  mud_ = alloy(scA->mud_, scB->mud_ , xa);
  ad_ = alloy(scA->ad_, scB->ad_ , xa);
  N0_ = alloy(scA->N0_, scB->N0_ , xa);
  an_ = alloy(scA->an_, scB->an_ , xa);
  a_ = alloy(scA->a_, scB->a_ , xa);
  aa_ = alloy(scA->aa_, scB->aa_ , xa);
  
  if (formula_ == 1)
  {
    assert(const_mob_ != NULL);
    assert(scA->const_mob_ != NULL);
    assert(scB->const_mob_ != NULL);
    const_mob_->build_alloy(scA->const_mob_, scB->const_mob_, xa);
  }
}

