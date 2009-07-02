// $Id$

#include "DopingDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"



TIBER_MODULE(DopingDependentMobility, doping_dependent)




void
DopingDependentMobility::read_database(void)
{
  Database& db = get_database();
  db.set_section("mobility/doping_dependent");

  formula_ = db.get("mobility_formula", formula_, true);

  std::vector<double> empty(2, 0);

  if (formula_ == 1)
  {
    // Model of Masetti et al.

    std::vector<double> data(empty);
    db.get("mumin1", data);
    mumin_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("mumin2", data);
    am_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("mu1", data);
    mud_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = std::vector<double>(2, 1);
    db.get("Cr", data);
    ad_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("Cs", data);
    N0_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("alpha", data);
    an_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("beta", data);
    a_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("Pc", data);
    aa_ =  get_carrier_type() == 'e' ? data[0] : data[1];
  }
  else
  {
    // Model of Arora et al.

    std::vector<double> data(empty);
    db.get("mumin", data, true);
    mumin_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("am", data);
    am_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("mud", data);
    mud_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("ad", data);
    ad_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = std::vector<double>(2, 1);
    db.get("N0", data);
    N0_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("aN", data);
    an_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("A", data);
    a_ =  get_carrier_type() == 'e' ? data[0] : data[1];

    data = empty;
    db.get("aA", data);
    aa_ =  get_carrier_type() == 'e' ? data[0] : data[1];
  }

}



void
DopingDependentMobility::create_submodels(void)
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
DopingDependentMobility::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  const DopingDependentMobility* scA =
    dynamic_cast<const DopingDependentMobility*>(comp_A);
  const DopingDependentMobility* scB =
    dynamic_cast<const DopingDependentMobility*>(comp_B);

  // formula should be the same, so we make some sanity check
  if (scA->formula_ != scB->formula_)
    throw InitFailedException("Doping dependent mobility has to use the same "
        "formula for both components of the alloy " + get_material()->get_name());
  formula_ = scA->formula_;
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
    assert(scA->const_mob_ != NULL);
    assert(scB->const_mob_ != NULL);
    destroy(const_mob_);
    const_mob_ = create_submodel_copy(scA->const_mob_);
    const_mob_->set_driftdiffusionproperties(&get_driftdiffusionproperties());
    const_mob_->init_alloy(scA->const_mob_, scB->const_mob_, xa);
  }
}

