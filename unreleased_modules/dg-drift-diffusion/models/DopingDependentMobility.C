// $Id: DopingDependentMobility.C 4268 2016-10-21 15:09:28Z maufder $

#include "DopingDependentMobility.h"
#include "DriftDiffusionProperties.h"

#include "Material.h"
#include "Database.h"

#include "TiberModule.h"





void
DopingDependentMobility::read_database(void)
{
  const Database& db = get_database();
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
DopingDependentMobility::prepare_submodels(void)
{
  if (formula_ == 1)
  {
    ModelOptions opts;
    opts.set_option("type", "constant");
    opts.set_option("particle", get_option("particle", "electron"));

    const_mob_ = MobilityModelInterface::create("constant", get_material(), opts);
    if (const_mob_ == NULL)
    {
      std::string msg("DopingDependentMobility: Could not ");
      msg += "create constant mobility model needed for formula of Masetti.";
      throw InitFailedException(msg);
    }

    add_submodel("const_mobility", const_mob_);
  }
}



void
DopingDependentMobility::do_init(void)
{
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
DopingDependentMobility::get_derivative_grad_potential(libMesh::RealGradient& dm)
{
  dm.zero();
}

void
DopingDependentMobility::get_derivative_grad_fermi(libMesh::RealGradient& dm)
{
  dm.zero();
}


