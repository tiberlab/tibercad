#ifndef _EXCPOLPROPS_H_
#define _EXCPOLPROPS_H_

#include "tiber_config.h"
#include "SimulationInterface.h"
#include "elem.h"

#include <string>
#include <map>

class ExcPolProps {
  public:
#ifdef HAVE_CONSTEXPR
    //constexpr static double a0 = 1e3;
    //constexpr static double b0 = 1e-9;
    constexpr static double pol_tau0 = 1e-11;
    constexpr static double my_KB = 1.3806488e-23;
    constexpr static double pn_lim = 1e-10;
#else
    //static const double a0 = 1e3;
    //static const double b0 = 1e-9;
    static const double pol_tau0 = 1e-11;
    static const double my_KB = 1.3806488e-23;
    static const double pn_lim = 1e-10;
#endif
    
    double b_final;
    double a_final;

    double pol_tau;
    double pol_energy;
    double density_renormalization;
    double ren_alpha;
    double int_g;
    double int_f;

    double V;
    double T;
    double dE_in_kT;
    double p;

    double lwell;

    double max_density;
    SimulationInterface* maxwell;
    ID hopID;

    std::map<ID, double> hopCache;

    ExcPolProps() {
      pol_tau = pol_tau0;
      pol_energy = 0;

      density_renormalization = 1;
      ren_alpha = 1;

      b_final = 0;
      a_final = 0;

      int_f = 0;
      int_g = 0;
      maxwell = NULL;
      V = 0;

      p = 0;
      max_density = 0;
    }

    void set(double temperature) {
      max_density = 0;
      T = temperature;

      maxwell = SimulationInterface::find_simulation("maxwell");
      hopID = maxwell->get_solution_id("XHopfield");
      if (!maxwell->is_solved()) {
        return;
      }

      ID id1 = maxwell->get_solution_id("EigenValue");
      ID id2 = maxwell->get_solution_id("WPolaritonImag");
      ID id3 = maxwell->get_solution_id("WPolariton");

      std::map<ID, std::vector<double> > map;
      map.insert(std::make_pair(id1, std::vector<double>(1)));
      map.insert(std::make_pair(id2, std::vector<double>(1)));
      map.insert(std::make_pair(id3, std::vector<double>(1)));
      maxwell->get_solution(map);

      pol_energy = map[id3][0] * Constants::hbar;
      // 1 - exp(dE)*1/(1-exp(-pn))
      double Ex0 = maxwell->get_options().get_option("Wexc", 1.0) * Constants::e;

      dE_in_kT = std::exp((pol_energy -  Ex0) / my_KB / T);
      double Me = Constants::me * 0.2;//TODO ???? 0.13 0.27 0.2

      lwell = maxwell->get_options().get_option("lwell", 3) * 1.0e-7;// in maxwell in nm, here in cm

      double Swell = 1e-6;
      double rootS = 1e-3;

      p = (2 * M_PI * Constants::hbar * Constants::hbar) /(Me * my_KB * T) * 1.0e4 * lwell; // m^-2 -> cm^-2

      pol_tau = std::abs(map[id2][0]) == 0 ? 1e-11 : 1.0 / std::abs(map[id2][0]);

      V = std::abs(map[id1][0] - map[id3][0]) * Constants::hbar;

      std::cout << "lwell " << lwell << "\n";
      std::cout << "Polariton life time: " << pol_tau << " V: " << V << " dE: " << dE_in_kT << " p (no lwell): " << p/lwell << " Ex0 " << Ex0 << " pol_energy " << pol_energy << "\n";
      flush(std::cout);

      hopCache.clear();

      double ae = -4.08 * Constants::e;
      double ah = 2.1 * Constants::e;
      double ro = 6.15e3;
      double u = 6.56e3;
      double a_b_GaN = 3e-9;
      double E_b_GaN = 45e-3; //in eV

      double MeeS = 6 * E_b_GaN * Constants::e * a_b_GaN * a_b_GaN;
      double as = (Ex0 - pol_energy + my_KB * T) * M_PI * (ae - ah) * (ae - ah) / ro / Constants::hbar / Constants::hbar / u / u / u * 1e4; // m2 -> cm2
      double bs = M_PI / 4 / Constants::hbar * MeeS * MeeS / my_KB / T * std::exp(-2 * (Ex0 - pol_energy) / my_KB / T) * 1e8; // m4 -> cm4

      std::cout << "a & b: " << as << " " << bs << " MeeS " << MeeS << " exp value:" << (-2 * (Ex0 - pol_energy) / my_KB / T)<< "\n";
      //double lwell = 4e-7; // in cm
      // 1d -> 1
      // 2d -> root(S)
      // 3d -> S

      b_final = bs * lwell / rootS;
      a_final = as / rootS;
    }

    // The aim is to calculate density_renormalization

/*    void calculate(const NumericVector<Number>& x, double one_integral, double density_integral, double density_square_integral, double T) {
      double Me = Constants::me * 0.2;//TODO ???? 0.13 0.27 0.2
      double myKB = 1.3806488e-23;

      double as = 1;
      double bs = 1e-9;//1e-9;
      //double p = 0.3; //   exp(-V/kT)*m*kT/(2pih^2)
      double p = 0;//std::exp(-V/myKB/T) * Me * myKB * T / (2 * M_PI * Constants::hbar * Constants::hbar) / 1.0e4; // m^-2 -> cm^-2

      double lwell = 4e-7; // in cm
      // 1d -> 1
      // 2d -> root(S)
      // 3d -> S
      double Swell = 1e-3; // in cm 10mkm x 10mkm

      b1_final = bs * lwell / Swell;
      a1_final_1 = as / Swell;

      // shit ?
      a1_final_2 = - bs / Swell / lwell * p;
      c1_final = - as / Swell / lwell * p;

      b2_final = b1_final;
      a2_final = as / Swell;

      //find DR and J
      //J = pol_tau * (b2_final / density_renormalization / density_renormalization * density_square_integral + a2_final / density_renormalization * density_integral);

      // Now we want to solve aa * x^2 + bb * x + cc = 0

      double aa = 1 - pol_tau * c1_final * one_integral;
      double bb = -1 - pol_tau * (a1_final_1 + a1_final_2) * density_integral;
      double cc = - pol_tau * b1_final * density_square_integral;


      //density_renormalization = (-bb + std::sqrt(bb*bb - 4 * aa * cc)) / 2 / aa;


      std::cout << "one_integral= " << one_integral << "\n";
      std::cout << "density_integral= " << density_integral << "\n";
      std::cout << "density_square_integral= " << density_square_integral << "\n";
      std::cout << "aa= " << aa << "\n";
      std::cout << "bb= " << bb << "\n";
      std::cout << "cc= " << cc << "\n";

      //std::cout << "DR= " << density_renormalization << "\n";
      //flush(std::cout);
    }*/

    inline double getXHopfield(const Elem* elem, const Point& point, bool local) {
      if (maxwell == NULL || !maxwell->is_solved()) {
        return 0;
      } else {
        ID subdomain = elem->subdomain_id();

        if (hopCache.count(subdomain)) {
          return hopCache[subdomain];
        } else {
          double x = 0;
          //Utils::Timer tt;
          maxwell->get_solution(elem, hopID, x, point, local);
          //std::cout << "Hop time: "  << tt.elapsed_string() << "\n";
          hopCache[subdomain] = x;
          return x;
        }
      }
    }

/*    double get_phonon_scattering(double density, const Elem* elem, const Point& point) {
      double XHop = getXHopfield(elem, point, false);

      return XHop * density_renormalization * ( (a1_final_1 * density + c1_final * density_renormalization) * int_f   + a2_final / density_renormalization * density);
    }

    double get_exc_exc_scattering(double density, const Elem* elem, const Point& point) {
      double XHop = getXHopfield(elem, point, false);

      return XHop * density_renormalization * ( (b1_final * density * density / density_renormalization + a1_final_2 * density) * int_f + b2_final * density * density / density_renormalization / density_renormalization);
    }

    double get_scattering_derivative(double density, const Elem* elem, const Point& point) {
      double XHop = getXHopfield(elem, point, false);

      return XHop * density_renormalization * ( (a1_final_1 + a1_final_2 + 2 * b1_final * density / density_renormalization) * int_f   +
                                                a2_final / density_renormalization + 2 * b2_final * density / density_renormalization / density_renormalization);


    }
*/
/*
    void findDR() {
      double DR = density_renormalization;

      // solve 1 - DR = F(DR*n)
      double F = f(DR);

      for (int i = 0; i < 30; i++) {
        if (F*F/(F*F + (1 - DR)*(1 - DR)) < )
      }
    }
*/

    // div on m_density
    inline double div_Gpol(double m_density, const Elem* elem, const Point& point) {
      return div_f(m_density, elem, point) / density_renormalization + div_g(m_density, elem, point) * pol_tau * int_f * ren_alpha;
    }

    inline double Gpol_exc(double m_density, const Elem* elem, const Point& point) {
      return f_exc(m_density, elem, point) + g_exc(m_density, elem, point) * pol_tau * int_f * density_renormalization * ren_alpha;
    }

    inline double Gpol_phon(double m_density, const Elem* elem, const Point& point) {
      return f_phon(m_density, elem, point) + g_phon(m_density, elem, point) * pol_tau * int_f * density_renormalization * ren_alpha;
    }

    inline double f(double m_density, const Elem* elem, const Point& point) {
      return (f_exc(m_density, elem, point) + f_phon(m_density, elem, point));
    }

    inline double g(double m_density, const Elem* elem, const Point& point) {
      return (g_exc(m_density, elem, point) + g_phon(m_density, elem, point));
    }

    inline double f_exc(double m_density, const Elem* elem, const Point& point) {
      double n_density = m_density / density_renormalization;
      double XHop = getXHopfield(elem, point, false);

      return XHop * b_final * n_density * n_density;
    }

    inline double f_phon(double m_density, const Elem* elem, const Point& point) {
      double n_density = m_density / density_renormalization;
      double XHop = getXHopfield(elem, point, false);

      return XHop * a_final * n_density;
    }

    inline double div_f(double m_density, const Elem* elem, const Point& point) {
      double n_density = m_density / density_renormalization;
      double XHop = getXHopfield(elem, point, false);

      return XHop * ( 2 * b_final * n_density + a_final);
    }

    inline double g_exc(double m_density, const Elem* elem, const Point& point) {
      double n_density = m_density / density_renormalization;
      double XHop = getXHopfield(elem, point, false);

      if (dE_in_kT < 1e-20) {
        return f_exc(m_density, elem, point);
      } else {
        if (p*n_density < pn_lim) {
          return f_exc(m_density, elem, point) * (1 - dE_in_kT / p / n_density);
        } else {
          return f_exc(m_density, elem, point) * (1 - dE_in_kT/(1 - std::exp(-p*n_density)));
        }
      }
    }

    inline double g_phon(double m_density, const Elem* elem, const Point& point) {
      double n_density = m_density / density_renormalization;
      double XHop = getXHopfield(elem, point, false);

      if (dE_in_kT < 1e-20) {
        return f_phon(m_density, elem, point);
      } else {
        if (p*n_density < pn_lim) {
          return f_phon(m_density, elem, point) * (1 - dE_in_kT/p/n_density);
        } else {
          return f_phon(m_density, elem, point) * (1 - dE_in_kT/(1 - std::exp(-p*n_density)));
        }
      }
    }

    inline double div_g(double m_density, const Elem* elem, const Point& point) {
      double n_density = m_density / density_renormalization;
      double XHop = getXHopfield(elem, point, false);

      if (dE_in_kT < 1e-20) {
        return div_f(m_density, elem, point);
      } else {
        if (p*n_density < pn_lim) {
          return div_f(m_density, elem, point) * (1 - dE_in_kT/p/n_density) +
		  f(m_density, elem, point) * dE_in_kT/p/n_density/n_density;
        } else {
          return div_f(m_density, elem, point) * (1 - dE_in_kT/(1 - std::exp(-p*n_density))) + f(m_density, elem, point) * dE_in_kT/(1 - std::exp(-p*n_density))/(1 - std::exp(-p*n_density)) * p * std::exp(- p * n_density);
        }
      }
    }
};





#endif /* _EXCPOLPROPS_H_*/
