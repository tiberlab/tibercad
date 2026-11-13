/*
 * This file is part of tiberCAD.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tiberCAD. If not, see <https://www.gnu.org/licenses/>.
 */

/*!
 * \file LatticeMatch.cpp
 * \brief tiberCAD API implementation.
 */

#include "tibercad/math/LatticeMatch.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <tuple>

using namespace LatticeMatch;

namespace
{
  //----------------------------------------------------------------------
  // small shared helpers
  //----------------------------------------------------------------------

  double spectral_norm(const Eigen::Matrix3d& M)
  {
    // induced 2-norm = largest singular value = sqrt of largest eigenvalue of M^T*M
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(M.transpose() * M, Eigen::EigenvaluesOnly);
    return std::sqrt(std::max(0.0, es.eigenvalues()(2)));
  }

  void polar(const Eigen::Matrix3d& F, Eigen::Matrix3d& R, Eigen::Matrix3d& U)
  {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(F.transpose() * F);
    Eigen::Vector3d w = es.eigenvalues().cwiseMax(0.0).cwiseSqrt();
    U = es.eigenvectors() * w.asDiagonal() * es.eigenvectors().transpose();
    R = F * U.inverse();
  }

  double rot_angle_deg(const Eigen::Matrix3d& R)
  {
    double c = (R.trace() - 1.0) / 2.0;
    c = std::min(1.0, std::max(-1.0, c));
    return std::acos(c) * 180.0 / M_PI;
  }

  //! Strain description of a deformation gradient F (superstrate -> substrate)
  Candidate describe(const Eigen::Matrix3d& F)
  {
    Candidate d;
    d.F = F;
    d.H_norm = spectral_norm(F - Eigen::Matrix3d::Identity());
    polar(F, d.R, d.U);
    d.biot = d.U - Eigen::Matrix3d::Identity();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(d.biot, Eigen::EigenvaluesOnly);
    d.principal = es.eigenvalues();

    d.rot_deg = rot_angle_deg(d.R);
    d.hydro = std::cbrt(F.determinant()) - 1.0;
    return d;
  }

  //----------------------------------------------------------------------
  // lattice reduction: LLL and Buerger
  //----------------------------------------------------------------------

  void gram_schmidt(const Eigen::Matrix3d& B, Eigen::Matrix3d& Bs, Eigen::Matrix3d& mu)
  {
    Bs = B;
    mu.setZero();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < i; ++j)
      {
        mu(i, j) = B.row(i).dot(Bs.row(j)) / Bs.row(j).dot(Bs.row(j));
        Bs.row(i) -= mu(i, j) * Bs.row(j);
      }
  }

  //! Lenstra-Lenstra-Lovasz lattice basis reduction (rows = basis vectors)
  Eigen::Matrix3d lll(Eigen::Matrix3d B, double delta = 0.99)
  {
    Eigen::Matrix3d Bs, mu;
    gram_schmidt(B, Bs, mu);

    int k = 1;
    while (k < 3)
    {
      for (int j = k - 1; j >= 0; --j)
      {
        double q = std::round(mu(k, j));
        if (q != 0.0)
        {
          B.row(k) -= q * B.row(j);
          gram_schmidt(B, Bs, mu);
        }
      }

      if (Bs.row(k).dot(Bs.row(k)) >=
          (delta - mu(k, k - 1) * mu(k, k - 1)) * Bs.row(k - 1).dot(Bs.row(k - 1)))
      {
        ++k;
      }
      else
      {
        B.row(k).swap(B.row(k - 1));
        gram_schmidt(B, Bs, mu);
        k = std::max(k - 1, 1);
      }
    }
    return B;
  }

  //! All integer vectors in {-3,...,3}^3 except the zero vector (342 of them)
  const std::vector<Eigen::Vector3i>& buerger_coefficients(void)
  {
    static const std::vector<Eigen::Vector3i> coef = []
    {
      std::vector<Eigen::Vector3i> v;
      v.reserve(342);
      for (int a = -3; a <= 3; ++a)
        for (int b = -3; b <= 3; ++b)
          for (int c = -3; c <= 3; ++c)
            if ((a != 0) || (b != 0) || (c != 0))
              v.emplace_back(a, b, c);
      return v;
    }();
    return coef;
  }

  //! Shortest right-handed basis (lengths = successive minima in 3D)
  Eigen::Matrix3d buerger(Eigen::Matrix3d B)
  {
    B = lll(B);
    double vol = std::abs(B.determinant());

    const auto& coef = buerger_coefficients();
    std::vector<Eigen::Vector3d> V(coef.size());
    for (std::size_t i = 0; i < coef.size(); ++i)
      V[i] = B.transpose() * coef[i].cast<double>();

    std::vector<std::size_t> order(V.size());
    for (std::size_t i = 0; i < order.size(); ++i)
      order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&V](std::size_t a, std::size_t b)
        { return V[a].squaredNorm() < V[b].squaredNorm(); });

    Eigen::Vector3d b0 = V[order[0]];

    Eigen::Vector3d b1 = Eigen::Vector3d::Zero();
    for (std::size_t k = 1; k < order.size(); ++k)
    {
      const Eigen::Vector3d& v = V[order[k]];
      if (b0.cross(v).norm() > 1e-8 * b0.dot(b0))
      {
        b1 = v;
        break;
      }
    }

    Eigen::Vector3d b2 = Eigen::Vector3d::Zero();
    for (std::size_t k = 1; k < order.size(); ++k)
    {
      const Eigen::Vector3d& v = V[order[k]];
      Eigen::Matrix3d M;
      M.row(0) = b0; M.row(1) = b1; M.row(2) = v;
      if (std::abs(std::abs(M.determinant()) - vol) < 1e-6 * vol)
      {
        b2 = v;
        break;
      }
    }

    Eigen::Matrix3d out;
    out.row(0) = b0; out.row(1) = b1; out.row(2) = b2;
    if (out.determinant() < 0)
      out.row(2) = -out.row(2);
    return out;
  }

  //----------------------------------------------------------------------
  // Hermite normal forms and sublattices
  //----------------------------------------------------------------------

  //! All index-N sublattices of Z^3 (row Hermite normal forms, each once)
  void hnf3(int N, std::vector<Eigen::Matrix3i>& out)
  {
    out.clear();
    for (int a = 1; a <= N; ++a)
    {
      if (N % a != 0)
        continue;
      int Na = N / a;
      for (int d = 1; d <= Na; ++d)
      {
        if (Na % d != 0)
          continue;
        int f = Na / d;
        for (int b = 0; b < d; ++b)
          for (int c = 0; c < f; ++c)
            for (int e = 0; e < f; ++e)
            {
              Eigen::Matrix3i H;
              H << a, b, c,
                   0, d, e,
                   0, 0, f;
              out.push_back(H);
            }
      }
    }
  }

  struct Sublattice
  {
    unsigned int N;
    Eigen::Matrix3d B;    // Buerger-reduced basis (rows)
    Eigen::Matrix3i M;    // integer coefficients: B ~= M * A
    Eigen::Vector3d mins; // normalized successive minima
  };

  std::vector<Sublattice> sublattices(const Eigen::Matrix3d& A, unsigned int Nmax)
  {
    Eigen::Matrix3d Ainv = A.inverse();
    double V0 = std::abs(A.determinant());

    std::vector<Sublattice> result;
    std::vector<Eigen::Matrix3i> Hs;
    for (unsigned int n = 1; n <= Nmax; ++n)
    {
      hnf3((int)n, Hs);
      for (const auto& H : Hs)
      {
        Eigen::Matrix3d b = buerger(H.cast<double>() * A);

        Eigen::Matrix3d Mf = b * Ainv;
        Eigen::Matrix3i M;
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            M(i, j) = (int)std::lround(Mf(i, j));

        Sublattice s;
        s.N = n;
        s.B = b;
        s.M = M;
        double denom = std::cbrt((double)n * V0);
        for (int i = 0; i < 3; ++i)
          s.mins(i) = b.row(i).norm() / denom;

        result.push_back(s);
      }
    }
    return result;
  }

  //! r's cells are sublattices of o's on both sides -> o is the finer match
  bool is_super_of(const Candidate& r, const Candidate& o)
  {
    Eigen::Matrix3d X = r.M1.cast<double>() * o.M1.cast<double>().inverse();
    Eigen::Matrix3d Y = r.M2.cast<double>() * o.M2.cast<double>().inverse();

    auto close_to_int = [](const Eigen::Matrix3d& M)
    {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if (std::abs(M(i, j) - std::round(M(i, j))) > 1e-6)
            return false;
      return true;
    };
    return close_to_int(X) && close_to_int(Y);
  }

  //----------------------------------------------------------------------
  // unimodular {-1,0,1} matrices, for axis-wise re-basing in stage 2
  //----------------------------------------------------------------------

  const std::vector<Eigen::Matrix3i>& unimodular_W(void)
  {
    static const std::vector<Eigen::Matrix3i> W = []
    {
      std::vector<Eigen::Matrix3i> w;
      const int vals[3] = {-1, 0, 1};
      for (int i0 = 0; i0 < 3; ++i0)
      for (int i1 = 0; i1 < 3; ++i1)
      for (int i2 = 0; i2 < 3; ++i2)
      for (int i3 = 0; i3 < 3; ++i3)
      for (int i4 = 0; i4 < 3; ++i4)
      for (int i5 = 0; i5 < 3; ++i5)
      for (int i6 = 0; i6 < 3; ++i6)
      for (int i7 = 0; i7 < 3; ++i7)
      for (int i8 = 0; i8 < 3; ++i8)
      {
        Eigen::Matrix3i M;
        M << vals[i0], vals[i1], vals[i2],
             vals[i3], vals[i4], vals[i5],
             vals[i6], vals[i7], vals[i8];
        if (std::lround(M.cast<double>().determinant()) == 1)
          w.push_back(M);
      }
      return w;
    }();
    return W;
  }

  //----------------------------------------------------------------------
  // stage 2 helpers
  //----------------------------------------------------------------------

  void lattice_vectors(const Eigen::Matrix3d& A, double Lmax,
      std::vector<Eigen::Vector3i>& U_out, std::vector<double>& L_out)
  {
    Eigen::Matrix3d Ainv = A.inverse();
    int rng[3];
    for (int k = 0; k < 3; ++k)
      rng[k] = (int)std::ceil(Lmax * Ainv.col(k).norm());

    U_out.clear();
    L_out.clear();
    for (int i = -rng[0]; i <= rng[0]; ++i)
      for (int j = -rng[1]; j <= rng[1]; ++j)
        for (int k = -rng[2]; k <= rng[2]; ++k)
        {
          if ((i == 0) && (j == 0) && (k == 0))
            continue;
          // keep only one of each +-pair: the one whose first nonzero
          // component is positive
          int first_nz = (i != 0) ? i : ((j != 0) ? j : k);
          if (first_nz < 0)
            continue;

          Eigen::Vector3i u(i, j, k);
          double L = (A.transpose() * u.cast<double>()).norm();
          if ((L > 1e-9) && (L <= Lmax))
          {
            U_out.push_back(u);
            L_out.push_back(L);
          }
        }
  }

  //! Exact nearest lattice point of A to each target
  void nearest_lattice_points(const Eigen::Matrix3d& A,
      const std::vector<Eigen::Vector3d>& targets,
      std::vector<Eigen::Vector3i>& U2_out, std::vector<double>& dist_out)
  {
    Eigen::Matrix3d Ainv = A.inverse();
    std::size_t n = targets.size();
    U2_out.assign(n, Eigen::Vector3i::Zero());
    dist_out.assign(n, std::numeric_limits<double>::max());

    for (std::size_t t = 0; t < n; ++t)
    {
      Eigen::Vector3d x = Ainv.transpose() * targets[t];
      Eigen::Vector3i base(static_cast<int>(std::floor(x(0))),
                            static_cast<int>(std::floor(x(1))),
                            static_cast<int>(std::floor(x(2))));

      for (int oi = -1; oi <= 2; ++oi)
        for (int oj = -1; oj <= 2; ++oj)
          for (int ok = -1; ok <= 2; ++ok)
          {
            Eigen::Vector3i u = base + Eigen::Vector3i(oi, oj, ok);
            Eigen::Vector3d real = A.transpose() * u.cast<double>();
            double d = (real - targets[t]).norm();
            if (d < dist_out[t])
            {
              dist_out[t] = d;
              U2_out[t] = u;
            }
          }
    }
  }

  Eigen::Vector3i primitive_direction(const Eigen::Vector3i& u)
  {
    int g = std::gcd(std::gcd(std::abs(u(0)), std::abs(u(1))), std::abs(u(2)));
    if (g == 0)
      g = 1;
    return Eigen::Vector3i(u(0) / g, u(1) / g, u(2) / g);
  }

  //! General search over matched lattice vectors (allows sheared/oblique cells)
  void general_candidates(const Eigen::Matrix3d& C1, const Eigen::Matrix3d& C2,
      double tol, double Lmax, unsigned int n_pairs, double max_N,
      std::vector<Eigen::Matrix3i>& Pm_out, std::vector<Eigen::Matrix3i>& Qm_out,
      std::vector<Eigen::Matrix3d>& F_out)
  {
    std::vector<Eigen::Vector3i> U1;
    std::vector<double> L;
    lattice_vectors(C1, Lmax, U1, L);

    std::vector<Eigen::Vector3d> targets(U1.size());
    for (std::size_t i = 0; i < U1.size(); ++i)
      targets[i] = C1.transpose() * U1[i].cast<double>();

    std::vector<Eigen::Vector3i> U2;
    std::vector<double> dist;
    nearest_lattice_points(C2, targets, U2, dist);

    struct Entry { Eigen::Vector3i u1, u2; double l; };
    std::vector<Entry> kept;
    for (std::size_t i = 0; i < U1.size(); ++i)
      if (dist[i] <= tol * L[i] / (1.0 - tol))
        kept.push_back({U1[i], U2[i], L[i]});

    // deduplicate by primitive direction of u1, keeping the shortest
    std::map<std::array<int, 3>, Entry> best;
    for (auto& e : kept)
    {
      Eigen::Vector3i p = primitive_direction(e.u1);
      std::array<int, 3> key{p(0), p(1), p(2)};
      auto it = best.find(key);
      if ((it == best.end()) || (e.l < it->second.l))
        best[key] = e;
    }

    std::vector<Entry> P;
    P.reserve(best.size());
    for (auto& kv : best)
      P.push_back(kv.second);
    std::sort(P.begin(), P.end(), [](const Entry& a, const Entry& b) { return a.l < b.l; });

    if (P.size() < 3)
      return;

    std::vector<Eigen::Vector3d> W1(P.size());
    for (std::size_t i = 0; i < P.size(); ++i)
      W1[i] = C1.transpose() * P[i].u1.cast<double>();
    double V1 = std::abs(C1.determinant());

    unsigned int n = std::min<unsigned int>(n_pairs, (unsigned int)P.size());
    for (unsigned int i = 0; i < n; ++i)
    {
      for (unsigned int j = i + 1; j < n; ++j)
      {
        Eigen::Vector3d cij = W1[i].cross(W1[j]);
        for (std::size_t k = 0; k < P.size(); ++k)
        {
          double h = std::abs(W1[k].dot(cij));
          if (!((h > 1e-6 * V1) && (h <= (max_N + 0.5) * V1)))
            continue;

          Eigen::Matrix3i Pm, Qm;
          Pm.row(0) = P[i].u1; Pm.row(1) = P[j].u1; Pm.row(2) = P[k].u1;
          Qm.row(0) = P[i].u2; Qm.row(1) = P[j].u2; Qm.row(2) = P[k].u2;

          Eigen::Matrix3d S1 = Pm.cast<double>() * C1;
          Eigen::Matrix3d S2 = Qm.cast<double>() * C2;

          double d2 = S2.determinant();
          if (std::abs(d2) <= 1e-6)
            continue;
          double d1 = S1.determinant();
          if ((d2 > 0) != (d1 > 0))
            continue;

          Eigen::Matrix3d F = S1.transpose() * S2.transpose().inverse();
          double hn = spectral_norm(F - Eigen::Matrix3d::Identity());
          if (hn <= tol)
          {
            Pm_out.push_back(Pm);
            Qm_out.push_back(Qm);
            F_out.push_back(F);
          }
        }
      }
    }
  }

  //! Axis-wise rational approximation of the residual dilation between two canonical cells
  std::vector<std::tuple<Eigen::Matrix3i, Eigen::Matrix3i, Eigen::Matrix3d>>
  axis_candidates(const Eigen::Matrix3d& C1, const Eigen::Matrix3d& C2, double tol,
      unsigned int qmax)
  {
    Eigen::Matrix3d T = C1 * C2.inverse();
    const auto& W = unimodular_W();

    double best_off = std::numeric_limits<double>::max();
    Eigen::Matrix3i best_V = Eigen::Matrix3i::Identity();
    Eigen::Matrix3d best_T = T;

    for (const auto& w : W)
    {
      Eigen::Matrix3d Wd = w.cast<double>();
      Eigen::Matrix3d TV = Wd * T * Wd.inverse();

      double off = 0.0;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          if (a != b)
            off = std::max(off, std::abs(TV(a, b)));
      off += 1e-9 * (double)w.cwiseAbs().sum();

      if (off < best_off)
      {
        best_off = off;
        best_V = w;
        best_T = TV;
      }
    }

    std::vector<std::vector<std::pair<int, int>>> opts(3);
    for (int i = 0; i < 3; ++i)
    {
      double t = best_T(i, i);
      double bestv = std::numeric_limits<double>::max();
      for (unsigned int q = 1; q <= qmax; ++q)
      {
        std::set<int> ps = { std::max(1, (int)std::floor(q * t)),
                              std::max(1, (int)std::ceil(q * t)) };
        for (int p : ps)
        {
          double e = std::abs((double)q * t / (double)p - 1.0);
          if (e < bestv - 1e-12)
          {
            bestv = e;
            opts[i].push_back({(int)q, p});
          }
        }
      }
    }

    std::vector<std::tuple<Eigen::Matrix3i, Eigen::Matrix3i, Eigen::Matrix3d>> out;
    for (auto& o0 : opts[0])
      for (auto& o1 : opts[1])
        for (auto& o2 : opts[2])
        {
          Eigen::Matrix3i Dp = Eigen::Vector3i(o0.first, o1.first, o2.first).asDiagonal();
          Eigen::Matrix3i Dq = Eigen::Vector3i(o0.second, o1.second, o2.second).asDiagonal();
          Eigen::Matrix3i P = Dp * best_V;
          Eigen::Matrix3i Q = Dq * best_V;

          Eigen::Matrix3d F = (P.cast<double>() * C1).transpose() *
                              (Q.cast<double>() * C2).transpose().inverse();
          double hn = spectral_norm(F - Eigen::Matrix3d::Identity());
          if (hn <= tol)
            out.push_back(std::make_tuple(P, Q, F));
        }
    return out;
  }
}


std::vector<Candidate>
LatticeMatch::canonical_cells(const Eigen::Matrix3d& A1, const Eigen::Matrix3d& A2,
    double tol, unsigned int Nmax1, unsigned int Nmax2, unsigned int top,
    double scale_min, double scale_max)
{
  auto subs1 = sublattices(A1, Nmax1);
  auto subs2 = sublattices(A2, Nmax2);

  double V1 = std::abs(A1.determinant());
  double V2 = std::abs(A2.determinant());

  std::vector<Eigen::Matrix3d> B2inv(subs2.size());
  for (std::size_t j = 0; j < subs2.size(); ++j)
    B2inv[j] = subs2[j].B.inverse();

  std::vector<Candidate> res;
  for (std::size_t i = 0; i < subs1.size(); ++i)
  {
    const Sublattice& s1 = subs1[i];
    for (std::size_t j = 0; j < subs2.size(); ++j)
    {
      const Sublattice& s2 = subs2[j];

      double lam = std::cbrt(((double)s1.N * V1) / ((double)s2.N * V2));
      if (!((lam > scale_min) && (lam < scale_max)))
        continue;

      bool mins_ok = true;
      for (int k = 0; k < 3; ++k)
        if (std::abs(s1.mins(k) / s2.mins(k) - 1.0) > 1.05 * tol + 1e-9)
        {
          mins_ok = false;
          break;
        }
      if (!mins_ok)
        continue;

      Eigen::Matrix3d W0 = (s1.B * B2inv[j]) / lam;
      Eigen::Matrix3i W;
      double maxdiff = 0.0;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
        {
          double w = std::round(W0(a, b));
          W(a, b) = (int)w;
          maxdiff = std::max(maxdiff, std::abs(W0(a, b) - w));
        }
      if (maxdiff >= 0.45)
        continue;
      if (std::lround(W.cast<double>().determinant()) != 1)
        continue;

      Eigen::Matrix3d S2 = lam * (W.cast<double>() * s2.B);
      Eigen::Matrix3d F = s1.B.transpose() * S2.transpose().inverse();

      Candidate d = describe(F);
      if (d.H_norm > tol)
        continue;

      d.N1 = s1.N;
      d.N2 = s2.N;
      d.M1 = s1.M;
      d.M2 = W * s2.M;
      res.push_back(d);
    }
  }

  std::sort(res.begin(), res.end(), [](const Candidate& a, const Candidate& b)
      {
        if (a.N1 != b.N1) return a.N1 < b.N1;
        if (a.N2 != b.N2) return a.N2 < b.N2;
        return a.H_norm < b.H_norm;
      });

  std::vector<Candidate> out;
  for (auto& r : res)
  {
    bool is_super = false;
    for (auto& o : out)
      if (is_super_of(r, o))
      {
        is_super = true;
        break;
      }
    if (is_super)
      continue;
    out.push_back(r);
    if (out.size() >= top)
      break;
  }
  return out;
}


Candidate
LatticeMatch::real_strain(const Eigen::Matrix3d& A1, const Eigen::Matrix3d& A2,
    const Candidate& canon)
{
  Eigen::Matrix3d C1 = canon.M1.cast<double>() * A1;
  Eigen::Matrix3d C2 = canon.M2.cast<double>() * A2;
  Eigen::Matrix3d F = C1.transpose() * C2.transpose().inverse();

  Candidate d = describe(F);
  d.N1 = canon.N1;
  d.N2 = canon.N2;
  d.M1 = canon.M1;
  d.M2 = canon.M2;
  return d;
}


void
LatticeMatch::preserve_orientation(Candidate& c, const Eigen::Matrix3d& A1)
{
  Eigen::Matrix3d C1 = c.M1.cast<double>() * A1;

  for (int i = 0; i < 3; ++i)
  {
    int jbest = 0;
    double best_cos = -1.0;
    for (int j = 0; j < 3; ++j)
    {
      double cos_ij = std::abs(C1.row(i).dot(A1.row(j))) /
          (C1.row(i).norm() * A1.row(j).norm());
      if (cos_ij > best_cos)
      {
        best_cos = cos_ij;
        jbest = j;
      }
    }

    if (C1.row(i).dot(A1.row(jbest)) < 0.0)
    {
      c.M1.row(i) = -c.M1.row(i);
      c.M2.row(i) = -c.M2.row(i);
    }
  }
}


std::vector<Candidate>
LatticeMatch::real_supercells(const Eigen::Matrix3d& A1, const Eigen::Matrix3d& A2,
    const Candidate& canon, double tol, double Lmax, unsigned int n_pairs,
    unsigned int max_N, unsigned int qmax, bool shear, unsigned int top)
{
  Eigen::Matrix3d C1 = canon.M1.cast<double>() * A1;
  Eigen::Matrix3d C2 = canon.M2.cast<double>() * A2;

  struct Row { long N1; double H; Eigen::Matrix3i R1, R2; Eigen::Matrix3d F; };
  std::vector<Row> rows;

  for (auto& t : axis_candidates(C1, C2, tol, qmax))
  {
    const Eigen::Matrix3i& P = std::get<0>(t);
    const Eigen::Matrix3i& Q = std::get<1>(t);
    const Eigen::Matrix3d& F = std::get<2>(t);

    Eigen::Matrix3i R1 = P * canon.M1;
    Eigen::Matrix3i R2 = Q * canon.M2;
    long N1 = std::lround(std::abs(R1.cast<double>().determinant()));
    double H = spectral_norm(F - Eigen::Matrix3d::Identity());
    rows.push_back({N1, H, R1, R2, F});
  }

  if (shear)
  {
    double max_N_val = (max_N == 0) ? 1e9 : (double)max_N;
    std::vector<Eigen::Matrix3i> Pm, Qm;
    std::vector<Eigen::Matrix3d> Fg;
    general_candidates(A1, A2, tol, Lmax, n_pairs, max_N_val, Pm, Qm, Fg);
    for (std::size_t i = 0; i < Pm.size(); ++i)
    {
      long N1 = std::lround(std::abs(Pm[i].cast<double>().determinant()));
      double H = spectral_norm(Fg[i] - Eigen::Matrix3d::Identity());
      rows.push_back({N1, H, Pm[i], Qm[i], Fg[i]});
    }
  }

  std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b)
      {
        if (a.N1 != b.N1) return a.N1 < b.N1;
        return a.H < b.H;
      });

  std::vector<Candidate> front;
  double best = std::numeric_limits<double>::max();
  for (auto& r : rows)
  {
    if (r.H < best - 1e-9)
    {
      best = r.H;
      Candidate d = describe(r.F);
      d.N1 = (unsigned int)r.N1;
      d.N2 = (unsigned int)std::lround(std::abs(r.R2.cast<double>().determinant()));
      d.M1 = r.R1;
      d.M2 = r.R2;
      front.push_back(d);
    }
  }
  if (front.size() > top)
    front.resize(top);
  return front;
}
