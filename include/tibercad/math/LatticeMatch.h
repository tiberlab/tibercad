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
 * \file LatticeMatch.h
 * \brief Public tiberCAD API header.
 *
 * \note This header is part of the public tiberCAD API.
 *       API version: 3.5
 */

#ifndef TC_LATTICEMATCH_H
#define TC_LATTICEMATCH_H

#include <Eigen/Dense>

#include <vector>


//! Lattice matching of two crystals at a fixed, given relative orientation
/*!
 * This namespace finds commensurate (coincident) supercells of two
 * lattices, given as 3x3 matrices whose ROWS are lattice vectors (in a
 * shared Cartesian frame, i.e. already in the desired relative
 * orientation). It does not use any plane/direction information.
 *
 * The approach is in two stages:
 *
 * -# canonical_cells(): find the smallest sublattices (by enumerating
 *    Hermite Normal Form sublattices of each lattice and Buerger-reducing
 *    them) that coincide in SHAPE, after the second lattice is given a
 *    free, compatible isotropic scale factor. Anisotropic stretch, shear
 *    and rotation beyond that are residual, and bounded by \c tol.
 * -# real_supercells(): starting from a canonical correspondence, search
 *    for actual supercells using the REAL lattice constants (no free
 *    rescaling), producing a Pareto front trading off supercell size
 *    against residual strain. Two search strategies are merged: cheap
 *    axis-wise rational approximation of the residual dilation along the
 *    canonical cell's own axes, and a general search over matched lattice
 *    vectors (allowing sheared/oblique supercells) up to a length cutoff.
 *
 * Both stages report their result as a deformation gradient \c F mapping
 * the second lattice (the "superstrate") onto the first (the
 * "substrate"), together with its polar decomposition and strain measures.
 */
namespace LatticeMatch
{
  //! A candidate coincident/commensurate supercell pair
  struct Candidate
  {
    //! Supercell size of the first lattice, in primitive cells
    unsigned int N1 = 0;

    //! Supercell size of the second lattice, in primitive cells
    unsigned int N2 = 0;

    //! Integer coefficients of the first lattice's supercell (rows), w.r.t. the input basis A1
    Eigen::Matrix3i M1 = Eigen::Matrix3i::Identity();

    //! Integer coefficients of the second lattice's supercell (rows), w.r.t. the input basis A2
    Eigen::Matrix3i M2 = Eigen::Matrix3i::Identity();

    //! Deformation gradient mapping the second lattice onto the first
    Eigen::Matrix3d F = Eigen::Matrix3d::Identity();

    //! Rotation part of the polar decomposition F = R*U
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();

    //! Stretch part of the polar decomposition F = R*U
    Eigen::Matrix3d U = Eigen::Matrix3d::Identity();

    //! Biot strain, U - I
    Eigen::Matrix3d biot = Eigen::Matrix3d::Zero();

    //! Principal strains, i.e. the eigenvalues of \c biot
    Eigen::Vector3d principal = Eigen::Vector3d::Zero();

    //! Spectral norm (induced 2-norm) of F - I
    double H_norm = 0.0;

    //! Rotation angle of \c R, in degrees
    double rot_deg = 0.0;

    //! Hydrostatic part, det(F)^(1/3) - 1
    double hydro = 0.0;
  };

  /*!
   * \brief Find the smallest coincident sublattices (shape match only)
   *
   * \param A1 the first (substrate) lattice vectors, as rows
   * \param A2 the second (superstrate) lattice vectors, as rows, in the
   *        same Cartesian frame (i.e. already in the desired relative
   *        orientation)
   * \param tol tolerance on the spectral norm of the residual deformation
   *        gradient (stretch + shear + rotation), after a free isotropic
   *        rescaling of \p A2
   * \param Nmax1 maximum supercell size to consider for \p A1, in
   *        primitive cells
   * \param Nmax2 maximum supercell size to consider for \p A2, in
   *        primitive cells
   * \param top maximum number of candidates to return
   * \param scale_min minimum acceptable compatible isotropic scale factor
   * \param scale_max maximum acceptable compatible isotropic scale factor
   *
   * \return candidates sorted by (N1, N2, H_norm), with candidates that
   *         are themselves supercells of an already-returned (smaller)
   *         candidate removed
   */
  std::vector<Candidate> canonical_cells(const Eigen::Matrix3d& A1, const Eigen::Matrix3d& A2,
      double tol = 0.03, unsigned int Nmax1 = 8, unsigned int Nmax2 = 12,
      unsigned int top = 5, double scale_min = 0.3, double scale_max = 3.0);

  /*!
   * \brief Find real commensurate supercells built on a canonical correspondence
   *
   * \param A1 the first (substrate) lattice vectors, as rows
   * \param A2 the second (superstrate) lattice vectors, as rows
   * \param canon a canonical correspondence, as returned by canonical_cells()
   * \param tol tolerance on the spectral norm of the residual deformation
   *        gradient, now with the REAL (not freely rescaled) lattice
   *        constants
   * \param Lmax cutoff length (in the same units as \p A1, \p A2) for the
   *        general oblique/sheared supercell search
   * \param n_pairs number of shortest matched-vector candidates to pair up
   *        in the general search
   * \param max_N cap on the first lattice's supercell size considered in
   *        the general search, in primitive cells (0 = unbounded)
   * \param qmax maximum denominator for the axis-wise rational
   *        approximation search
   * \param shear if true, also run the general oblique/sheared search (in
   *        addition to the cheap axis-wise rational approximation); if
   *        false, only axis-wise (shear-free) supercells are considered
   * \param top maximum number of Pareto-front points to return
   *
   * \return a Pareto front (size N1 vs. H_norm), sorted by increasing N1
   *         with strictly decreasing H_norm
   */
  std::vector<Candidate> real_supercells(const Eigen::Matrix3d& A1, const Eigen::Matrix3d& A2,
      const Candidate& canon, double tol = 0.01, double Lmax = 60.0,
      unsigned int n_pairs = 40, unsigned int max_N = 0, unsigned int qmax = 40,
      bool shear = true, unsigned int top = 10);

  /*!
   * \brief Deformation gradient of a canonical correspondence, using the
   *        REAL lattice vectors
   *
   * canonical_cells() matches lattice shape only, after giving \p A2 a
   * free isotropic rescaling -- its \c Candidate::F is therefore not the
   * actual strain needed to match the real lattices. This instead builds
   * the supercells \p canon.M1 * \p A1 and \p canon.M2 * \p A2 from the
   * real lattice vectors (no rescaling, no further search) and reports
   * the deformation gradient between them, with the same strain measures
   * as the other candidates.
   *
   * \param A1 the first (substrate) lattice vectors, as rows
   * \param A2 the second (superstrate) lattice vectors, as rows
   * \param canon a canonical correspondence, as returned by canonical_cells()
   *
   * \return the candidate \p canon, with \c F (and the strain measures
   *         derived from it) replaced by the real-lattice values
   */
  Candidate real_strain(const Eigen::Matrix3d& A1, const Eigen::Matrix3d& A2,
      const Candidate& canon);

  /*!
   * \brief Re-sign a candidate's supercell rows to match the substrate
   *        lattice's own axis directions
   *
   * A supercell row of \c M1 and the corresponding row of \c M2 can be
   * jointly negated without changing \c F at all (both lattice vectors
   * feeding that row of the linear map get negated together, and a linear
   * map applied to \c -v is just \c -F(v)) -- so this freedom carries no
   * physical content. Buerger reduction does not care about this sign
   * (every direction is tied in length with its mirror), so left alone it
   * can report a supercell basis that looks like an arbitrary axis flip of
   * the input lattice -- most noticeably when it falls on a direction the
   * caller cares about keeping, such as a growth direction. This picks,
   * for each row, whichever of the two signs has non-negative projection
   * onto the \p A1 row it is most closely aligned with, and negates both
   * \c M1 and \c M2's row if needed to match. \c F and the strain measures
   * derived from it are unaffected.
   *
   * \param c the candidate to re-sign, in place
   * \param A1 the first (substrate) lattice vectors, as rows -- the
   *        reference orientation to preserve
   */
  void preserve_orientation(Candidate& c, const Eigen::Matrix3d& A1);
}

#endif // TC_LATTICEMATCH_H
