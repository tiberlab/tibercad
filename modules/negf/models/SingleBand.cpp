/*  
 * This file is part of the tiberCAD module negf.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file SingleBand.C
 * \brief tiberCAD negf module implementation.
 *
 * \note This file is part of module negf.
 */


#include "SingleBand.h"
#include "tibercad/io/Database.h"
#include "tibercad/physics/Material.h"

#include "tibercad/module/TiberModule.h"



SingleBand::SingleBand(const ModelOptions& options)
   : HamiltonianModel(options)
{
   //_inv_mass_crys(0);
   //_inv_mass(0);
}


void
SingleBand::do_init(void)
{

  _model = "single_band"; 

  _num_bands = 1;

  std::string particle = get_option("particle","el");

  if (particle=="el")
  {
    _band_type = "Ec";
  }

  if (particle=="hl")
  {
    _band_type = "Ev";
  }

  read_database();

  set_invmass_tensor();

  //std::cerr<<"invmass"<<std::endl<<_inv_mass<<std::endl;
}

void SingleBand::read_database(void)
{
  const Database& db = get_database();

  double m_l = 0.0;
  double m_t = 0.0;

  if (_band_type == "Ec")
  {
    db.set_section("conductionband");

    std::string valley = get_option("valley","G");

    if (valley=="G")
    {
      m_l = db.get("m_G", 1.0);
      m_t = db.get("m_G", 1.0);
    }

    if (valley=="X1" || valley=="X2" || valley=="X3")
    {
      m_l = db.get("m_X_l", 1.0);
      m_t = db.get("m_X_t", 1.0);
    }

    if (valley=="L1" || valley=="L2" || valley=="L3")
    {
      m_l = db.get("m_L_l", 1.0);
      m_t = db.get("m_L_t", 1.0);
    }

  }

  if (_band_type == "Ev")
  {
    db.set_section("valenceband");

    std::string band = get_option("band","HH");

    if (band=="HH")
    {
      m_l = -1 * db.get("m_HH", 1.0);
      m_t = -1 * db.get("m_HH", 1.0);
    }

    if (band=="LH")
    {
      m_l = -1 * db.get("m_LH", 1.0);
      m_t = -1 * db.get("m_LH", 1.0);
    }
  }

  _inv_mass_crys(0,0) = 1.0/m_l;
  _inv_mass_crys(1,1) = 1.0/m_t;
  _inv_mass_crys(2,2) = 1.0/m_t;
}

void SingleBand::set_invmass_tensor(void)
{
  std::string band = get_option("band","G");
  std::string valley = get_option("valley",band);

  // vector u defines rotation with respect to e1 = (1,0,0)
  // e.g. u = (0 1 0)

  VectorValue<double> u(0.0);

  if (valley=="G"){ u(0) = 1.0; _degeneracy = 1.0; }

  if (valley=="X"){ u(0) = 1.0; _degeneracy = 6.0; }

  if (valley=="X1"){ u(0) = 1.0; _degeneracy = 2.0; }
  if (valley=="X2"){ u(1) = 1.0; _degeneracy = 2.0; }
  if (valley=="X3"){ u(2) = 1.0; _degeneracy = 2.0; }

  if (valley=="L")
  {
    u(0) = 1.0;  u(1) = 1.0;  u(2) = 1.0;
    _degeneracy = 6.0;
  }

  if (valley=="L1")
  {
     u(0) = 1.0;  u(1) = 1.0;  u(2) = 1.0;
     _degeneracy = 2.0;
  }
  if (valley=="L2")
  {
    u(0) = 1.0;  u(1) = 0.0;  u(2) = -1.0;
    _degeneracy = 2.0;
  }
  if (valley=="L3")
  {
    u(0) = 1.0;  u(1) = -2.0;  u(2) = 1.0;
    _degeneracy = 2.0;
  }

  if (valley=="HH"){ u(0) = 1.0; _degeneracy = 1.0; }

  if (valley=="LH"){ u(0) = 1.0; _degeneracy = 1.0; }


  _degeneracy = get_option("degeneracy", _degeneracy);

  // set rotation R of the valley (G, X, L)

  VectorValue<double> e1(3, 0.0);
  e1(0)=1.0;
  double norm = u*u;

  VectorValue<double> k = e1.cross(u)/norm;

  double cos = u*e1/norm;

  double sin = sqrt(1-cos*cos);

  TensorValue<double> K(0.0);

  K(0,1) = -k(2); K(0,2) = k(1);
  K(1,0) = k(2); K(1,2) = -k(0);
  K(2,0) = -k(1); K(2,1) = k(0);

  TensorValue<double> R;

  R = sin*K + (1-cos)*K*K;

  R(0,0) = R(0,0) + 1; R(1,1) = R(1,1) + 1; R(2,2) = R(2,2) + 1;

  // Set rotation from Crystal to Device
  const RealTensor& Rot = get_material()->get_rotation_matrix();


  //std::cerr<<"Rot"<<std::endl<<Rot<<std::endl;
  _inv_mass = Rot * R * _inv_mass_crys * R.transpose() * Rot.transpose();

}
