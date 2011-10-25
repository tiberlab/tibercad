// $Id$

#include "SBuserHamiltonian.h"
#include "Database.h"
#include "Constants.h"

using namespace std;

//======================================================================//
SBuserHamiltonian::~SBuserHamiltonian()
{
  
}

//=======================================================================//
SBuserHamiltonian::SBuserHamiltonian(const ModelOptions& options)
  : SBbulkHamiltonian(options)
{
  edge = 0.0;
 
  
  imass = Tensor2Sym(1);
}



//======================================================================//


void SBuserHamiltonian::do_init( )
{


  EFAbulkHamiltonian::do_init();

  // for now this is automatically a valence band
  if (get_option("particle", "") == "hl")
  {
    // band edge
    const Database& db = get_database();
    db.set_section("valenceband");
    edge = db.get("E_v", 0.0);

    // mass
    const string& band = get_option("band", "hh");
  }

  // override of band edge
  edge = get_option("band_edge", edge);
  edge /= Constants::Hartree;

  // one degenerate band
  kp_bands.resize(1,0);

  kp_bands_map.insert(std::make_pair(2,0));
  // NOTE: this duplicates the band, which makes OpticsKP use two
  // identical bands with opposite spin
  // BUT: maybe optics cannot do anything with this single band anyways?
  kp_bands_map.insert(std::make_pair(3,0));

  double mass = get_option("mass", 1.0);

  if (mass == 0.0) throw InitFailedException("User-defined Hamiltonian: zero mass");

  imass = -(1.0 / mass) * Tensor2Sym(1.0);
  

  calculate_Hamiltonian_gen();
  
}

//======================================================================//
