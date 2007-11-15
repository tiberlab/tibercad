#include "Atom.h"


Atom::Atom(){}
	
Atom::~Atom(){}
	
Atom::Atom(std::string& init_specie, Tensor1& init_position){
	position = init_position;
	specie = init_specie;
}


