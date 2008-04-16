#include "Atom.h"


Atom::Atom(){
        id = 0;
	flag = 0;
}
	
Atom::Atom(std::string& init_specie, Tensor1& init_position){
	position = init_position;
	specie = init_specie;
	id = 0;
	flag = 0;
}

Atom::~Atom(){
}
	



