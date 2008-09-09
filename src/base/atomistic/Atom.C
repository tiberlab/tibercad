#include "Atom.h"


Atom::Atom(){
	_specie = 'none';
        _region_id = 0;
	_flag = 0;
}

Atom::Atom(std::string& init_specie, Tensor1& init_position){
	_position = init_position;
	_specie = init_specie;
	_region_id = 0;
	_flag = 0;
}

Atom::~Atom(){
}




