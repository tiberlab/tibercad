#include "BondMap.h"




BondMap::BondMap(void)
{
	_bond_map = NULL;
};


BondMap::~BondMap(void)
{
_grid_cell.clear();

//	for (unsigned int i = 0; i < _size; i++)  {
//		delete _bond_map[i];
//	}
//delete _bond_map;

};


void
BondMap::set_cutoff()
{
	_cutoff["Si"] = 1.5;
	_cutoff["Ga"] = 1.2;
	_cutoff["N"] = 1.2;
	_cutoff["Al"] = 1.2;
}

void
BondMap::do_init(const std::vector<Atom> &basis, const Tensor2Gen& period)
{
	int x, y;

	std::cout << "BondMap::do_init " << std::endl;
std::cout << "basis size is " << basis.size() << std::endl;

std::cout << "In BondMap pointer is " << this << std::endl;
	std::cout << "setting cutoff" << std::endl;

	//!set cutoff distancies
	set_cutoff();

	std::cout << "Periodical vectors are " << period << std::endl;

	_period = period;

	//Define the minimum spacing of the grid. The smaller it is, the faster is bonds calculations
	//Cannot be smaller than the higher bond lenght. (In Amstrong)
	const double minimum_spacing = 10.0;

	// Define the addressing grid
	//------------------------------------------------
	Tensor1 edge_min, edge_max;
	define_edges(basis, edge_min, edge_max);

	_origin = edge_min;

	// Maximum edges must be slightly enlarged to avoid strict inclusion errors
	edge_max(1) = edge_max(1) + 0.01;
	edge_max(2) = edge_max(2) + 0.01;
	edge_max(3) = edge_max(3) + 0.01;


	define_grid(minimum_spacing, edge_min, edge_max);

	//------------------------------------------------

	//Resize grid_cell according to grid dimensions
	_grid_cell.resize(_n_x);

	std::cout << "Allocating _grid_cell " << std::endl;

	for (x = 0; x < _n_x; x++)
	{
		_grid_cell[x].resize(_n_y);

		for (y = 0; y < _n_y; y++)
		{
			_grid_cell[x][y].resize(_n_z);

		}
	}

std::cout << "grid is sized " << _n_x << " " << _n_y << " " << _n_z << std::endl;

	//! Allocate bond map
	_bond_map = new unsigned int* [basis.size()];
	for (unsigned int i = 0; i < basis.size(); i++)  {
		_bond_map[i] = new unsigned int [9];
		for (unsigned int j = 0; j < 9; j++) _bond_map[i][j] = 0;
	}

}


void
BondMap::do_solve(const std::vector<Atom>& basis, const Tensor2Gen& period)
{

	// Include atoms in grid cells
	include_atoms(basis);

	//Process cells themselves
	for (unsigned int x = 0; x < _n_x; x++)
	    {
	    	for (unsigned int y = 0; y < _n_y; y++)
	    	{
	    		for (unsigned int z = 0; z < _n_z; z++)
	    		{
	    			process_surrounding_cell(basis, x, y, z);
	    		}
	    	}
	    }

	std::cout << "Trying to process corners " << std::endl;
	//Process cells and neighbours

	//process_internal(basis);

//	for (unsigned int i = 0; i < basis.size(); i++)
//	{
//	std::cout << "bond_map["<<i<<"][8] is " << _bond_map[i][8] << std::endl;
//	}

}



void
BondMap::include_atoms(const std::vector<Atom>& basis)
{
	unsigned int x, y, z;
	for (unsigned int i = 0; i < basis.size(); i++)
	{
		x = ( (floor(basis[i].get_position()(1) - _origin(1)) / _x_spacing ));
		y = ( (floor(basis[i].get_position()(2) - _origin(2)) / _y_spacing ));
		z = ( (floor(basis[i].get_position()(3) - _origin(3)) / _z_spacing ));

		_grid_cell[x][y][z].push_back(i);
	}

}


//void
//BondMap::process_atoms(const std::vector<Atom>& basis, const unsigned int i,
//		const unsigned int j)
//{
//	Tensor1 period;
//	period(1) = 0.0; period(2) = 0.0; period(3) = 0.0;
//	process_atoms(basis, i, j, period);
//}



void
BondMap::process_atoms(const std::vector<Atom>& basis, const unsigned int i,
		const unsigned int j, const Tensor1& period)
{

	unsigned int put_here = 0;
	bool not_already_signed;
	double cutofftmp;
	Tensor1 position1, position2;

	if ((i != j)/*&&(i==8)&&(j==97)*/)
	{


		cutofftmp = _cutoff[basis[i].get_specie().c_str()] + _cutoff[basis[j].get_specie().c_str()];
		//std::cout << "CUTOFFTMP IS " << cutofftmp<< std::endl;
		position1 = basis[i].get_position();
		position2 = basis[j].get_position() + period;

		//std::cout << "Distance is " << norm( position1 - position2) << std::endl;

		if ( cutofftmp == 0.0 ) std::cout << "WARNING, A CUTOFF DISTANCE IS NOT DEFINED " <<
		basis[i].get_specie() << basis[j].get_specie() << " wow" << std::endl;





		if ( norm( position1 - position2) < cutofftmp ){

			not_already_signed = true;
			put_here = _bond_map[i][8];
			//std::cout << " put_here is " << put_here << std::endl;


			//						for (unsigned int j1 = 0; j1 < 9; j1++)
			//							std::cout << "_bond_map[:][:] " <<  _bond_map[i][j1];

			if (put_here != 0)
			{
				for (unsigned int n = 0; n < put_here; n++){
					//std::cout << "first atom, n is " << n << std::endl;
					if (_bond_map[i][n] == j)
					{//std::cout << " already signed! ";
						not_already_signed = false;
						break;
					}
				}
			}
			else
				not_already_signed = true;
			if (not_already_signed)
			{
				//std::cout << "PROCESSING i : PUTTING " << j<< " IN ATOM " << i<< std::endl;
				_bond_map[i][put_here] = j;
				_bond_map[i][8]++;
			}

			not_already_signed = true;
			put_here = _bond_map[j][8];


			if (put_here != 0)
			{
				for (unsigned int n = 0; n < put_here; n++){
					//std::cout << "second atom, n is " << n << std::endl;

					if (_bond_map[j][n] == i)
					{
						not_already_signed = false;
						break;
					}
				}
			}
			else
				not_already_signed = true;
			if (not_already_signed)
			{
				//std::cout << "PROCESSING j: PUTTING " << i<< " IN ATOM " << j<< std::endl;
				_bond_map[j][put_here] = i;
				_bond_map[j][8]++;
			}

		}



	}



}




//void
//BondMap::process_atom_with_cell(const std::vector<Atom>& basis, const unsigned int i,
//		const unsigned int x, const unsigned int y, const unsigned int z)
//{
//	Tensor1 period;
//	period(1) = 0.0; period(2) = 0.0; period(3) = 0.0;
//	process_atom_with_cell(basis, i, x, y, z, period);
//}


void
BondMap::process_atom_with_cell(const std::vector<Atom>& basis, const unsigned int i,
		const unsigned int x, const unsigned int y, const unsigned int z, const Tensor1& period)
{
	//std::cout << "process_atom_with_cell... " ;

	for (unsigned int j = 0; j < _grid_cell[x][y][z].size(); j++)
	{
		process_atoms(basis, i, _grid_cell[x][y][z][j], period);
	}

	//std::cout << "done" << std::endl;
}



void
BondMap::process_cell(const std::vector<Atom>& basis, const unsigned int x1,
		const unsigned int y1, const unsigned int z1,
		int x2, int y2, int z2)
{
Tensor1 shift(0);
Tensor1 period_x(0), period_y(0), period_z(0);

period_x(1) = _period(1,1); period_x(2) = _period(2,1); period_x(3) = _period(3,1);
period_y(1) = _period(1,2); period_y(2) = _period(2,2); period_y(3) = _period(3,2);
period_z(1) = _period(1,3); period_z(2) = _period(2,3); period_z(3) = _period(3,3);


	//This section check informations about periodicity and computes a right shift
	//vector "shift"
	//------------------------------------------
if (x2 == -1)
{
	x2 = _n_x - 1; shift = shift - period_x;
}
if (x2 == _n_x)
{
	x2 = 0; shift = shift + period_x;
}
if (y2 == -1)
{
	y2 = _n_y - 1; shift = shift - period_y;
}
if (y2 == _n_y)
{
	y2 = 0; shift = shift + period_y;
}
if (z2 == -1)
{
	z2 = _n_z - 1; shift = shift - period_z;
}
if (z2 == _n_z)
{
	z2 = 0; shift = shift + period_z;
}
	//------------------------------------------

	for (unsigned int i = 0; i < _grid_cell[x1][y1][z1].size(); i++)
	{
		process_atom_with_cell(basis, _grid_cell[x1][y1][z1][i], x2, y2, z2, shift);
	}

}



//void
//BondMap::process_single_cell(const std::vector<Atom>& basis, const unsigned int x,
//		const unsigned int y, const unsigned int z)
//{
//	for (unsigned int i = 0; i < _grid_cell[x][y][z].size(); i++)
//	{
//		process_atom_with_cell(basis, _grid_cell[x][y][z][i], x, y, z);
//	}
//}
//

//void
//BondMap::process_internal(const std::vector<Atom>& basis)
//{
//unsigned int x, y, z;
//
//for (x = 0; x < _n_x; x++)
//{
//	for (y = 0; y < _n_y; y++)
//	{
//		for (z = 0; z < _n_z; z++)
//		{
//			process_surrounding_cell(basis, x, y, z);
//		}
//	}
//}
//}


void
BondMap::process_surrounding_cell(const std::vector<Atom>& basis, const unsigned int x, const unsigned int y,
		const unsigned int z)
{
	//N.B. as this function is used to process internal cell, no periodicity information is used

	int x1, y1, z1;
    int tmp_x, tmp_y, tmp_z;

    tmp_x = x; tmp_y = y; tmp_z = z;

	for (x1 = tmp_x -1; x1 <= tmp_x + 1; x1++)
    {

    	for (y1 = tmp_y -1; y1 <= tmp_y + 1; y1++)
    	{

    		for (z1 = tmp_z -1; z1 <= tmp_z + 1; z1++)
    		{
    			process_cell(basis, x, y, z, x1, y1, z1);
    		}
    	}
    }

}



void
BondMap::define_edges(const std::vector<Atom>& basis, Tensor1& edge_min, Tensor1& edge_max)
{
	edge_min(1) = 0.0; edge_min(2) = 0.0; edge_min(3) = 0.0;
	edge_max(1) = 0.0; edge_max(2) = 0.0; edge_max(3) = 0.0;

	for (unsigned int i = 0; i < basis.size(); i++)
	{
		if (basis[i].get_position()(1) < edge_min(1)) edge_min(1) = basis[i].get_position()(1);
		if (basis[i].get_position()(1) > edge_max(1)) edge_max(1) = basis[i].get_position()(1);

		if (basis[i].get_position()(2) < edge_min(2)) edge_min(2) = basis[i].get_position()(2);
		if (basis[i].get_position()(2) > edge_max(2)) edge_max(2) = basis[i].get_position()(2);

		if (basis[i].get_position()(3) < edge_min(3)) edge_min(3) = basis[i].get_position()(3);
		if (basis[i].get_position()(3) > edge_max(3)) edge_max(3) = basis[i].get_position()(3);
	}


}


void
BondMap::define_grid(const double minimum_spacing, const Tensor1& edge_min, const Tensor1& edge_max)
{
	double lenght_x, lenght_y, lenght_z;

	lenght_x = edge_max(1) - edge_min(1);
	lenght_y = edge_max(2) - edge_min(2);
	lenght_z = edge_max(3) - edge_min(3);

	std::cout << "lenghts are " << lenght_x << lenght_y << lenght_z << std::endl;

	_n_x = (floor(lenght_x / minimum_spacing)) + 1; _x_spacing = lenght_x / _n_x;
	_n_y = (floor(lenght_y / minimum_spacing)) + 1; _y_spacing = lenght_y / _n_y;
	_n_z = (floor(lenght_z / minimum_spacing)) + 1; _z_spacing = lenght_z / _n_z;

//	if (_n_x == 0) _n_x = 1;
//	if (_n_y == 0) _n_y = 1;
//	if (_n_y == 0) _n_y = 1;


}
