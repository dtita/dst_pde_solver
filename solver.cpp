#include "solver.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

namespace dauphine
{
	mesh::mesh(double time_space, double space_space, double maturity, double spot_max)
		: m_time_space(time_space), m_space_space(space_space), m_maturity(maturity),m_spot_max(spot_max)
	{
	}
	mesh::~mesh()
	{
	}
	double mesh::get_mesh_maturity() const{
		return m_maturity;
	}
}

