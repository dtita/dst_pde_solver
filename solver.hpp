#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

namespace dauphine
{
	class mesh {
	public:
		mesh(double time_space, double space_space,double maturity, double spot_max);
		double get_mesh_time_space() const;
		double get_mesh_maturity() const;
		double get_mesh_space_space() const;
		double get_mesh_spot_max() const;
		~mesh();
	private:
		double m_time_space;
		double m_space_space;
		double m_maturity;
		double m_spot_max;
	};
}

#endif
