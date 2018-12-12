#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

namespace dauphine
{
	class mesh {
	public:
		//mesh();
		mesh(double dt, double dx,double maturity, std::vector<double> spot_boundaries);
		double get_mesh_dt() const;
		double get_mesh_maturity() const;
		double get_mesh_dx() const;
		std::vector<double> get_mesh_spot_boundaries() const;
		~mesh();
	private:
		double m_dt;
		double m_dx;
		double m_maturity;
		std::vector<double> m_spot_boundaries;
	};

	class payoff {
	public:
		payoff(double(*f)(std::vector<double>));
		double function_operator(std::vector<double> arguments);
	private:
		double (*m_f)(std::vector<double>);
	};


	class boundary_conditions {
	};
	class rate {
	};
	class volatility {
	};
	class theta {
	};
	class price_vector {
	};
}

#endif
