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


	//pour initialiser payoff, rate, vol, boundaries
	class initial_function {
	public:
		initial_function(double(*f)(std::vector<double>));
		double function_operator(std::vector<double> arguments);
		~initial_function();
	private:
		double(*m_f)(std::vector<double>);
	};


	class price_vector {
	public: 
		price_vector();
	private:
	};

}

#endif
