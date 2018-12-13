#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

namespace dauphine
{
	class mesh {
	public:
		//mesh();
		mesh(double dt, double dx,double maturity, double spot, double theta, std::vector<double> spot_boundaries);
		double get_mesh_dt() const;
		double get_mesh_maturity() const;
		double get_mesh_dx() const;
		double get_mesh_spot() const;
		double get_mesh_theta() const;
		std::vector<double> get_mesh_spot_boundaries() const;
		~mesh();
	private:
		double m_dt;
		double m_dx;
		double m_maturity;
		double m_spot;
		double m_theta;
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
	double diag_coeff(mesh m, initial_function rate, std::vector<double> arguments);
	double subdiag_coeff(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments);
	double updiag_coeff(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments);

	//class matrix_elements {
	//public:
	//	matrix_elements(initial_function rate);
	//};
	//class diag_coeff : public matrix_elements{
	//public:
	//	double diag_function();
	//};


	class price_vector {
	public: 
		price_vector(initial_function payoff, mesh m, std::vector<double> arguments);
	private:
		initial_function m_payoff;
		mesh m_m;
	};

}

#endif
