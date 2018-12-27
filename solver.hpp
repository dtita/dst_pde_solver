#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

namespace dauphine
{
	class mesh {
	public:
		//mesh();
		mesh(double dt, double dx,double maturity, double spot, std::vector<double> spot_boundaries);
		double get_mesh_dt() const;
		double get_mesh_maturity() const;
		double get_mesh_dx() const;
		double get_mesh_spot() const;
		std::vector<double> spot_vector();
        
		std::vector<double> get_mesh_spot_boundaries() const;
		~mesh();
	private:
		double m_dt;
		double m_dx;
		double m_maturity;
		double m_spot;
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
	double diag_coeff(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments);
	double subdiag_coeff(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments);
	double updiag_coeff(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments);
	std::vector<double> initial_price_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments,initial_function payoff);
	std::vector<double> column_up(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff,std::vector<double> up_price);
	std::vector<double> price_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff, std::vector<double> col_up);
	std::vector<double> price_today(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff);

  
    // Autre methode - TEST
    
    std::vector<double> up_vector(mesh m, initial_function rate,initial_function vol, std::vector<double> arguments);
    std::vector<double> sub_vector(mesh m, initial_function rate,initial_function vol, std::vector<double> arguments);
    std::vector<double> diag_vector(mesh m, initial_function rate,initial_function vol, std::vector<double> arguments);
    //std::vector<double>spot_vector_bis(mesh m, initial_function rate,initial_function vol, std::vector<double> arguments);
   // void tridiag_algorithm(const std::vector<double>& a,const std::vector<double>& b,const std::vector<double>& c,const std::vector<double>& d,std::vector<double>& f);
    std::vector<double> tridiagonal_solver(std::vector<double> & a, std::vector<double> & b, std::vector<double> & c, std::vector<double> & f);
}

#endif
