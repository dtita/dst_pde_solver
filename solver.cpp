#include "solver.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>
namespace dauphine
{

	mesh::mesh(double dt, double dx, double maturity, double spot, std::vector<double> boundaries)
		: m_dt(dt), m_dx(dx), m_maturity(maturity),m_spot(spot), m_spot_boundaries(boundaries)
	{
	}

	mesh::~mesh()
	{
	}
	double mesh::get_mesh_maturity() const {
		return m_maturity;
	}
	double mesh::get_mesh_dt() const {
		return m_dt/365;
	}
	double mesh::get_mesh_dx() const {
		return m_dx;
	}
	double mesh::get_mesh_spot() const {
		return m_spot;
	}

	std::vector<double> mesh::get_mesh_spot_boundaries() const {
		return m_spot_boundaries;
	}


	std::vector<double> mesh::spot_vector() {
		int size = floor((m_spot_boundaries[0] - m_spot_boundaries[1]) / m_dx)+1;
		std::vector<double> result(size);
		for (std::size_t i = 0; i < result.size(); ++i)
		{
			result[i] =log(m_spot_boundaries[1]+i*m_dx);
		}
		return result;
	}


	initial_function::initial_function(double(*f)(std::vector<double>))
		: m_f(f)
	{
	}
	double initial_function::function_operator(std::vector<double> arguments) {
		return m_f(arguments);
	}
	initial_function::~initial_function() {
	}

	double diag_coeff(mesh m, initial_function rate,initial_function vol, std::vector<double> arguments) {
		if (arguments[0] == arguments[2]|| arguments[0] == arguments[3])
		{
			return 1;
		}
		else
		{
			return 1 / m.get_mesh_dt() + rate.function_operator(arguments)*arguments[4]+1/(m.get_mesh_dx()*m.get_mesh_dx())*arguments[4]*pow(vol.function_operator(arguments),2);
		}
	}
	double subdiag_coeff(mesh m, initial_function rate,initial_function vol, std::vector<double> arguments) {
			return -1 /2* arguments[4] / (m.get_mesh_dx()*m.get_mesh_dx())*pow(vol.function_operator(arguments), 2) +1/(4*m.get_mesh_dx())*arguments[4] *(pow(vol.function_operator(arguments),2)-rate.function_operator(arguments));
	}
	double updiag_coeff(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments) {

			return -1 / 2 * arguments[4] / (m.get_mesh_dx()*m.get_mesh_dx())*pow(vol.function_operator(arguments), 2) + 1 / (4 * m.get_mesh_dx())*arguments[4] *(-pow(vol.function_operator(arguments), 2) +rate.function_operator(arguments));
	}

	std::vector<double> initial_price_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff) {
		std::vector<double> result=m.spot_vector();
		for (std::size_t i = 0; i < result.size(); i++) {
			arguments[0] = (exp(result[i]));
			if (payoff.function_operator(arguments) == 0) {
				result[i] = 0;
			}
			else {
				result[i] = (payoff.function_operator(arguments));
			}
		}
		return result;
	}
	std::vector<double> column_up(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff,std::vector<double> up_price) {
		std::vector<double> result = up_price;
		std::vector<double> result2 = m.spot_vector();
		arguments[4] = arguments[4] - 1;
		for (std::size_t i = 1; i < result.size()-1; i++) {
			arguments[0] = result2[i];
			result[i] = result[i] * diag_coeff(m, rate,vol ,arguments) + result[i - 1] * subdiag_coeff(m, rate, vol, arguments) + result[i + 1] * updiag_coeff(m, rate, vol, arguments);
		}
		return result;
	}

	//algo used : https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	std::vector<double> price_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff,std::vector<double> col_up) {
		std::vector<double> result(col_up.size());
		double W = 0;
		std::vector<double> arguments_up = arguments;
		std::vector<double> arguments_down = arguments;
		std::vector<double> B(col_up.size());
		std::vector<double> D(col_up.size());
		B[0] = 0;
		D[0] = col_up[0];
		for (std::size_t i = 1; i < col_up.size()-1; i++) {
			arguments_down[0] = arguments[3] + (i-1)*m.get_mesh_dx();
			arguments[0] = arguments[3] + i*m.get_mesh_dx();
			arguments_up[0] = arguments[3] + (i+1)*m.get_mesh_dx();
			B[i] = subdiag_coeff(m, rate, vol, arguments) / (diag_coeff(m, rate, vol, arguments_down) - B[i - 1] * updiag_coeff(m, rate, vol, arguments_down));
			D[i] = (col_up[i] - D[i - 1] * updiag_coeff(m, rate, vol, arguments_down)) / (diag_coeff(m, rate, vol, arguments_down) - B[i - 1] * updiag_coeff(m, rate, vol, arguments_down));
			//W = subdiag_coeff(m, rate, vol, arguments)/diag_coeff(m, rate,vol, arguments_down);
			//B[i] = diag_coeff(m, rate,vol, arguments)-W*updiag_coeff(m,rate,vol,arguments_down);
			//D[i] = col_up[i] - W*col_up[i - 1];
		}
		std::size_t i( col_up.size() - 1);
		B[i] = B[i - 1] * updiag_coeff(m, rate, vol, arguments_down);
		D[i] = (col_up[i] - D[i - 1] * updiag_coeff(m, rate, vol, arguments_down)) / (diag_coeff(m, rate, vol, arguments_down) - B[i - 1] * updiag_coeff(m, rate, vol, arguments_down));

		result[col_up.size() - 1] = (col_up[col_up.size()-1]);
		//result[i] = (D[i]);
		for (std::size_t i = col_up.size() - 2; i >0; i--) {
			arguments[0] = arguments[3] + i*m.get_mesh_dx();
			result[i] = D[i] - B[i] * result[i + 1];
			//result[i] =( (D[i] - updiag_coeff(m, rate, vol, arguments)*result[i + 1])/B[i]);
			i = i;
		}
		return result;
	}
	std::vector<double> price_today(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff) {
		std::vector<double> ini_price(initial_price_vector(m, rate, vol, arguments, payoff));
		std::vector<double> col_up(column_up(m, rate, vol, arguments, payoff, ini_price));
		double dt = m.get_mesh_dt();
		std::vector<double> result2 = col_up;
		int nb_step =floor( arguments[1] / dt);
		std::vector<double> result1(col_up.size());
		for (int i = 1; i < nb_step; i++) {
			arguments[1] = arguments[1] - m.get_mesh_dt();
			result1=(price_vector(m, rate, vol, arguments, payoff, result2));
			result2 = column_up(m,rate,vol,arguments,payoff,result1);
		}
		return result1;


	}
	

}

