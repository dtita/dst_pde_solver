#include "solver.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

namespace dauphine
{
	//mesh::mesh()
	//	: m_dt(1), m_dx(1), m_maturity(1), m_spot_boundaries([1,1])
	//{
	//}
	mesh::mesh(double dt, double dx, double maturity, double spot,double theta, std::vector<double> boundaries)
		: m_dt(dt), m_dx(dx), m_maturity(maturity),m_spot(spot),m_theta(theta), m_spot_boundaries(boundaries)
	{
	}

	mesh::~mesh()
	{
	}
	double mesh::get_mesh_maturity() const {
		return m_maturity;
	}
	double mesh::get_mesh_dt() const {
		return m_dt;
	}
	double mesh::get_mesh_dx() const {
		return m_dx;
	}
	double mesh::get_mesh_spot() const {
		return m_spot;
	}
	double mesh::get_mesh_theta() const {
		return m_theta;
	}
	std::vector<double> mesh::get_mesh_spot_boundaries() const {
		return m_spot_boundaries;
	}
	std::vector<double> mesh::spot_vector() {
		int size = floor((m_spot_boundaries[0] - m_spot_boundaries[1]) / m_dx)+1;
		std::vector<double> result(size);
		for (std::size_t i = 0; i < result.size(); ++i)
		{
			result[i] =( m_spot_boundaries[1]+i*m_dx);
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

	double diag_coeff(mesh m, initial_function rate, std::vector<double> arguments) {
		if (arguments[0] == arguments[2]|| arguments[0] == arguments[3])
		{
			return 1;
		}
		else
		{
			return 1 / m.get_mesh_dt() + rate.function_operator(arguments)*m.get_mesh_theta()+1/(m.get_mesh_dx()*m.get_mesh_dx())*m.get_mesh_theta();
		}
	}
	double subdiag_coeff(mesh m, initial_function rate,initial_function vol, std::vector<double> arguments) {
		if (arguments[0] == arguments[2])
		{
			return 0;
		}
		else
		{
			return -1 /2*m.get_mesh_theta()/ (m.get_mesh_dx()*m.get_mesh_dx()) +1/(4*m.get_mesh_dx())*m.get_mesh_theta()*(pow(vol.function_operator(arguments),2)-rate.function_operator(arguments));
		}
	}
	double updiag_coeff(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments) {
		if (arguments[0] == arguments[2])
		{
			return 0;
		}
		else
		{
			return -1 / 2 * m.get_mesh_theta() / (m.get_mesh_dx()*m.get_mesh_dx()) + 1 / (4 * m.get_mesh_dx())*m.get_mesh_theta()*(-pow(vol.function_operator(arguments), 2) +rate.function_operator(arguments));
		}
	}

	std::vector<double> initial_price_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff) {
		std::vector<double> result=m.spot_vector();
		for (std::size_t i = 0; i < result.size(); i++) {
			arguments[0] = result[i];
			if (payoff.function_operator(arguments) == 0) {
				result[i] = 0;
			}
			else {
				result[i] = log(payoff.function_operator(arguments));
			}
		}
		return result;
	}

}
//std::vector<double> bound = m.get_mesh_spot_boundaries();
//int size = (bound[0] - bound[1]) / m.get_mesh_dx();
//std::vector<double> result(size);
//	double rate::function_operator(std::vector<double> arguments) {
//		return m_f(arguments);
//	}
//	rate::~rate() {
//	}
//
//	price_vector::price_vector() {
//	}
//}

