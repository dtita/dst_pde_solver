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
	mesh::mesh(double dt, double dx, double maturity,std::vector<double> boundaries)
		: m_dt(dt), m_dx(dx), m_maturity(maturity),m_spot_boundaries(boundaries)
	{
	}

	mesh::~mesh()
	{
	}
	double mesh::get_mesh_maturity() const{
		return m_maturity;
	}
	double mesh::get_mesh_dt() const {
		return m_dt;
	}
	double mesh::get_mesh_dx() const {
		return m_dx;
	}
	std::vector<double> mesh::get_mesh_spot_boundaries() const {
		return m_spot_boundaries;
	}
	payoff::payoff(double(*f)(std::vector<double>))
		: m_f(f)
	{
	}
	double payoff::function_operator(std::vector<double> arguments) {
		return m_f(arguments);
	}
	payoff::~payoff() {
	}
	double rate::function_operator(std::vector<double> arguments) {
		return m_f(arguments);
	}
	rate::~rate() {
	}

	price_vector::price_vector() {
	}
}

