#include <iostream>
#include "closed_form.hpp"
#include "solver.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace dauphine
{	
	//Initial functions to be changed according to the user's needs
	//modify this function to modify the payoff, arguments=[spot,maturity,mesh_boundaries_up,mesh_boundaries_down]
	double payoff_function(std::vector<double> arguments) {
		return std::max(arguments[0]-100,0.);
	}
	double rate_function(std::vector<double> arguments) {
		return 0.01;
	}
	double volatility_function(std::vector<double> arguments) {
		return 0.20;
	}
	double boundaries_up(std::vector<double> arguments) {
		return  std::max(arguments[2] - 100, 0.);
	}
	double boundaries_down(std::vector<double> arguments) {
		return std::max(arguments[3] - 100, 0.);
	}

	void test()
	{
		// make sure the arguments of the payoff are properly defined
		int number_arguments(4);
		double spot=100;
		double strike=100;
		double mesh_up_boundaries=150;
		double mesh_down_boundaries=50;

		std::vector<double> arguments(number_arguments);
		arguments[0] = spot;
		arguments[1] = strike;
		arguments[2] = mesh_up_boundaries;
		arguments[3] = mesh_down_boundaries;
		initial_function payoff(payoff_function);
		initial_function rate(rate_function);
		initial_function volatility(volatility_function);
		initial_function up_boundaries(boundaries_up);
		initial_function down_boundaries(boundaries_down);

		std::vector<double> mesh_boundaries(2);
		mesh_boundaries[0] = 0;
		mesh_boundaries[1] = 15;
		mesh m(1,1,1,100,0.5,mesh_boundaries);
		std::cout << payoff.function_operator(arguments) << std::endl;
		std::cout << diag_coeff(m, rate,arguments) << std::endl;

	}
}
int main(int argc, char* argv[])
{
	std::cout << dauphine::bs_price(100,100,0.2,5,true) << std::endl;
	dauphine::test();
    return 0;
}
