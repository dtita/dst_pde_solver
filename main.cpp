#include <iostream>
#include "closed_form.hpp"
#include "solver.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace dauphine
{
	//modify this function to modify the payoff
	double payoff_function(std::vector<double> arguments) {
		return std::max(arguments[0]-arguments[1],0.);
	}
	double rate_function(std::vector<double> arguments) {
		return arguments[0];
	}
	double volatility_function(std::vector<double> arguments) {
		return arguments[0];
	}

	double boundaries_up(std::vector<double> arguments) {
		return arguments[0];
	}
	double boundaries_down(std::vector<double> arguments) {
		return arguments[0];
	}

	void test()
	{
		// make sure the arguments of the payoff are properly defined
		int number_arguments_payoff(2);
		double fwd=100;
		double strike=100;

		std::vector<double> arguments(number_arguments_payoff);
		arguments[0] = fwd;
		arguments[1] = strike;
		initial_function payoff(payoff_function);
		initial_function rate(rate_function);
		initial_function volatility(volatility_function);
		initial_function up_boundaries(boundaries_up);
		initial_function down_boundaries(boundaries_down);

		std::vector<double> mesh_boundaries(2);
		mesh_boundaries[0] = 0;
		mesh_boundaries[1] = 15;
		mesh m(1,1,1,mesh_boundaries);
		std::cout << payoff.function_operator(arguments) << std::endl;

	}
}
int main(int argc, char* argv[])
{
	std::cout << dauphine::bs_price(100,100,0.2,5,true) << std::endl;
	dauphine::test();
    return 0;
}
