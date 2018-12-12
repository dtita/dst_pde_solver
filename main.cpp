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
	double rates_function(std::vector<double> arguments) {
		return arguments[0];
	}
	double volatility_function(std::vector<double> arguments) {
		return arguments[0];
	}

	std::vector<double> boundaries_function(std::vector<double> arguments) {
		return arguments;
	}

	void test()
	{
		// make sure the arguments of the payoff are properly defined
		int number_arguments_payoff(2);
		std::vector<double> arguments(number_arguments_payoff);
		arguments[0] = 320;
		arguments[1] = 50;

		payoff call_payoff(payoff_function);

		std::vector<double> boundaries(2);
		boundaries[0] = 0;
		boundaries[1] = 15;
		mesh m(1,1,1,boundaries);
		std::cout << call_payoff.function_operator(arguments) << std::endl;

	}
}
int main(int argc, char* argv[])
{
	std::cout << dauphine::bs_price(100,100,0.2,5,true) << std::endl;
	dauphine::test();
    return 0;
}
