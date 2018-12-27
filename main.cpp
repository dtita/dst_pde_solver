#include <iostream>
#include "closed_form.hpp"
#include "solver.hpp"
#include "solverT.hpp"
#include "params.hpp"
#include "mesh.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace dauphine
{	
	//Initial functions to be changed according to the user's needs
	//modify this function to modify the payoff, arguments=[spot,maturity,mesh_boundaries_up,mesh_boundaries_down,theta]
	double payoff_function(std::vector<double> arguments) {
		return std::max(arguments[0]-100,0.);
	}
	double rate_function(std::vector<double> arguments) {
		return 0.05/365;
	}
	double volatility_function(std::vector<double> arguments) {
		return 0.;
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
		int number_arguments(5);
		double spot=100.;
		double maturity=1./365;
		double mesh_up_boundaries=150.;
		double mesh_down_boundaries=50.;
		double theta = 0.5;

		std::vector<double> arguments(number_arguments);
		arguments[0] = spot;
		arguments[1] = maturity;
		arguments[2] = mesh_up_boundaries;
		arguments[3] = mesh_down_boundaries;
		arguments[4] = theta;
		std::vector<double> mesh_boundaries(2);
		mesh_boundaries[0] = 150.0;
		mesh_boundaries[1] = 50.0;

		initial_function payoff(payoff_function);
		initial_function rate(rate_function);
		initial_function volatility(volatility_function);
		initial_function up_boundaries(boundaries_up);
		initial_function down_boundaries(boundaries_down);

        
        mesh m(1.,1,1./365,100.,mesh_boundaries);
        //mesh m(1.,1,1.,100.,mesh_boundaries);

		std::vector<double> result = price_today(m,rate,volatility,arguments,payoff);
        //
        //std::cout <<"Price: "<<result[80] << std::endl;
        //
        ////std::vector<double> col=a_vector(1,10);
        //
        //for (std::size_t i = 0; i <= result.size(); ++i)
        //{
        //    std::cout << i<<": "<<result[i] << std::endl;
        //}
        
	}
}


int main(int argc, char* argv[])
{
	std::cout <<"Price BS: " << dauphine::bs_price(100,100,0.2,1,true) << std::endl;
	dauphine::test();
    return 0;
}
