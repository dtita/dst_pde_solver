#include <iostream>
#include "closed_form.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "volatility.hpp"
#include "rates.hpp"
#include "payoff.hpp"
#include "boundaries.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace dauphine
{
	void test()
	{
		// make sure the arguments of the payoff are properly defined
		int number_arguments(5);
		double spot=100.;
		double maturity=1.;
		
		//the client should be able to specify the value of θ
		double theta = 0.5;

        //Nbr of points
        int nb_x = 1001;
    
        //Creation mesh
        
	//The client should be able to specify to specify the mesh
	//mesh(double dt, double dx,double maturity(in years), double spot,boundaries)
        
        std::vector<double> mesh_boundaries(2);
        mesh_boundaries[0] = std::exp(std::log(spot) - 1.);
        mesh_boundaries[1] = std::exp(std::log(spot) + 1.);

        mesh m(1./365.,1001,1.,100.,mesh_boundaries);
        
    
    //Creation vol
        volatility vol;
        
	//Creation rates
        rates rate;

	//Creation payoff
        payoff p;
        
    // Creation boundaries
        
        boundaries bnd;


    //Compute price
		std::vector<double> result = price_today(theta,m,rate,vol,p,bnd,false); // Use true if time-dependent
		
    // Spot index
        int i = (nb_x - 1) / 2;
        
    //Greeks
		double price = result[i];
		double delta = (result[i] - result[i-1]) / (m.spot_vect[i] - m.spot_vect[i-1]);
		double gamma = (result[i + 1] - 2 * result[i] + result[i - 1]) / (pow((m.spot_vect[i+1] - m.spot_vect[i-1])/2.0, 2));
		double theta_product = (result[0]-price);
		// for the vega, we should have a discretization w.r.t the volatility
        
    //Print results (Price and Greeks)
		std::cout << "Price: " <<price << std::endl;
		std::cout << "Delta: " << delta << std::endl;
		std::cout << "Gamma: " << gamma << std::endl;
		std::cout << "Theta: " << theta_product << std::endl;

	}
}


int main(int argc, char* argv[])
{
    std::cout <<"Price BS: " << dauphine::bs_price(100,100,0.20,1.0,true) << std::endl;
	dauphine::test();
    return 0;
}
