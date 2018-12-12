#include <iostream>
#include "closed_form.hpp"
#include "solver.hpp"

namespace dauphine
{
	void test()
	{
		std::vector<double> boundaries(2);
		boundaries[0] = 0;
		boundaries[1] = 15;
		mesh m(1,1,1,boundaries);
		std::cout << m.get_mesh_maturity() << std::endl;
	}
}
int main(int argc, char* argv[])
{
	std::cout << "hello" << std::endl;
	std::cout << dauphine::vanilla_payoff(120,100,true) << std::endl;
	std::cout << dauphine::bs_price(100,100,0.2,5,true) << std::endl;
	dauphine::test();
    return 0;
}
