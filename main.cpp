#include <iostream>
#include "closed_form.hpp"

int main(int argc, char* argv[])
{
	std::cout << "hello" << std::endl;
	std::cout << dauphine::vanilla_payoff(120,100,true) << std::endl;
	std::cout << dauphine::bs_price(100,100,0.2,5,true) << std::endl;
    return 0;
}
