#include "payoff.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>

namespace dauphine
{
payoff::payoff()
{
}

//Declare the fonction for the payoff
double payoff::get_payoff(const double& time, const double& spot) const
{
    return 0; //Call option with strike price 100
}

payoff::~payoff()
{
}

bs_call::bs_call()
{
}
    
//Declare the fonction for the payoff
double bs_call::get_payoff(const double& time, const double& spot) const
{        
    return std::max(spot-100,0.); //Call option with strike price 100
}

bs_call::~bs_call()
{
}
    
    
}
