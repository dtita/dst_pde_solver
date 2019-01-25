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

bs_call::bs_call(const double& strike)
    : m_strike(strike)
{
}
    
//Declare the fonction for the payoff
double bs_call::get_payoff(const double& time, const double& spot) const
{        
    return std::max(spot-m_strike,0.); //Call option with strike price m_strike
}

bs_call::~bs_call()
{
}
    
    
}
