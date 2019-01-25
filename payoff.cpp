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
double payoff::get_payoff(const double& fwd) const
{
    return 0; 
}

payoff::~payoff()
{
}

bs_call::bs_call(const double& strike)
    : m_strike(strike)
{
}
    
//Declare the fonction for the payoff
double bs_call::get_payoff(const double& fwd) const
{        
    return std::max(fwd-m_strike,0.); //Call option with strike price m_strike
}

bs_call::~bs_call()
{
}
    
    
}
