//
//  volatility.cpp
//  
//
//  Created by David Tiobang on 21/01/2019.
//

#include "volatility.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>
//#include "solver.hpp"

namespace dauphine
{
volatility::volatility()
{
}
    
double volatility::get_volatility(const double& time, const double& spot) const
{
    //spot: Use in return if the volatility is path-dependent (depend on spot S)
    //time: Use in return if the volatility is time-dependent (depend on time t)
        
        return 0; //Here in the exemple the vol is constant and equal 0%
}

volatility::~volatility()
{
}

vol_const::vol_const(const double& vol)
    :m_vol(vol)
{
}

//Declare the fonction for the volatility
double vol_const::get_volatility(const double& time, const double& spot) const
{
    return m_vol; 
}

vol_const::~vol_const()
{
}
    
    
    
}
