#include "rates.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>

namespace dauphine
{
rates::rates()
{
}

rates::~rates()
{
}

//Declare the fonction for the rates
double rates_const::get_rates(const double& time,const double& spot) const
{
    return 0.0; //Constant rate (equal to 0.0)
}
}

