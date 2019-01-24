




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
double rates::get_rates(const double& time, const double& spot) const

{
    return 0.0; //Constant rate (equal to 0.0)
}

rates::~rates()
{
}

rates_const::rates_const(const double& rate)
    :m_rate(rate)
{
}


double rates_const::get_rates(const double& time, const double& spot) const
{
    return m_rate; 
}

rates_const::~rates_const()
{
}


}

