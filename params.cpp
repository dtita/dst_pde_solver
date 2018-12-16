//
//  params.cpp
//  
//
//  Created by David Tiobang on 14/12/2018.
//

#include "params.hpp"
#include "solver.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

namespace dauphine
{
    
    params::params(double dt, double dx, double maturity, double spot, std::vector<double> boundaries)
    : m_dt(dt), m_dx(dx), m_maturity(maturity),m_spot(spot), m_spot_boundaries(boundaries)
    {
    }
    
    params::~params()
    {
    }
    double params::get_maturity() const {
        return m_maturity;
    }
    double params::get_dt() const {
        return m_dt;
    }
    double params::get_dx() const {
        return m_dx;
    }
    double params::get_spot() const {
        return m_spot;
    }
    
    
}
