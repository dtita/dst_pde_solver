//
//  boundaries.cpp
//  
//
//  Created by David Tiobang on 22/01/2019.
//

#include "boundaries.hpp"
#include "mesh.hpp"
#include "rates.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <vector>

namespace dauphine
{
    boundaries::boundaries()
    
    {
        
    }
    
    
    boundaries::~boundaries()
    {
    }
    
    //Declare the fonction for the implement the boundary conditions
    double bound_dirichlet::bound_up(const double& f, const double& time, const double& spot,const rates_const& rate,const mesh& m) const
    {
        return f*exp(-rate.get_rates(time,spot)*m.get_mesh_dt());
    }
    
    double bound_dirichlet::bound_down(const double& f) const
    {

        return f;
    }
}
