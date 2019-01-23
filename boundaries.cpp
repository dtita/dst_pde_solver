//
//  boundaries.cpp
//  
//
//  Created by David Tiobang on 22/01/2019.
//

#include "boundaries.hpp"
#include "mesh.hpp"
#include "rates.hpp"
//#include "solver.hpp"
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
    
    //Declare the fonction for the implement the boundary conditions
    double boundaries::bound_up(double f,std::vector<double> arguments,rates rate, mesh m) const
    {
        //arguments[0]: Use in return if the condition is path-dependent (depend on spot S)
        //arguments[1]: Use in return if the condition is time-dependent (depend on time t)
        
        return f*exp(-rate.get_rates(arguments)*m.get_mesh_dt());
    }
    
    double boundaries::bound_down(double f) const
    {
        //arguments[0]: Use in return if the condition is path-dependent (depend on spot S)
        //arguments[1]: Use in return if the condition is time-dependent (depend on time t)
        
        return f;
    }

    
    boundaries::~boundaries()
    {
    }
    
}
