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
	double boundaries::get_boundaries(const double& f, const double& time, const double& spot, const rates_const& rate, const mesh& m) const
	{
		//spot: Use in return if the boundarie volatility is path-dependent (depend on spot S)
		//time: Use in return if the boundarie is time-dependent (depend on time t)
		//rate: Use in return if the boundarie volatility is rate-dependent 
		//vol: Use in return if the boundarie is vol-dependent 

		return 0; //Here in the exemple the vol is constant and equal 0%
	}
    
    boundaries::~boundaries()
    {
    }
    
    //Declare the fonction for the implement the boundary conditions
	bound_up_dirichlet::bound_up_dirichlet()
	{
	}

	bound_up_dirichlet::~bound_up_dirichlet()
	{
	}

	bound_down_dirichlet::bound_down_dirichlet()
	{
	}
	bound_down_dirichlet::~bound_down_dirichlet()
	{
	}

    double bound_up_dirichlet::get_boundaries(const double& f, const double& time, const double& spot,const rates_const& rate,const mesh& m) const
    {
        return f*exp(-rate.get_rates(time,spot)*m.get_mesh_dt());
    }
    
    double bound_down_dirichlet::get_boundaries(const double& f, const double& time, const double& spot, const rates_const& rate, const mesh& m) const
    {
        return f;
    }
}
