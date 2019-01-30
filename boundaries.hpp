//
//  boundaries.hpp
//  
//
//  Created by David Tiobang on 22/01/2019.
//

#ifndef boundaries_hpp
#define boundaries_hpp
#include "rates.hpp"

#include "mesh.hpp"
#include <stdio.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <vector>

namespace dauphine
{
    class boundaries
    {
    public:
        boundaries();
        // Should be virtual pure.
		virtual double get_boundaries(const double& f, const double& time, const double& spot, const rates& rate, const mesh& m) const;
        virtual ~boundaries();

        // Missing entity semantic: explicitly delete copy and move semantic
    };
    
    class bound_up_dirichlet : public boundaries
    {
    public:
		explicit bound_up_dirichlet();
        virtual double get_boundaries(const double& f, const double& time, const double& spot, const rates& rate, const mesh& m) const;
		virtual ~bound_up_dirichlet(); 
    };
	class bound_down_dirichlet : public boundaries
	{
	public:
		explicit bound_down_dirichlet();
		virtual double get_boundaries(const double& f, const double& time, const double& spot, const rates& rate, const mesh& m) const;
		virtual ~bound_down_dirichlet();
	};
    
}

#endif /* boundaries_hpp */
