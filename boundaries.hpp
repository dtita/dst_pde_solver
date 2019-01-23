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
        virtual ~boundaries();

    };
    
    class bound_dirichlet : public boundaries
    {
    public:
        
        double bound_up(const double& f, const double& time, const double& spot, const rates_const& rate, const mesh& m) const;
        double bound_down(const double& f) const;
        
        
    };
    
}

#endif /* boundaries_hpp */
