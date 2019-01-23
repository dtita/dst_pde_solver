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
        
        double bound_up(double f,std::vector<double> arguments,rates_const rate, mesh m) const;
        double bound_down(double f) const;
        
        
    };
    
}

#endif /* boundaries_hpp */
