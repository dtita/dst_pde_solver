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
        double bound_up(double f,std::vector<double> arguments,rates rate, mesh m) const;
        double bound_down(double f) const;
        virtual ~boundaries();

    };
    
}

#endif /* boundaries_hpp */
