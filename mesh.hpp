//
//  mesh.hpp
//  
//
//  Created by David Tiobang on 16/12/2018.
//

#ifndef mesh_hpp
#define mesh_hpp
//#ifndef SOLVER_HPP
//#define SOLVER_HPP
#include <stdio.h>
#include "params.hpp"

#include <vector>

namespace dauphine
{
    class Mesh
    {
    public:
        
        Mesh(double dt, double dx,double maturity, double spot,double theta, std::vector<double> spot_boundaries, std::vector<double> boundaries);
        std::vector<double> spot_vector(params p);
        std::vector<double> get_mesh_spot_boundaries() const;
        virtual ~Mesh();
        
    private:
       
        std::vector<double> m_spot_boundaries;
        params p;
    };
    
}

#endif /* mesh_hpp */
