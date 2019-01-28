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
#include <vector>

namespace dauphine
{
    class mesh
    {
    public:
        //mesh();
        mesh(double dt, int dx,double maturity, double spot, std::vector<double> spot_boundaries);
        // No const required in the return type since you return values
        double const get_mesh_dt() const;
        double const get_mesh_maturity() const;
        double const get_mesh_dx() const;
        double const get_mesh_spot() const;
        // These should be private memebers with
        // constant getters
        std::vector<double> spot_vect;
        std::vector<double> t_vect;
        double d_x;
        ~mesh();
    private:
        // no need for const here
        const double m_dt;
        // Why storing m_dx AND d_x? This is confusing and m_dx is never used
        double m_dx;
        double m_maturity;
        double m_spot;
        std::vector<double> m_spot_boundaries;
 
    };
    
}

#endif /* mesh_hpp */
