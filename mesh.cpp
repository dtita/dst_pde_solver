//
//  mesh.cpp
//  
//
//  Created by David Tiobang on 16/12/2018.
//

#include "mesh.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>

namespace dauphine
{
    mesh::mesh( double dt, int dx, double maturity, double spot, std::vector<double> boundaries)
    : m_dt(dt), m_dx(dx), m_maturity(maturity), m_spot(spot), m_spot_boundaries(boundaries)
    {
        std::vector<double> result(dx,0.);
        int inter_t =floor( maturity / dt);
        std::vector<double> result2(inter_t);
        
        double log_spot_min = std::log(boundaries[0]);
        double log_spot_max = std::log(boundaries[1]);
        double dlog = (log_spot_max - log_spot_min) /( dx-1);
        
        for (int i = 0; i < dx; ++i)
        {
            result[i] = std::exp(log_spot_min + i * dlog);
        }
        // Would be more efficient to directly compute spot_vect
        spot_vect = result;
        d_x = dlog;
        
        for (int i = 0; i < inter_t; ++i)
        {
            result2[i] = maturity-i*dt;
        }
        // Would be more efficient to directly compute t_vect
        t_vect = result2;
    }
    mesh::~mesh()
    {
    }
    double const mesh::get_mesh_maturity() const {
        return m_maturity;
    }
    double const mesh::get_mesh_dt() const {
        return m_dt;
    }
    double const mesh::get_mesh_dx() const {
        return m_dx;
    }
    double const mesh::get_mesh_spot() const {
        return m_spot;
    }

}
