//
//  params.hpp
//  
//
//  Created by David Tiobang on 14/12/2018.
//
#ifndef params_hpp
#define params_hpp
//#ifndef SOLVER_HPP
//#define SOLVER_HPP

#include <vector>
#include <stdio.h>

namespace dauphine
{
    class params{
        
        public:
            params(double dt, double dx,double maturity, double spot, std::vector<double> spot_boundaries);
            double get_dt() const;
            double get_maturity() const;
            double get_dx() const;
            double get_spot() const;
            std::vector<double> get_spot_boundaries() const;
            virtual ~params();
        
        protected:
            double m_dt;
            double m_dx;
            double m_maturity;
            double m_spot;
            std::vector<double> m_spot_boundaries;
        
    };
    
}

#endif /* params_hpp */
