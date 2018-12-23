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

// Class regroupant les parametres necessaires a la construction de la mesh et resolution du systeme

namespace dauphine
{
    class params{
        
        public:
            params(double dt, double dx,double maturity, double spot,double theta, std::vector<double> spot_boundaries);
            double get_dt() const;
            double get_maturity() const;
            double get_dx() const;
            double get_spot() const;
            double get_theta() const;
            std::vector<double> get_spot_boundaries() const;
            virtual ~params();
        
        protected:
            double m_dt;
            double m_dx;
            double m_maturity;
            double m_spot;
            double m_theta;
            std::vector<double> m_spot_boundaries;
        
    };

// class permettant de generer des fonctions pour la vol, rate...
 //   class function{
 //       public:
 //           new_function(double(*f)(std::vector<double>));
 //           double function_operator(std::vector<double> arguments);
 //           virtual ~new_function();
 //       private:
//            double(*m_f)(std::vector<double>);
        
 //   };
}

#endif /* params_hpp */
