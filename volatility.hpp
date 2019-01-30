//
//  volatility.hpp
//  
//
//  Created by David Tiobang on 21/01/2019.
//

#ifndef volatility_hpp
#define volatility_hpp

#include <vector>
//#include "solver.hpp"
namespace dauphine
{
    class volatility
    {
        public:
            // No need for explicit since the constructor does not accept any parameter
            explicit volatility();
            // Should be pure virtual method
            virtual double get_volatility(const double& time, const double& spot) const;
            virtual ~volatility();

            // Missing entity semantic: explicitly delete copy and move semantic
       
    };
    
    class vol_const : public volatility
    {
    public:
         explicit vol_const(const double& vol);
         virtual double get_volatility(const double& time, const double& spot) const;
         virtual ~vol_const();
    private:
        double m_vol;
        
    };

}
#endif /* volatility_hpp */
