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
            explicit volatility();
            virtual double get_volatility(const double& time, const double& spot) const;
            virtual ~volatility();
       
    };
    
    class vol_const : public volatility
    {
    public:
         explicit vol_const(const double& vol);
         double get_volatility(const double& time, const double& spot) const;
         virtual ~vol_const();
    private:
        double m_vol;
        
    };

}
#endif /* volatility_hpp */
