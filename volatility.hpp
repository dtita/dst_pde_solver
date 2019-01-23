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
            virtual ~volatility();
       
    };
    
    class vol_const : public volatility
    {
    public:
        double get_volatility(std::vector<double> arguments) const;
        
    };

}
#endif /* volatility_hpp */
