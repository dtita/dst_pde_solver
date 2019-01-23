#ifndef rates_hpp
#define rates_hpp

#include <vector>

namespace dauphine
{
    class rates
    {
        public:
            explicit rates();
            virtual ~rates();
       
    };
    
    class rates_const : public rates
    {
    public:
        double get_rates(const double& time, const double& spot) const;
        
        
    };

}

#endif 
