#ifndef payoff_hpp
#define payoff_hpp

#include <vector>


namespace dauphine
{
    class payoff
    {
        public:
            explicit payoff();
            virtual ~payoff();
       
    };
    
    class bs_call: public payoff
    {
    public:
        
        double get_payoff(const double& time, const double& spot) const;
    };

}
#endif 
