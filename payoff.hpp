#ifndef payoff_hpp
#define payoff_hpp

#include <vector>
#include "rates.hpp"



namespace dauphine
{
    class payoff
    {
        public:
            explicit payoff();
            virtual double get_payoff(const double& fwd) const;
            virtual ~payoff();
       
    };
    
    class bs_call: public payoff
    {
    public:
        explicit bs_call(const double& strike);
        virtual double get_payoff(const double& fwd) const;
        virtual ~bs_call();
    private:
        double m_strike;
    };

}
#endif 
