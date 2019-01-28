#ifndef payoff_hpp
#define payoff_hpp

#include <vector>
#include "rates.hpp"



namespace dauphine
{
    class payoff
    {
        public:
            // No need for explicit since the constructor does not accept any parameter
            explicit payoff();
            // Should be pure virtual method
            virtual double get_payoff(const double& fwd) const;
            virtual ~payoff();

            // Missing entity semantic: explicitly delete copy and move semantic
       
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
