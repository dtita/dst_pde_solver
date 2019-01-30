#ifndef rates_hpp
#define rates_hpp

#include <vector>



namespace dauphine
{
    class rates
    {
        public:
            // No need for explicit since the constructor does not accept any parameter
            explicit rates();
            // Should be pure virtual method
            virtual double get_rates(const double& time, const double& spot) const;
            virtual ~rates();
       
            // Missing entity semantic: explicitly delete copy and move semantic
    };
    
    class rates_const : public rates
    {
    public:
        explicit rates_const(const double& rate);
        virtual double get_rates(const double& time, const double& spot) const;
        virtual ~rates_const();
     private:
        double m_rate;
        
    };

}

#endif 
