#ifndef rates_hpp
#define rates_hpp

#include <vector>



namespace dauphine
{
    class rates
    {
        public:
            explicit rates();
            virtual double get_rates(const double& time, const double& spot) const;
            virtual ~rates();
       
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
