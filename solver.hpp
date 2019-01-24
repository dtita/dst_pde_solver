#ifndef SOLVER_HPP
#define SOLVER_HPP
#include "mesh.hpp"
#include "volatility.hpp"
#include "rates.hpp"
#include "payoff.hpp"
#include "boundaries.hpp"

#include <vector>

namespace dauphine
{
    //Creation initial price vector
    std::vector<double> initial_price_vector(const mesh& m, const payoff& p);
    
    //Coeffs Matrix
    std::vector<double> up_vector(const mesh& m, const rates& rate, const volatility& vol, const double& time, double spot, const double& theta);
    std::vector<double> sub_vector(const mesh& m, const rates& rate, const volatility& vol, const double& time, double spot, const double& theta);
	std::vector<double> diag_vector(const mesh& m, const rates& rate, const volatility& vol, const double& time, double spot, const double& theta);
    
    //Tridiag solver
    std::vector<double> tridiagonal_solver(const std::vector<double>&  a, std::vector<double>  b,  const std::vector<double>&  c, std::vector<double>  f);

    // Result vector
    std::vector<double> price_today(const double& theta, const mesh& m, const rates& rate, const volatility& vol, const  payoff& p, const boundaries& bnd_down, const boundaries& bnd_up, const bool& time_S_dependent);
}

#endif
