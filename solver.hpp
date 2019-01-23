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
    std::vector<double> initial_price_vector(mesh m, payoff p);
    
    //Coeffs Matrix
    std::vector<double> up_vector(mesh m, rates rate,volatility vol, std::vector<double> arguments);
    std::vector<double> sub_vector(mesh m, rates rate,volatility vol, std::vector<double> arguments);
    std::vector<double> diag_vector(mesh m, rates rate,volatility vol, std::vector<double> arguments);
    
    //Tridiag solver
    std::vector<double> tridiagonal_solver(std::vector<double>  a, std::vector<double>  b, std::vector<double>  c, std::vector<double>  f,boundaries bnd);

    // Result vector
    std::vector<double> price_today(double theta, mesh m, rates rate, volatility vol, payoff p,boundaries bnd, bool time_S_dependent);
}

#endif
