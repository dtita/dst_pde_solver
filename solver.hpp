#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

namespace dauphine
{
	class mesh {
	public:
		//mesh();
		mesh(double dt, double dx,double maturity, std::vector<double> spot_boundaries);
		double get_mesh_dt() const;
		double get_mesh_maturity() const;
		double get_mesh_dx() const;
		std::vector<double> get_mesh_spot_boundaries() const;
		~mesh();
	private:
		double m_dt;
		double m_dx;
		double m_maturity;
		std::vector<double> m_spot_boundaries;
	};


	//essayer de voir si on peut faire une famille avec les classes payoff, vol, rates, boundaries
	class payoff {
	public:
		payoff(double(*f)(std::vector<double>));
		double function_operator(std::vector<double> arguments);
		~payoff();
	private:
		double (*m_f)(std::vector<double>);
	};

	class rate {
	public:
		rate(double(*f)(std::vector<double>));
		double function_operator(std::vector<double> arguments);
		~rate();
	private:
		double(*m_f)(std::vector<double>);
	};

	class volatility {
	public:
		volatility(double(*f)(std::vector<double>));
		double function_operator(std::vector<double> arguments);
		~volatility();
	private:
		double(*m_f)(std::vector<double>);
	};

	class boundaries {
	public:
		boundaries(double(*f)(std::vector<double>));
		std::vector<double> function_operator(std::vector<double> arguments);
		~boundaries();
	private:
		std::vector<double>(*m_f)(std::vector<double>);
	};


	class price_vector {
	public: 
		price_vector();
	private:
	};

}

#endif
