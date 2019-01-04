#include "solverT.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>
namespace dauphine
{

	mesh::mesh(double dt, double dx, double maturity, double spot, std::vector<double> boundaries)
		: m_dt(dt), m_dx(dx), m_maturity(maturity), m_spot(spot), m_spot_boundaries(boundaries)
	{
		double size = floor((m_spot_boundaries[0] - m_spot_boundaries[1]) / m_dx) + 1;
		std::vector<double> result(size);
		std::vector<double> result2(size - 1);

		for (std::size_t i = 0; i < size; ++i)
		{
			result[i] = m_spot_boundaries[1] + i*m_dx;

		}
		spot_vect = result;

		for (std::size_t i = 0; i < size - 1; ++i)
		{
			result2[i] = log(spot_vect[i + 1]) - log(spot_vect[i]);

		}
		d_x = result2;


	}

	mesh::~mesh()
	{
	}
	double mesh::get_mesh_maturity() const {
		return m_maturity;
	}
	double mesh::get_mesh_dt() const {
		return m_dt / 365;
	}
	double mesh::get_mesh_dx() const {
		return m_dx;
	}
	double mesh::get_mesh_spot() const {
		return m_spot;
	}

	std::vector<double> mesh::get_mesh_spot_boundaries() const {
		return m_spot_boundaries;
	}


	std::vector<double> mesh::spot_vector() {
		int size = floor((m_spot_boundaries[0] - m_spot_boundaries[1]) / m_dx) + 1;
		std::vector<double> result(size);
		for (std::size_t i = 0; i < result.size(); ++i)
		{
			result[i] = log(m_spot_boundaries[1] + i*m_dx);
		}
		return result;
	}


	initial_function::initial_function(double(*f)(std::vector<double>))
		: m_f(f)
	{
	}
	double initial_function::function_operator(std::vector<double> arguments)
	{
		return m_f(arguments);
	}

	initial_function::~initial_function()
	{
	}



	// les fonctions ici sont pour les calculs des coeffs de la matrice tridiagonale sub = diag dessous, diag = diagonale,
	// il y a peut être une erreur dans les coeffs donc il faudrait qu'on fasse tous les trois le calcul et on compare pour être sur
	std::vector<double> diag_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments)
	{
		//std::vector<double> a=m.spot_vector();
		std::vector<double> a = m.spot_vect;
		long size = a.size();
		std::vector<double> result(size, 0.0);
		result[0] =1.0;
		result[size-1] =1.0;
		for (std::size_t i = 1; i < size - 1; ++i)
		{
			//result[i]=1.0+arguments[4]*m.get_mesh_dt()*((pow(vol.function_operator(arguments),2)/pow(m.get_mesh_dx(),2))+rate.function_operator(arguments));

			result[i] = 1.0 + arguments[4] * m.get_mesh_dt()*((pow(vol.function_operator(arguments), 2) / pow(m.d_x[i], 2)) + rate.function_operator(arguments));

		}
		return result;
	}


	std::vector<double> sub_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments)
	{
		//std::vector<double> a=m.spot_vector();
		std::vector<double> a = m.spot_vect;
		long size = a.size();
		std::vector<double> result(size, 0.0);
		result[size-1] =0;
		result[0] = 0;
		for (std::size_t i = 1; i < size - 1; ++i)
		{
			//result[i]=-0.5*arguments[4]*m.get_mesh_dt()*((pow(vol.function_operator(arguments),2)/pow(m.get_mesh_dx(),2))+((pow(vol.function_operator(arguments),2)-rate.function_operator(arguments))/(2.0*m.get_mesh_dx())));

			result[i] = -0.5*arguments[4] * m.get_mesh_dt()*((pow(vol.function_operator(arguments), 2) / pow(m.d_x[i], 2)) + ((pow(vol.function_operator(arguments), 2) - rate.function_operator(arguments)) / (2.0*m.d_x[i])));

		}
		return result;
	}


	std::vector<double> up_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments)
	{
		//std::vector<double> a=m.spot_vector();
		std::vector<double> a = m.spot_vect;
		long size = a.size();
		std::vector<double> result(size, 0.0);

		result[0] =0.0;
		result[size-1] = 0.0;
		for (std::size_t i = 1; i < size - 1; ++i)
		{

			//result[i]=0.5*arguments[4]*m.get_mesh_dt()*((-pow(vol.function_operator(arguments),2)/pow(m.get_mesh_dx(),2))+((pow(vol.function_operator(arguments),2)-rate.function_operator(arguments))/(2.0*m.get_mesh_dx())));
			result[i] = 0.5*arguments[4] * m.get_mesh_dt()*((-pow(vol.function_operator(arguments), 2) / pow(m.d_x[i], 2)) + ((pow(vol.function_operator(arguments), 2) - rate.function_operator(arguments)) / (2.0*m.d_x[i])));
		}
		return result;
	}

	// ici je compute le vecteur à la maturité pour avoir le prix à matu et pouvoir faire backward, donc c'est juste appliqué le payoff pour le spot
	std::vector<double> initial_price_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff) {
		std::vector<double> result = m.spot_vect;
		for (std::size_t i = 0; i < result.size(); i++) {
			//arguments[0] = (exp(result[i]));
			arguments[0] = result[i];
			if (payoff.function_operator(arguments) == 0) {
				result[i] = 0;
			}
			else {
				result[i] = payoff.function_operator(arguments);
			}
		}
		return result;
	}

	// une fois que j'ai un vecteur de prix je voudrais le transformer en le mutlipliant par la matrice M(1-theta)
	// comme ça la partie de droite du problème ne sera qu'un vecteur
	std::vector<double> column_up(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff, std::vector<double> up_price) {
		std::vector<double> result = up_price;
		std::vector<double> result2 = m.spot_vect;
		std::vector<double> a = sub_vector(m, rate, vol, arguments); //a(theta)
		std::vector<double> b = diag_vector(m, rate, vol, arguments); //b(theta)
		std::vector<double> c = up_vector(m, rate, vol, arguments); //c(theta)
		arguments[4] = arguments[4] - 1;
		for (std::size_t i = 1; i < result.size() - 1; i++) {
			//arguments[0] = result2[i]; porque ? var de sigma et r dans les coeff maybe
			result[i] = up_price[i] * b[i] + up_price[i - 1] * a[i] + up_price[i + 1] *c[i];
		}
		return result;
	}

	//algo used : https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	// là j'utilise en gros l'algo du lien au dessus pour solver, mais je pense qu'il faudrait le refaire à la main pour être
	// sur que ça marche vraiment, dans le wiki y a la démo et ça montre comment le faire à la main
	std::vector<double> price_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff, std::vector<double> col_up) {
		int size = col_up.size();
		std::vector<double> result(size);
		double W = 0;
		//std::vector<double> arguments_up = arguments;
		//std::vector<double> arguments_down = arguments;
		std::vector<double> B(size);
		std::vector<double> D(size);
		std::vector<double> a = sub_vector(m, rate, vol, arguments); //a(theta)
		std::vector<double> b = diag_vector(m, rate, vol, arguments); //b(theta)
		std::vector<double> c = up_vector(m, rate, vol, arguments); //c(theta)
		result[0] = col_up[0];
		result[size - 1] = col_up[size - 1]*exp(-rate.function_operator(arguments)*m.get_mesh_dt());
		for (std::size_t i = 1; i < size - 1; i++) {
			//arguments_down[0] = arguments[3] + (i - 1)*m.get_mesh_dx();
			//arguments[0] = arguments[3] + i*m.get_mesh_dx();
			//arguments_up[0] = arguments[3] + (i + 1)*m.get_mesh_dx();
			W = a[i] / b[i-1];
			B[i] = b[i] - W*c[i-1];
			D[i] = col_up[i] - W*col_up[i - 1];
		}
		for (std::size_t i = col_up.size() - 2; i > 0; i--) {
			//arguments[0] = arguments[3] + i*m.get_mesh_dx();
			result[i] = (D[i] - c[i]*result[i + 1]) / B[i];
			i = i;
		}
		return result;
	}


	//une fois que j'ai l'algo pour solver au dessus je boucle jusqu'à arriver à t=0 et avoir le prix initial
	std::vector<double> price_today(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff) {
		std::vector<double> ini_price(initial_price_vector(m, rate, vol, arguments, payoff));
		std::vector<double> col_up(column_up(m, rate, vol, arguments, payoff, ini_price));
		double dt = m.get_mesh_dt();
		std::vector<double> result2 = col_up;
		int nb_step = floor(arguments[1] / dt);
		std::vector<double> result1(col_up.size());
		for (int i = 1; i <= nb_step; i++) {
			//arguments[1] = arguments[1] - m.get_mesh_dt();
			result1 = (price_vector(m, rate, vol, arguments, payoff, result2));
			result2 = column_up(m, rate, vol, arguments, payoff, result1);
		}
		return result1;
	}

}