#include "solver.hpp"
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
		std::vector<double> result2(size - 1,log(m_dx));

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
		return m_dt / 365.0;
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

	// AUTRE METHODE - TEST

	// ici je compute le vecteur ‡ la maturitÈ pour avoir le prix ‡ matu et pouvoir faire backward, donc c'est juste appliquÈ le payoff pour le spot

	std::vector<double> initial_price_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff) {
		//std::vector<double> result=m.spot_vector();
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

	std::vector<double> diag_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments)
	{
		//std::vector<double> a=m.spot_vector();
		std::vector<double> a = m.spot_vect;
		long size = a.size();
		std::vector<double> result(size, 0.0);
		result[0] = 1.0;
		result[size - 1] = 1.0;
		for (std::size_t i = 1; i < size - 1; ++i)
		{
			//result[i]=1.0+arguments[4]*m.get_mesh_dt()*((pow(vol.function_operator(arguments),2)/pow(m.get_mesh_dx(),2))+rate.function_operator(arguments));

			result[i] = 1.0 + arguments[4] * m.get_mesh_dt()*((pow(vol.function_operator(arguments), 2) / pow(m.d_x[i-1], 2)) + rate.function_operator(arguments));

		}
		return result;
	}


	std::vector<double> sub_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments)
	{
		//std::vector<double> a=m.spot_vector();
		std::vector<double> a = m.spot_vect;
		long size = a.size();
		std::vector<double> result(size, 0.0);
		result[size - 1] = 0;
		result[0] = 0;
		for (std::size_t i = 1; i < size - 1; ++i)
		{
			//result[i]=-0.5*arguments[4]*m.get_mesh_dt()*((pow(vol.function_operator(arguments),2)/pow(m.get_mesh_dx(),2))+((pow(vol.function_operator(arguments),2)-rate.function_operator(arguments))/(2.0*m.get_mesh_dx())));

			result[i] = -0.5*arguments[4] * m.get_mesh_dt()*((pow(vol.function_operator(arguments), 2) / pow(m.d_x[i-1], 2)) + ((pow(vol.function_operator(arguments), 2) - rate.function_operator(arguments)) / (2.0*m.d_x[i-1])));
			result[i] = (arguments[4] *m.get_mesh_dt()*(-(0.5*pow(vol.function_operator(arguments), 2) / pow(m.d_x[i-1], 2) + 0.25*(pow(vol.function_operator(arguments), 2) - rate.function_operator(arguments) / m.d_x[i-1]))));
		}
		return result;
	}


	std::vector<double> up_vector(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments)
	{
		//std::vector<double> a=m.spot_vector();
		std::vector<double> a = m.spot_vect;
		long size = a.size();
		std::vector<double> result(size, 0.0);

		result[0] = 0.0;
		result[size - 1] = 0.0;
		for (std::size_t i = 1; i < size - 1; ++i)
		{

			//result[i]=0.5*arguments[4]*m.get_mesh_dt()*((-pow(vol.function_operator(arguments),2)/pow(m.get_mesh_dx(),2))+((pow(vol.function_operator(arguments),2)-rate.function_operator(arguments))/(2.0*m.get_mesh_dx())));
			//result[i] = (arguments[4] *m.get_mesh_dt()*((-0.5*pow(vol.function_operator(arguments), 2) / pow(m.d_x[i-1], 2) + 0.25*(vol.function_operator(arguments), 2) - rate.function_operator(arguments) / m.d_x[i - 1])));
			//result[i] = 0.5*arguments[4] * m.get_mesh_dt()*((-pow(vol.function_operator(arguments), 2) / pow(m.d_x[i-1], 2)) + ((pow(vol.function_operator(arguments), 2) - rate.function_operator(arguments)) / (2.0*m.d_x[i-1])));
			result[i] = (arguments[4] * m.get_mesh_dt()*((-0.5*pow(vol.function_operator(arguments), 2) / pow(m.d_x[i - 1], 2) + 0.25*(pow(vol.function_operator(arguments), 2) - rate.function_operator(arguments) / m.d_x[i - 1]))));

		}
		return result;
	}


	// Triadiag algo qui fonctionne !
	std::vector<double> tridiagonal_solver(std::vector<double>  a, std::vector<double>  b, std::vector<double>  c, std::vector<double>  f)
	{

		long n = f.size();
		std::vector<double> x(n);
		x[0] = f[0];
		for (int i = 1; i < n; i++) {

			double m = a[i] / b[i - 1];
			b[i] -= m*c[i - 1];
			f[i] -= m*f[i - 1];
		}
		// solve for last x value
		x[n - 1] = f[n - 1] / b[n - 1];

		// solve for remaining x values by back substitution
		for (int i = n - 2; i >= 0; i--)
			x[i] = (f[i] - c[i] * x[i + 1]) / b[i];

		return x;

	}


	std::vector<double> price_today(mesh m, initial_function rate, initial_function vol, std::vector<double> arguments, initial_function payoff)
	{


		std::vector<double> f_old = initial_price_vector(m, rate, vol, arguments, payoff);


		long N = f_old.size();
		std::vector<double> d(N, 0.0);
		std::vector<double> f_new(N, 0.0);
		double dt = m.get_mesh_dt();

		int nb_step = floor(arguments[1] / dt);


		std::vector<double> arg = arguments;
		arg[4] = arguments[4] - 1; //theta-1

		//Coeffs de la Matrice
		std::vector<double> a_1 = sub_vector(m, rate, vol, arg); //a(theta-1)
		std::vector<double> b_1 = diag_vector(m, rate, vol, arg); //b(theta-1)
		std::vector<double> c_1 = up_vector(m, rate, vol, arg); //c(theta-1)

		std::vector<double> a = sub_vector(m, rate, vol, arguments); //a(theta)
		std::vector<double> b = diag_vector(m, rate, vol, arguments); //b(theta)
		std::vector<double> c = up_vector(m, rate, vol, arguments); //c(theta)

		//Condition aux bords (Test pour un call)
		d[N - 1] = f_old[N - 1]; //Smax
		d[0] = f_old[0];

		//return c_1;
		for (int j = 0; j < nb_step; j++)

		{
			// Creation 2nd membre
			//d[N - 1] = f_old[N - 1] * exp(-rate.function_operator(arguments)*m.get_mesh_dt());

			for (long i = 1; i < N - 1; i++)
			{
				//d[i] = c_1[i]*f_old[i+1]+b_1[i]*f_old[i]+a_1[i-1]*f_old[i-1];

				d[i] = c_1[i] * f_old[i + 1] + b_1[i] * f_old[i] + a_1[i] * f_old[i - 1];

			}

			/* //Condition aux bords (Test pour un call)
			 d[N-1]=f_old[N-1]; //Smax
			 d[0]=f_old[0]; */

			 // Now we solve the tridiagonal system
			f_new = tridiagonal_solver(a, b, c, d);
			f_old = f_new;

		}
		return f_old;
		return f_new;
	}
}
    
    



