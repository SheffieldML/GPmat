// This code copyright Manuel Sanchez

#include <iostream>
#include "Optimize.h"
#include "auto_array.hpp"



using namespace Optimize;


float SCG::iterate(
	float* w, 
	const ErrFunc* E, 
	unsigned int nits, float errTH)
{
unsigned int dim = E->dim();

bool success;
float alpha, beta, delta, Delta, lambda, hLambda, nu, sigma, rr;
std::auto_array<float> dE(new float[dim]), dEp(new float[dim]);
std::auto_array<float> p(new float[dim]), r(new float[dim]), rp(new float[dim]), s(new float[dim]);
std::auto_array<float> wp(new float[dim]);

	// 1
	lambda = m_reg; hLambda = 0.0f;
	E->grad(r.get(), w); float* r_ = r.get(); Vector::Ops_N::cpt(r_, r_, dim);
	memcpy(p.get(), r.get(), dim*sizeof(float));
	success = true;

unsigned int i, k;

	for(k = 1; k <= nits; ++k)
	{
		/*
		std::cout << "lambda_" << k << ": " << lambda << " ";

		std::cout << "w_" << k << ": ";

		for(i = 0; i < dim; ++i)
			std::cout << w[i] << " ";

		std::cout << "r_" << k << ": ";

		for(i = 0; i < dim; ++i)
			std::cout << r.get()[i] << " ";

		std::cout << "s_" << k << ": ";

		for(i = 0; i < dim; ++i)
			std::cout << s.get()[i] << " ";

		std::cout << "p_" << k << ": ";

		for(i = 0; i < dim; ++i)
			std::cout << p.get()[i] << " ";

		std::cout << std::endl;
		*/

	float pp;

		Vector::Ops_N::dot(pp, p.get(), p.get(), dim);

		// 2
		if(success)
		{
			sigma = m_step/sqrt(pp); 
			
			for(i = 0; i < dim; ++i)
				wp.get()[i] = w[i] + sigma*p.get()[i];

			E->grad(dE.get(), w);
			E->grad(dEp.get(), wp.get());

			for(i = 0; i < dim; ++i)
				s.get()[i] = (dEp.get()[i] - dE.get()[i])/sigma;

			Vector::Ops_N::dot(delta, p.get(), s.get(), dim);
		}

		// 3
		for(i = 0; i < dim; ++i)
			s.get()[i] += (lambda - hLambda)*p.get()[i];

		delta += (lambda - hLambda)*pp;

		// 4
		if(delta <= 0.0f)
		{
			for(i = 0; i < dim; ++i)
				s.get()[i] += (lambda - 2.0f*delta/pp)*p.get()[i];

			hLambda = 2.0f*(lambda - delta/pp);
			delta = -delta + lambda*pp; lambda = hLambda;
		}

		// 5
		Vector::Ops_N::dot(nu, p.get(), r.get(), dim);
		alpha = nu/delta;

		for(i = 0; i < dim; ++i)
			wp.get()[i] = w[i] + alpha*p.get()[i];

		// 6
		Delta = 2.0f*delta*(E->eval(w) - E->eval(wp.get()))/(nu*nu);

		// 7
		if(Delta >= 0.0f)
		{
			memcpy(w, wp.get(), dim*sizeof(float)); // w_k = w_kp = w_k + alpha_k*p_k;
			E->grad(rp.get(), w); float* rp_ = rp.get(); Vector::Ops_N::cpt(rp_, rp_, dim);
			hLambda = 0; success = true;

			// 7.a
			if(k % dim == 0)
				memcpy(p.get(), rp.get(), dim*sizeof(float));
			else
			{
			float rprp, rrp;

				Vector::Ops_N::dot(rprp, rp.get(), rp.get(), dim);
				Vector::Ops_N::dot(rrp, r.get(), rp.get(), dim);

				beta = (rprp - rrp)/nu;

				for(i = 0; i < dim; ++i)
					p.get()[i] = rp.get()[i] + beta*p.get()[i];
			}

			memcpy(r.get(), rp.get(), dim*sizeof(float));

			// 7.b
			if(Delta >= 0.75f) lambda *= 0.5f;
		}
		else
			{ hLambda = lambda; success = false; }

		// 8
		if(Delta < 0.25f) lambda *= 4.0f;

		// 9
		Vector::Ops_N::dot(rr, r.get(), r.get(), dim);

		if(rr < errTH*errTH) break;
	}

	return rr;
}
