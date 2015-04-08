#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__


#include "SF_LinAlg.h"



using namespace SF;


namespace Optimize
{


struct ErrFunc
{

	virtual unsigned int dim(void) const = 0;

	virtual float eval(const float* pnt) const = 0;
	virtual void grad(float* grad, const float* pnt) const = 0;

};


struct Method
{

	virtual ~Method() { } 

	virtual float iterate(
		float* pnt, 
		const ErrFunc* err, 
		unsigned int nits, float errTH = 0.0f) = 0;

};


class SCG: public Method
{

	SCG(const SCG&);
	SCG& operator=(const SCG&);

	float m_step, m_reg;

public:

	SCG(float step = 1e-4f, float reg = 1.0f): m_step(step), m_reg(reg) { }

	float iterate(
		float* pnt, 
		const ErrFunc* E, 
		unsigned int nits, float errTH = 0.0f);

};


template<size_t DIM>
class Quad: public Optimize::ErrFunc
{

	float m_A[DIM][DIM];
	float m_b[DIM];
	float m_c;

public:

	Quad(const float* A, const float* b, float c) 
		{ memcpy(m_A, A, DIM*DIM*sizeof(float)); memcpy(m_b, b, DIM*sizeof(float)); m_c = c; }

	~Quad() { };

	inline unsigned int dim(void) const 
		{ return DIM; }

	float eval(const float* pnt) const 
	{
	float res, pAp, bp;
	static float Ap[DIM];

		Matrix::Ops_N::cvp(Ap, m_A, pnt, DIM, DIM);
		Vector::Ops_N::dot(pAp, pnt, Ap, DIM);
		Vector::Ops_N::dot(bp, m_b, pnt, DIM);

		return 0.5f*pAp + bp + m_c;
	}


	void grad(float* grad, const float* pnt) const
	{
	static float Ap[DIM], pA[DIM];

		Matrix::Ops_N::cvp(Ap, m_A, pnt, DIM, DIM);
		Matrix::Ops_N::rvp(pA, pnt, m_A, DIM, DIM);
		Vector::Ops_N::add(grad, Ap, pA, DIM);
		Vector::Ops_N::scl(grad, grad, 0.5f, DIM);
		Vector::Ops_N::add(grad, grad, m_b, DIM);
	}


};


}



#endif