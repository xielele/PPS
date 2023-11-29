/*! \file	Fitting.h
	
	\brief	polynomial data fitting
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Jan.10 2005 create
	\date	Sep.26 2005	redesign and rewrite
*/

#ifndef	AZ_FITTING_H
#define	AZ_FITTING_H

#include <vector>

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{
	//!\brief C.T=X
	//! C=(c0,c1,c2,...,cn), X=(x1, x2, ..., xm)
	//!\param	pT independent variable
	//!\param	pX variable
	//!\param	size training data size
	//!\param	order maximum order
	//!\param	pC coefficent(final reuslt)
	//!\return	void
	void poly_fit(double* pT, double* pX, unsigned int size, unsigned int order, double* pC);

	//!\brief C.T=X
	//! C=(c0,c1,c2,...,cn), X=(x1, x2, ..., xm)
	//!\param	T independent variable
	//!\param	X variable
	//!\param	order maximum order
	//!\param	C coefficent(final reuslt)
	//!\return	void
	void poly_fit(std::vector<double>& T, std::vector<double>& X, unsigned int order, std::vector<double>& C);

	//!\brief C.T=X
	//! C=(c0,c1,c2,...,cn), X=(x1, x2, ..., xm)
	//!\param	pT1 independent variable
	//!\param	pT2 independent variable
	//!\param	pX variable
	//!\param	size training data size
	//!\param	order maximum order
	//!\param	pC coefficent(final reuslt)
	//!\return	void
	void poly_fit(double* pT1, double* pT2, double* pX, unsigned int size, unsigned int order, double* pC);

	//!\brief C.T=X
	//! C=(c0,c1,c2,...,cn), X=(x1, x2, ..., xm)
	//!\param	T1 independent variable
	//!\param	T2 independent variable
	//!\param	X variable
	//!\param	order maximum order
	//!\param	C coefficent(final reuslt)
	//!\return	void
	void poly_fit(std::vector< double >& T1, std::vector< double >& T2, std::vector< double >& X, unsigned int order, std::vector< double >& C);

} //namespace alg

} //namespace az

#endif //AM_FITTING_H
