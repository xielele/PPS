/*! \file	AR.h
	
	\brief	multivariant autoregression model
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Aug.07 2007 create
*/

#ifndef	AZ_AR_H
#define	AZ_AR_H

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{
	//!\brief univariate autoregression model
	//! Xn+1 = A1*Xn + A2*Xn-1 + ... + Ap*Xn-p+1 + epsilon
	//! input (output) matrix and vector are column-wise
	//!\param	pX training data
	//!\param	dim variable dimension
	//!\param	size training data size
	//!\param	order order of AR model
	//!\param	pA coefficent(final reuslt)
	//!\param	pV variance
	//!\return	bool
	bool aruv(double** pX, unsigned int dim, unsigned int size, unsigned int order, double** pA, double* pV);

	//!\brief univariate autoregression model
	//! Xn+1 = C + A1*Xn + A2*Xn-1 + ... + Ap*Xn-p+1 + epsilon
	//! input (output) matrix and vector are column-wise
	//!\param	pX training data
	//!\param	dim variable dimension
	//!\param	size training data size
	//!\param	order order of AR model
	//!\param	pA coefficent(final reuslt)
	//!\param	pC mean 
	//!\param	pV variance
	//!\return	bool
	bool aruv(double** pX, unsigned int dim, unsigned int size, unsigned int order, double** pA, double* pC, double* pV);

	//!\brief multivariate autoregression model
	//! Xn+1 = C + A1*Xn + A2*Xn-1 + ... + Ap*Xn-p+1 + epsilon
	//! input (output) matrix and vector are column-wise
	//!\param	pX training data
	//!\param	dim variable dimension
	//!\param	size training data size
	//!\param	order order of AR model
	//!\param	pA coefficent(final reuslt)
	//!\param	pC mean 
	//!\param	pV variance
	//!\return	bool
	bool armv(double** pX, unsigned int dim, unsigned int size, unsigned int order, double** pA, double* pC, double* pV);
} //namespace alg

} //namespace az

#endif //AZ_AR_H
