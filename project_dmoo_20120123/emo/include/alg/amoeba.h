/*! \file	amoeba.h
	
	\brief	modified from Numerical Recipes in C++ (2nd edition), Chapter 10.4 Downhill Simplex Method in Multidimensions
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Nov.01 2007 create
*/

#ifndef	AZ_AMOEBA_H
#define	AZ_AMOEBA_H

#include <vector>

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{
	typedef std::vector<double> VEC;
	typedef std::vector< VEC > MAT;

	class amoeba
	{
	public:
		//!\brief	main function of simplex method
		//!\param	p input matrix (row-wise), the vertices
		//!\param	y the objectives of the vertices
		//!\param	ftol objective tolerance
		//!\param	funk objective function
		//!\param	nfunk number of function evaluations
		//!\param	NMAX maximum number of function evaluations
		//!\return	void
		void opt(MAT &p, VEC &y, const double ftol, unsigned int &nfunk, const unsigned int NMAX);

		virtual double objective(VEC& x)=0;
	protected:
		//!\brief	get the row sum
		//!\param	p input matrix (row-wise)
		//!\param	psum output sum
		//!\return	void
		void get_psum(MAT &p, VEC &psum);
		
		//!\brief	one try of simplex method
		//!\param	p input matrix (row-wise), the vertices
		//!\param	y the objectives of the vertices
		//!\param	psum sum of the matrix
		//!\param	funk objective function
		//!\param	ihi index of the worst point
		//!\param	fac fraction of reflect
		double amotry(MAT &p, VEC &y, VEC &psum, const unsigned int ihi, const double fac);
	};
} //namespace alg

} //namespace az

#endif //AZ_AMOEBA_H
