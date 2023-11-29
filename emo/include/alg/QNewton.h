/*! \file	QNewton.h
	
	\brief	Quaci-Newton optimization method, minimize a function
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Jan.10 2005 create the file
	\date	Sep.25 2005 rewrite the file
*/

#ifndef	AZ_QNEWTON_H
#define	AZ_QNEWTON_H

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

//!\brief objective function type
typedef double	(*FTYPE)(double*);

//!\brief gradient function type
typedef void	(*GTYPE)(double*, double*); 

//!\brief Quaci-Newton method
class QNewton
{
protected:
	int		mDim,		//!< dimension of variable
			mIter,		//!< iterate number
			mEval;		//!< number of objective evaluation

	double	mTolG,		//!< tolerance of gradient		
			mTolX,		//!< tolerance of variable
			mStepMax,	//!< maximum step length
			mF,			//!< current objective value
			mFBest,		//!< best objective value
			*pBoundUpp,	//!< the upper bound of variables
			*pBoundLow,	//!< the lower bound of variables
			*pG,		//!< gradient
			*pDG,		//!< temporal variale
			*pHDG,		//!< temporal variale
			*pP,		//!< temporal variale
			*pScild,	//!< temporal variale
			*pX,		//!< current variables
			*pXBest,	//!< the best varialbes
			**pHessin;	//!< the hessin matrix

	FTYPE	pFun;		//!< pointer to objective function
	GTYPE	pGun;		//!< pointer to gradient function
public:
	//!\brief constractor
	//!\param dim	variable dimension
	//!\param fun	pointer to objective function
	//!\param gun	pointer to gradient function
	//!\param x0	initial variables
	//!\param low	lower bound of variables
	//!\param up	upper bound of variables
	//!\param tolx		tolerance of variable
	//!\param tolg		tolerance of gradient
	//!\return no
	QNewton(					
		int		dim,
		FTYPE	fun,
		GTYPE	gun,
		double	*x0,
		double	*low,
		double	*up,
		double	tolx	= 1.0e-16,
		double	tolg	= 1.0e-10);

	//!\brief destractor
	//!\return no
	~QNewton();

	//!\brief	get variable dimension
	//!\return	variable dimension		
	inline int	Dimension()		{return mDim;}

	//!\brief	get the objective evaluation number
	//!\return	objective evaluation number
	inline int	Evaluations()	{return mEval;}
	
	//!\brief	get the iteration number
	//!\return	iteration number
	inline int	Iteration()		{return mIter;}
	
	//!\brief	get the best objective value
	//!\return	best objective value
	inline double  F()	{return mFBest;}

	//!\brief	get the best objective value
	//!\return	best objective value
	inline double* X()  {return pXBest;}

	//!\brief	one step
	//!\return	whether the run is successful
	bool Step();

protected:
	//!\brief	evaluate a vector
	//!\param	px vector
	//!\return	objective value
	double FEva(double* px);

	//!\brief	calculate the gradient of a vector
	//!\param	px vector
	//!\param	pg gradient
	//!\return	void
	void GEva(double* px, double* pg);

	//!\brief	linear search
	//!\param	x start point
	//!\param   f objective
	//!\param	g gradient
	//!\param	p temporal variable
	//!\return	whether the search is successful
	bool lnsrch(
		double*& x, 
		double & f, 
		double*& g, 
		double*& p);

};//class QNewton

} //namespace alg

} //namespace az

#endif	//AZ_QNEWTON_H
