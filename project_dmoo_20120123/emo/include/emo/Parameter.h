/*! \file	Parameter.h
	
	\brief	EA parameters

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Sep.27 2005 rewrite & reorganize structure
*/
#ifndef AZ_PARAMETER_H
#define AZ_PARAMETER_H

#include <list>
#include <map>
#include <vector>
#include <string>

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
#define	MAXDOUBLE 1.7E+308	//!< const value

//!\brief	evaluator definition
//!\param	F objective vector
//!\param	E equality constraint vector
//!\param	I inequality constraint vector
//!\param	X variable vector
//!\return	void
typedef void (*EVALUATOR)(	std::vector< double >& F, 
							std::vector< double >& E, 
							std::vector< double >& I, 
							std::vector< double >& X);

//!\brief parameter class contains all parameters
class CParameter
{
protected:
	bool
		mXCoding;	//!< coding in decision space
	unsigned int 
		mFSize,		//!< objective number
		mESize,		//!< equality constraint number
		mISize;		//!< inequality constraint number
	double 
		mPm,		//!< probability of mutation
		mPc,		//!< probability of crossover
		mTolX,		//!< toleance for X
		mTolF,		//!< toleance for F
		mTolC,		//!< toleance for constraints
		mETA_SBX,	//!< eta for SBX(simulated binary crossover)
		mETA_PM,	//!< eta for PM(polynomial mutation)
		mETA_PCX;	//!< eta for PCX(parent-centric recombination)
	EVALUATOR
		pEvaluator;	//!< pointer to the evaluator
	std::vector<double>
		mVBoundupp,		//!< upper bound of variables
		mVBoundlow,		//!< lower boudn of variables
		mVXrealupp,		//!< upper bound of variables
		mVXreallow;		//!< lower boudn of variables
	std::string 
		mProblem;		//!< problem name
public:
	//!\brief	constructor
	//!\return	void
	CParameter() {mXCoding = false;mFSize=0;mESize=0;mISize=0;}

	//!\brief	get x-coding state
	inline bool& XCoding() {return mXCoding;}

	//!\brief	get probability of mutation
	//!\return	probability of mutation
	inline double& Pm() {return mPm;}

	//!\brief	get probability of crossover
	//!\return	probability of crossover
	inline double& Pc() {return mPc;}

	//!\brief	get toleance for X
	//!\return	toleance for X
	inline double& TolX() {return mTolX;}

	//!\brief	get toleance for F
	//!\return	toleance for F
	inline double& TolF() {return mTolF;}

	//!\brief	get toleance for constraints
	//!\return	toleance for constraints
	inline double& TolC() {return mTolC;}

	//!\brief	get eta for SBX
	//!\return	eta for SBX
	inline double& ETA_SBX() {return mETA_SBX;} 

	//!\brief	get eta for PM
	//!\return	eta for PM
	inline double& ETA_PM() {return mETA_PM;} 

	//!\brief	get eta for PCX
	//!\return	eta for PCX
	inline double& ETA_PCX() {return mETA_PCX;} 
	
	//!\brief	set objective number
	//!\param	s new objective number
	//!\return	objective number
	inline unsigned int FSize(unsigned int s) {mFSize=s;return s;}

	//!\brief	get objective number
	//!\return	objective number
	inline unsigned int FSize() {return mFSize;}
	
	//!\brief	set variable number
	//!\param	s new variable number
	//!\return	variable number
	inline unsigned int XSize(unsigned int s) {mVBoundupp.resize(s); mVBoundlow.resize(s); mVXrealupp.resize(s); mVXreallow.resize(s);return s;}
	
	//!\brief	get variable number
	//!\return	variable number		
	inline unsigned int XSize() {return (unsigned int)mVBoundupp.size();	}

	//!\brief	set equality constraint number
	//!\param	s new equality constraint number
	//!\return	equality constraint number		
	inline unsigned int ESize(unsigned int s) {mESize=s; return s;}

	//!\brief	get equality constraint number
	//!\return	equality constraint number		
	inline unsigned int ESize() {return mESize;}

	//!\brief	set inequality constraint number
	//!\param	s new inequality constraint number
	//!\return	inequality constraint number	
	inline unsigned int ISize(unsigned int s) {mISize=s; return s;}

	//!\brief	get inequality constraint number
	//!\return	inequality constraint number	
	inline unsigned int ISize() {return mISize;}

	//!\brief	get lower bound of i-th variable
	//!\param	i variable index
	//!\return	lower bound of i-th variable
	inline double& XLow(unsigned int i) {return mVBoundlow[i];	}

	//!\brief	get upper bound of i-th variable
	//!\param	i variable index
	//!\return	upper bound of i-th variable
	inline double& XUpp(unsigned int i) {return mVBoundupp[i];	}

	//!\brief	get lower bound vector
	//!\return	lower bound
	inline std::vector<double>& XLow() {return mVBoundlow;	}

	//!\brief	get upper bound
	//!\return	upper bound of i-th variable
	inline std::vector<double>& XUpp() {return mVBoundupp;	}

	//!\brief	get (real) upper bound of i-th variable
	//!\param	i variable index
	//!\return	upper bound of i-th variable
	inline double& XRealUpp(unsigned int i) {return mVXrealupp[i];	}

	//!\brief	get (real) lower bound of i-th variable
	//!\param	i variable index
	//!\return	lower bound of i-th variable
	inline double& XRealLow(unsigned int i) {return mVXreallow[i];	}

	//!\brief	set evaluator
	//!\param	peva pointer to new evaluator
	//!\return	pointer to new evaluator
	inline EVALUATOR Evaluator(EVALUATOR peva) {pEvaluator=peva;return peva;}

	//!\brief	get evaluator
	//!\return	pointer to new evaluator
	inline EVALUATOR Evaluator() {return pEvaluator;}

	//!\brief	set problem name
	//!\param	str new problem name
	//!\return	problem name
	inline std::string& Problem(std::string str) {mProblem=str;return mProblem;}

	//!\brief	get problem name
	//!\return	problem name
	inline std::string& Problem() {return mProblem;}

}; //class CParameter

} //namespace mea

} //namespace az

#endif //AZ_PARAMETER_H
