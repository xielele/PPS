/*! \file	AlgD.h
	
	\brief	Framwork for Dynamic MOEA with GTM
	\brief  The basic idea is to set: x = x1 + x2
	\brief  the point x1(t) (centre point, coordinate origin, etc.) is moving with time
	\brief	the distribution P(x2) is more stable, i.e, P(x2(t-1)) is similar to P(x2(t))
	\brief	thus, the moving of x1 is modeled by time series (AR model in the code)
	\brief  the distribution of x2 is modeled by GTM (other method like LPCA can also be applied)
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
		
	\date	Aug.15 2007 create
*/


#ifndef	AZ_ALGORITHM_DMOO_H
#define	AZ_ALGORITHM_DMOO_H

#include <ctime>
#include <vector>
#include "alg/Matrix.h"
#include "alg/GTM.h"
#include "emo/Parameter.h"
#include "emo/PopulationMO.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{

//!\brief namespace of dynamic evolutionary algoirhtm 
namespace dea
{

//!\brief	Dynamic MOEA with GTM
class DMOO
{
protected:
	// parameters for algorithm
	// prediction method
	enum PMETHOD {	INIRIS,		//!< randomly initialize population 
					INIFPS,		//!< predict several points
					INIPPS,		//!< predict the centre by AR model
					INIPPSBP,	//!< predict the centre by BP model
					INIPPSNR};  //!< nonlinear regression model
	PMETHOD mStrategy;			//!< predict method used in the algorithm
	unsigned int
		mPopSize,	//!< population size
		mMaxStep,	//!< maximum generation
		mStep	,	//!< current step
		mEvas	;	//!< the calculation number of the objectives
	CPopulationMO	mPop;	//!< current population
	CParameter*	pPar;		//!< pointer to the parameter object
	bool	mbToChange,		//!< indicator to change in the next step
			mbChange;		//!< whether change in current step
	std::string		mOptimizer; //!< optimizer inside the time windows

	// parameters for problem
	unsigned int
		mTaoT	;	//!< the length of each state
	double 
		mT0		,	//!< the initial time
		mDelT	;	//!< the step size of time
	
	// parameters for modeling (time series)
	double			mAlpha;			//!< parameter to control the variance
	unsigned int	mMaxOrder;		//!< maximum order of time series
	CPopulationMO	mBest0,			//!< best history population: mBest = mPop - mC
					mBest1;			//!< best history population: mBest = mPop - mC
	std::list< std::vector<double> > hC;	//!< list of history centre points
	std::vector<double>	mC,					//!< current centre point
						pC,					//!< the predicted centre
						mStdC;				//!< standard deviation of predicted center
	bool	mIsPre;							//!< whether there is a prediction in the generation
public:
	//!\brief	constractor
	//!\param	strategy	prediction strategy
	//!\param	popsize		population size
	//!\param	stepmax		maximum steps
	//!\param	taot		the length of each state in case of generations
	//!\param	nt			servety of changes
	//!\param	torder		order of time series
	//!\parame	t0			initial state of time T
	//!\param	alpha		variance control parameter
	//!\param	par			parameter object
	//!\return	void
	DMOO(
		unsigned int	strategy,
		std::string&	optimizer,
		unsigned int	popsize	,
		unsigned int	stepmax	,
		unsigned int	taot	,
		unsigned int	nt		,
		unsigned int	torder	,
		double			t0		,
		double			alpha   ,
		CParameter&		par		);

	~DMOO();

	//!\brief	get the pointer to the parameter object
	//!\return	pointer to the parameter object
	inline CParameter& P() {return *pPar;}

	//!\brief	get the population
	//1\return	reference to population
	inline CPopulationMO& Population() {return mPop;}

	//!\brief	get the objective evaluation times
	//!\return	objective evaluation times
	inline unsigned int& EvaTimes() {return mEvas;}

	//!\brief	check to see whether the terminate condition is satified
	//!\return	true if terminate
	inline bool IsTerminate() {return mStep >= mMaxStep;}
	
	//!\brief	get the current step
	//!\return	current step
	inline unsigned int CurStep() {return mStep;}

	//!\brief	whether to change in the next step
	//!\return	bool
	inline bool IsToChange() {return mbToChange;}

	//!\brief	save population
	//!\return	void
	inline void Write(std::string name) {Population().Write(name);}

	//!\brief	reset parameters
	//!\return	void
	void Reset();

	//!\brief	one step 
	//!\return	current step
	unsigned int Step();
protected:
	//!\brief	make all new solutions in feasible region
	//!\param	pop offspring population
	//!\return	void
	void Check(CPopulationMO& pop);

	//!\brief	make all new solutions are different from old solutions
	//!\param	popnew offspring population
	//!\param	pop old population
	//!\return	void
	void Check(CPopulationMO& popnew, CPopulationMO& pop);

	//!\prief	to predict the new weight matrix and new centre point
	//!\return	void
	void Predict();

	bool EnvironmentChange();

	void InitPPS();

	void InitFPS();

	void InitRIS(bool all);

	//!\brief	Principal Curve Generator
	//!\param	popnew offspring population
	//!\param	size offspring population size
	//!\return	offspring population
	CPopulationMO& Generate(CPopulationMO& popnew, unsigned int size);

}; //class DMOO

} //namespace dea

} //namespace mea

} //namespace az

#endif //AZ_ALGORITHM_DMOO_H
