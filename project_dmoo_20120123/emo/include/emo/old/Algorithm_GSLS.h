/*! \file	Algorithm_GSLS.h
	
	\brief	Framwork for Global Search (GS), Local Search (LS) hybrid MOEA
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
		
	\date	Mar.14 2006 create
*/


#ifndef	AZ_ALGORITHMGSLS_H
#define	AZ_ALGORITHMGSLS_H

#include "Parameter.h"
#include "PopulationMO.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief	EA framework
//!\param	COG	offspring generator
//!\param	CLS local search method
//!\param	CES environmental selection
template<class COG, class CLS, class CES>
	class CGSLSEA : public COG, public CLS, public CES
{
protected:
	unsigned int
		mPopSize,	//!< population size
		mStep	,	//!< current step
		mStepMax,	//!< maximum steps
		mEvas	,	//!< the calculation number of the objectives
		mLSize	,	//!< the local search point size
		mSSize	;	//!< the sample size for each local search point
	CPopulationMO mPop;	//!< population
	CParameter*	pPar;	//!< pointer to the parameter object
public:
	//!\brief	constractor
	//!\param	stepmax maximum steps
	//!\param	par parameter object
	//!\param	pop initial population
	//!\return	void
	CGSLSEA	(
		unsigned int	stepmax	,
		CParameter&		par		,
		CPopulationMO&	pop		)
		: mPop(pop) 
	{
		mStepMax 	= stepmax;
		pPar		= &par;
		mStep		= 0;
		mEvas		= 0;
		mPopSize 	= mPop.Size();
		mLSize		= 10;
		mSSize		= 5;
	}

	//!\brief	destructor
	//!\return	void
	~CGSLSEA() {}

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
	inline bool IsTerminate() {return mStep >= mStepMax;}
	
	//!\brief	get the current step
	//!\return	current step
	inline unsigned int CurStep() {return mStep;}

	//!\brief	one step 
	//!\return	current step
	unsigned int Step()
	{
		CPopulationMO popnew(P()), poplocal(P());

		//-------------GLOBAL SEARCH--------------------------
		//sample new solutions
		COG::Generate(popnew, Population());

		//check the range of new solutions
		Check(popnew);

		//remove these repeated solutions
		Check(popnew, Population());

		//evaluate new solutions	
		popnew.Evaluate();
		mEvas += popnew.Size();	

		Population().Combine(popnew);

		//Population().Write( "pop.txt" );
		
		//-------------LOCAL SEARCH---------------------------
		if(mStep % 2 == 0)
		{
			CLS::Generate(poplocal, Population());

			mEvas += poplocal.Size();

			poplocal.Write("Local.txt");

			Population().Combine(poplocal);
		}
		
		//-------------SELECTION------------------------------
		//environmental select 
		CES::Select(Population(), mPopSize);

		return ++mStep;
	}
protected:
	//!\brief	make all new solutions in feasible region
	//!\param	pop offspring population
	//!\return	void
	void Check(CPopulationMO& pop)
	{
		for(unsigned int i=0; i<pop.Size(); i++)
			pop[i].Check();
	}

	//!\brief	make all new solutions are different from old solutions
	//!\param	popnew offspring population
	//!\param	pop old population
	//!\return	void
	void Check(CPopulationMO& popnew, CPopulationMO& pop)
	{
		CPopulationMO tmp(P());
		for(unsigned int i=0; i<popnew.Size(); i++)
			if(!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
		popnew.Clear();
		popnew.Combine(tmp);
	}
}; //class CGSLSEA

} //namespace mea

} //namespace az

#endif //AZ_ALGORITHMGSLS_H
