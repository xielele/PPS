/*! \file	Alg.h
	
	\brief	virutal class for Framwork for MOEA
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
		
	\date	May.22 2007 create
	\date	Jun.30 2008 redesign
*/


#ifndef	AZ_ALGORITHM_BASE_H
#define	AZ_ALGORITHM_BASE_H

#include <ctime>
#include "LogFile.h"
#include "Parameter.h"
#include "PopulationMO.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
	//!\brief	base class of multi-objective evolutionary algorithm
	class EABasic
	{
	protected:
		unsigned int
			mPopSize,			//!< population size
			mStep	,			//!< current step
			mStepMax,			//!< maximum steps
			mEvas	;			//!< the calculation number of the objectives
		double 
			mTM,				//!< time used to generate offspring
			mTS;				//!< time used to selection
		CPopulationMO	mPop;	//!< population
		CParameter	mPar;		//!< the parameter object
	public:
		//!\brief	constractor
		//!\param	stepmax maximum steps
		//!\param	par parameter object
		//!\param	pop initial population
		//!\return	void
		EABasic(
			unsigned int	stepmax	,
			CParameter&		par		,
			CPopulationMO&	pop		)
			: mPop(pop)
		{
			mStepMax 	= stepmax;
			mPar		= par;
			mStep		= 0;
			mEvas		= 0;
			mPopSize 	= mPop.Size();
			mTM			= mTS	= 0.0;
		}

		//!\brief	destructor
		//!\return	void
		virtual ~EABasic() {}

		//!\brief	get the pointer to the parameter object
		//!\return	pointer to the parameter object
		inline CParameter& P() {return mPar;}

		//!\brief	get the objective evaluation times
		//!\return	objective evaluation times
		inline unsigned int& EvaTimes() {return mEvas;}

		//!\brief	check to see whether the terminate condition is satified
		//!\return	true if terminate
		inline bool IsTerminate() {return mStep >= mStepMax;}
		
		//!\brief	get the maximal step
		//!\return	maximal step
		inline unsigned int MaxStep() {return mStepMax;}

		//!\brief	get the current step
		//!\return	current step
		inline unsigned int CurStep() {return mStep;}

		//!\brief	get the run time
		//!\return	void
		inline void Time(double& tm, double& ts) {tm = mTM/CLOCKS_PER_SEC; ts = mTS/CLOCKS_PER_SEC;}

		//!\brief	get the population
		//!\return	reference to population
		CPopulationMO& Population() {return mPop;}

		//!\breif	save data
		//!\brief	path path name to save
		void Save(std::string& path)
		{
			std::string fname = path + std::string(".pop");
			Population().Write(fname.c_str());
		}

		//!\brief	one step 
		//!\return	current step
		virtual unsigned int Step() = 0;
	protected:
		//!\brief	make all new solutions in feasible region
		//!\param	pop offspring population
		//!\return	void
		void Check(CPopulationMO& pop)
		{
			for(unsigned int i=0; i<pop.Size(); i++) pop[i].Check();
		}

		//!\brief	make all new solutions are different from old solutions
		//!\param	popnew offspring population
		//!\param	pop old population
		//!\return	void
		void Check(CPopulationMO& popnew, CPopulationMO& pop)
		{
			CPopulationMO tmp(P());
			for(unsigned int i=0; i<popnew.Size(); i++) if(!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
			popnew.Clear();
			popnew.Combine(tmp);
		}
	}; //class EABasic

	//!\brief	EA framework with Single Population for Multi-objective Optimization
	//!\param	COG	offspring generator
	//!\param	CES environmental selection
	template<class COG, class CES>
		class EAEngine : public EABasic, public COG, public CES
	{
	public:
		//!\brief	constractor
		//!\param	stepmax maximum steps
		//!\param	par parameter object
		//!\param	pop initial population
		//!\return	void
		EAEngine(
			unsigned int	stepmax	,
			CParameter&		par		,
			CPopulationMO&	pop		)
			: EABasic(stepmax, par, pop)
		{
		}

		//!\brief	destructor
		//!\return	void
		virtual ~EAEngine() {}

		//!\brief	one step 
		//!\return	current step
		virtual unsigned int Step()
		{
			clock_t t1,t2,t3;
			CPopulationMO popnew(P());

			t1 = clock();

			//sample new solutions
			COG::Generate(Population().Size(), popnew, Population());
			
			//check the range of new solutions
			Check(popnew);

			//remove these repeated solutions
			Check(popnew, Population());

			//evaluate new solutions	
			popnew.Evaluate();
			mEvas += popnew.Size();	

			t2 = clock();

			//environmental select 
			Population().Combine(popnew);
			CES::Select(Population(), mPopSize);
			
			t3 = clock();

			mTM += double(t2-t1);
			mTS += double(t3-t2);

			return ++mStep;
		}
	}; //class EAEngine

} //namespace mea

} //namespace az

#endif //AZ_ALGORITHM_BASE_H
