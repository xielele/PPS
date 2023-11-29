/*! \file	AlgDualPop.h
	
	\brief	Framwork for MOEA with dual populations, one for focuses on F-space and one focuses on X-space
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
		
	\date	Feb.21 2007 create
*/


#ifndef	AZ_ALGORITHM_DUALPOP_H
#define	AZ_ALGORITHM_DUALPOP_H

#include <algorithm>
#include <ctime>
#include <vector>
#include "LogFile.h"
#include "Parameter.h"
#include "PopulationMO.h"
#include "Sel.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
	//!\brief	EA framework
	//!\param	COG	offspring generator
	//!\param	CES environmental selection
	template<class COG, class CES>
		class CEADUALPOP : public COG, public CES
	{
	protected:
		unsigned int
			mPopSize,	//!< population size
			mStep	,	//!< current step
			mStepMax,	//!< maximum steps
			mEvas	;	//!< the calculation number of the objectives
		double 
			mTM,		//!< time used to model
			mTS;		//!< time used to selection
		CPopulationMO	mPop,	//!< population, only used to generate new trial solutions
						mPopF,	//!< archive which focuses on F-space
						mPopX;	//!< archive which focuses on X-space
		CParameter*	pPar;	//!< pointer to the parameter object
		std::vector<unsigned int> mIndex;	//!< selection index
	public:
		//!\brief	constractor
		//!\param	stepmax maximum steps
		//!\param	par parameter object
		//!\param	pop initial population
		//!\return	void
		CEADUALPOP	(
			unsigned int	stepmax	,
			CParameter&		par		,
			CPopulationMO&	pop		)
			: mPop(pop), mPopF(pop), mPopX(pop)
		{
			mStepMax 	= stepmax;
			pPar		= &par;
			mStep		= 0;
			mEvas		= 0;
			mPopSize 	= mPop.Size();
			mTM			= mTS	= 0.0;
			mIndex.resize(mPopSize);
			for(unsigned int i=0; i<mPopSize; i++) mIndex[i] = i;
		}

		//!\brief	destructor
		//!\return	void
		~CEADUALPOP() {}

		//!\brief	get the pointer to the parameter object
		//!\return	pointer to the parameter object
		inline CParameter& P() {return *pPar;}

		//!\brief	get the population
		//1\return	reference to population
		inline CPopulationMO& Population() {return mPopF;}

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

		//!\breif	save data
		//!\brief	path data file path
		void Save(std::string& path)
		{
		   std::string fname = path + std::string(".popf");
		   mPopF.Write(fname.c_str());
		   fname = path + std::string(".popx");
		   mPopX.Write(fname.c_str());
		}

		//!\brief	one step 
		//!\return	current step
		unsigned int Step()
		{
			//unsigned int i;
			clock_t t1,t2,t3;
			CPopulationMO popnew(P());

			// select points from PopX and PopF to generate offspring

			//std::random_shuffle(mIndex.begin(), mIndex.end());
			//for(i=0; i<(unsigned int)(mPopSize*1.0/5.0); i++)			mPop[i] = mPopF[mIndex[i]];
			//for(i=(unsigned int)(mPopSize*1.0/5.0); i<mPopSize; i++)	mPop[i] = mPopX[mIndex[i]];

			mPop.Clear();
			mPop.Copy(mPopX);

			//mPop = mPopF;
			//mPop.Copy(mPopX);

			//mPop = mPopF;
			//==============================================
			// note: there might be some duplicate points 
			//==============================================

			t1 = clock();
			//sample new solutions
			COG::Generate(popnew, mPop);
			
			//check the range of new solutions
			Check(popnew);

			//remove these repeated solutions
			Check(popnew, mPop);

			//evaluate new solutions	
			popnew.Evaluate();
			mEvas += popnew.Size();	

			//popnew.Write("meda.new");

			t2 = clock();

			//environmental select 
			mPopF.Copy(popnew);
			az::mea::sel::SMaxMin selF;
			selF.Select(mPopF,mPopSize);

			mPopX.Copy(popnew);
			az::mea::sel::SMaxMinX selX;
			//popnew.Write("meda.popn");
			//mPopX.Write("meda.popx0");
			//mPopX.Read("meda.show");
			selX.Select(mPopX,mPopSize);

			//mPopF.Write("meda.popf");
			//mPopX.Write("meda.popx");
			
			t3 = clock();

			mTM += double(t2-t1);
			mTS += double(t3-t2);

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
	}; //class CEA

} //namespace mea

} //namespace az

#endif //AZ_ALGORITHM_DUALPOP_H
