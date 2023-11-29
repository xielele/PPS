/*! \file	Alg2Pop.h
	
	\brief	Framwork for MOEA with bi-populations, one for focuses on F-space and one focuses on X-space
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
		
	\date	Feb.21 2007 create
	\date	Apr.22 2007 redesign and rewrite
	\date	May.09 2007 rewrite Step 1 in Step() function
*/


#ifndef	AZ_ALGORITHM_BIPOP_H
#define	AZ_ALGORITHM_BIPOP_H

#include <algorithm>
#include <vector>
#include "Alg.h"
#include "Sel.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
	//!\brief	EA framework for bi-population version
	//!\param	COG	offspring generator
	//!\param	CES environmental selection
	template<class COG, class CES>
		class CBPMEA : public CBMEA, public COG, public CES
	{
	protected:
		CPopulationMO	mPopF,	//!< archive which focuses on F-space
						mPopX;	//!< archive which focuses on X-space
	public:
		//!\brief	constractor
		//!\param	stepmax maximum steps
		//!\param	par parameter object
		//!\param	pop initial population
		//!\return	void
		CBPMEA	(
			unsigned int	stepmax	,
			CParameter&		par		,
			CPopulationMO&	pop		)
			: CBMEA(stepmax, par, pop), mPopF(pop), mPopX(pop)
		{
		}

		//!\brief	destructor
		//!\return	void
		~CBPMEA() {}

		//!\brief	get the population
		//!\return	reference to population
		CPopulationMO& Population() {return mPopF;}

		//!\brief	get the references of the populations
		//!\param	popf reference to pop focusing on F space
		//!\param	popx reference to pop focusing on X space
		//!\return	void
		void Population(CPopulationMO*& popf, CPopulationMO*& popx) 
		{
			popf = &mPopF;
			popx = &mPopX;
		}

		//!\breif	save data
		//!\brief	path data file path
		void Save(std::string& path)
		{
		   std::string fname = path + std::string(".popf");
		   mPopF.Write(fname.c_str());
		   fname = path + std::string(".popx");
		   mPopX.Write(fname.c_str());
		   CPopulationMO tmp(mPopF);
		   tmp.Copy(mPopX);
		   fname = path + std::string(".pop");
		   tmp.Write(fname.c_str());
		}

		//!\brief	one step 
		//!\return	current step
		unsigned int Step()
		{
			clock_t t1,t2,t3;
			CPopulationMO popnew(P());

			t1 = clock();

			//==============================================
			//Step 1: select points from PopX and PopF to generate offspring
			//strategy, use the union of popF and popX
			mPop.Unite(mPopF,mPopX);
			//==============================================
			
			//==============================================
			//Step 2: sample new solutions
			//sample 
			COG::Generate(mPopF.Size(), popnew, mPop);
			
			//==============================================
			//Step 3: evaluate new solutions
			//check the range of new solutions, let them leagel
			Check(popnew);
			//remove these repeated solutions
			Check(popnew, mPopF);
			Check(popnew, mPopX);
			//evaluate new solutions	
			popnew.Evaluate(); mEvas += popnew.Size();	

			//popnew.Write("meda.new");
			//==============================================

			// t2-t1 is the time to generate new solutions
			t2 = clock(); 
			
			//==============================================
			//Step 4: environmental selection
			//update popF
			mPopF.Copy(popnew);
			az::mea::sel::SRegF2 selF;
			selF.Select(mPopF,mPopSize);
			//update popX
			mPopX.Copy(popnew);
			az::mea::sel::SMaxMinX selX;
			selX.Select(mPopX,mPopSize);
			//==============================================

			//mPop.Write("All.pop");
			//popnew.Write("New.pop");

			// t3-t2 is the time to select
			t3 = clock();

			mTM += double(t2-t1);
			mTS += double(t3-t2);

			return ++mStep;
		}
	}; //class CBPMEA

} //namespace mea

} //namespace az

#endif //AZ_ALGORITHM_BIPOP_H
