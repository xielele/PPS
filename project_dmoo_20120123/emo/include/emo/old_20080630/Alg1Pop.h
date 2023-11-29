/*! \file	Alg1Pop.h
	
	\brief	Framwork for MOEA with Single Population
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
		
	\date	May.19 2005 create
	\date	May.24 2005 modify
	\date	Sep.28 2005	redesign and rewrite
	\date	May.22 2007 redesign and rewrite
*/


#ifndef	AZ_ALGORITHM_SINGLEPOP_H
#define	AZ_ALGORITHM_SINGLEPOP_H

#include "Alg.h"

#include "Sel.h"
#include "SelG.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
	//!\brief	EA framework with Single Population for Multi-objective Optimization
	//!\param	COG	offspring generator
	//!\param	CES environmental selection
	template<class COG, class CES>
		class CSPMEA : public CBMEA, public COG, public CES
	{
	protected:
		az::mea::gsel::Order0 gsel;
	public:
		//!\brief	constractor
		//!\param	stepmax maximum steps
		//!\param	par parameter object
		//!\param	pop initial population
		//!\return	void
		CSPMEA(
			unsigned int	stepmax	,
			CParameter&		par		,
			CPopulationMO&	pop		)
			: CBMEA(stepmax, par, pop)
		{
			//gsel.Init(pop);
		}

		//!\brief	destructor
		//!\return	void
		~CSPMEA() {}

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
		unsigned int Step()
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
			//mPop = gsel.Select(popnew);
			
			t3 = clock();

			mTM += double(t2-t1);
			mTS += double(t3-t2);

			return ++mStep;
		}
	}; //class CSPMEA

} //namespace mea

} //namespace az

#endif //AZ_ALGORITHM_SINGLEPOP_H
