/*! \file	AlgED.h

	\brief	Framwork for Differential Evolution

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex,
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Oct.14 2005 create
	\date	Dec.02 2005 add COG to a DE-like framework
	\date	May.22 2007 redesign & rewrite
*/


#ifndef	AZ_ALGORITHMDE_H
#define	AZ_ALGORITHMDE_H

#include "Alg.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief	DE-like algorithm framework
//!\param	COG offspring generator
//!\param	CES environmental selection
template<class COG, class CES>
	class CMDE : public EABasic, public COG, public CES
{
public:
	//!\brief	constractor
	//!\param	stepmax maximum steps
	//!\param	par parameter object
	//!\param	pop initial population
	//!\return	void
	CMDE(
		unsigned int	stepmax	,
		CParameter&		par		,
		CPopulationMO&	pop		)
		:EABasic(stepmax, par, pop)
		{
		}

	//!\brief	destructor
	//!\return	no
	virtual ~CMDE() {}

	//!\brief	one step
	//!\return	current step
	virtual unsigned int Step()
	{
		clock_t t1,t2,t3;
		CPopulationMO popnew(EABasic::P());

		t1 = clock();
		//sample new solutions
		COG::Generate(Population().Size(), popnew, Population());

		EABasic::mEvas += EABasic::Population().Size();

		Check(EABasic::Population());
		Check(popnew, EABasic::Population());

		t2 = clock();
		//environmental select
		EABasic::Population().Combine(popnew);
		CES::Select(Population(), mPopSize);

		t3 = clock();

		EABasic::mTM += double(t2-t1);
		EABasic::mTS += double(t3-t2);

		return ++EABasic::mStep;
	}
protected:
	//!\brief	make sure there is no repeated solutions in a population
	//!\param	pop population
	//!\return	void
	void Check(CPopulationMO& pop)
	{
		CPopulationMO tmp(EABasic::P());
		tmp.Combine(pop.At(0));
		for(unsigned int i=1; i<pop.Size(); i++)
			if(!tmp.IsContain(pop[i])) tmp.Combine(pop.At(i));
		pop.Clear();
		pop.Combine(tmp);
	}

	//!\brief	make all new solutions are different from old solutions
	//!\param	popnew offspring population
	//!\param	pop old population
	//!\return	void
	void Check(CPopulationMO& popnew, CPopulationMO& pop)
	{
		CPopulationMO tmp(EABasic::P());
		for(unsigned int i=0; i<popnew.Size(); i++)
			if(!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
		popnew.Clear();
		popnew.Combine(tmp);
	}
}; //class CDE

} //namespace mea

} //namespace az

#endif //AZ_ALGORITHMDE_H
