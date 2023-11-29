/*! \file	Generator_GTX.cpp
	
	\brief	Evolutionary Aglorithm Generators
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Mar.30 2006 redesign
*/
#include <cmath>
#include <algorithm>
#include <vector>
#include "emo/Parameter.h"
#include "alg/Random.h"
#include "emo/Gen.h"

namespace az
{
namespace mea
{
namespace gen
{

/////////////////////////////////////////////////////////////////////////////////////
#define	MaxA	 1.21	//!< constant value
#define	MinA	-0.21	//!< constant value
//GTX	constructor
XGTX::XGTX()
{
}

//GTX-Generator
CPopulationMO& XGTX::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{
	double SumA,tmpMin,tmpMax;
	unsigned int i,k,s;

	mvCoefficient.resize(pop.P().FSize());
	mvMatePool.resize(pop.Size());
	for(i=0; i<pop.Size(); i++) mvMatePool[i] = i;

	popnew.Resize(sizenew);

	for(s=0; s<popnew.Size(); s++)
	{	
		for(i=0; i<(unsigned int)mvCoefficient.size(); i++)	if((unsigned int)(mvMatePool.size()) - i > 1) 
				std::swap(mvMatePool[i], mvMatePool[ i + rnd::rand((unsigned int)(1), (unsigned int)(mvMatePool.size()) - i) ]);

		SumA = 0;
		for(k = 0; k < (unsigned int)mvCoefficient.size() - 1; k++)
		{
			tmpMin	= 1 - SumA - MaxA * (mvCoefficient.size()  - k - 1);
			tmpMax	= 1 - SumA - MinA * (mvCoefficient.size()  - k - 1);

			if(tmpMin < MinA) tmpMin	= MinA ;
			if(tmpMax > MaxA) tmpMax	= MaxA ;

			mvCoefficient[k] = rnd::rand(tmpMin, tmpMax);

			SumA += mvCoefficient[k];
		}
		mvCoefficient[ mvCoefficient.size() - 1 ]	= 1.0 - SumA;

		//crossover
		for(k=0; k<pop.P().XSize(); k++)
		{
			popnew[s][k] = 0.0;
			for(i=0;i<(unsigned int)mvCoefficient.size() ;i++)
				popnew[s][k] += pop[mvMatePool[i]][k] * mvCoefficient[i];

			//border handeling
			if(popnew[s][k]>pop.P().XUpp(k))		popnew[s][k] = 0.5*(pop[s][k] + pop.P().XUpp(k));
			else if(popnew[s][k]<pop.P().XLow(k))	popnew[s][k] = 0.5*(pop[s][k] + pop.P().XLow(k));
		}
	}
	return popnew;
}

} //namespace gen
} //namespace mea
} //namespace az
