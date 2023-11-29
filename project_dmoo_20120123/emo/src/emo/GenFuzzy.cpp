/*! \file	Generator_Fuzzy.cpp
	
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
// Fuzzy operator

// constructor
XFuzzy::XFuzzy()
{
}

// Fuzzy Generator
CPopulationMO& XFuzzy::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{
	unsigned int i,j,r1,r2;

	popnew.Resize(sizenew);
	
	for(i=0; i<popnew.Size(); i++)
	{
		do{r1=rnd::rand((unsigned int)(0), pop.Size());}while(r1==i);
		do{r2=rnd::rand((unsigned int)(0), pop.Size());}while(r2==i||r2==r1);
		
		for(j=0; j<pop.P().XSize(); j++)
		{
			if(rnd::rand(0.0,1.0)<0.5)
				popnew[i][j] = rnd::triangular(	pop[r1][j]-0.5*fabs(pop[r1][j]-pop[r2][j]),
												pop[r1][j]+0.5*fabs(pop[r1][j]-pop[r2][j]),
												pop[r1][j]);
			else
				popnew[i][j] = rnd::triangular(	pop[r2][j]-0.5*fabs(pop[r1][j]-pop[r2][j]),
												pop[r2][j]+0.5*fabs(pop[r1][j]-pop[r2][j]),
												pop[r2][j]);
			
			if(popnew[i][j]>pop.P().XUpp(j)) 
				popnew[i][j] = 0.5*(pop[i][j] + pop.P().XUpp(j));
			else if(popnew[i][j]<pop.P().XLow(j)) 
				popnew[i][j] = 0.5*(pop[i][j] + pop.P().XLow(j));
		}
	}
	return popnew;
}

} //namespace gen
} //namespace mea
} //namespace az
