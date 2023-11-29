/*! \file	GenESimplex.cpp
	
	\brief	Extended Simplex Operator
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Dec.08 2008 design
*/

#include <cmath>
#include <algorithm>
#include <vector>
#include "alg//Random.h"
#include "alg//Matrix.h"
#include "emo/Parameter.h"
#include "emo/Gen.h"

namespace az
{
namespace mea
{
namespace gen
{

// constructor
ESimplex::ESimplex()
{
}

// Generator
CPopulationMO& ESimplex::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{	
	unsigned int i,j;

	popnew.Resize(sizenew);

	//mMatDis.resize(pop.Size());
	//for(i=0; i<pop.Size(); i++) mMatDis[i].resize(pop.Size());
	//for(i=0; i<pop.Size(); i++)
	//{
	//	mMatDis[i][i] = 1.0E100;
	//	for(j=i+1; j<pop.Size(); j++)
	//		mMatDis[j][i] = mMatDis[i][j] = Distance(pop[i], pop[j]);
	//}

	//switch(pop.P().FSize())
	//{
	//case 2:
	//	mExten = 0.5;
	//	Generate1D(popnew,pop);
	//	break;
	//case 3:
	//	mExten = 1.414;
	//	Generate2D(popnew,pop);
	//	break;
	//default:
	//	break;
	//}

	return popnew;
}

} //namespace gen
} //namespace mea
} //namespace az
