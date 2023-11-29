/*! \file	SelDualMaxiMin.cpp

	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex,
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Apr.25 2008 create
*/
#include <cmath>
#include <algorithm>
#include <fstream>
#include <vector>
#include <list>

#include <float.h>

#include "alg/amoeba.h"
#include "emo/Sel.h"

namespace az
{
namespace mea
{
namespace sel
{
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& DualMaxiMin::SelectSort(CPopulationMO& pop, unsigned int size, unsigned int method)
	{
		unsigned int i,k,index;
		double fmin;

		if(pop.Size()<=size) return pop;

		mDimF = pop.P().FSize();
		mDimX = pop.P().XSize();
		mvSel.resize(pop.Size());
		for(i=0; i<pop.Size(); i++) mvSel[i] = false;

		// initialization
		pop.RankSort();
		InitMatrix(pop);

		// choose extreme points
		mvDisToSetF.resize(pop.Size()); for(i=0; i<pop.Size(); i++) mvDisToSetF[i] = 1.0E100;
		mvDisToSetX.resize(pop.Size()); for(i=0; i<pop.Size(); i++) mvDisToSetX[i] = 1.0E100;
		for(k=0; k<pop.P().FSize(); k++)
		{
			fmin = 1.0E200; index = 0;
			for(i=0; i<pop.Size(); i++) 
			{
				if(!mvSel[i] && pop[i].F(k)<fmin)
				{
					fmin  = pop[i].F(k);
					index = i;
				}
				mvSel[index]   = true;
				UpdateDistance(index);
			}
		}

		// select points
		k = 0;
		for(i=pop.P().FSize(); i<size; i++)
		{
			switch(method)
			{
			case 0:
				index = SelF(pop);
				break;
			case 1:
				index = SelX(pop);
				break;
			case 2:
			default:
				if(i%2 == 0) 
					index = SelF(pop);
				else
					index = SelX(pop);
				break;
			}
			mvSel[index] = true;
			UpdateDistance(index);
		}

		//Step 5: sort the final population
		int a=0,b=(int)pop.Size()-1;
		while(a<b)
		{
			while(a<(int)pop.Size()	&&  mvSel[a]) a++;
			while(b>=0				&& !mvSel[b]) b--;
			if(a<b)
			{
				pop.Swap(a,b);
				a++; b--;
			}
		}

		return pop;
	}

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& DualMaxiMin::Select(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()>size) SelectSort(pop, size, 2).Erase(size);

		return pop;
	}

	void DualMaxiMin::Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp)
	{
		CPopulationMO poptmp(pop);
		SelectSort(poptmp, (unsigned int)cenp.size(), 1);

		unsigned int i;
		for(i=0; i<cenp.size(); i++)
			cenp[i] = poptmp[i].F();; 
	}

	// init distance matrix
	void DualMaxiMin::InitMatrix(CPopulationMO& pop)
	{
		unsigned int i,k,s;
		mvDisX.resize(pop.Size()); for(i=0; i<pop.Size(); i++) mvDisX[i].resize(pop.Size());
		mvDisF.resize(pop.Size()); for(i=0; i<pop.Size(); i++) mvDisF[i].resize(pop.Size());

		for(i=0; i<pop.Size(); i++)
		{
			mvDisX[i][i] = 1.0E100;
			mvDisF[i][i] = 1.0E100;
			for(k=i+1; k<pop.Size(); k++)
			{
				mvDisX[i][k] = 0.0; for(s=0; s<mDimX; s++) mvDisX[i][k] += (pop[i][s] - pop[k][s])*(pop[i][s] - pop[k][s]); mvDisX[k][i] = mvDisX[i][k];
				mvDisF[i][k] = 0.0; for(s=0; s<mDimF; s++) mvDisF[i][k] += (pop[i].F(s) - pop[k].F(s))*(pop[i].F(s) - pop[k].F(s)); mvDisF[k][i] = mvDisF[i][k];
			}
		}
	}

	void DualMaxiMin::UpdateDistance(unsigned int index)
	{
		for(unsigned int i=0; i<mvDisX.size(); i++) 
		{
			if(!mvSel[i] && mvDisX[i][index]<mvDisToSetX[i]) mvDisToSetX[i] = mvDisX[i][index];
			if(!mvSel[i] && mvDisF[i][index]<mvDisToSetF[i]) mvDisToSetF[i] = mvDisF[i][index];
		}
	}

	unsigned int DualMaxiMin::SelF(CPopulationMO& pop)
	{
		unsigned int i, index = 0, rank = (unsigned int)(mvSel.size()+1); double maxdis = -1.0E100;
		for(i=0; i<mvSel.size(); i++) if(!mvSel[i])
		{
			if(pop[i].Rank()<rank || (pop[i].Rank()==rank && mvDisToSetF[i]>maxdis))
			{
				rank   = pop[i].Rank();
				maxdis = mvDisToSetF[i];
				index  = i;
			}
		}
		return index;
	}

	unsigned int DualMaxiMin::SelX(CPopulationMO& pop)
	{
		unsigned int i, index = 0, rank = (unsigned int)(mvSel.size()+1); double maxdis = -1.0E100;
		for(i=0; i<mvSel.size(); i++) if(!mvSel[i])
		{
			if(pop[i].Rank()<rank || (pop[i].Rank()==rank && mvDisToSetX[i]>maxdis))
			{
				rank   = pop[i].Rank();
				maxdis = mvDisToSetX[i];
				index  = i;
			}
		}
		return index;
	}
}//namespace sel
} //namespace mea
} //namespace az
