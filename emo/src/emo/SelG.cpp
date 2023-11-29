// SelG.cpp
	
#include <algorithm>
#include <cmath>
#include "emo/SelG.h"

//!\brief	az namespace, the top namespace
namespace az
{
//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief select strategies
namespace gsel
{

Order0::Order0()
{
	pPop = 0;
}

Order0::~Order0()
{
	if(pPop != 0) delete pPop;
}

void Order0::Init(CPopulationMO& pop)
{
	unsigned int i, nobj = pop.P().FSize();

	mPopSize = pop.Size();

	if(pPop != 0) delete pPop;
	pPop = new CPopulationMO(pop.P());
	(*pPop).Resize(pop.Size());

	mIdea.resize(nobj);
	mWeight.resize(mPopSize);
	for(i=0; i<mPopSize; i++)
	{
		mWeight[i].resize(nobj);
		
		// currently the method only works for bi-objective problems
		if(nobj == 2)
		{
			mWeight[i][0] = double(i)/(mPopSize-1.0);
			mWeight[i][1] = 1.0 - mWeight[i][0];
		}
	}
	for(i=0; i<nobj; i++)
	{
		mIdea[i] = 1.0E100;
	}
	
	UpdateIdea(pop);
	UpdatePop(pop, false);
}

CPopulationMO& Order0::Select(CPopulationMO& off)
{
	UpdateIdea(off);
	UpdatePop(off, true);
	return *pPop;
}

double Order0::SopFit(CIndividualMO& ind, unsigned int ith)
{
	//_TCHE1		
	double fitness = -1.0e+30, diff, feval;
	for(unsigned int n=0; n<ind.P().FSize(); n++)
	{
		diff = fabs(ind.F(n) - mIdea[n]);
		feval= mWeight[ith][n] == 0 ? 0.0001*diff : diff*mWeight[ith][n];
		if(feval>fitness) fitness = feval;
	}

	return fitness;
}

void Order0::UpdateIdea(CPopulationMO& off)
{
	for(unsigned int i=0; i<off.Size(); i++)
		for(unsigned int j=0; j<off.P().FSize(); j++)
			if(off[i].F(j) < mIdea[j])
				mIdea[j] = off[i].F(j);
}

void Order0::UpdatePop(CPopulationMO& off, bool state)
{
	bool replace;
	unsigned int i, j, p0, p1, r;
	double fit0, fit1;
	std::vector<unsigned int> index0(mPopSize), index1(mPopSize), copy(off.Size());
	for(i=0; i<mPopSize; i++) index0[i] = i;
	for(i=0; i<off.Size(); i++) index1[i] = i;
	std::random_shuffle(index0.begin(), index0.end());
	std::random_shuffle(index1.begin(), index1.end());
	for(i=0; i<off.Size(); i++) copy[i] = 2;	// each point can be used to update at most two sub-problems

	for(i=0; i<mPopSize; i++)
	{
		p0 = index0[i];

		fit0 = state ? SopFit((*pPop)[p0], p0) : 1.0E100;

		r = 0; replace = false;
		for(j=0; j<off.Size(); j++) if(copy[index1[j]]>0)
		{
			p1   = index1[j];
			fit1 = SopFit(off[p1], p0);
			if(fit1<fit0)
			{
				r       = p1;
				fit0    = fit1;
				replace = true;
			}
		}
		if(replace)
		{
			copy[r]--;
			(*pPop)[p0] = off[r];
		}
	}
}	

} //namespace gsel
} //namespace mea
} //namespace az
