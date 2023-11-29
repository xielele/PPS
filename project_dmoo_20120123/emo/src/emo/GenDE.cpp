/*! \file	Generator_DE.cpp
	
	\brief	Evolutionary Aglorithm Generators
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Mar.30 2006 redesign
*/

#include <algorithm>
#include <float.h>
#include "emo/Gen.h"


#if defined(WIN32)
    #define wxFinite(n) _finite(n)
#elif defined(_LINUX)
    #define wxFinite(n) finite(n)
#else
    #define wxFinite(n) ((n)==(n))
#endif

namespace az
{
namespace mea
{
namespace gen
{

// constructor
XDE::XDE()
{
	mF		= 1.0;
	mCR		= 1.0;
}

// set parameters
void XDE::Set(double f, double cr)
{
	mF		= f;
	mCR		= cr;
}

// generate new trial solutions
CPopulationMO& XDE::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{
	int dom;
	unsigned int i,j,r1,r2,r3,jrnd,index;
	
	CIndividualMO ind(pop.P());

	popnew.Clear();

	//generate new solutions
	for(i=0; i<pop.Size(); i++) pop[i].Rank(0);
	pop.IsSort(false);
	for(index=0; index<sizenew; index++)
	{
		do{r1=rnd::rand((unsigned int)(0), pop.Size());}while(r1==index);
		do{r2=rnd::rand((unsigned int)(0), pop.Size());}while(r2==index||r2==r1);
		do{r3=rnd::rand((unsigned int)(0), pop.Size());}while(r3==index||r3==r2||r3==r1);

		// generate one solution
		jrnd = rnd::rand((unsigned int)(0), pop.P().XSize());
		for(j=0; j<pop.P().XSize(); j++)
		{
			if(rnd::rand(0.0,1.0)<mCR||j==jrnd)
				ind[j] = pop[r1][j] + mF*(pop[r2][j]-pop[r3][j]);
			else
				ind[j] = pop[index][j];

			if(ind[j]>pop.P().XUpp(j)) 		ind[j] = 0.5*(pop[index][j] + pop.P().XUpp(j));
			else if(ind[j]<pop.P().XLow(j))	ind[j] = 0.5*(pop[index][j] + pop.P().XLow(j));
		}
		
		ind.Evaluate();

		dom = ind.CDominate(pop[index]);
		// dominate 'parent'
		if(dom>0)	
		{
			pop[index] = ind;
		}
		// non comparable with 'parent'
		else if(dom==0) popnew.Combine(ind);
	}

	return popnew;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

XNSDE::XNSDE()
{
	mF		= 0.8;
	mCR		= 0.4;
}

void XNSDE::Set(double f, double cr)
{
	mF		= f;
	mCR		= cr;
}

CPopulationMO& XNSDE::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{
	unsigned int j,r,r1,r2,r3,jrnd,index;
	
	popnew.Resize(sizenew);

	//generate new solutions
	for(index=0; index<popnew.Size(); index++)
	{
		//do{r1=rnd::rand((unsigned int)(0), pop.Size());}while(r1==index);
		//do{r2=rnd::rand((unsigned int)(0), pop.Size());}while(r2==index||r2==r1);
		//do{r3=rnd::rand((unsigned int)(0), pop.Size());}while(r3==index||r3==r2||r3==r1);

		do{r1=rnd::rand((unsigned int)(0), pop.Size()); r=rnd::rand((unsigned int)(0), pop.Size()); if(pop[r].Rank() < pop[r1].Rank()) r1=r;}while(r1==index);
		do{r2=rnd::rand((unsigned int)(0), pop.Size()); r=rnd::rand((unsigned int)(0), pop.Size()); if(pop[r].Rank() < pop[r2].Rank()) r2=r;}while(r2==index||r2==r1);
		do{r3=rnd::rand((unsigned int)(0), pop.Size()); r=rnd::rand((unsigned int)(0), pop.Size()); if(pop[r].Rank() < pop[r3].Rank()) r3=r;}while(r3==index||r3==r2||r3==r1);

		// generate one solution
		jrnd = rnd::rand((unsigned int)(0), pop.P().XSize());
		for(j=0; j<pop.P().XSize(); j++)
		{
			if(rnd::rand(0.0,1.0)<mCR||j==jrnd)
				popnew[index][j] = pop[r1][j] + mF*(pop[r2][j]-pop[r3][j]);
			else
				popnew[index][j] = pop[index][j];

			// if not a leagle float
			if( wxFinite(popnew[index][j]) == 0 ) popnew[index][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));

			if(popnew[index][j]>pop.P().XUpp(j)) 		popnew[index][j] = 0.5*(pop[index][j] + pop.P().XUpp(j));
			else if(popnew[index][j]<pop.P().XLow(j))	popnew[index][j] = 0.5*(pop[index][j] + pop.P().XLow(j));

			//if(popnew[index][j]>pop.P().XUpp(j)) 		popnew[index][j] = pop.P().XUpp(j);
			//else if(popnew[index][j]<pop.P().XLow(j))	popnew[index][j] = pop.P().XLow(j);
		}	

		PM(popnew[index]);
	}
	return popnew;
}

} //namespace gen
} //namespace mea
} //namespace az
