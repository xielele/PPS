/*! \file	Generator_LS.cpp
	
	\brief	Local Search strategies
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Apr.24 2006 create
*/

#include "mea/Generator_LS.h"

//!\brief EA Evolutionary Algorithms
namespace EA
{
//!\brief offspring generate strategies
namespace GEN
{
//!\brief local search strategies
namespace LS
{

//!\brief	constractor
LSBase::LSBase()
{
	mNoL = mNoE = 1;
}

//!\brief	set parameters
//!\param	noL	number of points to do local search
//!\param	noE maximum evaluations on each search point	
//!\return	void
void LSBase::Set( unsigned int noL, unsigned int noE)
{
	mNoL = noL;
	mNoE = noE;
}

//!\brief	select some points for local search
//!\param	pop		reference population
//!\param	start	start location in the reference population
//!\param	end		end location in the reference population
//!\param	size	the number of points to be selected
//!\param	index	the index of the selected points
//!\return	void	
void LSBase::SelectStartPoint(CPopulationMO& pop, unsigned int start, unsigned int end, unsigned int size, std::vector<unsigned int>& index)
{
	unsigned int i,j,k,next;
	double dis1,dis2;
	if(end-start<size)
	{
		index.resize(end-start);
		for(i=start; i<end; i++) index[i-start] = i;
		return;
	}

	if(size<pop.P().FSize()) size = pop.P().FSize()+1;

	std::vector<bool>	exist(end-start);
	std::vector<double>	mindis(end-start);
	for(i=0; i<exist.size(); i++) {exist[i] = false; mindis[i] = 1.0E200;}
	
	for(k=0; k<pop.P().FSize(); k++)
	{
		dis1 = 1.0E200; next = 0;
		for(i=0; i<exist.size(); i++) if(!exist[i] && pop[i+start].F(k)<dis1){dis1  = pop[i+start].F(k); next  = i;}
		exist[next]  = true;
		index.push_back(next+start);

		for(i=0; i<exist.size(); i++) if(!exist[i])
		{
			dis2 = 0.0;
			for(j=0; j<pop.P().FSize(); j++) dis2 += (pop[next+start].F(j)-pop[i+start].F(j))*(pop[next+start].F(j)-pop[i+start].F(j));
			if(dis2 < mindis[i]) mindis[i] = dis2;
		}
	}

	for(i=pop.P().FSize(); i<size; i++)
	{
		dis1 = -1.0E200; next = 0;
		for(j=0; j<exist.size(); j++) if(!exist[j] && mindis[j]>dis1) {dis1=mindis[j]; next=j;}
		exist[next] = true;
		index.push_back(next+start);

		for(j=0; j<exist.size(); j++) if(!exist[j])
		{
			dis2 = 0.0;
			for(k=0; k<pop.P().FSize(); k++) dis2 += (pop[next+start].F(k)-pop[j+start].F(k))*(pop[next+start].F(k)-pop[j+start].F(k));
			if(dis2 < mindis[j]) mindis[j] = dis2;
		}
	}
}

} //namespace LS

} //namespace GEN

} //namespace EA
