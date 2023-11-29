/*! \file	Selection_MaxMin.cpp
	
	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Mar.13 2006 create
	\date	Mar.30 2006 redesign
	\date	Feb.21 2007 add SMaxMinX
*/

#include <algorithm>
#include <fstream>
#include "emo/Sel.h"

namespace az
{
namespace mea
{
namespace sel
{

CPopulationMO& SMaxMin::SelectSort(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()<=size) return pop;

	unsigned int i, j, k, index, start, end;
	double fmin;
	std::vector< std::vector<double> > dismat;			// distance matrix
	std::vector< double >	dismin;						// minimum distance
	std::vector< bool >		exist;						// flag to show whether the point is selected

	//Step 1: rank sort
	pop.RankSort();

	//Step 2: find the sub set which "cover" the cut point
	start = end = 0;
	while(end<size)
	{
		start = end;
		while(end<pop.Size() && pop[end].Rank() == pop[start].Rank()) end++;
	}

	//Step 3: sort the sub set
	if(!pop[start].IsFeasible()) return pop;

	//Step 4: calculate the distance matrix
	dismat.resize(end-start);
	dismin.resize(end-start);
	exist.resize(end-start);
	for(i=0; i<dismat.size(); i++) dismat[i].resize(end-start);
	for(i=0; i<dismat.size(); i++)
	{
		exist[i]	= false;
		dismin[i]	= 1.0E200;
		dismat[i][i]= 1.0E200;
		for(j=i+1; j<dismat.size(); j++)
		{
			dismat[i][j] = 0.0;
			for(k=0; k<pop.P().FSize(); k++) dismat[i][j] += (pop[i+start].F(k) - pop[j+start].F(k))*(pop[i+start].F(k) - pop[j+start].F(k));
			dismat[i][j] = sqrt(dismat[i][j]);
			dismat[j][i] = dismat[i][j];
		}
	}
	
	//Step 5: select the solutions with minimum Fi(x)
	for(k=0; k<pop.P().FSize(); k++)
	{
		fmin = 1.0E200; index = 0;
		for(i=0; i<dismat.size(); i++) 
			if(!exist[i] && pop[i+start].F(k)<fmin)
			{
				fmin  = pop[i+start].F(k);
				index = i;
			}
			exist[index]   = true;
		for(i=0; i<dismat.size(); i++) if(!exist[i] && dismat[i][index]<dismin[i]) dismin[i] = dismat[i][index];
	}

	//Step 5: sort point one by one
	for(i=pop.P().FSize(); i<size-start; i++)
	{
		//find the one to be selected
		index = 0; fmin = -1.0E200;
		for(j=0; j<dismat.size(); j++) if(!exist[j] && dismin[j] > fmin) {fmin=dismin[j]; index=j;} 
		exist[index] = true;

		//update the minimum distance of other points
		for(j=0; j<dismat.size(); j++) if(!exist[j] && dismin[j] > dismat[j][index]) dismin[j] = dismat[j][index];
	}

	//Step 6: resort the population
	int a=0,b=(int)dismat.size()-1;
	while(a<b)
	{
		while(a<(int)dismat.size()	&&  exist[a]) a++;
		while(b>=0					&& !exist[b]) b--;
		if(a<b)
		{
			pop.Swap(a+start,b+start);
			a++; b--;
		}
	}

	return pop;
}

CPopulationMO& SMaxMin::Select(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()>size) SelectSort(pop, size).Erase(size);

	return pop;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

CPopulationMO& SMaxMinX::SelectSort(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()<=size) return pop;

	unsigned int i, j, k, index, start, end, sel;
	double fmin, dtmp;
	std::vector< std::vector<double> > dismat;			// distance matrix
	std::vector< double >	dismin;						// minimum distance
	std::vector< bool >		exist;						// flag to show whether the point is selected

	//Step 1: rank sort
	pop.RankSort();

	//Step 2: find the sub set which "cover" the cut point
	start = end = 0;
	while(end<size)
	{
		start = end;
		while(end<pop.Size() && pop[end].Rank() == pop[start].Rank()) end++;
	}

	//Step 3: sort the sub set
	if(!pop[start].IsFeasible()) return pop;

	//Step 4: calculate the distance matrix
	i=0; while(pop[i].Rank()==1) i++; if(i<size) start = i;	// sort all dominated solutions

	dismat.resize(end-start);
	dismin.resize(end-start);
	exist.resize(end-start);
	for(i=0; i<dismat.size(); i++) dismat[i].resize(end-start);
	for(i=0; i<dismat.size(); i++)
	{
		exist[i]	= false;
		dismin[i]	= 1.0E200;
		dismat[i][i]= 1.0E200;
		for(j=i+1; j<dismat.size(); j++)
		{
			dismat[i][j] = 0.0;
			for(k=0; k<pop.P().XSize(); k++) dismat[i][j] += (pop[i+start][k] - pop[j+start][k])*(pop[i+start][k] - pop[j+start][k]);
			dismat[j][i] = dismat[i][j];
		}
	}

	//Step 5: initialize the mindis
	if(start>0) // some nondominated points selected
	{
		for(i=0; i<end-start; i++)
		{
			dtmp = 0.0;
			for(j=0; j<start; j++)
			{
				for(k=0; k<pop.P().XSize(); k++) dtmp += (pop[i+start][k] - pop[j][k])*(pop[i+start][k] - pop[j][k]);
				if(dtmp < dismin[i]) dismin[i] = dtmp;
			}
		}
		sel = 0;
	}
	else	// select from nondominated points
	{
		// extreme points are selected firstly
		sel  = 0;
		for(k=0; k<pop.P().FSize(); k++)
		{
			fmin = 1.0E200; index = 0;
			for(i=0; i<dismat.size(); i++) 
				if(pop[i+start].F(k)<fmin)
				{
					fmin  = pop[i+start].F(k);
					index = i;
				}
			if(!exist[index])	// it is possible that one point optimizes all objs
			{
				sel ++;
				exist[index]   = true;
				for(i=0; i<dismat.size(); i++) if(!exist[i] && dismat[i][index]<dismin[i]) dismin[i] = dismat[i][index];
			}
		}
	}

	//Step 6: sort point one by one
	for(i=sel; i<size-start; i++)
	{
		//find the one to be selected
		index = 0; fmin = -1.0E200;
		for(j=0; j<dismat.size(); j++) if(!exist[j] && dismin[j] > fmin) {fmin=dismin[j]; index=j;} 
		exist[index] = true;

		//update the minimum distance of other points
		for(j=0; j<dismat.size(); j++) if(!exist[j] && dismin[j] > dismat[j][index]) dismin[j] = dismat[j][index];
	}

	//Step 6: resort the population
	int a=0,b=(int)dismat.size()-1;
	while(a<b)
	{
		while(a<(int)dismat.size()	&&  exist[a]) a++;
		while(b>=0					&& !exist[b]) b--;
		if(a<b)
		{
			pop.Swap(a+start,b+start);
			a++; b--;
		}
	}

	return pop;
}

CPopulationMO& SMaxMinX::Select(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()>size) SelectSort(pop, size).Erase(size);

	return pop;
}

}//namespace sel
} //namespace mea
} //namespace az
