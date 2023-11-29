/*! \file	SelDualCrowd.cpp
	
	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Feb.20 2007 create
	\date	Mar.21 2008 redesign
*/

#include <algorithm>
#include <numeric>
#include <vector>
#include <list>
#include "emo/Sel.h"

namespace az
{
namespace mea
{
namespace sel
{
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& DualCrowd::SelectSort(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()<=size) return pop;

		unsigned int start, end;

		//Step 1: rank sort
		pop.IsSort(false);
		pop.ERankSort();
		//pop.RankSort();

		//Step 2: find the sub set which "cover" the cut point
		start = end = 0;
		while(end<size)
		{
			start = end;
			while(end<pop.Size() && pop[end].Rank() == pop[start].Rank()) end++;
		}

		//Step 3: sort the sub set
		if(pop[start].IsFeasible() && end > start + 2 && start < size-2)
		{
			unsigned int i,j,k; double interval;

			std::vector<double>			cdistx(end-start), cdistf(end-start), cdist(end-start);
			std::vector<unsigned int>	index(end-start);
			
			for(i=0; i<(unsigned int)cdistx.size(); i++) cdistf[i] = cdistx[i] = 0.0;
			
			//calculate the crowded distance values in the objective space
			for(i=0; i<pop.P().FSize(); i++)
			{
				for(j=start; j<end; j++) index[j-start] = j;

				for(j=0; j<(unsigned int)index.size()-1; j++)
					for(k=j+1; k<(unsigned int)index.size(); k++)
						if(pop[index[j]].F(i) > pop[index[k]].F(i))
							std::swap(index[j], index[k]);

				interval =  pop[index[index.size()-1]].F(i) - pop[index[0]].F(i) ;
				for(j=1; j<(unsigned int)index.size()-1; j++) cdistf[index[j]-start] += (pop[index[j+1]].F(i) - pop[index[j-1]].F(i))/interval;
				cdistf[index[0]-start] = MAXDOUBLE;
				cdistf[index[index.size()-1]-start] = MAXDOUBLE;
			}
			//calculate the crowded distance values in the decision space
			for(i=0; i<pop.P().XSize(); i++)
			{
				for(j=start; j<end; j++) index[j-start] = j;

				for(j=0; j<(unsigned int)index.size()-1; j++)
					for(k=j+1; k<(unsigned int)index.size(); k++)
						if(pop[index[j]][i] > pop[index[k]][i])
							std::swap(index[j], index[k]);

				interval =  pop[index[index.size()-1]][i] - pop[index[0]][i] ;
				for(j=1; j<(unsigned int)index.size()-1; j++) cdistx[index[j]-start] += (pop[index[j+1]][i] - pop[index[j-1]][i])/interval;
				cdistx[index[0]-start] += 2.0*(pop[index[1]][i] - pop[index[0]][i])/interval;
				cdistx[index[index.size()-1]-start] = 2.0*(pop[index[index.size()-1]][i] - pop[index[index.size()-2]][i])/interval;
			}

			// average distance
			for(i=0; i<cdistf.size(); i++)
			{
				cdistf[i] /= double(pop.P().FSize());
				cdistx[i] /= double(pop.P().XSize());
			}
			double adf = std::accumulate(cdistf.begin(), cdistf.end(), 0.0)/double(cdistf.size());
			double adx = std::accumulate(cdistx.begin(), cdistx.end(), 0.0)/double(cdistf.size());
			for(i=0; i<cdistf.size(); i++)
			{
				if(cdistf[i] > adf || cdistx[i] > adx)
					cdist[i] = std::max(cdistf[i], cdistx[i]);
				else
					cdist[i] = std::min(cdistf[i], cdistx[i]);
			}
			
			//sort the sub-population according to the share value
			for(i=start; i<end; i++)
				for(j=i+1; j<end; j++)
					if(cdist[i-start] < cdist[j-start])
						pop.Swap(i, j);
		}

		return pop;
	}

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& DualCrowd::Select(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()>size) SelectSort(pop, size).Erase(size);

		return pop;
	}

}//namespace sel
} //namespace mea
} //namespace az
