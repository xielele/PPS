/*! \file	Selection_Crowd.cpp
	
	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Sep.27 2005 create
	\date	Mar.30 2006 redesign
*/

#include <vector>
#include <list>
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
	CPopulationMO& SCrowd::SelectSort(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()<=size) return pop;

		unsigned int start, end;

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
		if(pop[start].IsFeasible() && end > start + 2 && start < size-2)
		{
			unsigned int i,j,k; double interval;

			std::vector<double>			share(end-start);
			std::vector<unsigned int>	index(end-start);
			
			for(i=0; i<(unsigned int)share.size(); i++) share[i] = 0.0;
			
			//calculate the share values
			for(i=0; i<pop.P().FSize(); i++)
			{
				for(j=start; j<end; j++) index[j-start] = j;

				for(j=0; j<(unsigned int)index.size()-1; j++)
					for(k=j+1; k<(unsigned int)index.size(); k++)
						if(pop[index[j]].F(i) > pop[index[k]].F(i))
							std::swap(index[j], index[k]);

				interval =  pop[index[index.size()-1]].F(i) - pop[index[0]].F(i) ;
				for(j=1; j<(unsigned int)index.size()-1; j++) share[index[j]-start] += (pop[index[j+1]].F(i) - pop[index[j-1]].F(i))/interval;
				share[index[0]-start] = MAXDOUBLE;
				share[index[index.size()-1]-start] = MAXDOUBLE;
			}
			
			//sort the sub-population according to the share value
			for(i=start; i<end; i++)
				for(j=i+1; j<end; j++)
					if(share[i-start] < share[j-start])
						pop.Swap(i, j);

			share.clear();
			index.clear();
		}

		return pop;
	}

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& SCrowd::Select(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()>size) SelectSort(pop, size).Erase(size);

		return pop;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SCrowd2::SelectSort(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()<=size) return pop;

		unsigned int start, end;

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
		if(pop[start].IsFeasible() && end > start + 2 && start < size-2)
		{
			//unsigned int i,j,k,de; double interval;

			//while(end>size)
			//{
			//	std::vector<double>			share(end-start);
			//	std::vector<unsigned int>	index(end-start);
			//	
			//	for(i=0; i<(unsigned int)share.size(); i++) share[i] = 0.0;
			//	
			//	//calculate the share values
			//	for(i=0; i<pop.P().FSize(); i++)
			//	{
			//		for(j=start; j<end; j++) index[j-start] = j;

			//		for(j=0; j<(unsigned int)index.size()-1; j++)
			//			for(k=j+1; k<(unsigned int)index.size(); k++)
			//				if(pop[index[j]].F(i) > pop[index[k]].F(i))
			//					std::swap(index[j], index[k]);

			//		interval =  pop[index[index.size()-1]].F(i) - pop[index[0]].F(i) ;
			//		for(j=1; j<(unsigned int)index.size()-1; j++) share[index[j]-start] += (pop[index[j+1]].F(i) - pop[index[j-1]].F(i))/interval;
			//		share[index[0]-start] = MAXDOUBLE;
			//		share[index[index.size()-1]-start] = MAXDOUBLE;
			//	}
			//	
			//	//find the one to delete 
			//	de = start;
			//	for(i=start+1; i<end; i++) if(share[i-start] < share[de-start]) de=i;

			//	//move the one to the end
			//	pop.Swap(de, end-1);

			//	share.clear();
			//	index.clear();

			//	end--;
			//}

			//2008.12.04
			unsigned int i, j, k, kk, df = pop.P().FSize();
			std::vector< std::vector< unsigned int > >	index(df);
			// sort in each dimension
			for(i=0; i<df; i++)
			{
				index[i].resize(end-start); 
				for(j=start; j<end; j++) index[i][j-start] = j;
				for(j=start; j<end; j++) 
					for(k=j+1; k<end; k++)
						if(pop[index[i][j-start]].F(i) > pop[index[i][k-start]].F(i))
								std::swap(index[i][j-start], index[i][k-start]);
			}
			std::vector< bool > state(end-start); for(i=0; i<end-start; i++) state[i] = true;
			std::vector< double > share(end-start); 
			
			// remove one solution each time
			int il, ii, ir, de;
			for(kk=size; kk<end; kk++)
			{
				for(i=0; i<end-start; i++) share[i] = 0.0;

				for(i=0; i<df; i++)
				{
					double interval =  pop[index[i][end-start-1]].F(i) - pop[index[i][0]].F(i);
					il = -1; ii=0;
					while(ii<end-start)
					{
						while(ii<end-start && !state[index[i][ii]-start]) ii++;
						if(ii<end-start)
						{
							ir = ii+1; while(ir<end-start && !state[index[i][ir]-start]) ir++;
							if(il<0 || ir>=end-start)
								share[index[i][ii]-start] += MAXDOUBLE;
							else
								share[index[i][ii]-start] += (pop[index[i][ir]].F(i) - pop[index[i][il]].F(i))/interval;
						}
						il = ii;
						ii = ir;
					}
				}
				de = -1;
				for(i=0; i<end-start; i++) if(state[i]) if(de < 0 || share[i] < share[de]) de=i;
				state[de] = false;
			}
			
			il = 0; ir = end-start-1;
			while(il<ir)
			{
				while(il<end-start && state[il]) il++;
				while(ir>=0 && !state[ir]) ir--;
				if(il<ir)
				{
					pop.Swap(start+il, start+ir);
					il++;
					ir--;
				}
			}
		}

		return pop;
	}

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& SCrowd2::Select(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()>size) SelectSort(pop, size).Erase(size);

		return pop;
	}

}//namespace sel
} //namespace mea
} //namespace az
