/*! \file	Selection_Neighbor.cpp
	
	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Mar.30 2006 redesign
*/

//Selection.cpp

#include <algorithm>
#include "emo/Sel.h"

namespace az
{
namespace mea
{
namespace sel
{
///////////////////////////////////////////////////////////////////////////////
// crowded sort(SPEA2) strategy
CPopulationMO& SNeighbor::SelectSort(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()<=size) return pop;

	int start, end;

	//Step 1: rank sort
	pop.RankSort();

	//Step 2: find the sub set which "cover" the cut point
	start = end = 0;
	while(end<int(size))
	{
		start = end;
		while( end<int(pop.Size()) && pop[end].Rank() == pop[start].Rank() ) end++;
	}

	//Step 3: sort the sub set
	mRange	= pop.P().FSize();
	if(pop[start].IsFeasible())	Sort(pop, start, end, size-start);

	return pop;
}

//distance selection(SPEA2) strategy
CPopulationMO& SNeighbor::Select( CPopulationMO& pop, unsigned int size )
{
	SelectSort(pop, size).Erase(size);
	
	return pop;
}

//sort the population [start, end)
//after the sort, [start, start+remainsize) will be kept and others will be removed
CPopulationMO& SNeighbor::Sort(CPopulationMO& pop, int start, int end, int remainsize)
{
	if(end <= start + 2 ||remainsize <= int(pop.P().FSize()) )	return pop;

	int i, j, k; double tmp;

	mSize	= end - start;

	//Step 1: calculate the distance matrix
	mDis.resize(mSize-1);
	for( i=0; i<mSize-1; i++ ) mDis[i].resize(mSize-1-i);
	for( i=0; i<mSize-1; i++ )
		for( j=i+1; j<mSize; j++ )
		{
			tmp = 0.0;
			for(k=0; k<int(pop.P().FSize()); k++) 
				tmp += (pop[i+start].F(k)-pop[j+start].F(k))*(pop[i+start].F(k)-pop[j+start].F(k));
			DisMat(i,j) = tmp;
		}
	//Step 2: initialize node vector
	SetPar();

	//Step 3: eliminate one at a time
	for( i=0; i<mSize-remainsize; i++ ) EraseOne();

	//Step 4: move all solutions to be deleted to the 'tail'
	i=0; j=mSize-1;
	while( i<j )
	{
		while( i<mSize && vNode[i].bKeep  ) i++;
		while( j>=0 && !vNode[j].bKeep    ) j--;
		if(i < j)
		{
			pop.Swap(start+i, start+j);
			i++; j--;
		}
	}

	vNode.clear();
	mDis.clear();

	return pop;
}

double& SNeighbor::DisMat( int i, int j) 
{ 
	return i<j ? mDis[i][j-i-1] : mDis[j][i-j-1]; 
}

void SNeighbor::SetPar()
{
	int i, j, s, k;

	vNode.resize(mSize);
	for( i=0; i<mSize; i++ )
	{
		vNode[i].bKeep = true;
		vNode[i].vMin.resize( mRange );
		//find the mRange minimum to i
		vNode[i].vMin[0] = (i==0) ? 1 : 0;
		s=1; j= vNode[i].vMin[0]+1;
		while(s<mRange)
		{
			if(j==i) j++;
			for(k=s-1; k>=0; k--)
			{
				if( DisMat(i,j)<DisMat(i, vNode[i].vMin[k]) ) vNode[i].vMin[k+1] = vNode[i].vMin[k];
				else break;
			}
			vNode[i].vMin[k+1] = j;
			j++;s++;
		}
		for(; j<mSize; j++)
		{
			if(j==i) j++;
			if( j<mSize && DisMat(i,j)<DisMat(i, vNode[i].vMin[mRange-1]) )
			{
				for(k=mRange-2; k>=0; k--)
				{
					if( DisMat(i,j)<DisMat(i, vNode[i].vMin[k]) ) vNode[i].vMin[k+1] = vNode[i].vMin[k];
					else break;
				}
				vNode[i].vMin[k+1] = j;
			}
		}
		for(j=0; j<mRange; j++) vNode[vNode[i].vMin[j]].vInfluence.push_back(i);
	}
}

void SNeighbor::EraseOne()
{
	int minindex=-1,i,j,s,k,tmpA, tmpB, tmp=-1;
	std::list<int>::iterator it;

	//find the one with minimum distance
	for( i=0; i<mSize; i++ )
		if( vNode[i].bKeep && ( minindex<0 || DisMat(i,vNode[i].vMin[0]) < DisMat(minindex,vNode[minindex].vMin[0]) ))
				minindex = i;
	//two candidate solutions to erase
	tmpA = minindex; tmpB = vNode[minindex].vMin[0];
	//find the one to erase
	for(i=1; i<mRange; i++)
	{
		if( DisMat(tmpA, vNode[tmpA].vMin[i]) < DisMat(tmpB, vNode[tmpB].vMin[i])-1.0E-10 )
		{
			tmp = tmpA; break;
		}
		if( DisMat(tmpA, vNode[tmpA].vMin[i]) > DisMat(tmpB, vNode[tmpB].vMin[i])+1.0E-10 )
		{
			tmp = tmpB; break;
		}
	}
	if(tmp<0) tmp = rnd::rand(0.0,1.0) < 0.5 ? tmpA : tmpB;
	//erase tmp
	vNode[tmp].bKeep = false;
	//modify corresponding solutions
	it = vNode[tmp].vInfluence.begin();
	while(it!=vNode[tmp].vInfluence.end())
	{
		s = *it++;
		k = 0;
		while(vNode[s].vMin[k] != tmp) k++;
		for(;k<mRange-1;k++) vNode[s].vMin[k] = vNode[s].vMin[k+1];
		vNode[s].vMin[mRange-1] = -1;
		for(j=0; j<mSize; j++)
			if( j!=s && vNode[j].bKeep &&
				DisMat(s,j) >= DisMat(s, vNode[s].vMin[mRange-2]) && j != vNode[s].vMin[mRange-2] &&
				( vNode[s].vMin[mRange-1] < 0 || DisMat(s, j)<DisMat(s,vNode[s].vMin[mRange-1]) ) )
				vNode[s].vMin[mRange-1] = j;
		vNode[vNode[s].vMin[mRange-1]].vInfluence.push_back( s );
	}
	for(i=0; i<mRange; i++)
	{
		it = std::find( vNode[vNode[tmp].vMin[i]].vInfluence.begin(), vNode[vNode[tmp].vMin[i]].vInfluence.end(), tmp );
		vNode[vNode[tmp].vMin[i]].vInfluence.erase( it );
	}
}

// crowded sort strategy
CPopulationMO& SNeighbor2::SelectSort(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()<=size) return pop;

	int start, end;

	//Step 1: rank sort
	pop.RankSort();

	//Step 2: find the sub set which "cover" the cut point
	start = end = 0;
	while(end<int(size))
	{
		start = end;
		while( end<int(pop.Size()) && pop[end].Rank() == pop[start].Rank() ) end++;
	}

	//Step 3: sort the sub set
	//the only difference to CDistance
	mRange	= 2;//pop.P().FSize();
	if(pop[start].IsFeasible())	Sort(pop, start, end, size-start);

	return pop;
}

void SNeighbor2::EraseOne()
{
	int minindex=-1,i,j,s,k;
	std::list<int>::iterator it;

	//find the one with minimum distance
	for( i=0; i<mSize; i++ )
		if( vNode[i].bKeep && ( minindex<0 || DisMat(i,vNode[i].vMin[1]) < DisMat(minindex,vNode[minindex].vMin[1]) ))
				minindex = i;

	//erase
	vNode[minindex].bKeep = false;
	//modify corresponding solutions
	it = vNode[minindex].vInfluence.begin();
	while(it!=vNode[minindex].vInfluence.end())
	{
		s = *it++;
		k = 0;
		while(vNode[s].vMin[k] != minindex) k++;
		for(;k<mRange-1;k++) vNode[s].vMin[k] = vNode[s].vMin[k+1];
		vNode[s].vMin[mRange-1] = -1;
		for(j=0; j<mSize; j++)
			if( j!=s && vNode[j].bKeep &&
				DisMat(s,j) >= DisMat(s, vNode[s].vMin[mRange-2]) && j != vNode[s].vMin[mRange-2] &&
				( vNode[s].vMin[mRange-1] < 0 || DisMat(s, j)<DisMat(s,vNode[s].vMin[mRange-1]) ) )
				vNode[s].vMin[mRange-1] = j;
		vNode[vNode[s].vMin[mRange-1]].vInfluence.push_back( s );
	}
	for(i=0; i<mRange; i++)
	{
		it = std::find( vNode[vNode[minindex].vMin[i]].vInfluence.begin(), vNode[vNode[minindex].vMin[i]].vInfluence.end(), minindex );
		vNode[vNode[minindex].vMin[i]].vInfluence.erase( it );
	}
}

// assign rank
CPopulationMO& SStrength::RankAssignment(CPopulationMO& pop)
{
	int s, t,better; unsigned int i, j, size; double maxrank;

	//Step 0: Move all infeasible solutions to the tail of the sequence
	s=0; t=pop.Size()-1;
	while(s<t)
	{
		while(s<int(pop.Size()) && pop[s].IsFeasible()) s++;
		while(t>=0 && !pop[t].IsFeasible()) t--;
		if(s<t)	{pop.Swap(s,t); s++; t--;}
	}
	size = s;

	std::vector< std::vector< int> >	vDom(size); 
	std::vector< double > vS(size);
	for(i=0; i<size; i++)
	{
		vS[i] = 0;
		vDom[i].resize(size); vDom[i][i] = 0;
	}
	
	mFitness.resize(pop.Size());
	for(i=0; i<pop.Size(); i++) mFitness[i] = 0.0;

	//Step 1: assign domination matrix
	for(i=0; i<size; i++)
	{
		for (j=i+1 ; j<size; j++)
		{
			better = pop[i].Dominate(pop[j]);
			//i is dominated by j
			if(better<0) { vS[j]++; vDom[i][j] = -1; vDom[j][i] = 1; }
			//j is dominated by i
			else if(better>0) { vS[i]++; vDom[i][j] = 1; vDom[j][i] = -1; }
			else vDom[i][j] = vDom[j][i] = 0;
		}	//end for
	}	// end for

	//Step 2: assign rank values
	maxrank = 0;
	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++) if(vDom[i][j]<0) mFitness[i] += vS[j];
		if(mFitness[i]>maxrank) maxrank = mFitness[i];
	}

	//Step 3: assign rank to infeasible solutions
	for(i=size; i<pop.Size(); i++) mFitness[i] = maxrank+1.0;

	//Step 4: Sort the population by rank
	for(s=0; s<int(pop.Size()-1); s++) 
		for(t=s+1; t<int(pop.Size()); t++) 
			if(mFitness[t]<mFitness[s]) 
			{
				pop.Swap(s,t);
				std::swap(mFitness[s], mFitness[t]);
			}
	return pop;
}

// sort population, the first 'size' are best ones
CPopulationMO& SStrength::SelectSort(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()<=size) return pop;

	int start, end;

	//Step 1: rank sort
	RankAssignment(pop);

	//Step 2: find the sub set which "cover" the cut point
	start = end = 0;
	while(end<int(size))
	{
		start = end;
		while( end<int(pop.Size()) && mFitness[end] == mFitness[start] ) end++;
	}

	//Step 3: sort the sub set
	//condition 1: start > 0(is not the first layer)
	if(start>0 || !pop[start].IsFeasible())
	{
		StrengthSort(pop, start, end);
	}
	//condtion 2: start = 0(to elimate some points in the first layer)
	else
	{
		mRange	= pop.P().FSize();
		Sort(pop, start, end, size-start);
	}

	return pop;
}

// sort the population [start, end) according to strength
CPopulationMO& SStrength::StrengthSort(CPopulationMO& pop, unsigned int start, unsigned int end)
{
	unsigned int i,j,k,index;
	std::vector< double > vdis(pop.Size());
	index = (unsigned int)(sqrt(pop.Size()+0.0));
	for(i=start; i<end; i++)
	{
		for(j=0; j<pop.Size(); j++) 
			if(j!=i)
			{
				vdis[j] = 0.0;
				for(k=0; k<pop.P().FSize(); k++) vdis[j] += (pop[i].F(k)-pop[j].F(k))*(pop[i].F(k)-pop[j].F(k));
			}
			else vdis[j] = 1.0E100;
		std::sort(vdis.begin(), vdis.end());
		mFitness[i] += 1.0/(sqrt(vdis[index])+2.0);
	}

	for(i=start; i<end; i++)
		for(j=i+1; j<end; j++)
			if(mFitness[j]<mFitness[i])
			{
				pop.Swap(i,j);
				std::swap(mFitness[i],mFitness[j]);
			}

	return pop;
}

///////////////////////////////////////////////////////////////////////////////
// NSGA-II + SPEA2 strategy
CPopulationMO& SCrowdNeighbor::SelectSort(CPopulationMO& pop, unsigned int size)
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
	if(!pop[start].IsFeasible()) return pop;

	//distance matrix
	unsigned int i,j,k;
	std::vector< bool > exist(end);
	std::vector< std::vector<double> > mdis(end);
	for(i=0; i<end; i++)
	{
		exist[i] = true;
		mdis[i].resize(end);
	}
	for(i=0; i<end; i++)
	{
		for(j=i+1; j<end; j++)
		{
			mdis[i][j] = 0.0;
			for(k=0; k<pop.P().FSize(); k++) mdis[i][j] += (pop[i].F(k)-pop[j].F(k))*(pop[i].F(k)-pop[j].F(k));
			mdis[j][i] = mdis[i][j];
		}
		mdis[i][i] = 1.0E200;
	}
	
	double mindensity,density1,density2; unsigned int index;
	for(k=0; k<end-size; k++)
	{
		mindensity = 1.0E200; index = start;
		for(i=start; i<end; i++) if(exist[i])
		{
			density1 = density2 = 1.0E200;
			for(j=0; j<end; j++) if(exist[j])
			{
				if(mdis[i][j]<density1) {density2=density1; density1=mdis[i][j];}
				else if(mdis[i][j]<density2) density2=mdis[i][j];
			}
			if(density2<mindensity) {mindensity=density2; index=i;}
			//if(density1+density2<mindensity) {mindensity=density1+density2; index=i;}
		}
		exist[index]=false;
	}

	i=start; j=end-1;
	while(i<j)
	{
		while(i<end && exist[i]) i++;
		while(j>=start && !exist[j]) j--;
		if(i<j)
		{
			pop.Swap(i,j);
			i++; j--;
		}
	}

	return pop;
}

CPopulationMO& SCrowdNeighbor::Select(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()>size) SelectSort(pop, size).Erase(size);

	return pop;
}

///////////////////////////////////////////////////////////////////////////////
// NSGA-II + SPEA2 strategy
CPopulationMO& SCrowdNeighbor2::SelectSort(CPopulationMO& pop, unsigned int size)
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
	if(!pop[start].IsFeasible()) return pop;

	//distance matrix
	unsigned int i,j,k;
	std::vector< std::vector<double> > mdis(end);
	for(i=0; i<end; i++)
	{
		mdis[i].resize(end);
	}
	for(i=0; i<end; i++)
	{
		for(j=i+1; j<end; j++)
		{
			mdis[i][j] = 0.0;
			for(k=0; k<pop.P().FSize(); k++) mdis[i][j] += (pop[i].F(k)-pop[j].F(k))*(pop[i].F(k)-pop[j].F(k));
			mdis[j][i] = mdis[i][j];
		}
		mdis[i][i] = 1.0E200;
	}

	std::vector<double> density(end-start);
	double density1,density2;
	for(i=start; i<end; i++)
	{
		density1 = density2 = 1.0E200;
		for(j=0; j<end; j++)
		{
			if(mdis[i][j]<density1) {density2=density1; density1=mdis[i][j];}
			else if(mdis[i][j]<density2) density2=mdis[i][j];
		}
		density[i-start]=density2;
	}

	for(i=start; i<end; i++)
		for(j=i+1; j<end; j++) 
			if(density[i-start]>density[j-start])
			{
				pop.Swap(i,j);
				std::swap(density[i-start], density[j-start]);
			}

	return pop;
}

CPopulationMO& SCrowdNeighbor2::Select(CPopulationMO& pop, unsigned int size)
{
	if(pop.Size()>size) SelectSort(pop, size).Erase(size);

	return pop;
}

}//namespace sel
} //namespace mea
} //namespace az
