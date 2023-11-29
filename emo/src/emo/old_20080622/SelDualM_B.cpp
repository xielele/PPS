/*! \file	SelDualMM.cpp

	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex,
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Mar.17 2008 create
*/
#include <cmath>
#include <algorithm>
#include <fstream>
#include <vector>
#include <list>

#include <float.h>

#include "alg/amoeba.h"
#include "emo/Sel.h"

//#define AZ_TMP_OUT 1

namespace az
{
namespace mea
{
namespace sel
{
	const double PUNISH = 5.0;	// punishment

	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& DualMM::SelectSort(CPopulationMO& pop, unsigned int size)
	{
		unsigned int i,k,index;
		double fmin;

		if(pop.Size()<=size) return pop;

		// initialization
		mvSel.resize(pop.Size()); for(i=0; i<pop.Size(); i++) mvSel[i] = false;
		pop.RankSort();
		BuildSimplex(pop, size/2);
		InitMatrix(pop);

		// choose extreme points
		mvDisToSet.resize(pop.Size()); for(i=0; i<pop.Size(); i++) mvDisToSet[i] = 1.0E100;
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
			if(i%2 == 0) 
				index = SelF(k++, pop);
			else
				index = SelX(pop);
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
	CPopulationMO& DualMM::Select(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()>size) SelectSort(pop, size).Erase(size);

		return pop;
	}

	void DualMM::Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp)
	{
		unsigned int i,j,k,size = (unsigned int)cenp.size();

		CPopulationMO poptmp(pop);
		poptmp.RankSort();
		BuildSimplex(poptmp, size);
		for(i=0; i<size; i++) for(j=0; j<mDimF; j++)
		{
			cenp[i][j] = 0.0;
			for(k=0; k<mDimF; k++) cenp[i][j] += mvSimplex[k][j]*mvWeight[i][k];
		}
	}

	void DualMM::BuildSimplex(CPopulationMO& pop, unsigned int size)
	{
		unsigned int i, j, k;

		mDimF = pop.P().FSize();
		mDimX = pop.P().XSize();
		mvSimplex.resize(mDimF); for(i=0; i<mDimF; i++) mvSimplex[i].resize(mDimF);
		mNormal.resize(mDimF);

		// Step 1: find the points to form a simplex
		std::vector<double> vmin(mDimF); for(k=0; k<mDimF; k++) vmin[k] = 1.0E100;
		CPopulationMO poptmp(pop); poptmp.RankSort();
		for(i=0; i<mDimF; i++) mvSimplex[i][i] = -1.0E100;
		for(i=0; i<poptmp.Size(); i++) if(poptmp[i].Rank() == poptmp[0].Rank())
		{
			for(j=0; j<mDimF; j++) if( poptmp[i].F(j) > mvSimplex[j][j])
			{
				mvSimplex[j] = poptmp[i].F();
				for(k=0; k<mDimF; k++) if(mvSimplex[j][k] < vmin[k]) vmin[k] = mvSimplex[j][k];
			}
		}

#ifdef AZ_TMP_OUT
		std::ofstream ss("tmp0.dat");
		for(i=0; i<mDimF; i++)
		{
			for(j=0; j<mDimF; j++) ss<<mvSimplex[i][j]<<"\t"; ss<<std::endl;
		}
		for(j=0; j<mDimF; j++) ss<<mvSimplex[0][j]<<"\t"; ss<<std::endl;
		ss.close();
#endif

		// Step 2: find the Quasi-normal direction
		double norm = 0.0;
		if(mDimF == 2)
		{
			mNormal[0] = mvSimplex[1][1] - vmin[1];
			mNormal[1] = mvSimplex[0][0] - vmin[0];
		}
		else if(mDimF == 3)
		{
			mNormal[0] = ( (mvSimplex[1][1]-mvSimplex[0][1])*(mvSimplex[2][2]-mvSimplex[0][2])-(mvSimplex[1][2]-mvSimplex[0][2])*(mvSimplex[2][1]-mvSimplex[0][1]) );
			mNormal[1] = ( (mvSimplex[1][2]-mvSimplex[0][2])*(mvSimplex[2][0]-mvSimplex[0][0])-(mvSimplex[1][0]-mvSimplex[0][0])*(mvSimplex[2][2]-mvSimplex[0][2]) );
			mNormal[2] = ( (mvSimplex[1][0]-mvSimplex[0][0])*(mvSimplex[2][1]-mvSimplex[0][1])-(mvSimplex[1][1]-mvSimplex[0][1])*(mvSimplex[2][0]-mvSimplex[0][0]) );
			double tmp = ( 0-(mNormal[0]*mvSimplex[0][0]+mNormal[1]*mvSimplex[0][1]+mNormal[2]*mvSimplex[0][2]) );
			double tmp1= mNormal[0]*vmin[0] + mNormal[1]*vmin[1] + mNormal[2]*vmin[2] + tmp;
			if(tmp1>0) for(i=0; i<3; i++) mNormal[i] *= -1;
		}
		else
		{
			for(i=0; i<mDimF; i++)
			{
				mNormal[i] = 0.0;
				for(k=0; k<mDimF; k++) mNormal[i] += mvSimplex[k][i] - vmin[i];
			}
		}
		for(i=0; i<mDimF; i++) norm += mNormal[i]*mNormal[i];
		for(i=0; i<mDimF; i++) mNormal[i] /= sqrt(norm);

		// Step 3: move the simplex to make it 'dominates' all current points
		Trans = 1.0E100; double pro;
		for(i=0; i<pop.Size(); i++)
		{
			pro = 0.0;
			for(j=0; j<mDimF; j++) pro += (pop[i].F(j)-mvSimplex[0][j])*mNormal[j];
			if(pro<Trans) Trans = pro;
		}
		for(i=0; i<mDimF; i++) for(j=0; j<mDimF; j++) mvSimplex[i][j] += mNormal[j]*Trans;

		// Step 4: extend the simplex
		std::vector<double> O(mDimF);
		for(i=0; i<mDimF; i++) {O[i] = 0.0; for(j=0; j<mDimF; j++) O[i] += mvSimplex[j][i]; O[i] /= double(mDimF);}
		for(i=0; i<mDimF; i++) for(j=0; j<mDimF; j++) mvSimplex[i][j] += 0.25*(mvSimplex[i][j] - O[j]);

#ifdef AZ_TMP_OUT
		ss.open("tmp1.dat");
		for(i=0; i<mDimF; i++)
		{
			for(j=0; j<mDimF; j++) ss<<mvSimplex[i][j]<<"\t"; ss<<std::endl;
		}
		for(j=0; j<mDimF; j++) ss<<mvSimplex[0][j]<<"\t"; ss<<std::endl;
		ss.close();	
#endif

		// Step 5: set the weights
		if(mDimF == 2)
		{
			mvWeight.resize(size);
			for(i=0; i<size; i++)
			{
				mvWeight[i].resize(2);
				mvWeight[i][0] = double(i)/double(size-1);
				mvWeight[i][1] = 1.0 - mvWeight[i][0];
			}
		}
		else if(mDimF == 3)
		{
			unsigned int num = (unsigned int)sqrt(double(2*size))-1;
			while((num+1)*(num+2)/2 < size) num++;
			mvWeight.resize((num+1)*(num+2)/2);
			k = 0;
			for(i=0; i<=num; i++)
				for(j=0; j<=num-i; j++)
				{
					mvWeight[k].resize(3);
					mvWeight[k][0] = double(i)/double(num);
					mvWeight[k][1] = double(j)/double(num);
					mvWeight[k][2] = double(num-i-j)/double(num);
					k++;
				}
		}
		std::random_shuffle(mvWeight.begin(), mvWeight.end());
	}

	// init distance matrix
	void DualMM::InitMatrix(CPopulationMO& pop)
	{
		unsigned int i,k,s;
		mvDisX.resize(pop.Size()); for(i=0; i<pop.Size(); i++) mvDisX[i].resize(pop.Size());

		for(i=0; i<pop.Size(); i++)
		{
			mvDisX[i][i] = 1.0E100;
			for(k=i+1; k<pop.Size(); k++)
			{
				mvDisX[i][k] = 0.0;
				for(s=0; s<mDimX; s++) mvDisX[i][k] += (pop[i][s] - pop[k][s])*(pop[i][s] - pop[k][s]);
				mvDisX[k][i] = mvDisX[i][k];
			}
		}
	}

	void DualMM::UpdateDistance(unsigned int index)
	{
		for(unsigned int i=0; i<mvDisX.size(); i++) 
			if(!mvSel[i] && mvDisX[i][index]<mvDisToSet[i]) mvDisToSet[i] = mvDisX[i][index];
	}

	unsigned int DualMM::SelF(unsigned int k, CPopulationMO& pop)
	{
		unsigned int i, j, index=0; double dis0, dis1, mindis = 1.0E100;

		std::vector<double> P(mDimF);
		for(i=0; i<mDimF; i++) 
		{
			P[i] = 0.0;
			for(j=0; j<mDimF; j++) P[i] += mvSimplex[j][i]*mvWeight[k][j];
		}
		for(i=0; i<mvSel.size(); i++) if(!mvSel[i])
		{
			dis0 = 0.0; for(j=0; j<mDimF; j++) dis0 += (pop[i].F(j)-P[j])*mNormal[j];
			dis1 = 0.0; for(j=0; j<mDimF; j++) dis1 += (pop[i].F(j)-P[j])*(pop[i].F(j)-P[j]);
			dis1 = sqrt( dis1 - dis0*dis0 );
			if(dis0 + PUNISH*dis1 < mindis)
			{
				mindis = dis0 + PUNISH*dis1;
				index  = i;
			}

			//dis1 = 0.0; for(j=0; j<mDimF; j++) dis1 += (pop[i].F(j)-P[j])*(pop[i].F(j)-P[j]);
			//if(dis1 < mindis)
			//{
			//	mindis = dis1;
			//	index  = i;
			//}
		}
		return index;
	}

	unsigned int DualMM::SelX(CPopulationMO& pop)
	{
		unsigned int i, index = 0, rank = (unsigned int)(mvSel.size()+1); double maxdis = -1.0E100;
		for(i=0; i<mvSel.size(); i++) if(!mvSel[i])
		{
			if(pop[i].Rank()<rank || (pop[i].Rank()==rank && mvDisToSet[i]>maxdis))
			{
				rank   = pop[i].Rank();
				maxdis = mvDisToSet[i];
				index  = i;
			}
		}
		return index;
	}
}//namespace sel
} //namespace mea
} //namespace az
