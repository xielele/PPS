/*! \file	SelDualMM_OLD_old.cpp

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

namespace az
{
namespace mea
{
namespace sel
{
	const double PUNISH = 2.0;	// punishment

	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& DualMM_OLD::SelectSort(CPopulationMO& pop, unsigned int size)
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
		if(pop.P().FSize() == 2)
			BuildModel2(pop, size/2);
		else //if(pop.P().FSize() == 3)
			BuildModel3(pop, size/2);

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
	CPopulationMO& DualMM_OLD::Select(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()>size) SelectSort(pop, size).Erase(size);

		return pop;
	}

	void DualMM_OLD::Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp)
	{
		unsigned int i,k,size = (unsigned int)cenp.size();
		mDimF = pop.P().FSize();
		mDimX = pop.P().XSize();

		if(pop.P().FSize() == 2)
		{
			BuildModel2(pop, size);
			for(i=0; i<size; i++) for(k=0; k<2; k++)
				cenp[i][k] = A[k] + (B[k] - A[k])*(i+rnd::rand(0.0,1.0))/double(size); 

			//std::ofstream ss("tmp.dat");
			//for(i=0; i<size; i++) ss<<cenp[i][0]<<"\t"<<cenp[i][1]<<std::endl;
			//ss.close();
		}
		else
		{
			BuildModel3(pop, size);
			for(i=0; i<size; i++) for(k=0; k<3; k++)
				cenp[i][k] = A[k]*mvWeight[i][0] + B[k]*mvWeight[i][1] + C[k]*mvWeight[i][2]; 
		}
	}

	// build reference points
	void DualMM_OLD::BuildModel2(CPopulationMO& pop, unsigned int size)
	{
		unsigned int i, k;
		// Step 1: set initial reference points (two extreme points)
		unsigned int iA, iB;
		A.resize(mDimF); B.resize(mDimF);
		iA = iB = 0;
		for(i=1; i<pop.Size(); i++)
		{
			if((0.95*pop[i].F(0)+0.05*pop[i].F(1)) < (0.95*pop[iA].F(0)+0.05*pop[iA].F(1)))	iA = i;
			if((0.95*pop[i].F(1)+0.05*pop[i].F(0)) < (0.95*pop[iB].F(1)+0.05*pop[iB].F(0)))	iB = i;
		}
		while(iB == iA) iB = rnd::rand((unsigned int)0, pop.Size());
		if(pop[iB].F(0) < pop[iA].F(0)) std::swap(iA, iB);

		A = pop[iA].F(); B = pop[iB].F();

		//Step 2: calculate the projections onto AB and OAB
		//directions AB and OAB
		std::vector<double> DAB(mDimF),			// direction vector AB 
							DOAB(mDimF);		// direction vector OAB

		mvProL.resize(pop.Size()),	// projection on AB
		mvProOL.resize(pop.Size());	// projection on OAB
		double	LAB		=  0.0,					// length of AB
				LLR		=  0.0,					// length of LR
				tranOAB =  1.0E100,				// the translation on OAB
				promaxAB= -1.0E100,
				prominAB=  1.0E100;				// projection on AB
		for(i=0; i<mDimF; i++) {DAB[i]  = B[i]-A[i]; LAB += DAB[i]*DAB[i];} LAB = sqrt(LAB); 
		for(i=0; i<mDimF; i++)  DAB[i] /= LAB;
		DOAB[0]		= -DAB[1]; 
		DOAB[1]		=  DAB[0]; 
		for(i=0; i<pop.Size(); i++)
		{
			mvProOL[i] = 0.0; mvProL[i] = 0.0;
			for(k=0; k<mDimF; k++) {mvProOL[i] += (pop[i].F(k)-A[k])*DOAB[k]; mvProL[i] += (pop[i].F(k)-A[k])*DAB[k];}
			if(mvProOL[i]<tranOAB)  tranOAB = mvProOL[i];
			if(mvProL[i]<prominAB)	prominAB= mvProL[i];
			if(mvProL[i]>promaxAB)	promaxAB= mvProL[i];
		}
		for(i=0; i<pop.Size(); i++) mvProOL[i] += tranOAB;

		//Step 4: transformation of reference line AB
		prominAB = -0.25*LAB;	//extension
		promaxAB =  1.25*LAB;
		for(i=0; i<2; i++)
		{
			double tmp = B[i]-A[i];
			A[i] = A[i] + tranOAB*DOAB[i] - 0.25*tmp;
			B[i] = B[i] + tranOAB*DOAB[i] + 0.25*tmp;
		}

		LLR		 = promaxAB-prominAB;
		mvTarLen.resize(size);	//target projection
		for(i=0; i<size; i++) mvTarLen[i] = LLR*double(i)/double(size-1)+prominAB;
		std::random_shuffle(mvTarLen.begin(), mvTarLen.end());
	}

	// build reference points
	void DualMM_OLD::BuildModel3(CPopulationMO& pop, unsigned int size)
	{
		unsigned int i, j, k;
		// Step 1: set initial reference points (three extreme points)
		unsigned int iA, iB, iC;
		std::vector<double> O(mDimF); A.resize(mDimF); B.resize(mDimF); C.resize(mDimF);
		//CPopulationMO poptmp(pop);
		//poptmp.RankSort();
		iA = iB = iC = 0;
		for(i=1; i<pop.Size(); i++) //if(poptmp[i].Rank() == poptmp[0].Rank())
		{
			if((0.05*pop[i].F(0)+0.475*pop[i].F(1)+0.475*pop[i].F(2)) < (0.05*pop[iA].F(0)+0.475*pop[iA].F(1)+0.475*pop[iA].F(2)))	iA = i;
			if((0.05*pop[i].F(1)+0.475*pop[i].F(0)+0.475*pop[i].F(2)) < (0.05*pop[iB].F(1)+0.475*pop[iB].F(0)+0.475*pop[iB].F(2)))	iB = i;
			if((0.05*pop[i].F(2)+0.475*pop[i].F(0)+0.475*pop[i].F(1)) < (0.05*pop[iC].F(2)+0.475*pop[iC].F(0)+0.475*pop[iC].F(1)))	iC = i;
			
			//if(poptmp[i].F(0) > poptmp[iA].F(0)) iA = i;
			//if(poptmp[i].F(1) > poptmp[iB].F(1)) iB = i;
			//if(poptmp[i].F(2) > poptmp[iC].F(2)) iC = i;

			//if((0.95*poptmp[i].F(0)+0.025*poptmp[i].F(1)+0.025*poptmp[i].F(2)) > (0.95*poptmp[iA].F(0)+0.025*poptmp[iA].F(1)+0.025*poptmp[iA].F(2)))	iA = i;
			//if((0.95*poptmp[i].F(1)+0.025*poptmp[i].F(0)+0.025*poptmp[i].F(2)) > (0.95*poptmp[iB].F(1)+0.025*poptmp[iB].F(0)+0.025*poptmp[iB].F(2)))	iB = i;
			//if((0.95*poptmp[i].F(2)+0.025*poptmp[i].F(0)+0.025*poptmp[i].F(1)) > (0.95*poptmp[iC].F(2)+0.025*poptmp[iC].F(0)+0.025*poptmp[iC].F(1)))	iC = i;

		}
		A = pop[iA].F(); B = pop[iB].F(); C = pop[iC].F();

		//unsigned int iA1, iB1, iC1;
		//double disa = -1.0E100,disb = -1.0E100,disc = -1.0E100; 
		//for(i=0; i<poptmp.Size(); i++) if(poptmp[i].Rank() == poptmp[0].Rank())
		//{
		//	if(poptmp[i].F(0) < poptmp[iA].F(0)) if(poptmp[i].F(0)>disa) {iA1 = i; disa = poptmp[i].F(0);}
		//	if(poptmp[i].F(1) < poptmp[iB].F(1)) if(poptmp[i].F(1)>disb) {iB1 = i; disb = poptmp[i].F(1);}
		//	if(poptmp[i].F(2) < poptmp[iC].F(2)) if(poptmp[i].F(2)>disc) {iC1 = i; disc = poptmp[i].F(2);}
		//}
		//A = poptmp[iA].F(); B = poptmp[iB].F(); C = poptmp[iC].F();

		std::ofstream ss("tmp0.dat");
		for(i=0; i<3; i++) ss<<A[i]<<"\t"; ss<<std::endl;
		for(i=0; i<3; i++) ss<<B[i]<<"\t"; ss<<std::endl;
		for(i=0; i<3; i++) ss<<C[i]<<"\t"; ss<<std::endl;
		for(i=0; i<3; i++) ss<<A[i]<<"\t"; ss<<std::endl;
		ss.close();

		// Step 2: plane function
		Fa = ( (B[1]-A[1])*(C[2]-A[2])-(B[2]-A[2])*(C[1]-A[1]) );
		Fb = ( (B[2]-A[2])*(C[0]-A[0])-(B[0]-A[0])*(C[2]-A[2]) );
		Fc = ( (B[0]-A[0])*(C[1]-A[1])-(B[1]-A[1])*(C[0]-A[0]) );
		Fd = ( 0-(Fa*A[0]+Fb*A[1]+Fc*A[2]) );
		//normalize
		double norm = sqrt(Fa*Fa+Fb*Fb+Fc*Fc);
		Fa /= norm; Fb /= norm; Fc /= norm; Fd /= norm;
		////orthogonal diretion
		//if(A[0]*Fa + B[1]*Fb + C[2]*Fc + Fd > 0)
		//{
		//	Fa = -Fa;
		//	Fb = -Fb;
		//	Fc = -Fc;
		//	Fd = -Fd;
		//}

		// Step 3: move the plane to make it 'dominates' all current points
		Trans = 1.0E100; double pro;
		for(i=0; i<pop.Size(); i++)
		{
			pro = pop[i].F(0)*Fa + pop[i].F(1)*Fb + pop[i].F(2)*Fc + Fd;
			if(pro<Trans) Trans = pro;
		}
		A[0] += Fa*Trans; A[1] += Fb*Trans; A[2] += Fc*Trans;
		B[0] += Fa*Trans; B[1] += Fb*Trans; B[2] += Fc*Trans;
		C[0] += Fa*Trans; C[1] += Fb*Trans; C[2] += Fc*Trans;

		//ss.open("tmpff.dat");
		//ss<<Fa<<"\t"<<Fb<<"\t"<<Fc<<"\t"<<Fd<<"\t"<<Trans<<std::endl;
		//ss.close();

		//ss.open("tmp1.dat");
		//for(i=0; i<3; i++) ss<<A[i]<<"\t"; ss<<std::endl;
		//for(i=0; i<3; i++) ss<<B[i]<<"\t"; ss<<std::endl;
		//for(i=0; i<3; i++) ss<<C[i]<<"\t"; ss<<std::endl;
		//for(i=0; i<3; i++) ss<<A[i]<<"\t"; ss<<std::endl;
		//ss.close();

		// Step 4: extend ABC
		for(i=0; i<3; i++) O[i] = (A[i] + B[i] + C[i])/3.0;
		for(i=0; i<3; i++)
		{
			A[i] += 0.25*(A[i] - O[i]);
			B[i] += 0.25*(B[i] - O[i]);
			C[i] += 0.25*(C[i] - O[i]);
		}
		ss.open("tmp2.dat");
		for(i=0; i<3; i++) ss<<A[i]<<"\t"; ss<<std::endl;
		for(i=0; i<3; i++) ss<<B[i]<<"\t"; ss<<std::endl;
		for(i=0; i<3; i++) ss<<C[i]<<"\t"; ss<<std::endl;
		for(i=0; i<3; i++) ss<<A[i]<<"\t"; ss<<std::endl;
		ss.close();		


		// Step 5: set the weights
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
		std::random_shuffle(mvWeight.begin(), mvWeight.end());

		//ss.open("tmpcen.dat");
		//for(k=0; k<size; k++)
		//{
		//	for(i=0; i<3; i++) ss<<A[i]*mvWeight[k][0] + B[i]*mvWeight[k][1] + C[i]*mvWeight[k][2]<<"\t"; ss<<std::endl;
		//}
		//ss.close();	
	}

	// init distance matrix
	void DualMM_OLD::InitMatrix(CPopulationMO& pop)
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

	void DualMM_OLD::UpdateDistance(unsigned int index)
	{
		for(unsigned int i=0; i<mvDisX.size(); i++) 
			if(!mvSel[i] && mvDisX[i][index]<mvDisToSet[i]) mvDisToSet[i] = mvDisX[i][index];
	}

	unsigned int DualMM_OLD::SelF(unsigned int k, CPopulationMO& pop)
	{
		unsigned int i, index=0; double dis0, dis1, mindis = 1.0E100;

		if(pop.P().FSize() == 2)
		{
			for(i=0; i<mvSel.size(); i++) if(!mvSel[i])
				if(mvProOL[i] + PUNISH*fabs(mvProL[i]-mvTarLen[k]) < mindis)
				{
					mindis = mvProOL[i] + PUNISH*fabs(mvProL[i]-mvTarLen[k]);
					index  = i;
				}
		}
		else //if(pop.P().FSize() == 3)
		{
			std::vector<double> P(3);
			for(i=0; i<3; i++) P[i] = A[i]*mvWeight[k][0] + B[i]*mvWeight[k][1] + C[i]*mvWeight[k][2];
			for(i=0; i<mvSel.size(); i++) if(!mvSel[i])
			{
				//dis0 = pop[i].F(0)*Fa + pop[i].F(1)*Fb + pop[i].F(2)*Fc + Fd + Trans;
				//dis1 = sqrt( (pop[i].F(0)-P[0])*(pop[i].F(0)-P[0]) + (pop[i].F(1)-P[1])*(pop[i].F(1)-P[1]) + (pop[i].F(2)-P[2])*(pop[i].F(2)-P[2]) - dis0*dis0 );
				//if(dis0 + PUNISH*dis1 < mindis)
				//{
				//	mindis = dis0 + PUNISH*dis1;
				//	index  = i;
				//}

				double tmpmax = -1.0E100;
				for(unsigned int s=0; s<3; s++)
				{
					if((pop[i].F(s)-P[s])*mvWeight[k][s] > tmpmax) tmpmax = (pop[i].F(s)-P[s])*mvWeight[k][s];
				}
				if(tmpmax < mindis)
				{
					mindis = tmpmax;
					index  = i;
				}
			}
		}
		return index;
	}

	unsigned int DualMM_OLD::SelX(CPopulationMO& pop)
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
