/*! \file	SelRegF2.cpp

	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex,
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Apr.30 2007 create
	\date	May.04 2007 redesign
	\date	Oct.08 2007 redesign
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
	CPopulationMO& SRegF2::SelectSort(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()<=size) return pop;

		std::vector<bool> vExist;
		//SelectLine(pop, size, vExist, true);
		Select(pop, size, vExist);

		//Step 5: sort the final population
		int a=0,b=(int)pop.Size()-1;
		while(a<b)
		{
			while(a<(int)pop.Size()	&& !vExist[a]) a++;
			while(b>=0				&&  vExist[b]) b--;
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
	CPopulationMO& SRegF2::Select(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()>size) SelectSort(pop, size).Erase(size);

		return pop;
	}

	// ============================================================
	// Oct.17, 2007

	// select some reference central points from the population
	void SRegF2::Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp)
	{
		unsigned int noM = (unsigned int)cenp.size();

		mNoRef = 2;

		mIsModel = false;

		ChooseRef(pop);

		unsigned int i,k,d;
		std::vector<double> DAB(mDim);

		k=0; double clen = 0;
		for(i=0; i<noM; i++)
		{
			double len = mLen*(i+rnd::rand(0.0,1.0))/double(noM);
			if(i==0||len>clen)
			{
				clen+=mvLen[k]; k++; while(len>clen) {clen+=mvLen[k]; k++;}
				for(d=0; d<mDim; d++) DAB[d]  = (mvRef[k][d] - mvRef[k-1][d])/mvLen[k-1];
			}
			if(cenp[i].size()==mDim)
			{
				for(d=0; d<mDim; d++) cenp[i][d] = mvRef[k-1][d] + (len-clen+mvLen[k-1])*DAB[d];
			}
			else
			{
				unsigned int s, ind = 0; double dis, mindis = 1.0E100;
				for(s=0; s<pop.Size(); s++) 
				{
					dis = 0.0; 
					for(d=0; d<mDim; d++) dis+=pow(pop[s].F(d)-(mvRef[k-1][d] + (len-clen+mvLen[k-1])*DAB[d]),2.0);  
					if(dis<mindis) {mindis=dis; ind=s;}
				}
				cenp[i] = pop[ind].X();
			}
		}

#ifdef AZ_MODEL_OUT
		std::ofstream os("Ref.out");
		os<<std::scientific<<std::setprecision(10);
		for(i=0; i<noM; i++)
			os<<cenp[i][0]<<"\t"<<cenp[i][1]<<std::endl;
		os<<std::endl;
#endif
	}

	CPopulationMO& SRegF2::Select(CPopulationMO& pop, unsigned int size, std::vector<bool>& flag)
	{
		ChooseRef(pop);

		mIsModel = false;

		unsigned int i,k,d,s;
		std::vector<double> DAB(mDim), DOAB(mDim), proAB(pop.Size()), proOAB(pop.Size());
		flag.resize(pop.Size()); for(i=0; i<pop.Size(); i++) flag[i] = true;

		k=0; double clen = 0.0; 
		for(i=0; i<size; i++)
		{
			double maxab = -1.0E100;
			double len = mLen*double(i)/double(size-1);
			if((i==0||len>clen)&&i!=(size-1))
			{
				clen+=mvLen[k]; k++; while(len>clen) {clen+=mvLen[k]; k++;}
				for(d=0; d<mDim; d++) DAB[d]  = (mvRef[k][d] - mvRef[k-1][d])/mvLen[k-1];
				DOAB[0]	= -DAB[1]; DOAB[1] =  DAB[0]; 
				for(s=0; s<pop.Size(); s++) if(flag[s])
				{
					proOAB[s] = proAB[s] = 0.0;
					for(d=0; d<mDim; d++) {proOAB[s] += (pop[s].F(d)-mvRef[k-1][d])*DOAB[d]; proAB[s] += (pop[s].F(d)-mvRef[k-1][d])*DAB[d];}
					if(fabs(proOAB[s])>maxab) maxab=fabs(proOAB[s]);
				}
				
			}
			unsigned int id=pop.Size(); double optobj = 1.0E100, obj;
					
			for(s=0; s<pop.Size(); s++) if(flag[s])
			{
				obj = proOAB[s] + PUNISH*fabs(proAB[s]-(len-clen+mvLen[k-1]));
				//obj = proOAB[s] + pow(proAB[s]-(len-clen+mvLen[k-1]),2.0)/(maxab-fabs(proOAB[s])+1.0E-10);
				//obj = proOAB[s] + fabs(proAB[s]-(len-clen+mvLen[k-1]))*(fabs(proOAB[s])+fabs(proAB[s]-(len-clen+mvLen[k-1]))+1);
				//obj = proOAB[s] + fabs(proAB[s]-(len-clen+mvLen[k-1]))*std::max(fabs(proOAB[s]),fabs(proAB[s]-(len-clen+mvLen[k-1])));
				//obj = proOAB[s] + fabs(proAB[s]-(len-clen+mvLen[k-1]))*(fabs(proAB[s]-(len-clen+mvLen[k-1])));
				//obj = sqrt(pow(proAB[s]-(len-clen+mvLen[k-1]),2.0)+pow(proOAB[s],2.0));//+fabs(proAB[s]-(len-clen+mvLen[k-1]));
				if(id==pop.Size() || obj<optobj) {id=s; optobj=obj;}
			}
			flag[id] = false;
		}
		return pop;
	}

	// build reference points
	void SRegF2::ChooseRef(CPopulationMO& pop)
	{
		if(mIsModel) return;

		// Strategy 2: the central points are in the reference line
		unsigned int i, k;
		mDim = pop.P().FSize();
		mvRef.resize(mNoRef);
		for(i=0; i<mNoRef; i++) mvRef[i].resize(mDim);
		mvLen.resize(mNoRef-1);

		// Step 1: set initial reference points (two extreme points)
		unsigned int iA, iB;
		std::vector<double> A(mDim),B(mDim);
		iA = iB = 0;
		for(i=1; i<pop.Size(); i++)
		{
			if((0.95*pop[i].F(0)+0.05*pop[i].F(1)) < (0.95*pop[iA].F(0)+0.05*pop[iA].F(1)))	iA = i;
			if((0.95*pop[i].F(1)+0.05*pop[i].F(0)) < (0.95*pop[iB].F(1)+0.05*pop[iB].F(0)))	iB = i;
		}
		while(iB == iA) iB = rnd::rand((unsigned int)0, pop.Size());
		if(pop[iB].F(0) < pop[iA].F(0)) std::swap(iA, iB);

		A = pop[iA].F(); B = pop[iB].F();

		//Step 3: calculate the projections onto AB and OAB
		//directions AB and OAB
		std::vector<double> DAB(mDim),			// direction vector AB 
							DOAB(mDim),			// direction vector OAB
							proAB(pop.Size()),	// projection on AB
							proOAB(pop.Size());	// projection on OAB
		double	LAB		=  0.0,					// length of AB
				LLR		=  0.0,					// length of LR
				tranOAB =  1.0E100,				// the translation on OAB
				promaxAB= -1.0E100,
				prominAB=  1.0E100;				// projection on AB
		for(i=0; i<mDim; i++) {DAB[i]  = B[i]-A[i]; LAB += DAB[i]*DAB[i];} LAB = sqrt(LAB); 
		for(i=0; i<mDim; i++)  DAB[i] /= LAB;
		DOAB[0]		= -DAB[1]; 
		DOAB[1]		=  DAB[0]; 
		for(i=0; i<pop.Size(); i++)
		{
			proOAB[i] = 0.0; proAB[i] = 0.0;
			for(k=0; k<mDim; k++) {proOAB[i] += (pop[i].F(k)-A[k])*DOAB[k]; proAB[i] += (pop[i].F(k)-A[k])*DAB[k];}
			if(proOAB[i]<tranOAB)  tranOAB = proOAB[i];
			if(proAB[i]<prominAB)	prominAB= proAB[i];
			if(proAB[i]>promaxAB)	promaxAB= proAB[i];
		}
		//Step 4: transformation of reference line AB
		prominAB = -0.25*LAB;
		promaxAB =  1.25*LAB;
		//double alpha = 1.25*std::min(fabs(prominAB), fabs(promaxAB-LAB));
		//alpha = std::max(alpha, 0.25*LAB);
		//prominAB = -alpha;
		//promaxAB =  LAB+alpha;
		//prominAB = std::min(prominAB*0.5, -0.25*LAB);
		//promaxAB = LAB+ std::max((promaxAB-LAB)*0.5,  0.25*LAB);

		LLR		 = promaxAB-prominAB;

		//Step 5: set ference points
		if(mNoRef>2)
		{
			for(k=0; k<mNoRef; k++)
			{
				double lr = LLR*double(k)/double(mNoRef-1)+prominAB;

				double mindis = 1.0E100; unsigned int mini = 0; double dis;
				for(i=0; i<pop.Size(); i++)
				{
					dis = (proAB[i]-lr)*(proAB[i]-lr) +
						  (proOAB[i]-tranOAB)*(proOAB[i]-tranOAB);
					if(dis<mindis){mindis=dis; mini=i;}
				}
				for(i=0; i<mDim; i++) mvRef[k][i] = A[i]+lr*DAB[i]+(proOAB[mini]+0.1*tranOAB)*DOAB[i];
			}
		}
		else
		{
			mvRef[0][0]			= A[0] + prominAB*DAB[0] + tranOAB*DOAB[0];
			mvRef[0][1]			= A[1] + prominAB*DAB[1] + tranOAB*DOAB[1];
			mvRef[mNoRef-1][0]	= A[0] + promaxAB*DAB[0] + tranOAB*DOAB[0];
			mvRef[mNoRef-1][1]	= A[1] + promaxAB*DAB[1] + tranOAB*DOAB[1];
		}
		mLen = 0.0;
		for(k=0; k<mNoRef-1; k++)
		{
			mvLen[k] = 0.0; for(i=0; i<mDim; i++) mvLen[k] += (mvRef[k+1][i]-mvRef[k][i])*(mvRef[k+1][i]-mvRef[k][i]); mvLen[k] = sqrt(mvLen[k]);
			mLen += mvLen[k];
		}

		// converge too fast
		if(mLen <= 0.5*mLen0 && mCount < 5)
		{
			mLen  = mLen0;
			mvLen = mvLen0;
			mvRef = mvRef0;
			mCount++;
		}
		else
		{
			mLen0  = mLen;
			mvLen0 = mvLen;
			mvRef0 = mvRef;
			mCount = 0;
		}

		mIsModel = true;

#ifdef AZ_MODEL_OUT
	std::ofstream os("ModF.out");
	os<<std::scientific<<std::setprecision(10);
	os<<mvRef[0][0]<<"\t"<<mvRef[0][1]<<"\t"<<mvRef[mNoRef-1][0]<<"\t"<<mvRef[mNoRef-1][1]<<std::endl;
	for(i=0; i<mNoRef-1; i++)
	{
		for(k=0; k<mDim; k++) os<<mvRef[i][k]<<"\t";
		for(k=0; k<mDim; k++) os<<mvRef[i+1][k]<<"\t";
		os<<std::endl;
	}
	os.close();
#endif
	}
	
	// Nov.01, 2007
	double SRegF2::objective(std::vector<double>& x)
	{
		CIndividualMO ind(pLPOP->P());
		ind.X() = x;
		ind.Check();
		x = ind.X();
		ind.Evaluate();
		pLPOP->Copy(ind);
		double proOAB = 0.0, proAB = 0.0;
		for(unsigned int d=0; d<mDim; d++) {proOAB += (ind.F(d)-LC[d])*LDOAB[d]; proAB += (ind.F(d)-LC[d])*LDAB[d];}
		return proOAB + PUNISH*fabs(proAB);
	}
	CPopulationMO& SRegF2::LocalSearch(CPopulationMO& popc, CPopulationMO& pop, unsigned int size)
	{
		ChooseRef(pop);

		mIsModel = false;

		popc.Clear();
		pLPOP = &popc;

		unsigned int i,k,d,s;
		std::vector<double> proAB(pop.Size()), proOAB(pop.Size());
		LDAB.resize(mDim); LDOAB.resize(mDim); LC.resize(mDim);

		k=0; double clen = 0.0;
		for(i=0; i<size; i++)
		{
			double len = mLen*double(i)/double(size-1);
			if((i==0||len>clen)&&i!=(size-1))
			{
				clen+=mvLen[k]; k++; while(len>clen) {clen+=mvLen[k]; k++;}
				for(d=0; d<mDim; d++) LDAB[d]  = (mvRef[k][d] - mvRef[k-1][d])/mvLen[k-1];
				LDOAB[0]	= -LDAB[1]; LDOAB[1] =  LDAB[0]; 
				for(s=0; s<pop.Size(); s++)
				{
					proOAB[s] = proAB[s] = 0.0;
					for(d=0; d<mDim; d++) {proOAB[s] += (pop[s].F(d)-mvRef[k-1][d])*LDOAB[d]; proAB[s] += (pop[s].F(d)-mvRef[k-1][d])*LDAB[d];}
				}
			}
			
			for(d=0; d<mDim; d++) LC[d] = mvRef[k-1][d] + (len-clen+mvLen[k-1])*LDAB[d];

			// find a group of solutions to do local search
			unsigned int NL = 31, id;
			std::vector< std::vector<double> > X(NL); std::vector<double> Y(NL);
			double optobj = 1.0E100, optobj1 = -1.0E100, obj;
			for(d=0; d<NL; d++)
			{
				for(s=0; s<pop.Size(); s++)
				{
					obj = proOAB[s] + PUNISH*fabs(proAB[s]-(len-clen+mvLen[k-1]));
					if(obj > optobj1 && obj < optobj) {id=s; optobj=obj;}
				}
				Y[d]   = optobj1 = optobj;
				X[d]   = pop[id].X();
				optobj = 1.0E100;
			}
			unsigned int nfunk;
			amoeba::opt(X, Y, 1.0E-3, nfunk, 2*NL);
		}
		return popc;
	}


	// ============================================================
	// Oct.09, 2007
	// select some reference central points from the population
	void SRegF2::SelectCentre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp, bool inmodel)
	{
		unsigned int noM = (unsigned int)cenp.size();
		
		if(!inmodel)
		{
			// Strategy 1: the central points are in the population
			unsigned int i, j;
			std::vector<bool> vExist;
			SelectLine(pop, noM, vExist, true);
			i=0; for(j=0; j<pop.Size(); j++) if(!vExist[j]) cenp[i++] = pop[j].F();
		}
		else
		{
			// Strategy 2: the central points are in the reference line
			unsigned int i, k, iDim = pop.P().FSize();			// dimension of the space (2: default)
			std::vector<double> refL(iDim), refR(iDim);			// reference points, left and right

			// Step 1: set initial reference points (two extreme points)
			unsigned int iA, iB;
			iA = iB = 0;
			for(i=1; i<pop.Size(); i++)
			{
				if((0.95*pop[i].F(0)+0.05*pop[i].F(1)) < (0.95*pop[iA].F(0)+0.05*pop[iA].F(1)))	iA = i;
				if((0.95*pop[i].F(1)+0.05*pop[i].F(0)) < (0.95*pop[iB].F(1)+0.05*pop[iB].F(0)))	iB = i;
			}
			while(iB == iA) iB = rnd::rand((unsigned int)0, pop.Size());

			//Step 3: calculate the projections onto AB and OAB
			//directions AB and OAB
			std::vector<double> DAB(iDim),			// direction vector AB 
								DOAB(iDim);			// direction vector OAB
			double	LAB		= 0.0,					// length of AB
					tranOAB = 1.0E100,				// the translation on OAB
					proOAB,							// projection on OAB according to point A
					promaxAB,prominAB,proAB;		// projection on AB
			for(i=0; i<iDim; i++) {DAB[i]  = pop[iB].F(i)-pop[iA].F(i); LAB += DAB[i]*DAB[i];} LAB = sqrt(LAB); 
			for(i=0; i<iDim; i++)  DAB[i] /= LAB;
			DOAB[0]		= -DAB[1]; 
			DOAB[1]		=  DAB[0]; 
			promaxAB	= -1.0E100; 
			prominAB	=  1.0E100;
			for(i=0; i<pop.Size(); i++)
			{
				proOAB = 0.0; proAB = 0.0;
				for(k=0; k<iDim; k++) {proOAB += (pop[i].F(k)-pop[iA].F(k))*DOAB[k]; proAB += (pop[i].F(k)-pop[iA].F(k))*DAB[k];}
				if(proOAB<tranOAB)  tranOAB = proOAB;
				if(pop[i].Rank()==1 && proAB<prominAB)	prominAB= proAB;
				if(pop[i].Rank()==1 && proAB>promaxAB)	promaxAB= proAB;
			}
			//Step 4: transformation of reference line AB
			prominAB = std::min(prominAB, -0.25*LAB);
			promaxAB = std::max(promaxAB,  1.25*LAB);
			refL[0]	= pop[iA].F(0) + prominAB*DAB[0] + tranOAB*DOAB[0];
			refL[1]	= pop[iA].F(1) + prominAB*DAB[1] + tranOAB*DOAB[1];
			refR[0]	= pop[iA].F(0) + promaxAB*DAB[0] + tranOAB*DOAB[0];
			refR[1]	= pop[iA].F(1) + promaxAB*DAB[1] + tranOAB*DOAB[1];
			//refL[0]	= pop[iA].F(0) - 0.25*(pop[iB].F(0)-pop[iA].F(0)) + tranOAB*DOAB[0];
			//refL[1]	= pop[iA].F(1) - 0.25*(pop[iB].F(1)-pop[iA].F(1)) + tranOAB*DOAB[1];
			//refR[0]	= pop[iB].F(0) + 0.25*(pop[iB].F(0)-pop[iA].F(0)) + tranOAB*DOAB[0];
			//refR[1]	= pop[iB].F(1) + 0.25*(pop[iB].F(1)-pop[iA].F(1)) + tranOAB*DOAB[1];

			//Step 5: set ference points
			for(i=0; i<noM; i++)
			{
				double dis = double(i)/double(noM-1);
				for(k=0; k<iDim; k++) cenp[i][k] = refL[k] + dis*(refR[k]-refL[k]);
			}
		}
	}

	// Basic Strategy: select some best individuals to the next generation by using line model
	CPopulationMO& SRegF2::SelectLine(CPopulationMO& pop, unsigned int size, std::vector<bool>& flag, bool fix)
	{
		unsigned int i, k, index, iDim = pop.P().FSize();	// dimension of the space (2: default)
		std::vector<double> refL(iDim), refR(iDim);			// reference points, left and right

		// Step 1: set initial reference points
		unsigned int iA, iB;
		iA = iB = 0;
		for(i=1; i<pop.Size(); i++)
		{
			if((0.95*pop[i].F(0)+0.05*pop[i].F(1)) < (0.95*pop[iA].F(0)+0.05*pop[iA].F(1)))	iA = i;
			if((0.95*pop[i].F(1)+0.05*pop[i].F(0)) < (0.95*pop[iB].F(1)+0.05*pop[iB].F(0)))	iB = i;
		}
		while(iB == iA) iB = rnd::rand((unsigned int)0, pop.Size());
		
		//Step 2: calculate the projections onto AB and OAB
		//directions AB and OAB
		std::vector<double> DAB(iDim),			// direction vector AB 
							DOAB(iDim),			// direction vector OAB
							proAB(pop.Size()),	// projection on AB according to point A
							proOAB(pop.Size());	// projection on OAB according to point A
		double	LAB		=  0.0,					// length of AB
				LLR		=  0.0,					// length of LR
				tranOAB =  1.0E100,				// the translation on OAB	
				promaxAB= -1.0E100,
				prominAB=  1.0E100;
		for(i=0; i<iDim; i++) {DAB[i]  = pop[iB].F(i)-pop[iA].F(i); LAB += DAB[i]*DAB[i];} LAB = sqrt(LAB); 
		for(i=0; i<iDim; i++)  DAB[i] /= LAB;
		DOAB[0]		= -DAB[1]; 
		DOAB[1]		=  DAB[0]; 
		for(i=0; i<pop.Size(); i++)
		{
			proAB[i] = proOAB[i] = 0.0;
			for(k=0; k<iDim; k++) {proAB[i]+=(pop[i].F(k)-pop[iA].F(k))*DAB[k]; proOAB[i]+=(pop[i].F(k)-pop[iA].F(k))*DOAB[k];}
			if(proOAB[i]<tranOAB)	tranOAB  = proOAB[i];
			if(pop[i].Rank()==1 && proAB[i]<prominAB)	prominAB = proAB[i];
			if(pop[i].Rank()==1 && proAB[i]>promaxAB)	promaxAB = proAB[i];
		}

		//Step 3: calculate the reference line LR
		prominAB = std::min(prominAB, -0.25*LAB);
		promaxAB = std::max(promaxAB,  1.25*LAB);
		refL[0]	= pop[iA].F(0) + prominAB*DAB[0] + tranOAB*DOAB[0];
		refL[1]	= pop[iA].F(1) + prominAB*DAB[1] + tranOAB*DOAB[1];
		refR[0]	= pop[iA].F(0) + promaxAB*DAB[0] + tranOAB*DOAB[0];
		refR[1]	= pop[iA].F(1) + promaxAB*DAB[1] + tranOAB*DOAB[1];
		//refL[0]	= pop[iA].F(0) - 0.25*(pop[iB].F(0)-pop[iA].F(0)) + tranOAB*DOAB[0];
		//refL[1]	= pop[iA].F(1) - 0.25*(pop[iB].F(1)-pop[iA].F(1)) + tranOAB*DOAB[1];
		//refR[0]	= pop[iB].F(0) + 0.25*(pop[iB].F(0)-pop[iA].F(0)) + tranOAB*DOAB[0];
		//refR[1]	= pop[iB].F(1) + 0.25*(pop[iB].F(1)-pop[iA].F(1)) + tranOAB*DOAB[1];

		LLR = promaxAB-prominAB; //1.5*LAB;//

#ifdef AZ_MODEL_OUT
	// save reference line AB
	std::ofstream os("ModF.out");
	os<<std::scientific<<std::setprecision(10);
	os<<refL[0]<<"\t"<<refL[1]<<"\t"<<refR[0]<<"\t"<<refR[1]<<std::endl;
	os.close();
#endif

		//Step 5: select according to the reference line (reference points)
		flag.resize(pop.Size()); for(i=0; i<pop.Size(); i++) flag[i] = true;
		std::vector<unsigned int> vind(size); for(i=0; i<size; i++) vind[i]=i; std::random_shuffle(vind.begin(), vind.end());
		double proab, obj, optobj;
		unsigned int id;
		for(k=0; k<size; k++)
		{
			index = vind[k];	// the order of the subproblems is random
			if(fix)
				proab = double(index)*LLR/double(size)-0.25*LAB;
			else
				proab = rnd::rand(std::max(0.0,index-0.35), std::min(index+0.35,size-1.0))*LLR/double(size)-0.25*LAB;
			
			id = pop.Size(); optobj = 1.0E200;
			for(i=0; i<pop.Size(); i++) if(flag[i])
			{
				obj = (proOAB[i] + tranOAB) + PUNISH*fabs(proAB[i]-proab);
				if(id==pop.Size() || obj<optobj) {id=i; optobj=obj;}
			}
			flag[id] = false;
		}

		return pop;
	}
//
//	// ====================================
//	// Oct.09, 2007
//	// Complicated Strategy: approximate the PF by several line segments
//	CPopulationMO& SRegF2::SelectLineSegment(CPopulationMO& pop, unsigned int size, std::vector<bool>& flag, bool fix)
//	{
//		unsigned int i,k, iDim = pop.P().FSize();				// dimension of the space (2: default)
//		unsigned iRef = 7;										// iRef reference points
//		std::vector< std::vector<double> >	vRef(iRef);			// reference points
//		for(i=0; i<iRef; i++) vRef[i].resize(iDim);
//
//		// Step 1: set initial reference points (two extreme points)
//		unsigned int iA, iB;
//		iA = iB = 0;
//		for(i=1; i<pop.Size(); i++)
//		{
//			if((0.95*pop[i].F(0)+0.05*pop[i].F(1)) < (0.95*pop[iA].F(0)+0.05*pop[iA].F(1)))	iA = i;
//			if((0.95*pop[i].F(1)+0.05*pop[i].F(0)) < (0.95*pop[iB].F(1)+0.05*pop[iB].F(0)))	iB = i;
//			//if(pop[i].F(0)<pop[iA].F(0))	iA = i;
//			//if(pop[i].F(1)<pop[iB].F(1))	iB = i;
//		}
//		while(iB == iA) iB = rnd::rand((unsigned int)0, pop.Size());
//		// set two extreme reference points
//		// !!! here we use the default dimension 2
//		vRef[0][0]		= pop[iA].F(0) - 0.1*fabs(pop[iB].F(0)-pop[iA].F(0));
//		vRef[0][1]		= pop[iA].F(1) + 0.1*fabs(pop[iB].F(1)-pop[iA].F(1));
//		vRef[iRef-1][0]	= pop[iB].F(0) + 0.1*fabs(pop[iB].F(0)-pop[iA].F(0));
//		vRef[iRef-1][1]	= pop[iB].F(1) - 0.1*fabs(pop[iB].F(1)-pop[iA].F(1));
//
//		//Step 2: set other reference points
//		//ChooseRef( (iRef-1)/2, 0, iRef-1, vRef, pop);
//		ChooseRef(vRef, pop);
//
//#ifdef AZ_TMPOUT
//	std::ofstream os("ModF.out");
//	os<<std::scientific<<std::setprecision(10);
//	os<<vRef[0][0]<<"\t"<<vRef[0][1]<<"\t"<<vRef[iRef-1][0]<<"\t"<<vRef[iRef-1][1]<<std::endl;
//	for(i=0; i<iRef-1; i++)
//	{
//		for(k=0; k<iDim; k++) os<<vRef[i][k]<<"\t";
//		for(k=0; k<iDim; k++) os<<vRef[i+1][k]<<"\t";
//		os<<std::endl;
//	}
//	os.close();
//#endif
//
//		//Step 4: select
//		Select(pop, flag, vRef, size, fix);
//
//		return pop;
//	}
//
//	//!\brief	set reference points
//	//!\param	iC current reference point
//	//!\param	iL left reference point
//	//!\param	iR right reference point
//	//!\param	vRef set of reference points
//	//!\param	pop population
//	//!\return	void
//	void SRegF2::ChooseRef(unsigned int iC, unsigned int iL, unsigned int iR, std::vector< std::vector<double> >& vRef, CPopulationMO& pop)
//	{
//		unsigned int i, k, iDim = (unsigned int)vRef[0].size();
//
//		if(vRef[iL][0] < vRef[iR][0] && vRef[iL][1] < vRef[iR][1] ||
//		   vRef[iR][0] < vRef[iL][0] && vRef[iR][1] < vRef[iL][1])
//		{
//			for(i=0; i<iDim; i++) vRef[iC][i] = 0.5*(vRef[iL][i]+vRef[iR][i]);
//		}
//		else
//		{
//			//direction
//			double LenR=0.0;
//			std::vector<double> LR(iDim); for(i=0; i<iDim; i++) {LR[i] = vRef[iR][i]-vRef[iL][i]; LenR+=LR[i]*LR[i];}
//			LenR = sqrt(LenR); for(i=0; i<iDim; i++) LR[i] /= LenR;
//			std::vector<double> OLR(iDim); OLR[0]=-LR[1]; OLR[1]=LR[0];
//
//			//reference center point
//			std::vector<double> C(iDim); for(i=0; i<iDim; i++) C[i] = 0.5*(vRef[iL][i]+vRef[iR][i])-LenR*OLR[i];
//
//			//find the reference point from the population
//			double prox, proy, dis, mindis=1.0E100, t=0.0;
//			std::vector<unsigned int> ipop(pop.Size()); for(i=0; i<pop.Size(); i++) ipop[i] = i;std::random_shuffle(ipop.begin(),ipop.end());
//			for(i=0; i<pop.Size(); i++)
//			{
//				prox = proy = 0.0;
//				for(k=0; k<iDim; k++) {prox+=(pop[ipop[i]].F(k)-C[k])*LR[k]; proy+=(pop[ipop[i]].F(k)-C[k])*OLR[k];}
//				dis = proy + PUNISH*fabs(prox);
//				if(dis<mindis) {mindis=dis; t=proy;}
//			}
//			for(i=0; i<iDim; i++) vRef[iC][i] = C[i] + t*OLR[i];
//		}
//		//find the other reference points
//		if(iC>iL+1)
//		{
//			ChooseRef( (iC+iL)/2, iL, iC, vRef, pop);	// left
//			ChooseRef( (iC+iR)/2, iC, iR, vRef, pop);	// right
//		}
//	}
//
//	//!\brief	set reference points
//	//!\param	vRef set of reference points
//	//!\param	pop population
//	//!\return	void
//	void SRegF2::ChooseRef(std::vector< std::vector<double> >& vRef, CPopulationMO& pop)
//	{
//		unsigned int i, k, n, iDim = (unsigned int)vRef[0].size(), iRef = (unsigned int)vRef.size();
//		
//		//direction
//		double LenR=0.0;
//		std::vector<double> LR(iDim); for(i=0; i<iDim; i++) {LR[i] = vRef[iRef-1][i]-vRef[0][i]; LenR+=LR[i]*LR[i];}
//		LenR = sqrt(LenR); for(i=0; i<iDim; i++) LR[i] /= LenR;
//		std::vector<double> OLR(iDim); OLR[0]=-LR[1]; OLR[1]=LR[0];
//		
//		//projections along OLR (y) and LR (x)
//		std::vector<double> prox(pop.Size()), proy(pop.Size());
//		double miny = 0.0, absy = -1.0;
//		for(i=0; i<pop.Size(); i++)
//		{
//			prox[i] = proy[i] = 0.0;
//			for(k=0; k<iDim; k++) {prox[i]+=(pop[i].F(k)-vRef[0][k])*LR[k]; proy[i]+=(pop[i].F(k)-vRef[0][k])*OLR[k];}
//			if(proy[i]<miny)		miny = proy[i];
//		}
//
//		//select the reference point one by one
//		for(n=1; n<iRef-1; n++)
//		{
//			double x = LenR*double(n)/double(iRef-1);
//
//			double mindis = 1.0E100; unsigned int mini = 0; double dis;
//			for(i=0; i<pop.Size(); i++)
//			{
//				dis = (prox[i]-x)*(prox[i]-x)+(miny-proy[i])*(miny-proy[i]);
//				if(dis<mindis){mindis=dis; mini=i;}
//			}
//			for(i=0; i<iDim; i++) vRef[n][i] = vRef[0][i]+x*LR[i]+proy[mini]*OLR[i];
//			if(fabs(proy[mini])>absy)	absy = fabs(proy[mini]);
//		}
//
//		//transformation algong OLR
//		for(n=0; n<iRef; n++) for(i=0; i<iDim; i++) vRef[n][i] -= absy*OLR[i];
//	}
//
//	// select points to reference points
//	void SRegF2::Select(CPopulationMO& pop, std::vector<bool>& vExist, std::vector< std::vector<double> >& vRef, unsigned int size, bool fix)
//	{
//		unsigned int i, k, index, iRef=(unsigned int)vRef.size(), iDim=pop.P().FSize();
//
//		vExist.resize(pop.Size()); for(i=0; i<pop.Size(); i++) vExist[i]=true;
//
//		//calculate the length of each segment
//		std::vector<double> vLenS(iRef-1); double LenS=0.0, DetS=0.0;
//		for(i=0; i<iRef-1; i++) {vLenS[i] = 0.0; for(k=0; k<iDim; k++) vLenS[i] += (vRef[i+1][k]-vRef[i][k])*(vRef[i+1][k]-vRef[i][k]); vLenS[i] = sqrt(vLenS[i]); LenS += vLenS[i];}
//		DetS = LenS/double(size-1.0);
//
//		//directions
//		std::vector<double> C(iDim),LR(iDim),OLR(iDim);
//
//		double lenL=-vLenS[0], lenR=0.0, lenC; unsigned int iS=0;
//		double prox, proy, dis, mindis = 1.0E200; unsigned int id;
//
//		for(index=0; index<size; index++)
//		{
//			if(fix)
//				lenC = double(index)*DetS;
//			else
//				lenC = rnd::rand(std::max(0.0,index-0.35), std::min(index+0.35,size-1.0))*DetS;
//			//recalculate the directions
//			if(lenC>=lenR && iS<iRef-1)
//			{
//				for(k=0; k<iDim; k++) LR[k] = (vRef[iS+1][k]-vRef[iS][k])/vLenS[iS];
//				OLR[0]=-LR[1]; OLR[1]=LR[0];
//				lenL += vLenS[iS];
//				lenR += vLenS[iS++];
//			}
//
//			for(k=0; k<iDim; k++) C[k] = vRef[iS-1][k] + (lenC-lenL)*LR[k];
//
//			id = pop.Size(); mindis = 1.0E200;
//			for(i=0; i<pop.Size(); i++) if(vExist[i])
//			{
//				prox= proy = 0.0;
//				for(k=0; k<iDim; k++) {prox+=(pop[i].F(k)-C[k])*LR[k]; proy+=(pop[i].F(k)-C[k])*OLR[k];}
//				//dis	= proy + PUNISH*prox*prox/fabs(proy);
//				dis = proy + PUNISH*fabs(prox);
//				if(id==pop.Size() || dis<mindis) {id=i; mindis=dis;}
//			}
//			vExist[id] = false;
//		}
//	}
}//namespace sel
} //namespace mea
} //namespace az
