/*! \file	GenModRM.cpp
	
	\brief	Evolutionary Aglorithm Generator with Local PCA (RM-MEDA)
	\brief  "RM-MEDA: A Regularity Model-Based Multiobjective Estimation of Distribution Algorithm". IEEE Transaction on Evolutionary Computation
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Nov.29 2005 make great changes: noise, border checking
	\date	Apr.10 2006 redesign
	\date	Jul.18 2006 add quadratic models
	\date	Nov.12 2006 modify to uniform version
	\date	Jun.26 2006 rename and change Generate()
	\date	Sep.03 2007 modify the boundary checking procedure, it plays an important role in the algorithm
*/
#include <ctime>
#include <list>
#include <vector>
#include <cmath>
#include <float.h>
#include "alg/Fitting.h"
#include "alg/LocalPCA.h"
#include "alg/Kmeans.h"
#include "alg/HCSampler.h"
#include "emo/GenMod.h"

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
namespace mod
{

//Local PCA based EDA generator
//constructor
RM::RM()
{
	mLatentDim	= 0;
	mMaxCluster	= 0;
}

//initialize the LPCA
void RM::Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension)
{
	mLatentDim	= latent;
	mMaxCluster	= cluster;
	mTrainSteps	= trainsteps;
	mExtension	= extension;
}

// model-based generator
CPopulationMO& RM::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int c,i,j,k,noNew;

	popnew.Resize(sizenew);

	noNew = popnew.Size();
	if(noNew<1) return popnew;
	
	//Step 1: assign new data
	if( popref.Size() != mDataSize ) 
	{
		Clear();	
		mDataSize	= popref.Size();
		mDataDim	= popref.P().XSize();
		pData		= new double*[mDataSize];
		for( i=0; i<mDataSize; i++ ) pData[i] = new double[mDataDim];
	}
	for( i=0; i<mDataSize; i++ ) for( j=0; j<mDataDim; j++ ) pData[i][j] = popref[i][j];
	
	//Step 2: train with Local PCA
	alg::LocalPCA lpca;
	//alg::Kmeans lpca;
	lpca.Set(mMaxCluster, mDataSize, mDataDim, mLatentDim, mTrainSteps);
	for( i=0; i<mDataSize; i++ ) for( j=0; j<mDataDim; j++ ) lpca.mvX[i][j] = popref[i][j];
	lpca.Train();

	//Step 3: calculate the probability of each cluster, i.e the size to create in each cluster
	// 2 at least in each cluster
	unsigned int nt = 0;
	std::vector<unsigned int> nc(mMaxCluster);
	for(i=0; i<mMaxCluster; i++) { nc[i] = (lpca.mvNo[i] >= 1) ? 2:0; nt += nc[i];}
	double vt = 0.0;
	std::vector<double> vc(mMaxCluster);
	for(i=0; i<mMaxCluster; i++) 
	{ 
		vc[i] = 0.0;
		if(lpca.mvNo[i] >  1)
		{
			vc[i] = 1.0;
			for(j=0; j<mLatentDim; j++) vc[i] *= lpca.mvProMax[i][j] - lpca.mvProMin[i][j];
 		}
		vt += vc[i];
	}
	double ns = noNew - nt + 0.0;
	for(i=0; i<mMaxCluster; i++) { c = (unsigned int)(ns*vc[i]/vt); nc[i] += c; nt += c;}
	while(nt<noNew)
	{
		i = rnd::rand((unsigned int)0, (unsigned int)mMaxCluster);
		if(lpca.mvNo[i]>1) {nc[i]++; nt++;}
	}

#ifdef AZ_MODEL_OUT
	std::ofstream fhand("model.set");
	fhand<<"PCA"<<std::endl;
	fhand<<mMaxCluster<<std::endl;
#endif
	//Step 4: create new trial solutions in each cluster
	unsigned int np=0; 
	double sd;
	for(c=0; c<mMaxCluster; c++)
	{
		std::vector<unsigned int> dataindex;
		if(lpca.mvNo[c] > 0)
		{
			dataindex.resize(lpca.mvNo[c]);j=0;
			for(i=0; i<(unsigned int)lpca.mvIndex.size(); i++) if(lpca.mvIndex[i] == c) dataindex[j++] = i;
		}

		// only one point in the cluster
		if(lpca.mvNo[c] == 1)
		{
			sd  = 0.0;
			for(i=0; i<mDataDim; i++) sd += 0.1*(popref.P().XUpp(i)-popref.P().XLow(i));
			sd /= mDataDim;
			for(i=0; i<nc[c]; i++) 
			{
				for( j=0; j<mDataDim; j++ ) 
				{
					popnew[np+i][j] = lpca.mvMean[c][j] + sd*rnd::gaussian();

					if( wxFinite(popnew[np+i][j]) == 0 )			popnew[np+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
					else if(popnew[np+i][j] < popnew.P().XLow(j))	popnew[np+i][j] = 0.5*(popnew.P().XLow(j)+popref[dataindex[0]][j]);
					else if(popnew[np+i][j] > popnew.P().XUpp(j))	popnew[np+i][j] = 0.5*(popnew.P().XUpp(j)+popref[dataindex[0]][j]);
				}
			}
		}
		// more than one point in the cluster
		else if(lpca.mvNo[c] > 1)
		{
			std::vector< std::vector<double> > t(mLatentDim);
			for(i=0; i<mLatentDim; i++) t[i].resize(nc[c]);
			std::vector<double> low(mLatentDim),upp(mLatentDim);
			for(i=0; i<mLatentDim; i++) 
			{
				low[i] = lpca.mvProMin[c][i]-mExtension*(lpca.mvProMax[c][i]-lpca.mvProMin[c][i]);
				upp[i] = lpca.mvProMax[c][i]+mExtension*(lpca.mvProMax[c][i]-lpca.mvProMin[c][i]);
			}
			alg::LHC(t, low, upp);

			sd = 0.0;
			for(i=mLatentDim; i<mDataDim; i++) sd += fabs(lpca.mvEigenvalue[c][i]);
			sd = sqrt(sd / double(mDataDim-mLatentDim));
			//if(sd<0.05*fabs(lpca.mvEigenvalue[c][mLatentDim-1])) sd=0.05*fabs(lpca.mvEigenvalue[c][mLatentDim-1]);
			
			for(i=0; i<nc[c]; i++) for( j=0; j<mDataDim; j++ ) 
			{
				popnew[np+i][j] = lpca.mvMean[c][j] + sd*rnd::gaussian();
				for(k=0; k<mLatentDim; k++) popnew[np+i][j] += t[k][i]*lpca.mvEigenvector[c][k][j];

				if( wxFinite(popnew[np+i][j]) == 0 )			popnew[np+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
				else if(popnew[np+i][j] < popnew.P().XLow(j))	popnew[np+i][j] = 0.5*(popnew.P().XLow(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
				else if(popnew[np+i][j] > popnew.P().XUpp(j))	popnew[np+i][j] = 0.5*(popnew.P().XUpp(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
			}
		}
		np += nc[c];
#ifdef AZ_MODEL_OUT
		fhand<<lpca.mvProMin[c][0] * lpca.mvEigenvector[c][0][0] + lpca.mvMean[c][0] <<"\t"
			 <<lpca.mvProMin[c][0] * lpca.mvEigenvector[c][0][1] + lpca.mvMean[c][1] <<"\t"
			 <<lpca.mvProMax[c][0] * lpca.mvEigenvector[c][0][0] + lpca.mvMean[c][0] <<"\t"
			 <<lpca.mvProMax[c][0] * lpca.mvEigenvector[c][0][1] + lpca.mvMean[c][1] <<"\t"
			 <<sd<<std::endl;
#endif
	}

#ifdef AZ_MODEL_OUT
	fhand.close();
#endif	

	//// guide xover
	//GuidedXOver xover;
	//xover.XOver(popnew, popref);
	////for(i=0; i<popnew.Size(); i++) PM(popnew[i]);

	return popnew;
}

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az
