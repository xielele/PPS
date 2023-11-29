/*! \file	GenModRM_O.cpp
	
	\brief	Evolutionary Aglorithm Generator with Local PCA
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Nov.29 2005 make great changes: noise, border checking
	\date	Apr.10 2006 redesign
	\date	Jul.18 2006 add quadratic models
	\date	Sep.03 2007 rewrite to check the difference between RM and RM_OLD
*/
#include <ctime>
#include <list>
#include <vector>
#include <cmath>
#include <float.h>
#include "alg/HCSampler.h"
#include "emo/LogFile.h"
#include "emo/GenMod.h"

//#define SD_DIS

#ifdef AZ_MODEL_OUT
#include <fstream>
#endif

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
RM_OLD::RM_OLD()
{
	mLatentDim	= 0;
	mMaxCluster	= 0;
}

//initialize the LPCA
void RM_OLD::Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension)
{
	mLatentDim	= latent;
	mMaxCluster	= cluster;
	mTrainSteps	= trainsteps;
	mExtension	= extension;
}
	
//model-based generator
CPopulationMO& RM_OLD::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int c,i,j,k,noNew;

	//Step 1: clear the return population
	popnew.Resize(sizenew);
	noNew = popnew.Size();
	if(noNew<1) return popnew;
	
	//Step 2: assign new data
	if( popref.Size() != mDataSize ) 
	{
		Clear();	
		mDataSize	= popref.Size();
		mDataDim	= popref.P().XSize();
		pData		= new double*[mDataSize];
		for( i=0; i<mDataSize; i++ ) pData[i] = new double[mDataDim];
	}
	for( i=0; i<mDataSize; i++ ) for( j=0; j<mDataDim; j++ ) pData[i][j] = popref[i][j];
	

	//Step 3: train Local PCA model
	mLocalPCA.Initialize( mDataDim, mLatentDim, mMaxCluster );
	mLocalPCA.Train( mTrainSteps, mDataSize, pData );

	//Step 4: calculate the probability of each cluster, i.e the size to create in each cluster
	// 2 at least in each cluster
	unsigned int nt = 0;
	std::vector<unsigned int> nc(mMaxCluster);
	for(i=0; i<mMaxCluster; i++) { nc[i] = (mLocalPCA.DataIndex(i).size() >= 1) ? 2:0; nt += nc[i];}
	double vt = 0.0;
	std::vector<double> vc(mMaxCluster);
	std::vector< std::vector<double> > pmin,pmax;
	pmin.resize(mMaxCluster);pmax.resize(mMaxCluster);
	for(i=0; i<mMaxCluster; i++)
	{
		pmin[i].resize(mLatentDim);pmax[i].resize(mLatentDim);
		for(j=0; j<mLatentDim; j++) pmin[i][j] = pmax[i][j] = 0.0;
		if(mLocalPCA.DataIndex(i).size() > 1)
		{
			std::vector<unsigned int> dataindex; double p;
			dataindex.resize(mLocalPCA.DataIndex(i).size());
			std::copy(mLocalPCA.DataIndex(i).begin(), mLocalPCA.DataIndex(i).end(),dataindex.begin());
			for(k=0; k<mLatentDim; k++)
			{
				pmin[i][k] = 1.0E100; pmax[i][k] =-1.0E100;
				for(j=0; j<dataindex.size(); j++)
				{
					p = Project(k,dataindex[j],i);	
					if(p<pmin[i][k]) pmin[i][k] = p;
					if(p>pmax[i][k]) pmax[i][k] = p;
				}
			}
		}
	}
	for(i=0; i<mMaxCluster; i++) 
	{ 
		vc[i] = 0.0;
		if(mLocalPCA.DataIndex(i).size() >  1)
		{
			vc[i] = 1.0;
			for(j=0; j<mLatentDim; j++) vc[i] *= pmax[i][j] - pmin[i][j];
 		}
		vt += vc[i];
	}
	double ns = noNew - nt + 0.0;
	for(i=0; i<mMaxCluster; i++) {c = (unsigned int)(ns*vc[i]/vt); nc[i] += c; nt += c;}
	while(nt<noNew)
	{
		i = rnd::rand((unsigned int)0, (unsigned int)mMaxCluster);
		if(mLocalPCA.DataIndex(i).size()>1) {nc[i]++; nt++;}
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
		if(mLocalPCA.DataIndex(c).size() > 0)
		{
			dataindex.resize(mLocalPCA.DataIndex(c).size());
			std::copy(mLocalPCA.DataIndex(c).begin(), mLocalPCA.DataIndex(c).end(),dataindex.begin());
		}

		// only one point in the cluster
		if(mLocalPCA.DataIndex(c).size() == 1)
		{
			sd  = 0.0;
			for(i=0; i<mDataDim; i++) sd += 0.1*(popref.P().XUpp(i)-popref.P().XLow(i));
			sd /= mDataDim;
			for(i=0; i<nc[c]; i++) 
			{
				for( j=0; j<mDataDim; j++ ) 
				{
					popnew[np+i][j] = mLocalPCA.Mean(c)[j] + sd*rnd::gaussian();

					if( wxFinite(popnew[np+i][j]) == 0 )			popnew[np+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
					else if(popnew[np+i][j] < popnew.P().XLow(j))	popnew[np+i][j] = 0.5*(popnew.P().XLow(j)+popref[dataindex[0]][j]);
					else if(popnew[np+i][j] > popnew.P().XUpp(j))	popnew[np+i][j] = 0.5*(popnew.P().XUpp(j)+popref[dataindex[0]][j]);
				}
			}
		}
		// more than one point in the cluster
		else if(mLocalPCA.DataIndex(c).size() > 1)
		{
			std::vector< std::vector<double> > t(mLatentDim);
			for(i=0; i<mLatentDim; i++) t[i].resize(nc[c]);
			std::vector<double> low(mLatentDim),upp(mLatentDim);
			for(i=0; i<mLatentDim; i++) 
			{
				low[i] = pmin[c][i]-mExtension*(pmax[c][i]-pmin[c][i]);
				upp[i] = pmax[c][i]+mExtension*(pmax[c][i]-pmin[c][i]);
			}
			alg::LHC(t, low, upp);

			sd = 0.0;
			for(i=mLatentDim; i<mDataDim; i++) sd += fabs(mLocalPCA.Eigenvalue(c)[i]);
			sd = sqrt(sd / double(mDataDim-mLatentDim));
			
			for(i=0; i<nc[c]; i++) for( j=0; j<mDataDim; j++ ) 
			{
				popnew[np+i][j] = mLocalPCA.Mean(c)[j] + sd*rnd::gaussian();
				for(k=0; k<mLatentDim; k++) popnew[np+i][j] += t[k][i]*mLocalPCA.Eigenvector(c)(j,k);

				if( wxFinite(popnew[np+i][j]) == 0 )			popnew[np+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
				else if(popnew[np+i][j] < popnew.P().XLow(j))	popnew[np+i][j] = 0.5*(popnew.P().XLow(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
				else if(popnew[np+i][j] > popnew.P().XUpp(j))	popnew[np+i][j] = 0.5*(popnew.P().XUpp(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
			}
		}
		np += nc[c];
#ifdef AZ_MODEL_OUT
		fhand<<pmin[c][0] * mLocalPCA.Eigenvector(c)(0,0) + mLocalPCA.Mean(c)[0] <<"\t"
			 <<pmin[c][0] * mLocalPCA.Eigenvector(c)(1,0) + mLocalPCA.Mean(c)[1] <<"\t"
			 <<pmax[c][0] * mLocalPCA.Eigenvector(c)(0,0) + mLocalPCA.Mean(c)[0] <<"\t"
			 <<pmax[c][0] * mLocalPCA.Eigenvector(c)(1,0) + mLocalPCA.Mean(c)[1] <<"\t"
			 <<sd<<std::endl;
#endif
	}
#ifdef AZ_MODEL_OUT
	fhand.close();
#endif	
	return popnew;
}

//project point into a dimension
double RM_OLD::Project(unsigned int dim, unsigned int index, unsigned int clu)
{
	unsigned int i; double len;
	len = 0.0;
	//projection
	for( i=0; i<mDataDim; i++ ) len += (pData[index][i]-mLocalPCA.Mean(clu)[i]) * mLocalPCA.Eigenvector(clu)(i,dim);
	return len;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az
