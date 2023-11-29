// LocalPCAO.cpp

#include <cmath>
#include <fstream>
#include "alg/Random.h"
#include "alg/LocalPCAO.h"

namespace az
{

namespace alg
{

//current active cluster
CLUSTER* pCluster=0;

//objective function for Quaci-Newton method
double F(double* r)
{
	unsigned int i,k;
	unsigned int dim = (*pCluster).mData.RowSize();
	double ff = 0.0;
	alg::FVECTOR tmpx(dim),tmp;
	for(k=0; k<(*pCluster).mData.ColSize(); k++)
	{
		for(i=0; i<dim; i++) tmpx[i] = (*pCluster).mData(i,k)-r[i];
		(*pCluster).mPi.LeftMultiply(tmpx, tmp);
		for(i=0; i<dim; i++) ff += tmp[i]*tmpx[i];
	}
	return ff;
}

//gradient function for Quaci-Newton method
void G(double* px, double* pg)
{
	unsigned int i, j, k;
	unsigned int dim = (*pCluster).mData.RowSize();
	alg::FVECTOR tmpx(dim),tmpg(dim),tmp;
	for(i=0; i<dim; i++) pg[i] = 0.0;
	for(k=0; k<(*pCluster).mData.ColSize(); k++)
	{
		for(i=0; i<dim; i++) tmpx[i] = (*pCluster).mData(i,k)-px[i];
		for(i=0; i<dim; i++)
			for(j=0; j<dim; j++)
				pg[i] += -tmpx[j]*((*pCluster).mPi(i, j) + (*pCluster).mPi(j, i));
	}
}

//construct function
LocalPCAO::LocalPCAO()
{
	mDataSize	= 0;
	mDataDim	= 0;
	mLatentDim	= 0;
	mClusterNo	= 0;
	pData		= 0;
	mbReset		= true;
}

//destruct function
LocalPCAO::~LocalPCAO()
{
	unsigned int i;
	for(i=0; i<mClusters.size(); i++) delete mClusters[i];
	mClusters.clear();
	mMinDis.clear();
}

//initialize 
void LocalPCAO::Initialize(unsigned int datadim, unsigned int latentdim, unsigned int maxcluster)
{
	unsigned int i;
	mDataDim	= datadim;
	mLatentDim	= latentdim;
	mbReset		= true;
	mUpp.resize(datadim);
	mLow.resize(datadim);
	for(i=0; i<datadim; i++)
	{
		mUpp[i] = 2000.0;
		mLow[i] = -2000.0;
	}
	if(mClusters.size() != maxcluster)
	{
		for(i=0; i<mClusters.size(); i++) delete mClusters[i];
		mClusters.resize(maxcluster);
		for(i=0; i<maxcluster; i++) mClusters[i] = new CLUSTER;
	}
}

//train
void LocalPCAO::Train(unsigned int steps, unsigned int size, double **pdata)
{
	unsigned int i,j,k,iter,clu,unchanged;
	//double tmp;
	alg::Matrix tmpPI1, tmpPI2;

	mDataSize	= size;
	pData		= pdata;

	//for(i=0; i<mDataDim; i++) {mUpp[i] = mLow[i] = pData[0][i];}
	//for(k=1; k<mDataSize; k++)
	//{
	//	for(i=0; i<mDataDim; i++) 
	//	{
	//		if(pData[k][i] > mUpp[i] ) mUpp[i] = pData[k][i];
	//		if(pData[k][i] < mLow[i] ) mLow[i] = pData[k][i];
	//	}
	//}
	//for(i=0; i<mDataDim; i++) 
	//{
	//	tmp = mUpp[i] - mLow[i];
	//	mUpp[i] += 0.1*tmp;
	//	mLow[i] -= 0.1*tmp;
	//}

	mMinDis.resize(mDataSize);

	//Step1: reset the reference vector
	if(mbReset)
	{
		mClusterNo	= (unsigned int)mClusters.size();
		mbReset		= false;

		for(i=0; i<mClusterNo; i++)
		{
			(*mClusters[i]).mPi.Identity(mDataDim);
			(*mClusters[i]).mPCA.Initialize(mDataDim);
		}
		
		alg::VINDEX tmpindex(mDataSize);
		for(j=0; j<mDataSize; j++) tmpindex[j]=j;
		for(i=0; i<mClusterNo; i++)
		{
			std::swap(tmpindex[i], tmpindex[rnd::rand(i, mDataSize)]);
			(*mClusters[i]).mCore.resize(mDataDim);
			for(j=0; j<mDataDim; j++) (*mClusters[i]).mCore[j] = pData[tmpindex[i]][j];
		}
	}

	//Step2: train iterations
	iter = 0;
	while(iter++<steps)
	{
		//Step 2.1: partition
		Partition();

		//Step 2.2: update reference vectors
		unchanged = 0;
		//Calcluate the mean, covariance matrix, eigenvalue and eigenvector
		for(clu=0; clu<mClusterNo; clu++)
		{
			//train pca
			(*mClusters[clu]).mPCA.Train();
			//calculate matrix PI
			if(mDataDim > mLatentDim)
			{
				tmpPI1.Resize(mDataDim, mDataDim - mLatentDim);
				for(j=mLatentDim; j<mDataDim; j++)
					for(k=0; k<mDataDim; k++)
						tmpPI1(k, j-mLatentDim) = (*mClusters[clu]).mPCA.Eigenvector()(k, j);
				tmpPI2 = tmpPI1;
				tmpPI2.Trans();
				tmpPI1.Multiply(tmpPI2, (*mClusters[clu]).mPi);
			}
			//Update reference vectors
			unchanged += (mClusterNo==1 || mDataDim==mLatentDim || !OptimizeRefe(clu)) ? 1:0;
		}
		//no reference vectors can be updated
		if(unchanged>=mClusterNo) break;
	}

//std::ofstream out("time.txt",std::ios::app);
//out<<iter<<"\t";
//out.close();

	CalDistance();
}

//partition data according to current reference vectors
void LocalPCAO::Partition()
{
	unsigned int cluster=0,index,clu,i;
	double dis,dismin;

	//Step1: clear current index list
	for(clu=0; clu<mClusterNo; clu++) (*mClusters[clu]).mIndex.clear();

	//Step2: find a closest cluster for each point
	for(index=0; index<mDataSize; index++)
	{
		dismin=1.0E100;
		for(clu=0; clu<mClusterNo; clu++)
		{
			dis = Distance(index,clu);
			if(dis<dismin) { dismin=dis; cluster=clu; }
		}
		(*mClusters[cluster]).mIndex.push_back(index);
	}

	//Step3: delete the empty clusters
	clu = 0;
	while(clu<mClusterNo)
	{
		while(clu<mClusterNo && (*mClusters[clu]).mIndex.size()>0) clu++;
		if(clu<mClusterNo)
		{
			mClusterNo--;
			std::swap(mClusters[clu], mClusters[mClusterNo]);
		}
	}

	//Step 4: assign data to pca
	for(clu=0; clu<mClusterNo; clu++)
	{
		(*mClusters[clu]).mData.Resize(mDataDim, (unsigned int)(*mClusters[clu]).mIndex.size());
		alg::LINDEX::iterator it = (*mClusters[clu]).mIndex.begin();
		index = 0;
		while(it != (*mClusters[clu]).mIndex.end())
		{
			for(i=0; i<mDataDim; i++) (*mClusters[clu]).mData(i,index) = pData[*it][i];
			index++;
			it++;
		}
		(*mClusters[clu]).mPCA.Data((*mClusters[clu]).mData);
	}
}

//calculate the distance from a point to a cluster
double LocalPCAO::Distance(unsigned int data, unsigned int clu)
{
	return Distance_Recdis(data, clu);
}

double LocalPCAO::Distance_Eucdis(unsigned int data, unsigned int clu)
{
	unsigned int i;
	double tmp = 0.0;
	for(i=0; i<mDataDim; i++)
		tmp += (pData[data][i] - (*mClusters[clu]).mCore[i])*(pData[data][i] - (*mClusters[clu]).mCore[i]);
	return tmp;
}

double LocalPCAO::Distance_Recdis(unsigned int data, unsigned int clu)
{
	unsigned int i;
	double tmp = 0.0;
	alg::FVECTOR tmpVec1, tmpVec2;
	tmpVec1.resize(mDataDim);
	for(i=0; i<mDataDim; i++) tmpVec1[i] = (pData[data][i] - (*mClusters[clu]).mCore[i]);
	(*mClusters[clu]).mPi.LeftMultiply(tmpVec1, tmpVec2);
	for(i=0; i<mDataDim; i++) tmp += tmpVec1[i]*tmpVec2[i];
	return tmp;
}

//update reference vector
bool LocalPCAO::Update(unsigned int clu)
{
	return OptimizeRefe(clu);
}

bool LocalPCAO::OptimizeRefe(unsigned int clu)
{
	unsigned int i;
	alg::FVECTOR tmpr	= (*mClusters[clu]).mCore;
	alg::FVECTOR tmpx	= (*mClusters[clu]).mPCA.Mean();

	double* px	= &(*(tmpx.begin()));
	double* plow= &(*(mLow.begin()));
	double* pup	= &(*(mUpp.begin()));

	pCluster = mClusters[clu];
	alg::QNewton mini(mDataDim, F, G, px, plow, pup);
	while(mini.Iteration() < 500 && mini.Evaluations() < 2000 && mini.Step());

	double error = 0.0;
	double*const ppx = mini.X();
	for(i=0; i<(unsigned int)tmpr.size(); i++)
	{
		tmpr[i] = ppx[i];
		error += (tmpr[i]-(*mClusters[clu]).mCore[i])*(tmpr[i]-(*mClusters[clu]).mCore[i]);
	}

	//tmpr = (*mClusters[clu]).mPCA.Mean();
	//double error = 0.0;
	//for(i=0; i<(unsigned int)tmpr.size(); i++)	error += (tmpr[i]-(*mClusters[clu]).mCore[i])*(tmpr[i]-(*mClusters[clu]).mCore[i]);
	//error = sqrt(error);

	(*mClusters[clu]).mCore = tmpr;

	return error>1.0e-5;
}

bool LocalPCAO::RebuildCenter(unsigned int clu)
{
	unsigned int i;
	alg::FVECTOR tmp((*mClusters[clu]).mCore.size());
	for(i=0; i<(*mClusters[clu]).mCore.size(); i++)
	{
		tmp[i] = 0.0;
		alg::LINDEX::iterator it = (*mClusters[clu]).mIndex.begin();
		while(it != (*mClusters[clu]).mIndex.end())
		{
			tmp[i] += pData[*it][i];
			it++;
		}
		tmp[i] /= double((*mClusters[clu]).mIndex.size());
	}

	double error = 0.0;
	for(i=0; i<(unsigned int)tmp.size(); i++)
	{
		error += (tmp[i]-(*mClusters[clu]).mCore[i])*(tmp[i]-(*mClusters[clu]).mCore[i]);
	}
	(*mClusters[clu]).mCore = tmp;

	return error>1.0e-5;
}

void LocalPCAO::CalDistance()
{
	unsigned int clu;
	alg::LINDEX::iterator it;
	for(clu=0; clu<mMinDis.size(); clu++) mMinDis[clu] = 0.0;
	for(clu=0; clu<mClusterNo; clu++)
	{
		it = (*mClusters[clu]).mIndex.begin();
		while(it != (*mClusters[clu]).mIndex.end())
		{
			mMinDis[*it] = sqrt(Distance(*it, clu));
			it++;
		}
	}
}

} //namespace alg

} //namespace az
