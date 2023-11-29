// Som.cpp

#include <algorithm>
#include <cmath>
#include "alg/Random.h"
#include "alg/SOM.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

CBatchSOM::CBatchSOM()
{
	mRow	= 0;
	mCol	= 0;
	mDim	= 0;
	mDataNo= 0;
	pData	= 0;
	pWeight= 0;
	pNear	= 0;
	pDis	= 0;
	
}

CBatchSOM::~CBatchSOM()
{
	Clear();
}

void CBatchSOM::Initialize(int row, int col, int dim)
{
	int i;

	mThetaMin	= 0.5; 
	mThetaMax	= 4.0;

	if(mRow != row || mCol != col || mDim != dim)
	{
		//delete the old space if possible
		Clear();

		//create new space
		mRow		= row;
		mCol		= col;
		mDim		= dim;
		mWeightNo	= mRow*mCol;
		pWeight		= new double*[mWeightNo];
		for(i=0; i<mWeightNo; i++) pWeight[i] = new double[mDim];
	}

	bInit = true;
}


void CBatchSOM::Train(int steps, int size, double **pdata)
{
	int i, j, step;

	if(size != mDataNo && pNear != 0) 
	{
		delete []pNear; pNear	= 0;
		delete []pDis;	pDis	= 0;
	}
	if(pNear == 0)pNear = new int[size];
	if(pDis == 0)	pDis  = new double[size];	
	mDataNo = size;
	pData	= pdata;

	//need to be initizalized
	if(bInit)
	{
		for(i=0; i<mDataNo; i++)	 pNear[i] = i;
		for(i=0; i<mWeightNo; i++) std::swap(pNear[i], pNear[i+rnd::rand(0,mWeightNo-i)]);
		for(i=0; i<mWeightNo; i++) for(j=0; j<mDim; j++) pWeight[i][j] = pData[pNear[i]][j];
		bInit = false;
	}

	for(step=0; step<steps; step++)
	{
		mTheta = mThetaMax - (mThetaMax - mThetaMin)*double(step) / double(steps);

		//step 1: for all data, find appropriate nearest weight index 
		for(i=0; i<mDataNo; i++) NearWeight(i);

		//step 2: update all weights
		for(i=0; i<mWeightNo; i++) UpdateWeight(i);
	}
}

void CBatchSOM::Clear()
{
	int i;
	if(pNear != 0) 
	{
		delete []pNear;
		pNear = 0;
	}
	
	if(pDis != 0)
	{
		delete []pDis;
		pDis = 0;
	}

	if(pWeight != 0)
	{
		for(i=0; i<mRow*mCol; i++) delete []pWeight[i];
		delete []pWeight;
		pWeight = 0;
	}
}

void CBatchSOM::NearWeight(int k)
{
	int i,j,min=0;
	double d, dis, minDis=1.0E20;
	for(i=0; i<mWeightNo; i++)
	{
		dis = 0.0;
		for(j=0; j<mDim; j++) 
		{
			d	= pWeight[i][j] - pData[k][j];
			dis+= d*d;
		}
		if(dis < minDis)
		{
			minDis	= dis;
			min		= i;
		}
	}
	pNear[k] = min;
	pDis[k]  = sqrt(minDis);
}

void CBatchSOM::UpdateWeight(int k)
{
	int d, i;
	double dis, dist=0;
	for(d=0; d<mDim; d++) pWeight[k][d] = 0.0;
	for(i=0; i<mDataNo; i++)
	{
		dis	 = Neighbourhood(k, pNear[i]);
		dist+= dis;
		for(d=0; d<mDim; d++)	pWeight[k][d] += dis*pData[i][d];	
	}
	for(d=0; d<mDim; d++)	pWeight[k][d] /= dist;
}

double CBatchSOM::Neighbourhood(int k, int n)
{
	int kr = k / mCol; int kc = k % mCol;
	int nr = n / mCol; int nc = n % mCol;
	double dis2 = (kr - nr)*(kr - nr) + (kc - nc)*(kc - nc);
	return exp(- dis2 / (mTheta*mTheta));
}

} //namespace alg

} //namespace az

