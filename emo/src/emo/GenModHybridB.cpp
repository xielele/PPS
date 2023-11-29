/*! \file	Generator_Model_HybridB.cpp
	
	\brief	Evolutionary Aglorithm Generator with Local PCA and SBX
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date   Nov.30 2005 change convergence test & new offspring size in each cluster
	\date	Apr.10 2006 redesign
*/
//#define SD_DIS

#include <list>
#include <vector>
#include <cmath>
#include <fstream>
#include "alg/Fitting.h"
#include "emo/GenMod.h"

#define TEST_CENTER 1

#ifdef TEST_CENTER
#include <fstream>
#endif

namespace az
{
namespace mea
{
namespace gen
{
namespace mod
{

//EDA+GA hybrid algoirithm B, EAD and GA runs in a single step
//constructor
ModelHybridB::ModelHybridB()
{
	mLatentDim	= 0;
	mMaxCluster	= 0;
	mbReset		= true;
}

//initialize the LocalPCA
void ModelHybridB::Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension, double rho)
{
	mLatentDim	= latent;
	mMaxCluster	= cluster;
	mTrainSteps	= trainsteps;
	mExtension	= extension;
	mRho		= rho;
	mbReset		= true;
}
	
//reset Local PCA
void ModelHybridB::Reset()
{
	mbReset = true;
}

//modle-based generator
CPopulationMO& ModelHybridB::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j;

	Reset();
	
	//Step 1: clear the return population
	popnew.Resize(sizenew);
	
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
	
	//Step 3: reset Local PCA if necessary
	if( mbReset )
	{
		mLocalPCA.Initialize( mDataDim, mLatentDim, mMaxCluster );
		mbReset = false;
	}

	//Step 4: train Local PCA model
	mLocalPCA.Train( mTrainSteps, mDataSize, pData );

	//Step 5: generate new solutions
	//1D
	if(mLatentDim == 1)			ModelGen1D(popnew, popref);
	//2D
	else if(mLatentDim == 2)	ModelGen2D(popnew, popref);

	return popnew;
}

//modle-based generator for 1D structure
CPopulationMO& ModelHybridB::ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j,clu,size,popcur,clusize,lasm,index;
	double sd,minT,maxT,start,step,totallength;
	std::vector<double> length(mLocalPCA.ClusterSize()),
						C( 3 ), 
						T(popref.Size()), 
						X(popref.Size()),
						t(popref.Size());
	std::vector<unsigned int> dataindex,
							  newsize(mLocalPCA.ClusterSize()),
							  clustertype(mLocalPCA.ClusterSize());
	CPopulationMO par(popref.P()),son(popref.P());

	//Step 1: check to determine which model to use in each cluster
	totallength=0.0; clusize = 0; lasm=0; index = 0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)
	{
		//EDA: modelling
		//=====================================================================================
		//strategy 1: convergence decision in X space 
		if(	mLocalPCA.DataIndex(clu).size()>2 && 
			fabs(mLocalPCA.Eigenvalue(clu)[1]) < mRho*mRho*fabs(mLocalPCA.Eigenvalue(clu)[0]))
		{
			lasm				= clu;
			clustertype[clu]	= 1;
			length[clu]			= sqrt(fabs(mLocalPCA.Eigenvalue(clu)[0]));
			totallength		   += length[clu];
			clusize			   += (unsigned int)mLocalPCA.DataIndex(clu).size();
			index			   += 1;
		}
		//=====================================================================================
		//strategy 2: convergence decision in F space 
		//if(	mLocalPCA.DataIndex(clu).size()>2 )
		//{ 
		//	dataindex.resize(mLocalPCA.DataIndex(clu).size());
		//	std::copy(mLocalPCA.DataIndex(clu).begin(), mLocalPCA.DataIndex(clu).end(),dataindex.begin());
		//	LEARN::Matrix mat(2, (unsigned int)mLocalPCA.DataIndex(clu).size());
		//	for(j=0; j<mat.ColSize(); j++)
		//	{
		//		mat(0,j) = popref[dataindex[j]].F(0);
		//		mat(1,j) = popref[dataindex[j]].F(1);
		//	}
		//	LEARN::PCA pca(mat);
		//	pca.Initialize(2);
		//	pca.Train();
		//	if(pca.Eigenvalue()[1] < mRho*mRho*pca.Eigenvalue()[0])
		//	{
		//		lasm				= clu;
		//		clustertype[clu]	= 1;
		//		length[clu]			= sqrt(fabs(mLocalPCA.Eigenvalue(clu)[0]));
		//		totallength		   += length[clu];
		//		clusize			   += (unsigned int)mLocalPCA.DataIndex(clu).size();
		//		index			   += 1;
		//	}
		//	else
		//	{
		//		clustertype[clu]	= 0;
		//	}
		//}
		//======================================================================================
		//GA: crossover & mutation
		else
		{
			clustertype[clu]	= 0;
		}
	}

#ifdef TEST_CENTER
	std::ofstream fhand("LPCA_model.dat");
	fhand<<index<<std::endl;
#endif

	//Step 2: allocate offspring size to each cluster
	size = 0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)
	{
		if(clustertype[clu]==0)
		{
			newsize[clu] = (unsigned int)mLocalPCA.DataIndex(clu).size();
		}
		else
		{
			//=============================================================================================================
			//strategy 1: according to the "length" of each cluter
			newsize[clu] = (clu < lasm) ? (unsigned int)(clusize*length[clu]/totallength) : (clusize-size);
			//=============================================================================================================
			////strategy 2: use the cluster point size to be offspring size
			//newsize[clu] = (clu < lasm) ? (unsigned int)mLocalPCA.DataIndex(clu).size() : (clusize-size);
			//=============================================================================================================
			size += newsize[clu];
		}
	}

	//Step 3: modelling and sampling
	popcur = 0; minT   = 1.0E100; maxT   =-1.0E100;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)
	{
		size = newsize[clu];

		if(size<1) continue;

		dataindex.resize(mLocalPCA.DataIndex(clu).size());
		std::copy(mLocalPCA.DataIndex(clu).begin(), mLocalPCA.DataIndex(clu).end(),dataindex.begin());

		switch(clustertype[clu])
		{
			//GA: crossover & mutation
			case 0:
				par.Clear(); son.Clear();

				for(i=0; i<dataindex.size(); i++) par.Combine(popref[dataindex[i]]);

				XSBX::Generate(par.Size(), son, par);

				for(i=0; i<size; i++) for(j=0; j<mDataDim; j++) popnew[popcur+i][j] = son[i][j];
				break;
			//EDA: linear model or quadratic model
			case 1:
			case 2:
				//=============================================================================
				//standard deviation strategy 1: by average distance
			#ifdef SD_DIS
				sd = 0.0;
				for(i=0; i<dataindex.size(); i++) sd += mLocalPCA.DisToCore(dataindex[i]);
				sd /= double(size)*sqrt(double(mDataDim));
				//=============================================================================
				//standard deviation strategy 2: by eigenvalues
			#else
				sd = 0.0;
				//for(i=1; i<mDataDim; i++) sd += sqrt(fabs(mLocalPCA.Eigenvalue(clu)[i]));
				//sd /= double(mDataDim-1.0);
				for(i=1; i<mDataDim; i++) sd += fabs(mLocalPCA.Eigenvalue(clu)[i]);
				sd = sqrt(sd / double(mDataDim-1.0));
			#endif
				//=============================================================================

				//calculate the projection into the 1st principal vector
				minT   = 1.0E100; maxT   =-1.0E100;
				for(i=0; i<dataindex.size(); i++)
				{
					T[i]				 = Project(0,dataindex[i],clu);	
					if(T[i]<minT)	minT = T[i];
					if(T[i]>maxT)	maxT = T[i];
				}

				//t in [-mExtension, 1+mExtension]
				start = minT-(maxT-minT)*mExtension;
				step  = (maxT-minT)*(1.0+2*mExtension)/double(size);
				for(i=0; i<size; i++)
				{
					t[i]  = rnd::rand(start,start+step);
					start+= step;
				}
				
				//linear model
				if(clustertype[clu]==1)
				{
					for( i=0; i<size; i++ )
					{
						for( j=0; j<mDataDim; j++ )
							popnew[popcur+i][j] =   t[i] * mLocalPCA.Eigenvector(clu)(j,0)
													+ mLocalPCA.Mean(clu)[j]
													+ sd*rnd::gaussian();
					}
				}
				//quadratic model
				else
				{
					//sample new solution
					for(j=0; j<mDataDim; j++)
					{
						for(i=0; i<dataindex.size(); i++) X[i] = popref[dataindex[i]][j];
						
						//build polynormial model
						alg::poly_fit( T, X, 2, C );

						for(i=0; i<size; i++) popnew[popcur+i][j] = C[0] + C[1]*t[i] + C[2]*t[i]*t[i] + sd*rnd::gaussian();
					}
				}
				break;
			default:
				break;
		}
		popcur += size;
#ifdef TEST_CENTER
		fhand<<minT * mLocalPCA.Eigenvector(clu)(0,0) + mLocalPCA.Mean(clu)[0] <<"\t"
			 <<minT * mLocalPCA.Eigenvector(clu)(1,0) + mLocalPCA.Mean(clu)[1] <<"\t"
			 <<maxT * mLocalPCA.Eigenvector(clu)(0,0) + mLocalPCA.Mean(clu)[0] <<"\t"
			 <<maxT * mLocalPCA.Eigenvector(clu)(1,0) + mLocalPCA.Mean(clu)[1] <<std::endl;
#endif
	}
#ifdef TEST_CENTER
	fhand.close();
#endif
	return popnew;
}

//modle-based generator for 2D structure
CPopulationMO& ModelHybridB::ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j,clu,size,size1,size2,popcur,clusize,lasm;
	double sd,minT1,maxT1,minT2,maxT2,step1,step2,totalarea;
	std::vector<double> area(mLocalPCA.ClusterSize()),
						t1(popref.Size()),
						t2(popref.Size());
	std::vector<unsigned int> dataindex,
							  newsize(mLocalPCA.ClusterSize()),
							  clustertype(mLocalPCA.ClusterSize());
	CPopulationMO par(popref.P()),son(popref.P());

	//Step 1: check to determine which model to use in each cluster
	totalarea=0.0; clusize = 0; lasm = 0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)
	{
		//EDA: modelling
		//=====================================================================================
		//strategy 1: convergence decision in X space 
		if(	mLocalPCA.DataIndex(clu).size()>2 && 
			fabs(mLocalPCA.Eigenvalue(clu)[2]) < mRho*mRho*fabs(mLocalPCA.Eigenvalue(clu)[1]))
		{
			lasm				= clu;
			clustertype[clu]	= 1;
			area[clu]			= sqrt(fabs(mLocalPCA.Eigenvalue(clu)[0]*mLocalPCA.Eigenvalue(clu)[1]));
			totalarea		   += area[clu];
			clusize			   += (unsigned int)mLocalPCA.DataIndex(clu).size();
		}
		//=====================================================================================
		//strategy 2: convergence decision in F space 
		//if(	mLocalPCA.DataIndex(clu).size()>2 )
		//{ 
		//	dataindex.resize(mLocalPCA.DataIndex(clu).size());
		//	std::copy(mLocalPCA.DataIndex(clu).begin(), mLocalPCA.DataIndex(clu).end(),dataindex.begin());
		//	LEARN::Matrix mat(3, (unsigned int)mLocalPCA.DataIndex(clu).size());
		//	for(j=0; j<mat.ColSize(); j++)
		//	{
		//		mat(0,j) = popref[dataindex[j]].F(0);
		//		mat(1,j) = popref[dataindex[j]].F(1);
		//		mat(2,j) = popref[dataindex[j]].F(2);
		//	}
		//	LEARN::PCA pca(mat);
		//	pca.Initialize(3);
		//	pca.Train();
		//	if(pca.Eigenvalue()[2] < mRho*mRho*pca.Eigenvalue()[1])
		//	{
		//		lasm				= clu;
		//		clustertype[clu]	= 1;
		//		area[clu]			= sqrt(fabs(mLocalPCA.Eigenvalue(clu)[0]*mLocalPCA.Eigenvalue(clu)[1]));
		//		totalarea		   += area[clu];
		//		clusize			   += (unsigned int)mLocalPCA.DataIndex(clu).size();
		//	}
		//	else
		//	{
		//		clustertype[clu]	= 0;
		//	}
		//}
		//======================================================================================
		//GA: crossover & mutation
		else
		{
			clustertype[clu]	= 0;
		}
	}

	std::ofstream file("cluster.dat",std::ios::app);
	file<<popref.Size()-clusize<<"\t";
	file.close();

	//Step 2: allocate offspring size to each cluster
	size = 0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)
	{
		if(clustertype[clu]==0)
		{
			newsize[clu] = (unsigned int)mLocalPCA.DataIndex(clu).size();
		}
		else
		{
			//=============================================================================================================
			//strategy 1: according to the "length" of each cluter
			newsize[clu] = (clu < lasm) ? (unsigned int)(clusize*area[clu]/totalarea) : (clusize-size);
			//=============================================================================================================
			////strategy 2: use the cluster point size to be offspring size
			//newsize[clu] = (clu < lasm) ? (unsigned int)mLocalPCA.DataIndex(clu).size() : (clusize-size);
			//=============================================================================================================
			size += newsize[clu];
		}
	}

	//Step 3: modelling and sampling
	popcur = 0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)
	{
		size = newsize[clu];

		if(size<1) continue;

		dataindex.resize(mLocalPCA.DataIndex(clu).size());
		std::copy(mLocalPCA.DataIndex(clu).begin(), mLocalPCA.DataIndex(clu).end(),dataindex.begin());

		switch(clustertype[clu])
		{
			//GA: crossover & mutation
			case 0:
				par.Clear(); son.Clear();

				for(i=0; i<dataindex.size(); i++) par.Combine(popref[dataindex[i]]);

				XSBX::Generate(par.Size(), son, par);

				for(i=0; i<size; i++) for(j=0; j<mDataDim; j++) popnew[popcur+i][j] = son[i][j];
				break;
			//EDA: linear model
			case 1:
				//=============================================================================
				//standard deviation strategy 1: by average distance
			#ifdef SD_DIS
				sd = 0.0;
				for(i=0; i<dataindex.size(); i++) sd += mLocalPCA.DisToCore(dataindex[i]);
				sd /= double(size)*sqrt(double(mDataDim));
				//=============================================================================
				//standard deviation strategy 2: by eigenvalues
			#else
				sd = 0.0;
				//for(i=1; i<mDataDim; i++) sd += sqrt(fabs(mLocalPCA.Eigenvalue(clu)[i]));
				//sd /= double(mDataDim-1.0);
				for(i=2; i<mDataDim; i++) sd += fabs(mLocalPCA.Eigenvalue(clu)[i]);
				sd = sqrt(sd / double(mDataDim-2.0));
			#endif
				//=============================================================================

				//calculate ranges of projections into the 1st and 2nd principal vectors
				minT1 = minT2 = 1.0E100; maxT1 = maxT2 = -1.0E100;
				for(i=0; i<dataindex.size(); i++)
				{
					step1 = Project(0,dataindex[i],clu);				
					if(step1<minT1)		minT1 = step1;
					else if(step1>maxT1)maxT1 = step1;
					
					step1 = Project(1,dataindex[i],clu);				
					if(step1<minT2)		minT2 = step1;
					else if(step1>maxT2)maxT2 = step1;
				}

				size1 = size2 = (unsigned int)( sqrt( size-1.0 ) );
				while( size1*size2 < size )
				{
					size1++;
					size2 = ROUND( double(size1)*(maxT2-minT2)/(maxT1-minT1) );
				}

				//uniform grid
				step1 = (maxT1-minT1)*(1.0+2*mExtension)/double(size1);
				step2 = (maxT2-minT2)*(1.0+2*mExtension)/double(size2);
				for( i=0; i<size1; i++ )
					for( j=0; j<size2; j++ )
						if(i*size2+j < size)
						{
							if(size1<2)
								t1[i*size2+j] = rnd::rand(minT1,maxT1);
							else
								t1[i*size2+j] = ( 1.0 + mExtension ) * minT1 - mExtension * maxT1 + rnd::rand( i*step1, (i+1.0)*step1 );
							if(size2<2)
								t2[i*size2+j] = rnd::rand(minT2,maxT2);
							else
								t2[i*size2+j] = ( 1.0 + mExtension ) * minT2 - mExtension * maxT2 + rnd::rand( j*step2, (j+1.0)*step2 );
						}

				//sample new solution
				for( i=0; i<size; i++ )
				{
					for( j=0; j<mDataDim; j++ )
					{
						popnew[popcur+i][j] =	mLocalPCA.Mean(clu)[j] +
												t1[i]*mLocalPCA.Eigenvector(clu)(j,0) + //dimension 1
												t2[i]*mLocalPCA.Eigenvector(clu)(j,1) + //dimension 2
												sd*rnd::gaussian();						//noise
					}
				}
				break;
			default:
				break;
		}
		popcur += size;
	}

	return popnew;
}

//project point into a dimension
double ModelHybridB::Project(unsigned int dim, unsigned int index, unsigned int clu)
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
