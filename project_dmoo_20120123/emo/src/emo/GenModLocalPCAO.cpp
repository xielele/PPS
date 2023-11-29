/*! \file	Generator_Model_LocalPCA.cpp
	
	\brief	Evolutionary Aglorithm Generator with Local PCA
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Nov.29 2005 make great changes: noise, border checking
	\date	Apr.10 2006 redesign
	\date	Jul.18 2006 add quadratic models
*/
#include <ctime>
#include <list>
#include <vector>
#include <cmath>
#include <float.h>
#include "alg/Fitting.h"
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
ModelLocalPCA::ModelLocalPCA()
{
	mLatentDim	= 0;
	mMaxCluster	= 0;
	mbReset		= true;
}

//initialize the LPCA
void ModelLocalPCA::Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension)
{
	mLatentDim	= latent;
	mMaxCluster	= cluster;
	mTrainSteps	= trainsteps;
	mExtension	= extension;
	mbReset		= true;
}
	
//reset Local PCA
void ModelLocalPCA::Reset()
{
	mbReset = true;
}

//model-based generator
CPopulationMO& ModelLocalPCA::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
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
	
	//Step 3: reset the weights of SOM if necessary
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

	//Step 6: guided crossover
	//GuidedXOver::XOver(popnew, popref);

	return popnew;
}

//modle-based generator for 1D structure
//linear model
CPopulationMO& ModelLocalPCA::ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j,clu,size,popcur;
	double sd,Tmin,Tmax;
	double totallength; unsigned int index;
	std::vector<double> length(mLocalPCA.ClusterSize()),
						T,
						t;
	std::vector<unsigned int> newsize(mLocalPCA.ClusterSize()),dataindex(mLocalPCA.ClusterSize());

	//Step 1: allocate offspring size to each cluster
	totallength=0.0; index=0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)if(mLocalPCA.DataIndex(clu).size()>=2)
	{
		length[index]		= sqrt(fabs(mLocalPCA.Eigenvalue(clu)[0]));
		//dataindex[index]	= clu;
		totallength		   += length[index++];
	}
	size = 0;
	for(clu=0; clu<index; clu++)
	{
		//=============================================================================================================
		//strategy 1: according to the "length" of each cluter
		newsize[clu] = (clu < (index-1)) ? (unsigned int)(popnew.Size()*length[clu]/totallength) : (popnew.Size()-size);
		//=============================================================================================================
		////strategy 2: use the cluster point size to be offspring size
		//newsize[clu] = (clu < (index-1)) ? mLocalPCA.DataIndex(dataindex[clu]).size(): (popref.Size()-size);
		//=============================================================================================================
		size += newsize[clu];
	}

#ifdef AZ_MODEL_OUT
	std::ofstream fhand("model.set");
	fhand<<"PCA"<<std::endl;
	fhand<<index<<std::endl;
#endif

	//Step2: build linear models, the trainning size > 2
	popnew.Resize(size); index=0; popcur=0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)if(mLocalPCA.DataIndex(clu).size()>=2)
	{
		//size = (unsigned int)mLocalPCA.DataIndex(clu).size();
		size = newsize[index++];

		if(size < 1) continue;

		dataindex.resize(mLocalPCA.DataIndex(clu).size());
		std::copy(mLocalPCA.DataIndex(clu).begin(), mLocalPCA.DataIndex(clu).end(),dataindex.begin());

		//=============================================================================
		//standard deviation strategy 1: by average distance
#ifdef SD_DIS
		sd = 0.0;
		for(i=0; i<dataindex.size(); i++) sd += mLocalPCA.DisToCore(dataindex[i])*mLocalPCA.DisToCore(dataindex[i]);
		sd = sqrt(sd/(double(dataindex.size())*(mDataDim-1)));
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
		T.resize(dataindex.size());
		Tmin   = 1.0E100; Tmax   =-1.0E100;
		for(i=0; i<dataindex.size(); i++)
		{
			T[i] = Project(0,dataindex[i],clu);	
			if(T[i]<Tmin)	Tmin = T[i];
			if(T[i]>Tmax)	Tmax = T[i];
		}
		for(i=0; i<dataindex.size(); i++)
		{
			T[i] = (T[i]-Tmin)/(Tmax-Tmin);
		}

		t.resize(size);
		for(i=0; i<size; i++)
		{
			t[i] = rnd::rand(double(i),double(i+1))*(1.0+2*mExtension)/double(size)-mExtension;
		}

		//linear model
		if(true)
		{
			for(i=0; i<size; i++)
			{
				for( j=0; j<mDataDim; j++ )
				{
					popnew[popcur+i][j] =   (t[i] * (Tmax-Tmin) + Tmin) * mLocalPCA.Eigenvector(clu)(j,0)
											+ mLocalPCA.Mean(clu)[j]
											+ sd*rnd::gaussian();

					if( wxFinite(popnew[popcur+i][j]) == 0 )
					{
						popnew[popcur+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
					}

					// border check strategy 1: DE strategy
					if(popnew[popcur+i][j] < popnew.P().XLow(j))		popnew[popcur+i][j] = 0.5*(popnew.P().XLow(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
					else if(popnew[popcur+i][j] > popnew.P().XUpp(j))	popnew[popcur+i][j] = 0.5*(popnew.P().XUpp(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
				}
			}
		}
		//quadratic model
		else
		{
			std::vector<double> X(dataindex.size()), C(3);
			for(j=0; j<mDataDim; j++ )
			{
				for(i=0; i<dataindex.size(); i++) X[i] = pData[dataindex[i]][j];

				//build polynormial model
				alg::poly_fit( T, X, 2, C );

				for(i=0; i<size; i++)
				{
					popnew[popcur+i][j] =	C[0] + 
											C[1]*t[i] + 
											C[2]*t[i]*t[i] + 
											sd*rnd::gaussian();		

					if( wxFinite(popnew[popcur+i][j]) == 0 )
					{
						popnew[popcur+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
					}

					// border check strategy 1: DE strategy
					if(popnew[popcur+i][j] < popnew.P().XLow(j))		popnew[popcur+i][j] = 0.5*(popnew.P().XLow(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
					else if(popnew[popcur+i][j] > popnew.P().XUpp(j))	popnew[popcur+i][j] = 0.5*(popnew.P().XUpp(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
				}
			}
		}
		popcur += size;

#ifdef AZ_MODEL_OUT
		fhand<<Tmin * mLocalPCA.Eigenvector(clu)(0,0) + mLocalPCA.Mean(clu)[0] <<"\t"
			 <<Tmin * mLocalPCA.Eigenvector(clu)(1,0) + mLocalPCA.Mean(clu)[1] <<"\t"
			 <<Tmax * mLocalPCA.Eigenvector(clu)(0,0) + mLocalPCA.Mean(clu)[0] <<"\t"
			 <<Tmax * mLocalPCA.Eigenvector(clu)(1,0) + mLocalPCA.Mean(clu)[1] <<"\t"
			 <<sd<<std::endl;
#endif
	}
#ifdef AZ_MODEL_OUT
	fhand.close();
#endif
	return popnew;
}

//modle-based generator for 2D structure
CPopulationMO& ModelLocalPCA::ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j,k,clu,size,gridc,gricr,popcur,index;
	double sd,Tmin0,Tmax0,Tmin1,Tmax1,totalarea;
	std::vector<double> area(mLocalPCA.ClusterSize()),
						T0,T1,
						t0,t1;
	std::vector<unsigned int> newsize(mLocalPCA.ClusterSize()), dataindex(mLocalPCA.ClusterSize());

	//Step 1: allocate offspring size to each cluster
	totalarea=0.0; index=0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)if(mLocalPCA.DataIndex(clu).size()>2)
	{
		area[index] = sqrt(fabs(mLocalPCA.Eigenvalue(clu)[0]*mLocalPCA.Eigenvalue(clu)[1]));
		//dataindex[index]	= clu;
		totalarea  += area[index++];
	}
	size = 0;
	for(clu=0; clu<index; clu++)
	{
		//=============================================================================================================
		//strategy 1: according to the "area" of each cluter
		newsize[clu] = (clu < (index-1)) ? (unsigned int)(popnew.Size()*area[clu]/totalarea) : (popnew.Size()-size);
		//=============================================================================================================
		////strategy 2: use the cluster point size to be offspring size
		//newsize[clu] = (clu < (index-1)) ? mLocalPCA.DataIndex(dataindex[clu]).size(): (popref.Size()-size);
		//=============================================================================================================
		size += newsize[clu];
	}

	//Step2: build models, the trainning size > 2
	popnew.Resize(size); index=0; popcur=0;
	for(clu=0; clu<mLocalPCA.ClusterSize(); clu++)if(mLocalPCA.DataIndex(clu).size()>1)
	{
		//size = (unsigned int)mLocalPCA.DataIndex(clu).size();
		size = newsize[index++];

		if(size < 1) continue;
		
		dataindex.resize(mLocalPCA.DataIndex(clu).size());
		std::copy(mLocalPCA.DataIndex(clu).begin(), mLocalPCA.DataIndex(clu).end(),dataindex.begin());

		//=============================================================================
		//standard deviation strategy 1: by average distance
#ifdef	SD_DIS
		sd = 0.0;
		for(i=0; i<dataindex.size(); i++) sd += mLocalPCA.DisToCore(dataindex[i])*mLocalPCA.DisToCore(dataindex[i]);
		sd = sqrt(sd/(double(dataindex.size())*(mDataDim-2.0)));
		//=============================================================================
		//standard deviation strategy 2: by eigenvalues
#else
		sd = 0.0;
		//for(i=2; i<mDataDim; i++) sd += sqrt(fabs(mLocalPCA.Eigenvalue(clu)[i]));
		//sd /= double(mDataDim-2.0);
		if(mDataDim>2)
		{
			for(i=2; i<mDataDim; i++) sd += fabs(mLocalPCA.Eigenvalue(clu)[i]);
			sd = sqrt(sd / double(mDataDim-2.0));
		}
#endif
		//=============================================================================

		//calculate ranges of projections into the 1st and 2nd principal vectors
		T0.resize(dataindex.size()); T1.resize(dataindex.size());
		Tmin0 = Tmin1 = 1.0E100; Tmax0 = Tmax1 = -1.0E100;
		for(i=0; i<dataindex.size(); i++)
		{
			T0[i] = Project(0,dataindex[i],clu);				
			if(T0[i]<Tmin0)	Tmin0 = T0[i];
			if(T0[i]>Tmax0)	Tmax0 = T0[i];
			
			T1[i] = Project(1,dataindex[i],clu);				
			if(T1[i]<Tmin1)	Tmin1 = T1[i];
			if(T1[i]>Tmax1)	Tmax1 = T1[i];
		}
		for(i=0; i<dataindex.size(); i++)
		{
			T0[i] = (T0[i] - Tmin0)/(Tmax0-Tmin0);
			T1[i] = (T1[i] - Tmin1)/(Tmax1-Tmin1);
		}

		//offspring grid
		gridc = gricr = (unsigned int)( sqrt( size-1.0 ) );
		while( gridc*gricr < size )
		{
			gridc++;
			gricr = ROUND( double(gridc)*(Tmax1-Tmin1)/(Tmax0-Tmin0) );
		}

		//assign index to offspring
		t0.resize(size); t1.resize(size);
		for(j=0; j<gridc; j++)
		{
			for(i=0; i<gricr; i++)
			{
				k=j*gricr+i;
				if(k<size)
				{
					t0[k] = rnd::rand(double(j),double(j+1))*(1.0+2*mExtension)/double(gridc)-mExtension;
					t1[k] = rnd::rand(double(i),double(i+1))*(1.0+2*mExtension)/double(gricr)-mExtension;
				}
			}
		}

		//sample new solution
		//linear model
		if(true)
		{
			for( i=0; i<size; i++ )
			{
				for( j=0; j<mDataDim; j++ )
				{
					popnew[popcur+i][j] =	mLocalPCA.Mean(clu)[j] +
											(t0[i]*(Tmax0-Tmin0) + Tmin0)*mLocalPCA.Eigenvector(clu)(j,0) + //dimension 1
											(t1[i]*(Tmax1-Tmin1) + Tmin1)*mLocalPCA.Eigenvector(clu)(j,1) + //dimension 2
											sd*rnd::gaussian();						//gaussian noise

					if( wxFinite(popnew[popcur+i][j]) == 0 )
					{
						popnew[popcur+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
					}

					// border check
					if(popnew[popcur+i][j] < popnew.P().XLow(j))		popnew[popcur+i][j] = 0.5*(popnew.P().XLow(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
					else if(popnew[popcur+i][j] > popnew.P().XUpp(j))	popnew[popcur+i][j] = 0.5*(popnew.P().XUpp(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
				}
			}
		}
		else
		{
			std::vector<double> X(dataindex.size()), C(6);
			for(j=0; j<mDataDim; j++ )
			{
				for(i=0; i<dataindex.size(); i++) X[i] = pData[dataindex[i]][j];

				//build polynormial model
				alg::poly_fit( T0, T1, X, 2, C );

				for(i=0; i<t0.size(); i++)
				{
					popnew[popcur+i][j] = C[0] + 
										 C[1]*t0[i] +
										 C[2]*t1[i] +
										 C[3]*t0[i]*t0[i] +
										 C[4]*t0[i]*t1[i] +
										 C[5]*t1[i]*t1[i] +
										 sd*rnd::gaussian();	

					if( wxFinite(popnew[popcur+i][j]) == 0 )
					{
						popnew[popcur+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
					}

					// border check
					if(popnew[popcur+i][j] < popnew.P().XLow(j))		popnew[popcur+i][j] = 0.5*(popnew.P().XLow(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
					else if(popnew[popcur+i][j] > popnew.P().XUpp(j))	popnew[popcur+i][j] = 0.5*(popnew.P().XUpp(j)+popref[dataindex[rnd::rand((unsigned int)0,(unsigned int)dataindex.size())]][j]);
				}
			}
		}
		popcur += size;
	}
	return popnew;
}

//project point into a dimension
double ModelLocalPCA::Project(unsigned int dim, unsigned int index, unsigned int clu)
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
