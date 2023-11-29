/*! \file	GenModRM3.cpp

	\brief	An Enhanced version of RM-MEDA (GenModRM3.cpp)

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex,
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Jul.06 2008 design
*/

#include <algorithm>
#include <list>
#include <vector>
#include <cmath>
#include <fstream>
#include <float.h>
#include <numeric>
#include "alg/HCSampler.h"
#include "alg/Kmeans.h"
#include "emo/GenMod.h"
#include "emo/LogFile.h"

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
//constructor
RM3::RM3()
{
	//set default parameters
	mExtension		= 0.5;
	mUPFExtension	= 0.5;
	mPCAThreshold	= 0.8;
	mNoM			= 20;
	mMaxM			= 30;
}

//deconstructor
RM3::~RM3()
{
}

//set parameters
void RM3::Set(double extension, double threshold, unsigned int maxmodel, unsigned int maxgen)
{
	mExtension		= extension;
	mPCAThreshold	= threshold;
	mMaxM			= maxmodel;
	mMaxG			= maxgen;
	mGen			= 0;

#ifdef AZ_MODEL_OUT
	LogFile log('C');
#endif
}

//model-based generator
CPopulationMO& RM3::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int index, newI;

	popnew.Resize(sizenew);

	//mNoM = mGen < 0.8*mMaxG ? rnd::rand((unsigned int)1, 4*mMaxM/5) : rnd::rand(4*mMaxM/5, mMaxM);
	//mNoM = mGen < 0.5*mMaxG ? rnd::rand((unsigned int)1, mMaxM/2) : rnd::rand(mMaxM/2, mMaxM);
	if(mGen >= (10 *(mMaxM-5))) 
		mNoM = mMaxM;
	else	
		mNoM = 5 + ((mGen/10) % (mMaxM-5));


	mNoNew = sizenew; if(mNoNew<1) return popnew;

	mvPops.resize(mNoM);

	//Step 1: clustering
	Partition(popref);

	//Step 2: build model
	// assign data to each model
	// dimension of the manifold in each cluster
	// std. (noise) of each model
	// post probability of each model
	// spend more cost on bad clusters
	BuildModel(popref, mNoNew);

	//Step 3: sample new solutions
#ifdef AZ_MODEL_OUT
	std::ofstream fhand("ModX.out");
	fhand<<mNoM<<std::endl;
	fhand.close();

	LogFile log;
	log<<mGen<<"\t";
#endif
	newI = 0;
	for(index=0; index<mNoM; index++)
	{
		if(mvPops[index].method == MODEL)
			ALModelGen(popnew, newI, popref, mvPops[index].pca, mvPops[index].std, mvPops[index].dim, mvPops[index].no);
#ifdef AZ_MODEL_OUT
	log<<mvDim[index]<<"\t";
#endif
	}
	popnew.Erase(newI);
#ifdef AZ_MODEL_OUT
	log<<"\n";
#endif

	mGen++;
	return popnew;
}

//modle-based generator
//linear model
CPopulationMO& RM3::ALModelGen(CPopulationMO& popnew, unsigned int& index, CPopulationMO& popref, alg::PCA& pca, double std, unsigned int dim, unsigned int size)
{
	unsigned int i,j,k;
	std::vector< double > Tmax, Tmin;
	std::vector< std::vector<double> > T;

	//Step 1: the model is not a single point, set parameters
	if(dim>0)
	{
		double pro, pro1, pro2;
		// projection into primary components
		Tmax.resize(dim); Tmin.resize(dim);
		for(i=0; i<dim; i++) { Tmin[i] = 1.0E200; Tmax[i] = -1.0E200; }
		for(i=0; i<pca.Data().ColSize(); i++)
		{
			for(j=0; j<dim; j++)
			{
				pro = Project(j,i,pca);
				if(pro > Tmax[j]) Tmax[j] = pro;
				if(pro < Tmin[j]) Tmin[j] = pro;
			}
		}

		// sample latent points in the manifold
		double ext = (pow(1.0+mExtension, 1.0/(dim+0.0))-1.0)*0.5;
		for(i=0; i<dim; i++)
		{
			pro1	= Tmin[i];
			pro2	= Tmax[i];
			Tmin[i] = pro1-ext*(pro2-pro1);
			Tmax[i] = pro2+ext*(pro2-pro1);
		}

		T.resize(dim);
		for(i=0; i<dim; i++) T[i].resize(size);
		alg::LHC(T, Tmin, Tmax);
	}

	//Step 2: sample new trial solutions in decision space
	unsigned int par1;
	for(i=0; i<size; i++)
	{
		par1 = rnd::rand((unsigned int)0, popref.Size());
		for( j=0; j<popref.P().XSize(); j++ )
		{
			popnew[index+i][j] = pca.Mean()[j] + std*rnd::gaussian();
			// !!! if dim == 0, it will be a Gaussian model and the following sentence will be ignored
			for(k=0; k<dim; k++) popnew[index+i][j] += pca.Eigenvector()(j,k)*T[k][i];
			// if it is not a leagle float
			if( wxFinite(popnew[index+i][j]) == 0 ) popnew[index+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));

			// boundary checking
			if(popnew[index+i][j] < popnew.P().XLow(j))			popnew[index+i][j] = 0.5*(popnew.P().XLow(j)+popref[par1][j]);
			else if(popnew[index+i][j] > popnew.P().XUpp(j))	popnew[index+i][j] = 0.5*(popnew.P().XUpp(j)+popref[par1][j]);
		}
	}
	index += size;

#ifdef AZ_MODEL_OUT
	std::ofstream fhand("ModX.out",std::ios::app);
	fhand<<std::scientific<<std::setprecision(5);
	switch(dim)
	{
	case 0:
		fhand<<"SPHERE"<<std::endl
			 <<pca.Mean()[0] <<"\t"
			 <<pca.Mean()[1] <<"\t"
			 <<pca.Mean()[2] <<"\t"
			 <<std;
		break;
	case 1:
		fhand<<"LINE"<<std::endl;
		break;
	case 2:
		fhand<<"RECT"<<std::endl;
		break;
	case 3:
		fhand<<"CUBE"<<std::endl;
		break;
	default:
		//fhand<<"ERROR"<<std::endl;
		fhand<<"CUBE"<<std::endl;
		dim = 3;
		break;
	}
	for(i=0; i<dim; i++)
	{
		 fhand<<Tmin[i] * pca.Eigenvector()(0,i) + pca.Mean()[0] <<"\t"
			  <<Tmin[i] * pca.Eigenvector()(1,i) + pca.Mean()[1] <<"\t"
			  <<Tmin[i] * pca.Eigenvector()(2,i) + pca.Mean()[2] <<"\t"
			  <<Tmax[i] * pca.Eigenvector()(0,i) + pca.Mean()[0] <<"\t"
			  <<Tmax[i] * pca.Eigenvector()(1,i) + pca.Mean()[1] <<"\t"
			  <<Tmax[i] * pca.Eigenvector()(2,i) + pca.Mean()[2] <<"\t";
	}
	fhand<<std<<std::endl;
	fhand.close();
#endif

	return popnew;
}

//project point into a dimension
double RM3::Project(unsigned int dim, unsigned int index, alg::PCA& pca)
{
	unsigned int i; double len;
	len = 0.0;
	//projection
	for( i=0; i<pca.Data().RowSize(); i++ ) len += (pca.Data()(i,index) - pca.Mean()[i]) * pca.Eigenvector()(i,dim);
	return len;
}

// choose the neighborhood of each reference point
void RM3::Partition(CPopulationMO& pop)
{
	unsigned int i, j, k, size = pop.Size();

	std::vector<unsigned int> index(size), location(mNoM+2);
	if(pop[0].ID() == pop[1].ID())
	{
		for(i=0; i<size; i++) index[i] = i;
	}
	else
	{
		for(i=0; i<size; i++) index[pop[i].ID()] = i;
	}

	unsigned int interval = size/(mNoM+1), left = size-interval*(mNoM+1);
	if(left>0) left = rnd::rand((unsigned int)0,left+1);
	location[0] = 0; location[mNoM+1] = size-1; for(i=1; i<mNoM+1; i++) location[i] = i*interval + left; 

	for(k=0; k<mNoM; k++)
	{
		mvPops[k].ratio = 0.0;
		mvPops[k].data.Resize(pop.P().XSize(), location[k+2]-location[k]+1);
		for(i=location[k]; i<=location[k+2]; i++)
		{
			for(j=0; j<pop.P().XSize(); j++) mvPops[k].data(j,i-location[k]) = pop[index[i]][j];
			if(pop[index[i]].Rank()>1) mvPops[k].ratio += 1.0;
		}
		mvPops[k].ratio /= double(location[k+2]-location[k]+1);
	}
}

// set model parameters
void RM3::BuildModel(CPopulationMO& pop, unsigned int sizenew)
{
	unsigned int i, index, dim = pop.P().XSize();

	// Step 0: select subpopulations to model and sample
	for(i=0; i<mNoM; i++) mvPops[i].method = MODEL;
	// 'good' subpopulations have a probability of 0.5 to be selected for modeling and sampling
	std::vector<unsigned int> vgood;
	for(i=0; i<mNoM; i++) if(mvPops[i].ratio < 0.5) vgood.push_back(i);
	std::random_shuffle(vgood.begin(), vgood.end());
	for(i=0; i<vgood.size()/2; i++) mvPops[vgood[i]].method = NOMETHOD;

	// Step 1: analyze data (pca)
	for(index=0; index<mNoM; index++) if(mvPops[index].method == MODEL)
	{
		mvPops[index].pca.Data(mvPops[index].data);
		mvPops[index].pca.Initialize(dim);
		mvPops[index].pca.Train();
	}

	//Step 2: set the dimension and std. of each subset
	for(index=0; index<mNoM; index++) {mvPops[index].dim = dim; mvPops[index].std = 0.0;}
	for(index=0; index<mNoM; index++) if(mvPops[index].method == MODEL)
	{
		//Step 2.1: calculate the standard deviation which represents the noise and the dimension of manifold of PS
		double ad	= 0.0;
		for(i=0; i<dim; i++) ad += fabs(mvPops[index].pca.Eigenvalue()[i]);
		for(i=0; i<dim; i++)
		{
			mvPops[index].std += fabs(mvPops[index].pca.Eigenvalue()[i]);
			if( mvPops[index].std >= mPCAThreshold*ad )
			{
				mvPops[index].dim = i+1;
				break;
			}
		}

//#define REPAIR_DIM 1
#ifdef REPAIR_DIM
		//Step 2.2: repair standard deviation and dimension
		if(dim > 10 && mvDim[index] >=4)
		{
			mvDim[index]	= 4;
			mvStd[index]	= 0.0;
			for(i=0; i<mvDim[index]; i++) mvStd[index] += fabs(mvPCA[index].Eigenvalue()[i]);
		}
#endif
		mvPops[index].std = sqrt((ad - mvPops[index].std) / double(dim-mvPops[index].dim));
	}

	// Step 3: post probability of each model (how many offspring to generate for each model)
	// probability = 1.0/mNoM
	unsigned int work	= 0;
	unsigned int kk		= 0;
	unsigned int tolNo	= 0;
	std::vector<unsigned int> iid; 
	for(index=0; index<mNoM; index++) if(mvPops[index].method != NOMETHOD) work++; iid.resize(work);
	for(index=0; index<mNoM; index++) if(mvPops[index].method != NOMETHOD)
	{
		mvPops[index].no = (unsigned int)(double(sizenew)/double(work));
		tolNo			+= mvPops[index].no;
		iid[kk++]		 = index;
	}
	for(index=0; index<sizenew-tolNo; index++) mvPops[iid[rnd::rand((unsigned int)0,work)]].no++;
}
} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az

////modle-based generator
////linear model
//CPopulationMO& RM3::GAGen(CPopulationMO& popnew, unsigned int& index, unsigned int clu, unsigned int size)
//{
//	unsigned int i,j,num,r1,r2,r3;
//	num = mvDatMat[clu].ColSize();
//	for(i=0; i<size; i++)
//	{
//		do{r1=rnd::rand((unsigned int)(0), num);}while(false);
//		do{r2=rnd::rand((unsigned int)(0), num);}while(r2==r1);
//		do{r3=rnd::rand((unsigned int)(0), num);}while(r3==r2||r3==r1);
//		for( j=0; j<popnew.P().XSize(); j++ )
//		{
//			popnew[index+i][j] = mvDatMat[clu](j,r1) + 0.5*(mvDatMat[clu](j,r2)-mvDatMat[clu](j,r3));
//			// if it is not a leagle float
//			if( wxFinite(popnew[index+i][j]) == 0 ) popnew[index+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
//
//			// boundary checking
//			if(popnew[index+i][j] < popnew.P().XLow(j))			popnew[index+i][j] = 0.5*(popnew.P().XLow(j)+mvDatMat[clu](j,r1));
//			else if(popnew[index+i][j] > popnew.P().XUpp(j))	popnew[index+i][j] = 0.5*(popnew.P().XUpp(j)+mvDatMat[clu](j,r1));
//		}
//		PM(popnew[index+i]);
//		popnew[index+i].OPT(1);
//	}
//	index += size;
//	return popnew;
//}