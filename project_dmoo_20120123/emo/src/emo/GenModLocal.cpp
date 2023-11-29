/*! \file	Generator_ModelPCA.cpp
	
	\brief	Evolutionary Aglorithm Generator: Model/PCA based generator
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Apr.07 2006 design
	\date	Jul.17 2006 modify, add six mating restriction strategies
*/

#include <algorithm>
#include <list>
#include <vector>
#include <cmath>
#include <fstream>
#include <float.h>
#include "alg//Fitting.h"
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
ModelPCA::ModelPCA()
{
	//set default parameters
	mLatentDim	= 1;
	mMaxGen		= 100;
		
	mModelSize	= 10;
	mDeltaT		= 10;

	mMatingStrategy	= 1;
	mCurGen			= 0;
}

//deconstructor
ModelPCA::~ModelPCA()
{
	//clear history pool
	if(pvHis.size() >  0) 
	{
		for(unsigned int i=0; i<pvHis.size(); i++) 
			if(pvHis[i] != 0) delete pvHis[i]; 
		pvHis.clear();
	}
}

//!\brief	set parameters
//!\param	extension	extension ratio
//!\param	maxgen		maximum generation
//!\param	modsize		the number of model points
//!\param	mate		mating restriction strategy
//!\param	threshold	convergence threshold
//!\return	void
void ModelPCA::Set(double extension, double threshold, unsigned int maxgen, unsigned int modsize, unsigned int mate)
{
	mExtension	= extension;
	mMaxGen		= maxgen;
	mModelSize	= modsize;
	mCurGen		= 0;

	mMatingStrategy	= mate;

	//set history pool data structure
	mDeltaT		= 10;
	if(pvHis.size() >  0) {for(unsigned int i=0; i<pvHis.size(); i++) delete pvHis[i]; pvHis.clear();}
	pvHis.resize(mDeltaT);
	if(pvHis.size() >  0) {for(unsigned int i=0; i<pvHis.size(); i++) pvHis[i] = 0;}
	mHisIndex	= 0;
	
	//convergence threshold
	mConThreshold	= threshold;

	mbConverged		= false;
	mSigmoidBeta	= 7.5;

	mCandidateSize	= mLastCandidateSize = mModelSize;
}

//model-based generator
CPopulationMO& ModelPCA::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i, j, index, newI;
	
	//Step 0: set some parameters
	//latent space dimension
	mLatentDim  = popref.P().FSize()-1;
	popnew.Resize(sizenew);
	//number of mvIndex
	index = (unsigned int)(2*double(popref.Size())/double(mModelSize));
	index = std::min(index, popref.Size()/2);
	mvIndex.resize(index);

    //==================================================================
	//Step 1: select the index points
	const unsigned int IndexStrategy = 1;
	std::vector<unsigned int> index_tmp;
	switch(IndexStrategy)
	{
	//Strategy 1: MaxiMin selection
	case 1:
		MaxiMin(popref, mvIndex, (unsigned int)mvIndex.size());
		break;
	//Strategy 2: random selection
	case 2:
	default:
		index_tmp.resize(popref.Size());
		for(i=0; i<popref.Size(); i++) index_tmp[i] = i;
		std::random_shuffle(index_tmp.begin(), index_tmp.end());
		for(i=0; i<mvIndex.size(); i++) mvIndex[i]	= index_tmp[i];
		break;
	}

	//===================================================================
	//Step 2: mating restriction
	double ratio,convergence;
	switch(mMatingStrategy)
	{
	// strategy 1: adaptive mating restriction: small -> large with sigmoid function
	// ^
	// |				____________
	// |			   /
	// |			  /	
	// |			 /	
	// |___________/
	// |
	// |---------------------------> t
	case 1:
		convergence = ConvergenceTest(popref);
		//find the converge point
		if(!mbConverged && mCurGen > 0.2*mMaxGen && convergence <= mConThreshold)
		{
			mbConverged	= true;
			mSigmoidBeta= (0.25 + 0.75*(mCurGen+0.0)/double(mMaxGen+0.0))*15.0;
		}
		
		if(!mbConverged)
		{
			mCandidateSize = mLastCandidateSize;
		}
		else
		{
			ratio			= 1.0/(1.0+exp( -15.0*(mCurGen)/double(mMaxGen) + mSigmoidBeta));
			mCandidateSize	= mLastCandidateSize + (unsigned int)((popref.Size()-mModelSize-mLastCandidateSize)* ratio);
		}
		break;
	// strategy 2: adaptive mating restriction: small -> large with sigmoid function
	// ^
	// |				____________
	// |			   /
	// |		 _____/	
	// |	    /	
	// |_______/
	// |
	// |---------------------------> t
	case 2:
		convergence = ConvergenceTest(popref);
		//test converge or not
		//not converge -> converge
		if(!mbConverged && mCurGen > 0.2*mMaxGen && convergence <= mConThreshold)
		{
			mbConverged			= true;
			mSigmoidBeta		= (0.25 + 0.75*(mCurGen+0.0)/double(mMaxGen+0.0))*15.0;
			mLastCandidateSize	= mCandidateSize;
		}
		else if(mCurGen <= 0.2*mMaxGen || convergence > mConThreshold)
		{
			mbConverged			= false;
			mLastCandidateSize	= mCandidateSize;
		}
		
		if(!mbConverged)
		{
			mCandidateSize = mLastCandidateSize;
		}
		else
		{
			ratio			= 1.0/(1.0+exp( -15.0*(mCurGen)/double(mMaxGen) + mSigmoidBeta));
			mCandidateSize	= mLastCandidateSize + (unsigned int)((popref.Size()-mModelSize-mLastCandidateSize)* ratio);
		}
		break;
	// strategy 3: adaptive mating restriction: small -> large
	// ^
	// |		    ____________
	// |		   |
	// |		   |
	// |		   |
	// |___________|
	// |
	// |---------------------------> t
	case 3:
		convergence = ConvergenceTest(popref);
		//find the converge point
		if(!mbConverged && mCurGen > 0.2*mMaxGen && convergence <= mConThreshold)
		{
			mbConverged	= true;
		}
		
		if(!mbConverged)
		{
			mCandidateSize	= mModelSize; //== mLastCandidateSize;
		}
		else
		{
			mCandidateSize	= popref.Size()-mModelSize;
		}
		break;
	// strategy 4: adaptive mating restriction: small->large->small->large...
	// ^
	// |		____	  ________
	// |	   |	|	 |
	// |	   |	|	 |
	// |	   |	|	 |	
	// |_______|	|____|
	// |
	// |---------------------------> t
	case 4:
		convergence = ConvergenceTest(popref);
		//test converge or not
		if(mCurGen > 0.2*mMaxGen && convergence <= mConThreshold)
		{
			mbConverged	= true;
		}
		else
		{
			mbConverged	= false;
		}
		
		if(!mbConverged)
		{
			mCandidateSize	= mModelSize; //==mLastCandidateSize;
		}
		else
		{
			mCandidateSize	= popref.Size()-mModelSize;
		}
		break;
	//strategy 5: minimal mating restriction
	// ^
	// |
	// |
	// |
	// |
	// |___________________________
	// |
	// |---------------------------> t
	case 5:
	default:
		mCandidateSize	= mModelSize; //==mLastCandidateSize;		
		break;
	//strategy 6: maximal mating restriction
	// ^
	// |___________________________
	// |
	// |
	// |
	// |
	// |
	// |---------------------------> t
	case 6:
		mCandidateSize	= popref.Size()-mModelSize;
		break;
	}
	//if(mCurGen<1)
	//{
	//	LOG::LogFile log;//('C');
	//	log<<"\t"<<std::endl<<mCandidateSize<<"\t";
	//}
	//else
	//{
	//	LOG::LogFile log;
	//	log<<mCandidateSize<<"\t";
	//}

	//===================================================================
	//Step 3: find candidate solutions for each index point
	{
		double dis;
		std::list<double>			ldis;						//neighborhood distance 
		std::list<unsigned int>		idis;						//neighborhood index
		std::list<double>::iterator			it1;
		std::list<unsigned int>::iterator	it2;
		std::vector<unsigned int>	candidate(popref.Size());	//candidate
		for(i=0; i<popref.Size(); i++) candidate[i] = i;
		
		mvModelIndex.resize((unsigned int)mvIndex.size());
		for(index=0; index<(unsigned int)mvIndex.size(); index++)
		{
			std::random_shuffle(candidate.begin(), candidate.end());
			mvModelIndex[index].resize(mModelSize);
			//selecte model points from candidate points
			if(mCandidateSize>mModelSize)
			{
				ldis.clear(); idis.clear();
				ldis.push_back(0.0); idis.push_back(mvIndex[index]);
				for(i=0; i<mCandidateSize; i++) if(mvIndex[index] != candidate[i])
				{
					//dis = 0; for(j=0; j<popref.P().XSize(); j++) dis += (popref[candidate[i]][j] - popref[mvIndex[index]][j])*(popref[candidate[i]][j] - popref[mvIndex[index]][j]);
					dis = 0; for(j=0; j<popref.P().FSize(); j++) dis += (popref[candidate[i]].F(j) - popref[mvIndex[index]].F(j))*(popref[candidate[i]].F(j) - popref[mvIndex[index]].F(j));
					it1 = ldis.begin(); it2 = idis.begin();
					while(it1 != ldis.end() && *it1 < dis) {it1++; it2++;}
					ldis.insert(it1, 1, dis); idis.insert(it2, 1, candidate[i]);
				}
				it2 = idis.begin();
				for(i=0; i<mModelSize; i++) mvModelIndex[index][i] = *it2++; 
			}
			else
			{	
				i=0; j=0;
				while(i<mModelSize) 
				{
					if(mvIndex[index] != candidate[j])
						mvModelIndex[index][i++] = candidate[j];
					j++;
				}
			}//end if
		}//end for
	}//end Step 3.

	//===================================================================
	//Step 4: assign offspring size to each group 
	mvOffSize.resize(mvIndex.size());
	const unsigned int OffStrategy = 1;
	switch(OffStrategy)
	{
	//Strategy 1: assign equally
	case 1:
		Assign_Naive(popref, popnew.Size(), mvOffSize);
		break;
	//Strategy 2: assign accroding to density
	case 2:
	default:
		Assign_Density(popref, popnew.Size(), mvOffSize);
		break;
	}

	//===================================================================
	//Step 5: create new solutions
#ifdef AZ_MODEL_OUT
	if(mLatentDim == 1)	
	{
		std::ofstream fhand("model.set");
		fhand<<"PCA"<<std::endl;
		fhand<<mvIndex.size()<<std::endl;
		fhand.close();
	}
#endif
	newI = 0;
	for(index=0; index<(unsigned int)mvIndex.size(); index++)
	{
		//Step 5.1: doing PCA in this cluster
		mMat.Resize(popref.P().XSize(), (unsigned int)mvModelIndex[index].size());
		for(i=0; i<mvModelIndex[index].size(); i++)
		{
			for(j=0; j<popref.P().XSize(); j++) mMat(j,i) = popref[mvModelIndex[index][i]][j];
		}

		alg::PCA pca(mMat);
		pca.Initialize(popref.P().XSize());
		pca.Train();

		//Step 3.3: create new solutions
		//1D structure
		if(mLatentDim == 1)			ModelGen1D(popnew, popref, pca, newI, mvOffSize[index]);
		//2D structure
		else if(mLatentDim == 2)	ModelGen2D(popnew, popref, pca, newI, mvOffSize[index]);
	}

	popnew.Erase(newI);

	mCurGen++;

	return popnew;
}

//modle-based generator for 1D structure
//linear model
CPopulationMO& ModelPCA::ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref, alg::PCA& pca, unsigned int& index, unsigned int size)
{
	unsigned int i,j;

	//standard deviation
	double sd = 0.0;
	if(popref.P().XSize()>1)
	{
		for(i=1; i<popref.P().XSize(); i++) sd += fabs(pca.Eigenvalue()[i]);
		sd = sqrt(sd / double(popref.P().XSize()-1.0));
	}

	//projection
	std::vector<double> t(size), T(pca.Data().ColSize());
	double Tmin = 1.0E200, Tmax = -1.0E200;
	for(i=0; i<T.size(); i++) 
	{
		T[i] = Project(0,i,pca);
		if(T[i] > Tmax) Tmax = T[i];
		if(T[i] < Tmin) Tmin = T[i];
	}
	for(i=0; i<T.size(); i++) T[i] = (T[i] - Tmin)/(Tmax-Tmin);

#ifdef AZ_MODEL_OUT
	std::ofstream fhand("model.set",std::ios::app);
	fhand<<Tmin * pca.Eigenvector()(0,0) + pca.Mean()[0] <<"\t"
		 <<Tmin * pca.Eigenvector()(1,0) + pca.Mean()[1] <<"\t"
		 <<Tmax * pca.Eigenvector()(0,0) + pca.Mean()[0] <<"\t"
		 <<Tmax * pca.Eigenvector()(1,0) + pca.Mean()[1] <<"\t"
		 <<sd<<std::endl;
	fhand.close();
#endif

	//index for new trail solutions t \in [-mExtension, 1.0+mExtension]
	for(i=0; i<t.size(); i++) t[i] = rnd::rand(double(i),double(i+1))*(1.0+2*mExtension)/double(t.size())-mExtension;

	//linear model
	if(true)
	{
		for(i=0; i<t.size(); i++)
		{
			for( j=0; j<popref.P().XSize(); j++ )
			{
				popnew[index+i][j] = (t[i]*(Tmax-Tmin)+Tmin) * pca.Eigenvector()(j,0)
									+ pca.Mean()[j]
									+ sd*rnd::gaussian();

				if( wxFinite(popnew[index+i][j]) == 0 )
				{
					popnew[index+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
				}

				// border check strategy
				if(popnew[index+i][j] < popnew.P().XLow(j))			popnew[index+i][j] = 0.5*(popnew.P().XLow(j)+pca.Data()(j, rnd::rand((unsigned int)0,(unsigned int)pca.Data().ColSize())));
				else if(popnew[index+i][j] > popnew.P().XUpp(j))	popnew[index+i][j] = 0.5*(popnew.P().XUpp(j)+pca.Data()(j, rnd::rand((unsigned int)0,(unsigned int)pca.Data().ColSize())));
			}
		}
	}
	//quadratic model
	else
	{
		std::vector<double> X(pca.Data().ColSize()), C(3);
		for(j=0; j<popref.P().XSize(); j++ )
		{
			for(i=0; i<X.size(); i++) X[i] = pca.Data()(j,i);

			//build polynormial model
			alg::poly_fit( T, X, 2, C );

			for(i=0; i<t.size(); i++)
			{
				popnew[index+i][j] = C[0] + 
									 C[1]*t[i] + 
									 C[2]*t[i]*t[i] + 
									 sd*rnd::gaussian();	

				if( wxFinite(popnew[index+i][j]) == 0 )
				{
					popnew[index+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
				}

				// border check strategy
				if(popnew[index+i][j] < popnew.P().XLow(j))			popnew[index+i][j] = 0.5*(popnew.P().XLow(j)+pca.Data()(j, rnd::rand((unsigned int)0,(unsigned int)pca.Data().ColSize())));
				else if(popnew[index+i][j] > popnew.P().XUpp(j))	popnew[index+i][j] = 0.5*(popnew.P().XUpp(j)+pca.Data()(j, rnd::rand((unsigned int)0,(unsigned int)pca.Data().ColSize())));
			}
		}
	}
	index += (unsigned int)(t.size());

	return popnew;
}

//modle-based generator for 2D structure
CPopulationMO& ModelPCA::ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref, alg::PCA& pca, unsigned int& index, unsigned int size)
{
	unsigned int i,j,k,gridr,gridc;

	//standard deviation
	double sd = 0.0;
	if(popref.P().XSize()>2)
	{
		for(i=2; i<popref.P().XSize(); i++) sd += fabs(pca.Eigenvalue()[i]);
		sd = sqrt(sd / double(popref.P().XSize()-2.0));
	}

	//projection in primary eigenvector
	std::vector<double> t0(size), t1(size), T0(pca.Data().ColSize()), T1(pca.Data().ColSize());
	double Tmin0,Tmin1,Tmax0,Tmax1;
	Tmin0 = Tmin1 = 1.0E200; Tmax0 = Tmax1 = -1.0E200;
	for(i=0; i<T0.size(); i++) 
	{
		T0[i] = Project(0,i,pca);
		if(T0[i] > Tmax0) Tmax0 = T0[i];
		if(T0[i] < Tmin0) Tmin0 = T0[i];
		T1[i] = Project(1,i,pca);
		if(T1[i] > Tmax1) Tmax1 = T1[i];
		if(T1[i] < Tmin1) Tmin1 = T1[i];
	}
	for(i=0; i<T0.size(); i++)
	{
		T0[i] = (T0[i] - Tmin0)/(Tmax0-Tmin0);
		T1[i] = (T1[i] - Tmin1)/(Tmax1-Tmin1);
	}

	//offspring grid
	gridc = gridr = (unsigned int)( sqrt( size-1.0 ) );
	while( gridc*gridr < size )
	{
		gridc++;
		gridr = ROUND( double(gridc)*(Tmax1-Tmax1)/(Tmax0-Tmin0) );
	}
	
	//assign index to offspring
	for(j=0; j<gridc; j++)
	{
		for(i=0; i<gridr; i++)
		{
			k=j*gridr+i;
			if(k<size)
			{
				t0[k] = rnd::rand(double(j),double(j+1))*(1.0+2*mExtension)/double(gridc)-mExtension;
				t1[k] = rnd::rand(double(i),double(i+1))*(1.0+2*mExtension)/double(gridr)-mExtension;
			}
		}
	}

	//linear model
	if(true)
	{
		for(i=0; i<t0.size(); i++)
		{
			for( j=0; j<popref.P().XSize(); j++ )
			{
				popnew[index+i][j] = (t0[i]*(Tmax0-Tmin0)+Tmin0) * pca.Eigenvector()(j,0)
									+(t1[i]*(Tmax1-Tmin1)+Tmin1) * pca.Eigenvector()(j,1)
									+ pca.Mean()[j]
									+ sd*rnd::gaussian();

				if( wxFinite(popnew[index+i][j]) == 0 )
				{
					popnew[index+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
				}

				// border check strategy
				if(popnew[index+i][j] < popnew.P().XLow(j))			popnew[index+i][j] = 0.5*(popnew.P().XLow(j)+pca.Data()(j, rnd::rand((unsigned int)0,(unsigned int)pca.Data().ColSize())));
				else if(popnew[index+i][j] > popnew.P().XUpp(j))	popnew[index+i][j] = 0.5*(popnew.P().XUpp(j)+pca.Data()(j, rnd::rand((unsigned int)0,(unsigned int)pca.Data().ColSize())));
			}
		}
	}
	//quadratic model
	else
	{
		std::vector<double> X(pca.Data().ColSize()), C(6);
		for(j=0; j<popref.P().XSize(); j++ )
		{
			for(i=0; i<X.size(); i++) X[i] = pca.Data()(j,i);

			//build polynormial model
			alg::poly_fit( T0, T1, X, 2, C );

			for(i=0; i<t0.size(); i++)
			{
				popnew[index+i][j] = C[0] + 
									 C[1]*t0[i] +
									 C[2]*t1[i] +
									 C[3]*t0[i]*t0[i] +
									 C[4]*t0[i]*t1[i] +
									 C[5]*t1[i]*t1[i] +
									 sd*rnd::gaussian();	

				if( wxFinite(popnew[index+i][j]) == 0 )
				{
					popnew[index+i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));
				}

				// border check strategy
				if(popnew[index+i][j] < popnew.P().XLow(j))			popnew[index+i][j] = 0.5*(popnew.P().XLow(j)+pca.Data()(j, rnd::rand((unsigned int)0,(unsigned int)pca.Data().ColSize())));
				else if(popnew[index+i][j] > popnew.P().XUpp(j))	popnew[index+i][j] = 0.5*(popnew.P().XUpp(j)+pca.Data()(j, rnd::rand((unsigned int)0,(unsigned int)pca.Data().ColSize())));
			}
		}
	}
	index += (unsigned int)(t0.size());

	return popnew;
}

//project point into a dimension
double ModelPCA::Project(unsigned int dim, unsigned int index, alg::PCA& pca)
{
	unsigned int i; double len;
	len = 0.0;
	//projection
	for( i=0; i<pca.Data().RowSize(); i++ ) len += (pca.Data()(i,index) - pca.Mean()[i]) * pca.Eigenvector()(i,dim);
	return len;
}

// assign offspring to each group with equal probability
void ModelPCA::Assign_Naive(CPopulationMO& pop, unsigned int tolsize, std::vector<unsigned int>& number)
{
	unsigned int i, tolNo=0;
	for(i=0; i<(unsigned int)number.size(); i++)
	{
		number[i] = (unsigned int)(double(tolsize) / double(number.size()));
		tolNo    += number[i];
	}
	for(i=0; i<tolsize-tolNo; i++) number[rnd::rand((unsigned int)0,(unsigned int)number.size())]++;
}

// assign offspring to each group according to density
void ModelPCA::Assign_Density(CPopulationMO& pop, unsigned int tolsize, std::vector<unsigned int>& number)
{
	unsigned int index, i, j, baseNo, leftNo, tolNo=0;
	double dis, T = 0.0;
	
	//Step 1: density estimation for each point
	std::vector<double>	density(pop.Size()), ndis(2);
	for(index=0; index<pop.Size(); index++)
	{
		ndis[0] = 1.0E20; ndis[1] = 1.0E30;
		for(i=0; i<pop.Size(); i++) if(i != index)
		{
			dis = 0; for(j=0; j<pop.P().FSize(); j++) dis += (pop[i].F(j) - pop[index].F(j))*(pop[i].F(j) - pop[index].F(j));
			if(dis<ndis[0]) {ndis[1]=ndis[0]; ndis[0]=dis;}
			else if(dis<ndis[1]) {ndis[1]=dis;}
		}
		density[index] = sqrt(ndis[1]);
	}
	
	//Step 2: density estimation for each group
	std::vector<double> dos(number.size());
	for(i=0; i<(unsigned int)number.size(); i++)
	{
		dos[i] = 0.0;
		for(j=0; j<mvModelIndex[i].size(); j++) dos[i] += density[mvModelIndex[i][j]];
		dos[i] /= double(mvModelIndex[i].size());
		T   += dos[i];
	}
	if(T==0.0) T = 1.0;
		
	//Step 3: assign
	baseNo = 1;
	leftNo = tolsize - baseNo*mvIndex.size();
	for(i=0; i<(unsigned int)mvIndex.size(); i++)
	{
		number[i] = baseNo + (unsigned int)(leftNo * dos[i]/T);
		tolNo	 += number[i];
	}
	for(i=0; i<tolsize-tolNo; i++) number[rnd::rand((unsigned int)0,(unsigned int)number.size())]++;
}

// MaxiMin select 
void ModelPCA::MaxiMin(CPopulationMO& pop, std::vector<unsigned int>& index, unsigned int size)
{
	unsigned int i, j, k, iindex;
	double fmin;
	std::vector< std::vector<double> > dismat;	// distance matrix
	std::vector< double >	dismin;				// minimum distance
	std::vector< bool >		exist;				// flag to show whether the point is selected

	if(size<pop.P().FSize()) size = pop.P().FSize();
	index.resize(size);

	//Step 1: calculate the distance matrix
	dismat.resize(pop.Size());
	dismin.resize(pop.Size());
	exist.resize(pop.Size());
	for(i=0; i<dismat.size(); i++) dismat[i].resize(pop.Size());
	for(i=0; i<dismat.size(); i++)
	{
		exist[i]	= false;
		dismin[i]	= 1.0E200;
		dismat[i][i]= 1.0E200;
		for(j=i+1; j<dismat.size(); j++)
		{
			dismat[i][j] = 0.0;
			for(k=0; k<pop.P().FSize(); k++) dismat[i][j] += (pop[i].F(k) - pop[j].F(k))*(pop[i].F(k) - pop[j].F(k));
			dismat[i][j] = sqrt(dismat[i][j]);
			dismat[j][i] = dismat[i][j];
		}
	}
	
	//Step 3: select the extream points, i.e the solutions with minimum Fi(x)
	for(k=0; k<pop.P().FSize(); k++)
	{
		fmin = 1.0E200; iindex = 0;
		for(i=0; i<dismat.size(); i++) 
			if(!exist[i] && pop[i].F(k)<fmin)
			{
				fmin  = pop[i].F(k);
				iindex= i;
			}
		exist[iindex] = true; index[k] = iindex; 
		for(i=0; i<dismat.size(); i++) if(!exist[i] && dismat[i][iindex]<dismin[i]) dismin[i] = dismat[i][iindex];
	}

	//Step 4: select point one by one
	for(i=pop.P().FSize(); i<size; i++)
	{
		//find the next one to be selected
		iindex = 0; fmin = -1.0E200;
		for(j=0; j<dismat.size(); j++) if(!exist[j] && dismin[j] > fmin) {fmin=dismin[j]; iindex=j;} 
		exist[iindex] = true; index[i] = iindex;

		//update the minimum distance of other points
		for(j=0; j<dismat.size(); j++) if(!exist[j] && dismin[j] > dismat[j][iindex]) dismin[j] = dismat[j][iindex];
	}
}

//convergence test
double ModelPCA::ConvergenceTest(CPopulationMO& pop)
{
	unsigned int i, j, k;
	double dis, mindis, convergence;
	
	//current position is empty, so the history pool is not full
	if(pvHis[mHisIndex]==0) 
	{
		//store current population into history pool
		pvHis[mHisIndex]	= new CPopulationMO(pop);
		//set a default test value
		convergence = 10.0*mConThreshold;
	}
	//the history pool is full, test the convergence
	else
	{
		convergence = 0.0;
		//onvergence test : average distance from current population to history population
		for(i=0; i<pop.Size(); i++)
		{
			mindis	= 1.0E200;
			for(j=0; j<(*pvHis[mHisIndex]).Size(); j++)
			{
				dis = 0.0; 
				//test in objective sapce
				for(k=0; k<pop.P().FSize(); k++) dis += (pop[i].F(k) - (*pvHis[mHisIndex])[j].F(k))*(pop[i].F(k) - (*pvHis[mHisIndex])[j].F(k));
				//test in decision space
				//for(k=0; k<pop.P().XSize(); k++) dis += (pop[i][k] - (*pvHis[mHisIndex])[j][k])*(pop[i][k] - (*pvHis[mHisIndex])[j][k]);
				if(dis < mindis) mindis = dis;
			}
			convergence += sqrt(mindis);
		}
		convergence /= double(pop.Size()*mDeltaT);
		
		//update the history pool
		(*pvHis[mHisIndex]) = pop;
	}

	//increase the pool index
	mHisIndex = (mHisIndex + 1) % pvHis.size();
	
	return convergence;
}

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az
