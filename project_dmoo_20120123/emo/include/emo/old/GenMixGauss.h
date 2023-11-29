/*! \file	GenMixGauss.h
	
	\brief	Mixture Gaussian Model-based generator
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Nov.18 2006 create
*/

#ifndef	AZ_GENMIXGAUSS_H
#define	AZ_GENMIXGAUSS_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <float.h>
#include "algorithm/Matrix.h"
#include "algorithm/Random.h"
#include "algorithm/MixGauss.h"
#include "PopulationMO.h"

#if defined(WIN32)
    #define wxFinite(n) _finite(n)
#elif defined(_LINUX)
    #define wxFinite(n) finite(n)
#else
    #define wxFinite(n) ((n)==(n))
#endif

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief	gen namespace, generators
namespace gen
{

//!\brief MixGauss, mixture Gaussian model-based generator
class MixGauss
{
protected:
	double			mErrTol;		//!< error tolerance
	unsigned int	mNClu,			//!< number of clusters
					mNX,			//!< number of trainning data
					mDX,			//!< dimension of trainning data
					mMaxIter;		//!< maximal trainning steps
public:
	//!\brief	set parameters
	//!\param	nclu	number of cluster
	//!\param	maxiter	maximal trainning steps
	//!\param	errtol	error tolerance
	//!\return	void
	void Set(unsigned int nclu, unsigned int maxiter=100, double errtol=1.0E-10)
	{
		mErrTol	= errtol;
		mNClu	= nclu;
		mMaxIter= maxiter;
		mNX		= 0;
		mDX		= 0;
	}

	//!\brief	Generator
	//!\param	popnew offspring population
	//!\param	popref parent population
	//!\return	offspring population
	mea::CPopulationMO& Generate(mea::CPopulationMO& popnew, mea::CPopulationMO& popref)
	{
		unsigned int c,m,n,k;
		mNX = popref.Size();
		mDX = popref.P().XSize();
		
		// model building
		alg::MixGauss mixGauss;
		mixGauss.Set(mNClu, mNX, mDX, mMaxIter, mErrTol);
		for(m=0; m<mNX; m++) mixGauss.mvX[m] = popref[m].X();
		mixGauss.Train();
			
		// set the point number of each cluster
		unsigned int nt = 0;
		std::vector<unsigned int> nc(mNClu);
		for(c=0; c<mNClu; c++) { nc[c] = 2+(unsigned int)((mNX-2*mNClu)*mixGauss.mvPc[c]); nt += nc[c];}
		while(nt<mNX)
		{
			c = rnd::rand((unsigned int)0, (unsigned int)mNClu);
			nc[c]++; nt++;
		}

		// sample new trial solutions
		popnew.Resize(mNX);
		nt = 0;
		std::vector<double> tmp(mDX);
		for(c=0; c<mNClu; c++)
		{
			//=================================================================================================
			// this part could be rewroten if using other method
			alg::Matrix cov(mDX,mDX), R(mDX,mDX);
			for(m=0; m<mDX; m++) for(n=0; n<mDX; n++) cov(m,n)	= mixGauss.mvSigma[c][m][n];
			alg::Cholesky(R, cov);
			
			std::ofstream file("tmp.txt");
			file<<cov<<std::endl<<R;
			file.close();

			for(m=nt; m<nt+nc[c]; m++)
			{
				for(n=0; n<mDX; n++) tmp[n] = rnd::gaussian();
				for(n=0; n<mDX; n++)
				{
					popnew[m][n] = mixGauss.mvMean[c][n];
					for(k=0; k<mDX; k++) popnew[m][n] += tmp[k]*R(k,n); 

					if( wxFinite(popnew[m][n]) == 0 ) popnew[m][n] = rnd::rand(popnew.P().XLow(n),popnew.P().XUpp(n));

					// border check strategy 1: DE strategy
					if(popnew[m][n] < popnew.P().XLow(n))		popnew[m][n] = 0.5*(popnew.P().XLow(n)+popref[m][n]);
					else if(popnew[m][n] > popnew.P().XUpp(n))	popnew[m][n] = 0.5*(popnew.P().XUpp(n)+popref[m][n]);
				}
			}
			//=================================================================================================
			nt += nc[c];
		}

		popnew.Write("popnew.txt");
		return popnew;
	}
};//class MixGaussian

} //namespace gen

} //namespace mea

} //namespace az

#endif //AZ_GENMIXGAUSS_H
