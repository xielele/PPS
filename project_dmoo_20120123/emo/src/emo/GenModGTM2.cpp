/*! \file	Generator_Model_GTM1.cpp
	
	\brief	Evolutionary Aglorithm Generator with GTM (mapping from latent space to parameter space)
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Apr.10 2006 redesign
*/

#include <cmath>
#include <fstream>
#include <float.h>
#include "alg/Matrix.h"
#include "emo/GenMod.h"
//#include "emo/LogFile.h"

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

//constractor
ModelGTM2::ModelGTM2()
{
	mbWeight	= true;
}

//destractor
ModelGTM2::~ModelGTM2()
{
}

//initialize the GTM
void ModelGTM2::Set(unsigned int noLatent, unsigned int noBaseFun, unsigned int dimLatent, unsigned int trainsteps, double extension)
{
	mTrainSteps	= trainsteps;
	mExtension	= extension;
	mLatentDim	= dimLatent;
	mNoLatent	= noLatent;
	mNoBaseFun	= noBaseFun;
	if(mLatentDim==1 && mNoLatent % 2 == 0) mNoLatent++;
	if(mLatentDim==2 && mNoBaseFun < 4) mNoBaseFun = 4;
	mR			= 1.0E100;
}
	
//reset the initial weights of GTM
void ModelGTM2::Reset()
{
	mbWeight = true;
}

//build GTM model and sample new solutions
CPopulationMO& ModelGTM2::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j;
	
	Reset();

	//Step 1: clear the return population
	popnew.Resize(sizenew);
	
	//Step 2: assign new data
	if( popref.Size() != mDataSize ) 
	{
		
		mDataSize	= popref.Size();
		mDataDim	= popref.P().XSize();
		mT.Resize(mDataSize,mDataDim);
	}
	for( i=0; i<mDataSize; i++ ) for( j=0; j<mDataDim; j++ ) mT(i,j) = popref[i][j];
	
	//Step 3: reset the weights of GTM if necessary
	if( mbWeight )
	{
		if(mLatentDim==1)
			mGTM.Initialize1(mFI,mW,mBeta,mT,mNoLatent,mNoBaseFun,2.0);
		else
			mGTM.Initialize2(mFI,mW,mBeta,mT,mNoLatent,mNoBaseFun,2.0);
		mbWeight = false;
	}

	//Step 4: train GTM model
	mGTM.Train(mW,mBeta,mT,mFI,mTrainSteps);

	//Step 5: build Principal Curve(Surface) model and sample new solutions
    //1D
	if(mLatentDim==1)
		ModelGen1D(popnew, popref);
	//2D
	else ModelGen2D(popnew, popref);

	// Step 7: check the range
	for(i=0; i<mDataSize; i++)
		for(j=0; j<mDataDim; j++)
		{
			if( wxFinite(popnew[i][j]) == 0 ) popnew[i][j] = rnd::rand(popnew.P().XLow(j),popnew.P().XUpp(j));

			// border check strategy 1: DE strategy
			if(popnew[i][j] < popnew.P().XLow(j))		popnew[i][j] = 0.5*(popnew.P().XLow(j)+popref[i][j]);
			else if(popnew[i][j] > popnew.P().XUpp(j))	popnew[i][j] = 0.5*(popnew.P().XUpp(j)+popref[i][j]);
	}

	return popnew;
}

//mapping solutions from latent space to decision variable space
//Mar.07 2006
CPopulationMO& ModelGTM2::ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int n,k,d;
	double dist;

	alg::Matrix X(popref.Size(), 1), MU(mNoBaseFun, 1);
	
	//assign latent variables
	for(n=0; n<popref.Size(); n++)
		X(n,0) = rnd::rand(double(n),double(n+1))*(2.0+2*mExtension)/double(popref.Size())-1.0-mExtension;
	//center points of base functions
	for(n=0; n<mNoBaseFun; n++)		
		MU(n, 0)= (-1.0 + double(n)*2.0/double(mNoBaseFun-1.0))*double(mNoBaseFun)/double(mNoBaseFun-1.0);
	//variance of base functions
	double sigma = 2.0 * (MU(1,0) - MU(0,0));

	alg::Matrix MDIS(X.RowSize(), MU.RowSize());

  	for (n = 0; n < MU.RowSize(); n++)
		for (k = 0; k < X.RowSize(); k++)
		{
			MDIS(k, n) = 0.0;
			for (d = 0; d<X.ColSize(); d++)
      		{
				dist = MU(n, d)-X(k, d);
				MDIS(k,n) += dist*dist;
			}
		}

	alg::Matrix FI(MDIS.RowSize(), MDIS.ColSize()+1);
    double tmp = -1.0/(2*sigma*sigma);
	for(n=0; n<MDIS.RowSize(); n++)
		for(k=0; k<MDIS.ColSize(); k++)
			FI(n,k) = exp(MDIS(n,k)*tmp);

	//Add bias basis function
	for(n=0; n<FI.RowSize(); n++) FI(n,FI.ColSize()-1) = 1.0;

	//create new solutions
	//Gaussian noise variance
	double radius  = 1.0 / sqrt(mBeta); if(radius>mR) radius = mR; else mR = radius;
	for(n=0; n<popnew.Size(); n++)
	{
		for(d=0; d<popnew.P().XSize(); d++)
		{
			popnew[n][d] = radius*rnd::gaussian();
			for(k=0; k<FI.ColSize(); k++) popnew[n][d] += FI(n,k)*mW(k,d);
		}
	}

	return popnew;
}

CPopulationMO& ModelGTM2::ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int n,k,d,t;
	double dist;

	unsigned int Xdim = (unsigned int)(sqrt(double(popnew.Size()+0.01)));
	unsigned int FIdim= (unsigned int)(sqrt(double(mNoBaseFun+0.01)));

	alg::Matrix X(popnew.Size(), 2), MU(FIdim*FIdim, 2);

	//assign latent variables
	t = 0;
	for(n=0; n<Xdim; n++)
		for(k=0; k<Xdim; k++)
		{
			X(t,0) = rnd::rand(double(n),double(n+1))*(2.0+2*mExtension)/double(Xdim)-1.0-mExtension;
			X(t,1) = rnd::rand(double(k),double(k+1))*(2.0+2*mExtension)/double(Xdim)-1.0-mExtension;
			t++;
		}
	while(t<popnew.Size())
	{
		X(t,0) = rnd::rand(-1.0-mExtension,1.0+mExtension);
		X(t,1) = rnd::rand(-1.0-mExtension,1.0+mExtension);
		t++;
	}

	//center points of base functions
	t = 0;
	for(n=0; n<FIdim; n++)
		for(k=0; k<FIdim; k++)
		{
			MU(t,0) = (double(n)*2.0/double(FIdim-1.0)-1.0)*double(FIdim)/double(FIdim-1.0);
			MU(t,1) = (double(FIdim-1.0-k)*2.0/double(FIdim-1.0)-1.0)*double(FIdim)/double(FIdim-1.0);
			t++;
		}
	//variance of base functions
	double sigma = 2.0 * (MU(1,1) - MU(2,1));

	alg::Matrix MDIS(X.RowSize(), MU.RowSize());

  	for (n = 0; n < MU.RowSize(); n++)
		for (k = 0; k < X.RowSize(); k++)
		{
			MDIS(k, n) = 0.0;
			for (d = 0; d<X.ColSize(); d++)
      		{
				dist = MU(n, d)-X(k, d);
				MDIS(k,n) += dist*dist;
			}
		}

	alg::Matrix FI(MDIS.RowSize(), MDIS.ColSize()+1);
    double tmp = -1.0/(2*sigma*sigma);
	for(n=0; n<MDIS.RowSize(); n++)
		for(k=0; k<MDIS.ColSize(); k++)
			FI(n,k) = exp(MDIS(n,k)*tmp);

	//Add bias basis function
	for(n=0; n<FI.RowSize(); n++) FI(n,FI.ColSize()-1) = 1.0;

	//create new solutions
	//Gaussian noise variance
	double radius  = 1.0 / sqrt(mBeta); if(radius>mR) radius = mR; else mR = radius;
	for(n=0; n<popnew.Size(); n++)
	{
		for(d=0; d<popnew.P().XSize(); d++)
		{
			popnew[n][d] = radius*rnd::gaussian();
			for(k=0; k<FI.ColSize(); k++) popnew[n][d] += FI(n,k)*mW(k,d);
		}
	}

	return popnew;
}

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az
