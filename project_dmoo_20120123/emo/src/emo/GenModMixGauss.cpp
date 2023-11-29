// GenModMixGauss.cpp

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <float.h>
#include "alg/Matrix.h"
#include "alg/Random.h"
#include "alg/MixGauss.h"
#include "emo/PopulationMO.h"
#include "emo/GenMod.h"

#if defined(WIN32)
    #define wxFinite(n) _finite(n)
#elif defined(_LINUX)
    #define wxFinite(n) finite(n)
#else
    #define wxFinite(n) ((n)==(n))
#endif

template< class T >
	void LOGG( unsigned int it, T x, unsigned int s )
{
	std::ofstream fhand;									
	if(s > 0 )	fhand.open("gen.txt",std::ios::app);		
	else		fhand.open("gen.txt");					
	fhand<<it<<"\t"<<x<<std::endl;							
	fhand.close();											
}

namespace az
{
namespace mea
{
namespace gen
{
namespace mod
{

void MixGauss::Set(unsigned int type, unsigned int dlat, unsigned int nclu, double exten, unsigned int maxiter, double errtol)
{
	mType		= type;
	mDLat		= dlat;
	mErrTol		= errtol;
	mNClu		= nclu;
	mExtension	= exten;
	mMaxIter	= maxiter;
	mNX			= 0;
	mDX			= 0;
}


mea::CPopulationMO& MixGauss::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int c,m,n,k;

	popnew.Resize(sizenew);

	mNX = popref.Size();
	mDX = popref.P().XSize();
	
	// model building
	alg::MixGauss mixGauss;
	mixGauss.Set(mType, mNClu, mNX, mDX, mDLat, mMaxIter, mErrTol);
	for(m=0; m<mNX; m++) mixGauss.mvX[m] = popref[m].X();
	mixGauss.Train();
	
	// set the point number of each cluster
	unsigned int nt = 0;
	std::vector<unsigned int> nc(mNClu);
	std::vector<double> pc(mNClu);double pct = 0.0;
	for(c=0; c<mNClu; c++) { pc[c] = mixGauss.mvPc[c]; pct += pc[c];}
	for(c=0; c<mNClu; c++) { nc[c] = 2+(unsigned int)((sizenew-2*mNClu)*pc[c]/pct); nt += nc[c];}
	while(nt>sizenew)
	{
		c = rnd::rand((unsigned int)0, (unsigned int)mNClu);
		if(nc[c]>2){nc[c]++; nt--;}
	}
	while(nt<sizenew)
	{
		c = rnd::rand((unsigned int)0, (unsigned int)mNClu);
		nc[c]++; nt++;
	}

#ifdef AZ_MODEL_OUT
	std::ofstream fhand("model.set");
	switch(mType)
	{
		case 0: //UMG
			fhand<<"UMG"<<std::endl;
			break;
		case 1:
		case 2:
			fhand<<"MMG"<<std::endl;
			break;
		default:
			break;
	}
	fhand<<mNClu<<std::endl;
#endif

	// sample new trial solutions
	nt = 0;
	std::vector<double> tmp(mDX);
	std::vector< std::vector<double> > cholR, cov;
	for(c=0; c<mNClu; c++)
	{
		switch(mType)
		{
		case 0: //UMG
			for(n=0; n<mDX; n++) tmp[n] = sqrt(fabs(mixGauss.mvSigma[c][n][0]+1.0E-50));
			for(m=nt; m<nt+nc[c]; m++)
			{
				for(n=0; n<mDX; n++)
				{
					popnew[m][n] = mixGauss.mvMean[c][n] + (1.0+mExtension)*tmp[n]*rnd::gaussian();

					if( wxFinite(popnew[m][n]) == 0 ) popnew[m][n] = rnd::rand(popnew.P().XLow(n),popnew.P().XUpp(n));

					// border check strategy 1: DE strategy
					if(popnew[m][n] < popnew.P().XLow(n))		popnew[m][n] = 0.5*(popnew.P().XLow(n)+popref[m][n]);
					else if(popnew[m][n] > popnew.P().XUpp(n))	popnew[m][n] = 0.5*(popnew.P().XUpp(n)+popref[m][n]);
				}
			}
#ifdef AZ_MODEL_OUT
			fhand<<mixGauss.mvMean[c][0]<<"\t"<<mixGauss.mvMean[c][1]<<"\t"<<tmp[0]<<"\t"<<tmp[1]<<std::endl;
#endif
			break;
		case 1: //MMG
			//Cholesky decomposition
			alg::Chol(cholR, mixGauss.mvSigma[c]);

			for(m=nt; m<nt+nc[c]; m++)
			{
				for(n=0; n<mDX; n++) tmp[n] = (1.0+mExtension)*rnd::gaussian();
				for(n=0; n<mDX; n++)
				{
					popnew[m][n] = mixGauss.mvMean[c][n];
					for(k=0; k<mDX; k++) popnew[m][n] += tmp[k]*cholR[n][k]; 

					if( wxFinite(popnew[m][n]) == 0 ) popnew[m][n] = rnd::rand(popnew.P().XLow(n),popnew.P().XUpp(n));

					// border check strategy 1: DE strategy
					if(popnew[m][n] < popnew.P().XLow(n))		popnew[m][n] = 0.5*(popnew.P().XLow(n)+popref[m][n]);
					else if(popnew[m][n] > popnew.P().XUpp(n))	popnew[m][n] = 0.5*(popnew.P().XUpp(n)+popref[m][n]);
				}
			}
#ifdef AZ_MODEL_OUT
			fhand<<mixGauss.mvMean[c][0]<<"\t"<<mixGauss.mvMean[c][1]<<"\t"<<cholR[0][0]<<"\t"<<0<<"\t"<<cholR[1][0]<<"\t"<<cholR[1][1]<<std::endl;
#endif
			break;
		case 2: //PPCA
			//find the covariance matrix
			cov.resize(mDX); for(m=0; m<mDX; m++) cov[m].resize(mDX);
			for(m=0; m<mDX;m++)
			{
				for(n=m; n<mDX; n++)
				{
					cov[m][n] = (m==n) ? mixGauss.mvCov[c] : 0.0;
					for(k=0; k<mDLat; k++) cov[m][n] += mixGauss.mvU[c][k][m]*mixGauss.mvU[c][k][n];
					cov[n][m] = cov[m][n];
				}
			}

			//Cholesky decomposition
			alg::Chol(cholR, cov);

			for(m=nt; m<nt+nc[c]; m++)
			{
				for(n=0; n<mDX; n++) tmp[n] = (1.0+mExtension)*rnd::gaussian();
				for(n=0; n<mDX; n++)
				{
					popnew[m][n] = mixGauss.mvMean[c][n];
					for(k=0; k<mDX; k++) popnew[m][n] += tmp[k]*cholR[n][k]; 

					if( wxFinite(popnew[m][n]) == 0 ) popnew[m][n] = rnd::rand(popnew.P().XLow(n),popnew.P().XUpp(n));

					// border check strategy 1: DE strategy
					if(popnew[m][n] < popnew.P().XLow(n))		popnew[m][n] = 0.5*(popnew.P().XLow(n)+popref[m][n]);
					else if(popnew[m][n] > popnew.P().XUpp(n))	popnew[m][n] = 0.5*(popnew.P().XUpp(n)+popref[m][n]);
				}
			}
#ifdef AZ_MODEL_OUT
			fhand<<mixGauss.mvMean[c][0]<<"\t"<<mixGauss.mvMean[c][1]<<"\t"<<cholR[0][0]<<"\t"<<0<<"\t"<<cholR[1][0]<<"\t"<<cholR[1][1]<<std::endl;
#endif
			break;
		default:
			break;
		}
		nt += nc[c];
	}

#ifdef AZ_MODEL_OUT
	fhand.close();
#endif

	return popnew;
}

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az
