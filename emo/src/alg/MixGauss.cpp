// MixGauss.cpp

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <float.h>

// which library is used to calculate the eigens

#ifdef AZ_GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_permutation.h>
#else
#include "alg/Matrix.h"
#endif

#include "alg/Random.h"
#include "alg/Kmeans.h"
#include "alg/MixGauss.h"

#if defined(WIN32)
    #define wxFinite(n) _finite(n)
#elif defined(_LINUX)
    #define wxFinite(n) finite(n)
#else
    #define wxFinite(n) ((n)==(n))
#endif

template< class T >
	void LOGM( unsigned int it, T x, unsigned int s )
{
	std::ofstream fhand;									
	if(s > 0 )	fhand.open("mix.txt",std::ios::app);		
	else		fhand.open("mix.txt");					
	fhand<<it<<"\t"<<x<<std::endl;							
	fhand.close();											
}																

namespace az
{
namespace alg
{

void MixGauss::Set(unsigned type, unsigned int nclu, unsigned int nx, unsigned int dx, unsigned int dlat, unsigned int maxiter, double errtol)
{
	unsigned int i,j;
	
	mType	= type;
	mNClu	= nclu;
	mNX		= nx;
	mDX		= dx; 
	mDLat   = dlat;
	mMaxIter= maxiter;
	mErrTol	= errtol;
	mIter	= 0;
	mLik	= 1.0E-50;

	mvX.resize(mNX);
	for(i=0; i<mNX; i++) mvX[i].resize(mDX);

	mvMean.resize(mNClu);
	for(i=0; i<mNClu; i++) mvMean[i].resize(mDX);

	mvPc.resize(mNClu);

	mvPxc.resize(mNClu); mvPcx.resize(mNClu);
	for(i=0; i<mNClu; i++) 
	{
		mvPxc[i].resize(mNX);
		mvPcx[i].resize(mNX);
	}
	
	switch(mType)
	{
	case 0:// UMG
		mvSigma.resize(mNClu);
		for(i=0; i<mNClu; i++) 
		{
			mvSigma[i].resize(mDX);
			for(j=0; j<mDX; j++) mvSigma[i][j].resize(1);
		}
		break;
	case 1:// MMG
		mvSigma.resize(mNClu);
		for(i=0; i<mNClu; i++) 
		{
			mvSigma[i].resize(mDX);
			for(j=0; j<mDX; j++) mvSigma[i][j].resize(mDX);
		}
		break;
	case 2://PPCA
		mvCov.resize(mNClu);
		mvLambda.resize(mNClu);
		mvU.resize(mNClu);
		for(i=0; i<mNClu; i++)
		{
			mvLambda[i].resize(mDLat);
			mvU[i].resize(mDLat);
			for(j=0; j<mDLat; j++) mvU[i][j].resize(mDX);
		}
		break;
	default:
		break;
	}
}

void MixGauss::Train()
{
	unsigned int	c,	// cluster index
					m,	// row index
					n;	// col index
	unsigned int	k;

	// Step 1: initialize by Kmeans
	Kmeans kmeans;
	kmeans.Set(mNClu, mNX, mDX, 1, mMaxIter);
	kmeans.mvX	= mvX;
	kmeans.Train();
	mvMean		= kmeans.mvMean;

	for(c=0; c<mNClu; c++) 
	{
		// compute probability P(c)
		mvPc[c] = kmeans.mvNo[c]/double(mNX);

		std::vector<unsigned int> index(kmeans.mvNo[c]);
		n=0;
		for(m=0; m<mNX; m++) if(kmeans.mvIndex[m]==c) index[n++] = m;

		switch(mType)
		{
		case 0: //UMG
			for(n=0; n<mDX; n++)
			{
				mvSigma[c][n][0] = 0.0;
				for(k=0; k<index.size(); k++) mvSigma[c][n][0] += (mvX[index[k]][n] - mvMean[c][n])*(mvX[index[k]][n] - mvMean[c][n]);
				mvSigma[c][n][0] /= double(kmeans.mvNo[c]-0.0);
			}
			break;
		case 1: //MMG
			// compute the covariance
			for(m=0; m<mDX; m++)
				for(n=m; n<mDX; n++)
				{
					mvSigma[c][m][n] = 0.0;
					for(k=0; k<index.size(); k++) mvSigma[c][m][n] += (mvX[index[k]][m] - mvMean[c][m])*(mvX[index[k]][n] - mvMean[c][n]);
					mvSigma[c][m][n] /= double(kmeans.mvNo[c]-0.0);
					mvSigma[c][n][m]  = mvSigma[c][m][n];
				}
			break;
		case 2: //PPCA
			mvCov[c] = 0.0;
			for(m=mDLat; m<mDX; m++)	mvCov[c]		+= fabs(kmeans.mvEigenvalue[c][m]); mvCov[c] /= double(mDX-mDLat);
			for(m=0; m<mDLat; m++)		mvLambda[c][m]   = kmeans.mvEigenvalue[c][m];
			for(m=0; m<mDLat; m++) for(n=0; n<mDX; n++) mvU[c][m][n] = kmeans.mvEigenvector[c][m][n];
			break;
		default:
			break;
		}
	}

	//no furter training for MMG
	if(mType == 1) return;

	// Step 2: trainning
	std::vector<double> pc(mNClu); std::vector< std::vector<double> > cov;
	if(mType == 2)
	{
		cov.resize(mDX);
		for(m=0; m<mDX; m++) cov[m].resize(mDX);
	}

	while(mIter++ < mMaxIter)
	{
		// compute probability P(x|c) and P(c|x)
		GaussPDF();

		// compute the likelihood
		double liktmp1=0.0,liktmp2;
		for(m=0; m<mNX; m++)
		{
			liktmp2  = 0.0;
			for(c=0; c<mNClu; c++) liktmp2 += mvPxc[c][m]*mvPc[c];
			if(!wxFinite(liktmp2) || liktmp2<1.0E-50)	liktmp2 = 1.0E-50;
			liktmp1 += log(liktmp2);
		}
		liktmp1 /= double(mNX);

		if(fabs(liktmp1/mLik-1.0) < mErrTol) break;

		mLik = liktmp1;
		
		// update Pc			
		double tpc = 0.0;
		for(c=0; c<mNClu; c++)
		{
			pc[c] = 0.0;
			for(m=0; m<mNX; m++) pc[c] += mvPcx[c][m];
			if(!wxFinite(pc[c]) || pc[c]<1.0E-50) pc[c] = 1.0E-50;
			tpc += pc[c];
		}
		for(c=0; c<mNClu; c++) { mvPc[c] = pc[c]/tpc; }// pc[c] *= mNX; }

		// update centers
		for(c=0; c<mNClu; c++)
		{
			for(n=0; n<mDX; n++)
			{
				mvMean[c][n]	= 0.0;
				for(m=0; m<mNX; m++) mvMean[c][n] += mvPcx[c][m] * mvX[m][n];
				mvMean[c][n]  /= pc[c];
			}
		}

		// update covariance
		switch(mType)
		{
		case 0: //UMG
			for(c=0; c<mNClu; c++)
			{
				for(n=0; n<mDX; n++)
				{
					mvSigma[c][n][0] = 0.0;
					for(m=0; m<mNX; m++) mvSigma[c][n][0] += (mvX[m][n] - mvMean[c][n])*(mvX[m][n] - mvMean[c][n])*mvPcx[c][m];
					mvSigma[c][n][0] /= pc[c];
					if(!wxFinite(mvSigma[c][n][0]) || mvSigma[c][n][0]<1.0E-50) mvSigma[c][n][0] = 1.0E-50;
				}
			}
			break;
		case 1: //MMG
			break;
		case 2: //PPCA
			for(c=0; c<mNClu; c++)
			{
				for(m=0; m<mDX; m++)
				{
					for(n=m; n<mDX; n++)
					{
						cov[m][n] = 0.0;
						for(k=0; k<mNX; k++) cov[m][n] += (mvX[k][m] - mvMean[c][m])*(mvX[k][n] - mvMean[c][n])*mvPcx[c][k];
						cov[m][n] /= pc[c];
						cov[n][m]  = cov[m][n];
					}
				}
				alg::Eigen(mvLambda[c], mvU[c], mDLat, cov);
			}
			break;
		default:
			break;
		}
	}// while(mIter++ < mMaxIter)
	Write("tmp.txt");
}

void MixGauss::Write(std::string file)
{
	unsigned int c,m,n;
	std::ofstream out(file.c_str());
	out<<std::scientific<<std::setprecision(20);
	
	out<<"Train Steps "<<mIter<<std::endl;
	
	out<<"Data"<<std::endl;
	for(m=0; m<mNX; m++) 
	{
		for(n=0; n<mDX; n++) out<<mvX[m][n]<<"\t";
		out<<std::endl;
	}

	out<<"==================P(c)================="<<std::endl;
	for(c=0; c<mNClu; c++) out<<mvPc[c]<<"\t"; out<<std::endl;

	out<<"=================P(x|c)================="<<std::endl;
	for(m=0; m<mNX; m++)
	{
		for(c=0; c<mNClu; c++) out<<mvPxc[c][m]<<"\t"; out<<std::endl;
	}

	out<<"=================P(c|x)================="<<std::endl;
	for(m=0; m<mNX; m++)
	{
		for(c=0; c<mNClu; c++) out<<mvPcx[c][m]<<"\t"; out<<std::endl;
	}

	for(c=0; c<mNClu; c++)
	{
		out<<std::endl<<"===========cluster "<<c<<"==========="<<std::endl;
		out<<"mean"<<std::endl;	for(n=0; n<mDX; n++) out<<mvMean[c][n]<<"\t";out<<std::endl;
		if(mType == 0 || mType == 1)
		{
			out<<"covariance"<<std::endl;	
			for(m=0; m<mDX; m++) { for(n=0; n<mvSigma[c][m].size(); n++) out<<mvSigma[c][m][n]<<"\t";out<<std::endl; }
		}
	}

	out.close();
}

void MixGauss::GaussPDF()
{
	const double PI = 3.1415926535897932385;
	unsigned int c, i,j,k,s;
	double tmp, cont;

	for(c=0; c<mNClu; c++)
	{
		switch(mType)
		{
		case 0:
			cont = 1.0;
			for(k=0; k<mDX; k++) cont *= mvSigma[c][k][0];
			cont = sqrt(pow(2*PI, (double)mDX)*(cont+1.0E-50));
			for(i=0; i<mNX; i++)
			{
				mvPxc[c][i] = 0.0;
				for(j=0; j<mDX; j++) mvPxc[c][i] += (mvX[i][j]-mvMean[c][j])*(mvX[i][j]-mvMean[c][j])/mvSigma[c][j][0];
				mvPxc[c][i] = exp(-0.5*mvPxc[c][i])/cont;	// Pxc
				if(!wxFinite(mvPxc[c][i])) mvPxc[c][i] = 0.0;
				mvPcx[c][i] = mvPxc[c][i] * mvPc[c];		// Pcx
			}
		case 1:
			break;
		case 2:	
			cont = mDX * log(2*PI) + mDX * log(mvCov[c]);
			for(i=0; i<mDLat; i++) cont -= log(1.0 - mvCov[c]/mvLambda[c][i]);
			for(i=0; i<mNX; i++)
			{
				mvPxc[c][i] = 0.0;
				for(j=0; j<mDX; j++)	mvPxc[c][i] += (mvX[i][j]-mvMean[c][j])*(mvX[i][j]-mvMean[c][j]);
				for(j=0; j<mDLat; j++)	
				{
					tmp = 0.0;
					for(k=0; k<mDX; k++) tmp += (mvX[i][j]-mvMean[c][j])*mvU[c][j][k];
					mvPxc[c][i] -= tmp*tmp*(1.0-mvCov[c]/mvLambda[c][j]);
				}
				mvPxc[c][i] /= mvCov[c];

				mvPxc[c][i] = exp(-0.5*(mvPxc[c][i]+cont));	// Pxc
				if(!wxFinite(mvPxc[c][i]))	mvPxc[c][i] = 0.0;
				if( mvPxc[c][i] > 1.0 )		mvPxc[c][i] = 1.0;
				mvPcx[c][i] = mvPxc[c][i] * mvPc[c];		// Pcx
			}
			break;
		default:
			break;
		}
	}

	// repire Pxc
	for(i=0; i<mNX; i++)
	{
		tmp = 0.0;
		for(c=0; c<mNClu; c++) tmp += mvPxc[c][i];
		if( tmp < 1.0E-10)
		{ 
			cont = 1.0E100; s = 0;
			switch(mType)
			{
				case 0:
					for(c=0; c<mNClu; c++) 
					{
						tmp = 0.0;
						for(j=0; j<mDX; j++) if(mvSigma[c][j][0] >0) tmp += (mvX[i][j]-mvMean[c][j])*(mvX[i][j]-mvMean[c][j])/mvSigma[c][j][0];
						if(tmp < cont) {cont = tmp; s = c;}
					}
					break;
				case 1:
					break;
				case 2:
					for(c=0; c<mNClu; c++) 
					{
						tmp = 0.0;
						for(j=0; j<mDLat; j++) for(k=0; k<mDX; k++) tmp += (mvX[i][j]-mvMean[c][j])*mvU[c][j][k];
						if(tmp < cont) {cont = tmp; s = c;}
					}
					break;
			}
			for(c=0; c<mNClu; c++) 
			{
				if(c==s) mvPxc[c][i] = 1.0 - (c-1.0)*1.0E-20; else mvPxc[c][i] = 1.0E-20;
				mvPcx[c][i] = mvPxc[c][i] * mvPc[c];
			}
		}
	}

	// repire and uniform Pcx
	for(c=0; c<mNClu; c++)
	{
		tmp = 0.0;
		for(i=0; i<mNX; i++) tmp += mvPcx[c][i];
		if( !wxFinite(tmp) || tmp <= 1.0E-50) for(i=0; i<mNX; i++) mvPcx[c][i] = 1.0/double(mNX); 
	}
	for(i=0; i<mNX; i++)
	{
		tmp = 0.0;
		for(c=0; c<mNClu; c++) tmp += mvPcx[c][i];
		if( !wxFinite(tmp) || tmp <= 1.0E-50) 
			for(c=0; c<mNClu; c++) mvPcx[c][i]  = 1.0/double(mNClu);
		else 
			for(c=0; c<mNClu; c++) mvPcx[c][i] /= tmp;
	}
}

//R*R' = Sigma
void Chol( std::vector< std::vector<double> >& R, std::vector< std::vector<double> >& sigma )
{
	unsigned int m,n;
	unsigned int mDX 	= (unsigned int)sigma.size();
	R.resize(mDX);
	for(m=0; m<mDX; m++) {R[m].resize(mDX);for(n=0; n<mDX; n++) R[m][n] = 0.0;}
#ifdef AZ_GSL
	gsl_matrix* ps = gsl_matrix_alloc(mDX,mDX);
	for(m=0; m<mDX; m++) for(n=0; n<mDX; n++) gsl_matrix_set(ps, m, n, sigma[m][n]);
	// close the default error handler
	gsl_set_error_handler_off();
	int err = gsl_linalg_cholesky_decomp(ps);
	if(err != 0)
	{
		for(m=0; m<mDX; m++) R[m][m] = sqrt(sigma[m][m]);
	}
	else
	{
		for(m=0; m<mDX; m++) for(n=0; n<=m; n++) 
		{
			R[m][n] = gsl_matrix_get(ps,m,n);
			if(!wxFinite(R[m][n])) R[m][n] = 0.0;
		}
	}
	gsl_matrix_free(ps);
#else
	az::alg::Matrix cov(mDX,mDX), r(mDX,mDX);
	for(m=0; m<mDX; m++) for(n=0; n<mDX; n++) cov(m,n)	= sigma[m][n];
	
	bool succes = alg::Cholesky(r, cov);

	if(succes)
	{
		for(m=0; m<mDX; m++) for(n=0; n<=m; n++) 
		{
			R[m][n] = r(m,n);
			if(!wxFinite(R[m][n])) R[m][n] = sqrt(fabs(sigma[m][n]));
		}
	}
	else
	{
		for(m=0; m<mDX; m++) R[m][m] = sqrt(sigma[m][m]);
	}
#endif
} 

} //namespace alg

} //namespace az
