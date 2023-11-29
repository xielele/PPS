// GTM.cpp

//#include "matlab.h"
//libmatlb.lib libmx.lib libmmfile.lib
//Matlab math library can solve A X = B easily, but the library can only be used by .EXE 

#include <iostream>
#include <cmath>
#include "alg/GTM.h"

namespace az
{

namespace alg
{

const double REALMAX= 1.7976931348623157e+308;
const double PI		= 3.1415926535898;

GTM::GTM()
{
	//mlfEnterNewContext(0, 0);
}

GTM::~GTM()
{
	//mlfRestorePreviousContext(0, 0);
}

//initialize 1D latent structure
//GTM_STP1.m
void GTM::Initialize1(Matrix& FI, Matrix& W, double& Beta, Matrix& T, unsigned int noLatentVar, unsigned int noBaseF, double s)
{
	unsigned int n;

	Matrix X(noLatentVar, 1), MU(noBaseF, 1);
	
	//assign latent variables
	for(n=0; n<noLatentVar; n++)
		X(n,0)	= -1.0 + double(n)*2.0/double(noLatentVar-1.0);
	//center points of base functions
	for(n=0; n<noBaseF; n++)		
		MU(n, 0)= (-1.0 + double(n)*2.0/double(noBaseF-1.0))*double(noBaseF)/double(noBaseF-1.0);
	//variance of base functions
	double sigma = s * (MU(1,0) - MU(0,0));
	
	//calculate FI
	GaussianGrid(FI, MU, X, sigma);

	//estimate the initial W and beta
	EstimateWB(W, Beta, T, X, FI);
}

//initialize 2D latent structure
//GTM_STP2.m
void GTM::Initialize2(Matrix& FI, Matrix& W, double& Beta, Matrix& T, unsigned int noLatentVar, unsigned int noBaseF, double s)
{
	unsigned int n,k,t;

	unsigned int Xdim = (unsigned int)(sqrt(double(noLatentVar+0.01)));
	unsigned int FIdim= (unsigned int)(sqrt(double(noBaseF+0.01)));

	Matrix X(Xdim*Xdim, 2), MU(FIdim*FIdim, 2);

	//assign latent variables
	t = 0;
	for(n=0; n<Xdim; n++)
		for(k=0; k<Xdim; k++)
		{
			X(t,0) = double(n)*2.0/double(Xdim-1.0)-1.0;
			X(t,1) = double(Xdim-1.0-k)*2.0/double(Xdim-1.0)-1.0;
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
	double sigma = s * (MU(1,1) - MU(2,1));
	
	//calculate FI
	GaussianGrid(FI, MU, X, sigma);

	//estimate the initial W and beta
	EstimateWB(W, Beta, T, X, FI);
}

//GTM_GBF.m
void GTM::GaussianGrid(Matrix& FI, Matrix& MU, Matrix& X, double& sigma)
{
	unsigned int n,k;
	double tmp;
	
	//calculate the squared distances between X and MU
	Matrix DISX_MU;
	DisGrid(DISX_MU, MU, X);
	
	//Calculate outputs of the Gaussian basis functions
	//\FI = exp\{-\frac{\|X-MU\|^2}{2\sigma^2}\}
	FI.Resize(DISX_MU.RowSize(), DISX_MU.ColSize()+1);
    tmp = -1.0/(2*sigma*sigma);
	for(n=0; n<DISX_MU.RowSize(); n++)
		for(k=0; k<DISX_MU.ColSize(); k++)
			FI(n,k) = exp(DISX_MU(n,k)*tmp);

	//Add bias basis function
	for( n=0; n<FI.RowSize(); n++ ) FI(n,FI.ColSize()-1) = 1.0;
}

//GTM_DIST
void GTM::DisGrid(Matrix& MDIS, Matrix& T, Matrix& Y)
{
	unsigned int n,k,d;
	double dist;

	MDIS.Resize(Y.RowSize(), T.RowSize());

  	for (n = 0; n < T.RowSize(); n++)
		for (k = 0; k < Y.RowSize(); k++)
		{
			MDIS(k, n) = 0.0;
			for (d = 0; d<Y.ColSize(); d++)
      		{
				dist = T(n, d)-Y(k, d);
				MDIS(k,n) += dist*dist;
			}
		}
}

void GTM::EstimateWB(Matrix& W, double& Beta, Matrix& T, Matrix& X, Matrix& FI)
{
	unsigned int n, k;
	double tmps;
	Matrix eVts,eVls,A,normX,B,trFI,inFI,tmp;
	
	//calculate the eigenvalues and eigenvectors
	PCA(eVts, eVls, T);

	//A = eVts(:,1:L)*diag(sqrt(eVls(1:L)));
	A.Resize(eVts.RowSize(), X.ColSize());
	for(k=0; k<A.ColSize(); k++)
	{
		tmps = sqrt(eVls(k,0));
		for(n=0; n<A.RowSize(); n++)
			A(n,k) = eVts(n,k)*tmps;
	}

	//normX = (X - ones(size(X))*diag(mean(X)))*diag(1./std(X));
	std::vector<double> mean,std;
	normX.Resize(X.RowSize(), X.ColSize());
	X.RowMean(mean); X.RowStd(std);
	for(n=0; n<normX.RowSize(); n++)
		for(k=0; k<normX.ColSize(); k++)
			normX(n,k) = (X(n,k)-mean[k])/std[k];

	//FI*W = normX*A'
	A.Trans();
	normX.Multiply(A,B);

	//Cholesky factorization
	Matrix A1,B1;
	trFI = FI; trFI.Trans();
	trFI.Multiply(FI,A1);
	trFI.Multiply(B,B1);
	//get W by Cholesky method
	CholeskySolve(W,A1,B1);
	if(W.RowSize()<1)//Cholesky failed
	{
		Matrix A2;
		pinv(A2,A1);
		A2.Multiply(B1, W);
	}

	//W(Mplus1,:) = mean(T);
	T.RowMean(mean);
	for(n=0; n<W.ColSize(); n++) W(W.RowSize()-1,n)=mean[n];
	
	//  interDistBeta = gtm_bi(FI*W);
	//  if (L < length(T(1,:))) 
	//    beta = min(interDistBeta,(1/eVls(L+1)));
	//  else
	//    beta = interDistBeta;

	FI.Multiply(W,Pro);
	Beta = EstimateB(Pro);

	if(X.ColSize() < T.ColSize()) 
	{
		tmps	= 1.0 / eVls(X.ColSize(), 0);
		if(tmps < Beta) Beta = tmps; 
	}
}

void GTM::PCA(Matrix& eVts, Matrix& eVls, Matrix& T)
{
	unsigned int i,dim = T.ColSize();

	Matrix cov;
	std::vector<double> vls(dim),mean;
	cov.Identity(dim);
	eVts.Identity(dim);
	eVls.Resize(dim,1);

	//Step 1: calculate the mean
	T.RowMean(mean);

	//Step 2: get the covariance matrix
	Matrix datatmp(T);
	datatmp.RowSub(mean);
	Matrix datatmp1(datatmp);
	datatmp.Trans();
	datatmp.Multiply(datatmp1, cov);
	cov.Divide(double( T.RowSize()>1 ? T.RowSize()-1.0 : T.RowSize() ));

	//Step 3: calculate the eigenvalue and eigenvector
	cov.Eig(vls, eVts);

	for(i=0; i<dim; i++) eVls(i,0) = vls[i];
}

//GTM_BI.m
double GTM::EstimateB(Matrix& Y)
{
	unsigned int n, k;
	double meanNN, tmpMin;
	Matrix YInterDist;

	//yInterDist  = gtm_dist(Y, Y) + diag(ones(length(Y(:,1)),1)*realmax);
	DisGrid(YInterDist, Y, Y);
	for(n=0; n<YInterDist.RowSize(); n++) YInterDist(n,n) += REALMAX;

	//meanNN = mean(min(yInterDist));
	meanNN 	= 0.0;
	n		= 0;
	for( k=0; k<YInterDist.ColSize(); k++ )
	{
		tmpMin = REALMAX;
		for( n=0; n<YInterDist.RowSize(); n++ ) 
			if(YInterDist(n,k) < tmpMin) tmpMin = YInterDist(n,k);
		meanNN += tmpMin;
	}
	meanNN /= (YInterDist.RowSize()+0.0);

	//beta        = (2/meanNN);
	return 2.0 / meanNN;
}

//gtm_trn.m
void GTM::Train(Matrix& W, double& Beta, Matrix& T,	Matrix& FI, unsigned int Cycles)
{
	unsigned int K,  N,  n, k, cycle;
	
	K = FI.RowSize();
	N = T.RowSize();

	GDist.Resize(K,N);
	GR.Resize(K,N);
	GMin.resize(N);
	GMax.resize(N);

	//FI_T = FI';
	Matrix FI_T(FI); FI_T.Trans();
	
	//calucluate the distance between Gaussian centers and sample points
	Matrix tmp;
	FI.Multiply(W,tmp);
	DisGrid(GDist,T,tmp);
	for(n=0; n<N; n++)
	{
		GMin[n] = REALMAX; GMax[n] = -REALMAX;
		for(k=0; k<K; k++)
		{
			if(GDist(k,n)<GMin[n]) GMin[n] = GDist(k,n);
			if(GDist(k,n)>GMax[n]) GMax[n] = GDist(k,n);
		}
	}

	Matrix sum(K,K),A,B;
	//training cycles
	for(cycle=0; cycle<Cycles; cycle++)
	{
		//calculate the log-likelihood, and update the likelihood matrix
		//double llh = LogLikeliHood(Beta, T.ColSize());
		LogLikeliHood(Beta, T.ColSize());
		
		//Equ(2.12) in "GTM:The Generative Topographic Mapping" to update W
		//A*W = B
		//A = full(FI'*spdiags(sum(gtmGlobalR')', 0, K, K)*FI);
		//B = FI'*R*FI
		//get R
		for(k=0; k<K; k++)
		{
			sum(k,k) = 0.0;
			for(n=0; n<N; n++) sum(k,k) += GR(k,n);
		}
		//get A
		FI_T.Multiply(sum,tmp);
		tmp.Multiply(FI,A);
		//get B
		FI_T.Multiply(GR,tmp);
		tmp.Multiply(T,B);
		//get W by Cholesky method
		CholeskySolve(W,A,B);
		if(W.RowSize()<1)//Cholesky failed
		{
			Matrix A1;
			pinv(A1,A);
			A1.Multiply(B, W);
		}

		//calucluate the distance between Gaussian centers and sample points
		FI.Multiply(W,Pro);
		DisGrid(GDist,T,Pro);
		for(n=0; n<N; n++)
		{
			GMin[n] = REALMAX; GMax[n] = -REALMAX;
			for(k=0; k<K; k++)
			{
				if(GDist(k,n)<GMin[n]) GMin[n] = GDist(k,n);
				if(GDist(k,n)>GMax[n]) GMax[n] = GDist(k,n);
			}
		}

		//update the variance
		//beta = ND / sum(sum(gtmGlobalDIST.*gtmGlobalR));
		Beta = 0.0;
		for( n=0; n<GR.RowSize(); n++ )
			for( k=0; k<GR.ColSize(); k++ ) Beta += GDist(n,k) * GR(n,k);
		Beta = double(T.RowSize()*T.ColSize())/Beta;
	}
}

//gtm_rspg.m
//calculate the log-likelihood
double GTM::LogLikeliHood(double beta, double D)	
{
	unsigned int N, K, n, k;
	double	result, tmp1, tmp2;

	K	= GDist.RowSize();	//Gaussian centre point number
	N	= GDist.ColSize();	//data number

	//In calculation mode > 0, the distances between Gaussian centres
	//and data points are shifted towards being centred around zero,  
	//w.r.t. the extreme (min- and max-) values.
	//Since the difference between any two distances is still the same, 
	//the responsabilities will be the same. The advantage is that 
	//the risk of responsabilities dropping below zero (in finite precision) 
	//in the exponentiation below, due to large distances, is decreased.
	//However, while we CAN calculate reliably with zero (0), we CAN'T 
	//calculate reliably with infinity (Inf). Hence, the shifting of distances 
	//must not be so large that the exponentiation yields infinity as result.
	std::vector<double> distCorr(N);
	for(n=0; n<N; n++)
	{
		tmp1	= (GMin[n] + GMax[n])/2.0;
		tmp2	= GMin[n] + 1400.0/beta;
		distCorr[n]= tmp1 < tmp2 ? tmp1:tmp2;
		for(k=0; k<K; k++) GDist(k, n) -= distCorr[n];
	}

	//Since the normalisation factor of the Gaussians is cancelled out
	//when normalising the responsabilities below (R = R*diag(1 ./ rSum))
	//it is omitted here. This, however, is corrected for when calculating
	//the log-likelihood further below.
	std::vector<double> rSum(N);
	for(n=0; n<N; n++)
	{
		rSum[n] = 0.0;
        for(k=0; k<K; k++) 
		{
			GR(k,n)  = exp(-beta*GDist(k,n)/2.0);
			rSum[n] += GR(k,n);
		}
		for(k=0; k<K; k++) GR(k,n) /= rSum[n];
	}
	
	//Equ(2.3) in "GTM:The Generative Topographic Mapping"
	result = 0.0;
	for(n=0; n<N; n++) result += log(fabs(rSum[n])) + distCorr[n]*(-beta/2.0);
	result += (N+0.0)*((D/2.0)*log(beta/(2.0*PI)) - log(K+0.0));

	return result;
}

} //namespace alg

} //namespace az
