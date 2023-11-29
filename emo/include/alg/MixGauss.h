/*! \file	MixGauss.h
	
	\brief	Mixture Gaussian model
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Nov.18 2006 create
*/

#ifndef	AZ_MIXGAUSS_H
#define	AZ_MIXGAUSS_H

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

//!\brief	Chol decomposition for a covariance matrix, used to sample Gaussian distribution random numbers, R*R' = Sigma
//!\param	R the output matrix
//!\param	sigma the covariance matrix
//!\return	void
void Chol( std::vector< std::vector<double> >& R, std::vector< std::vector<double> >& sigma );

//!\brief MixGaussian, mixture Gaussian model
class MixGauss
{
public:
	double			mErrTol,		//!< error tolerance
					mLik;			//!< likelihood
	unsigned int	mNClu,			//!< number of clusters
					mNX,			//!< number of trainning data
					mDX,			//!< dimension of trainning data
				    mDLat,			//!< dimension of latent sapce
					mMaxIter,		//!< maximal trainning steps
					mIter,			//!< real trainning steps
					mType;			//!< type of mixture Gaussian model: univariate or multivariate model
	std::vector< double >				mvPc;							//!< pripor probability P(c)
	std::vector< std::vector<double> >	mvX,							//!< trainning data, each row is a data vector
										mvMean,							//!< center of each cluster											
										mvPxc,							//!< probability P(x|c)
										mvPcx,							//!< probability P(c|x)
										mvLambda;						//!< eigenvalues for PPCA model
	std::vector< std::vector< std::vector<double> > >	mvSigma;		//!< covariance for Multivariant Model
	std::vector< double > 								mvCov;		    //!< covariance for Univariate model and PPCA model														
	std::vector< std::vector< std::vector<double> > > 	mvU;            //! primary eigenvectors in PPCA model

public:
	//!\brief	set parameters
	//!\param	type	model type
	//!\param	nclu	number of cluster
	//!\param	nx		number of trainning data
	//!\param	dx		dimension of trainning data
	//!\param   dlat    dimension of latent sapce
	//!\param	maxiter	maximal trainning steps
	//!\param	errtol	error tolerance
	//!\return	void
	void Set(unsigned type, unsigned int nclu, unsigned int nx, unsigned int dx, unsigned int dlat, unsigned int maxiter=100, double errtol=1.0E-5);

	//!\brief	train process
	//!\return	void
	void Train();

	//!\brief	write model into file
	//!\param	file file name
	//!\return	void
	void Write(std::string file);
protected:
	//!\brief	calcualte P(x|c) and P(c|x)
	//!\return	void
	void GaussPDF();
};//class MixGauss

} //namespace alg

} //namespace az

#endif //AZ_MIXGAUSS_H
