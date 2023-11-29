/*! \file	Model.h
	
	\brief	Estimation Distribution Models
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Nov.18 2006 create
*/

#ifndef	AZ_MODEL_H
#define	AZ_MODEL_H

#include <fstream>
#include <string>
#include <vector>

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

//!\brief	calcualte the eigenvalues and eigenvectors of a covariance matrix
//!\param	eva		sorted eigenvalues
//!\param	eve		sorted eigenvectors
//!\param	no		the main no-eigens are stored
//!\param	cov		the covariance matrix
//!\return	void
void Eigen(std::vector<double>& eva, std::vector< std::vector<double> >& eve, unsigned int no, std::vector< std::vector<double> >& cov);

//!\brief Model structure, used to store data
class Model
{
public:
	double			mErrTol;		//!< error tolerance
	unsigned int	mNClu,			//!< number of clusters
					mNX,			//!< number of trainning data
					mDX,			//!< dimension of trainning data
					mDLat,			//!< dimension of latent space
					mMaxIter,		//!< maximal trainning steps
					mIter;			//!< real trainning steps
	std::vector< unsigned int >			mvIndex,						//!< cluster index of each data
										mvNo;							//!< the number of points assigned to each cluster
	std::vector< std::vector<double> >	mvX,							//!< trainning data, each row is a data vector
										mvMean,							//!< center of each cluster
										mvEigenvalue,					//!< eigenvalue of ecah cluster
										mvProMin,						//!< minimal projection in primary dimensions
										mvProMax;						//!< maximal projection in primary dimensions
	std::vector< std::vector< std::vector<double> > >	mvEigenvector;	//!< eigenvector of each cluster

public:
	//!\brief	set parameters
	//!\param	nclu	number of cluster
	//!\param	nx		number of trainning data
	//!\param	dx		dimension of trainning data
	//!\param	dlat	dimension of latent space
	//!\param	maxiter	maximal trainning steps
	//!\param	errtol	error tolerance
	//!\return	void
	void Set(unsigned int nclu, unsigned int nx, unsigned int dx, unsigned int dlat, unsigned int maxiter=100, double errtol=1.0E-5);
	//!\brief	write model into file
	//!\param	file file name
	//!\return	void
	void Write(std::string file);

protected:
	//!\brief	calcualte the mean, eigenvalue and eigenvector of a given set of data
	//!\param	mean	mean value of the data
	//!\param	eva		sorted eigenvalues
	//!\param	eve		sorted eigenvectors
	//!\param	index	index of given set of data
	//!\return	void
	void Eigen(std::vector<double>& mean, std::vector<double>& eva, std::vector< std::vector<double> >& eve, std::vector< unsigned int >& index);
};// class Model

} //namespace alg

} //namespace az

#endif //AZ_MODEL_H
