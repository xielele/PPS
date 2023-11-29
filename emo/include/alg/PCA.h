/*! \file	PCA.h
	
	\brief	Principal Component Analysis(PCA)
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Dec.30 2004  create
	\date	Sep.25 2005 rewrite & reorganize structure
*/

#ifndef	AZ_PCA_H
#define	AZ_PCA_H

#include <iostream>
#include <string>
#include "alg/Matrix.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

	//!\brief PCA class
	//!\warning the data must be stored in a column matrix
	class PCA
	{
	protected:
		alg::Matrix*	pData,		//!< data matrix 
						mCov,			//!< colvariance matrix
						mEigenvector;	//!< eigenvectors
		alg::FVECTOR	mMean,			//!< data mean
						mEigenvalue;	//!< eigenvalue
	public:
		//!\brief	constructor
		//!\return	void
		PCA();

		//!\brief	constructor
		//!\param	data data matrix
		//!\return	void
		PCA(alg::Matrix& data);
		
		//!\brief	destructor
		//!\return	void
		~PCA();

		//!\brief	set data
		//!\param	data data matrix
		//!\return	data matrix reference
		inline alg::Matrix& Data(alg::Matrix& data) {pData=&data;return *pData;}

		//!\brief	get data matrix reference
		//!\return	data matrix reference
		inline alg::Matrix& Data() {return *pData;}

		//!\brief	get data mean vector
		//!\return	mean vector
		inline alg::FVECTOR& Mean() {return mMean;}
		
		//!\brief	get eigenvectors
		//!\return	eigenvector matrix
		inline alg::Matrix& Eigenvector() {return mEigenvector;}

		//!\brief	get eigenvalue vector
		//!\return	eigenvalue vector
		inline alg::FVECTOR& Eigenvalue() {return mEigenvalue;}	

		//!\brief	initialize colvariance matrix to be an identity
		//!\param	dim data dimension
		//!\return	void
		void Initialize(unsigned int dim);

		//!\brief	get the range of the projections in a dimension 
		//!\param	dim data dimension
		//!\param	min minimum projection
		//!\param	max maximum projection
		//!\return	void
		void PrimaryRange(unsigned int dim, double& min, double& max);
		
		//!\brief	calculate the mean, eigenvalue and eigenvectors 
		//!\return	void		
		void Train();

		//!\brief	write results into a stream
		//!\param	os output stream
		//!\param	pca PCA object
		//!\return	output stream
		friend std::ostream& operator<<(std::ostream& os, PCA& pca);

		//!\brief	write results into a stream
		//!\return	void
		void Write( std::ostream& os );

		//!\brief	write results into a file
		//!\param	name filename
		//!\return	void
		void Write( std::string& name );

		//!\brief	write results into a file
		//!\param	name filename
		//!\return	void
		void Write( const char* name );
	};

} //namespace alg

} //namespace az

#endif //AZ_PCA_H
