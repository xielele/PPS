/*! \file	GTM.h
	
	\brief	Generative Topgraphic Mapping
	the algorithm needs Matlab math libraries and it is modified from 
	the GTM Matlab toolbox(by Markus Svensen, 1996) 
	and the paper "GTM:The Generative Topographic Mapping"(C.M Bishop, M Svensen and C.K.I Williams)

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Mar.10 2005 create
	\date	Oct.19 2005 modify
*/

#ifndef	AZ_GTM_H
#define	AZ_GTM_H

#include <cmath>
#include "alg/Matrix.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{
	typedef alg::Matrix Matrix;
	
	//!\brief Generative Topographic Mapping(GTM)
	class GTM
	{
	protected:
		Matrix	GR,		//!< probability matrix 
				GDist,	//!< distance(between Gasussian centers and sample points) matrix
				Pro;	//!< project Pro = FI*W
		std::vector<double> GMin,	//!< minimum distance vector 
							GMax;	//!< maximum distance vector
	public:
		//!\brief	constructor
		//!\return	no
		GTM();

		//!\brief	destructor
		//!\return	no
		~GTM();

		//!\brief	initialize FI, W and Beta if latent structure is 1D
        //!\param	FI a matrix of the basis functions(output)
		//!\param	W initial weight matrix(output)
		//!\param	Beta initial value for the inverse variance(output)
		//!\param	T data matrix(input)
		//!\param	noLatentVar number of latent variable(input)
		//!\param	noBaseF number of base functions(input)
		//!\param	s the width of basis functions relative to the distance between two neighbouring basis function centres(input)
		//!\return	no
		void Initialize1(Matrix& FI, Matrix& W, double& Beta, Matrix& T, unsigned int noLatentVar, unsigned int noBaseF, double s);

		//!\brief	initialize FI, W and Beta if latent structure is 2D
        //!\param	FI a matrix of the basis functions(output)
		//!\param	W initial weight matrix(output)
		//!\param	Beta initial value for the inverse variance(output)
		//!\param	T data matrix(input)
		//!\param	noLatentVar number of latent variable(input)
		//!\param	noBaseF number of base functions(input)
		//!\param	s the width of basis functions relative to the distance between two neighbouring basis function centres(input)
		//!\return	no
		void Initialize2(Matrix& FI, Matrix& W, double& Beta, Matrix& T, unsigned int noLatentVar, unsigned int noBaseF, double s);


		//!\brief	a EM train step
		//!\param	W weight matrix(input-output)
		//!\param	Beta value for the inverse variance(input-output)
		//!\param	T data matrix(input)
		//!\param	FI matrix of the basis functions(input)
		//!\param	Cycles train cycles(input)
		//!\return	no
		void Train(Matrix& W, double& Beta, Matrix& T,	Matrix& FI, unsigned int Cycles);

		//!\brief	get project Pro = FI*W
		//!\return	project
		inline Matrix& Project(){return Pro;}
	protected:
		//!\brief	calculate the distance between T and Y
		//!\param	MIDS distance matrix(output)
		//!\param	T data matrix(input)
		//!\param	Y data matrix(input)
		//!\return	no
		void DisGrid(Matrix& MDIS, Matrix& T, Matrix& Y);

		//!\brief	calculates the output of Gaussian basis functions for a given set of input
		//!\param	FI matrix of the basis functions(output)
		//!\param	MU matrix containing the centers of the basis functions(input)
		//!\param	X latent variable matrix(input)
		//!\param	sigma a scalar giving the standard deviation of the radii-symmetric Gaussian basis functions(input)
		//!\return	no
		void GaussianGrid(Matrix& FI, Matrix& MU, Matrix& X, double& sigma);

		//!\brief	calculate the eigenvalues and eigenvectors of a matrix T'*T
		//!\param	eVts eigenvectors(output)
		//!\param	eVls eigenvalues(output)
		//!\param	T data matrix(input)
		//!\return	no
		void PCA(Matrix& eVts, Matrix& eVls, Matrix& T);

		//!\brief	estimate the inverse variance by a matrix
		//!\param	Y reference matrix(input)
		//!\return	inverse variance 
		double EstimateB(Matrix& Y);

		//!\brief	estimate weight and inverse variance 
		//!\param	W weight matrix(output)
		//!\param	Beta inverse variance(output)
		//!\param	T sample data points(input)
		//!\param	X latent variables(input)
		//!\param	FI matrix of the basis functions(input)
		void EstimateWB(Matrix& W, double& Beta, Matrix& T, Matrix& X, Matrix& FI);

		//!\brief	calculate the log-likelihood
		//!\param	beta inverse variance(input)
		//!\param	D data dimension(input)
		//!\return	the log-likelihood
		double LogLikeliHood(double beta, double D);	

		//!\breif	invert a matrix
		//!\param	mat matrix(input-output)
		//!\return	the inverted matrix
		Matrix& Invert(Matrix& mat);
	};

} //namespace alg

} //namespace az

#endif	//AZ_GTM_H
