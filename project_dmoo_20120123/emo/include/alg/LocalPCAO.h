/*! \file	LocalPCAO.h
	
	\brief	Local Principal Component Analysis(Local PCA)
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Jan.05 2005 create
	\date	Mar.30 2005 modify
	\date	Apr.15 2005 add VolumeDim()
	\date	Aug.09 2005 rewrite
	\date	Sep.25 2005 rewrite & reorganize structure
*/

#ifndef	AZ_LOCALPCAO_H
#define	AZ_LOCALPCAO_H

#include <vector>
#include "alg/PCA.h"
#include "alg/Matrix.h"
#include "alg/QNewton.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

//!\brief a cluster structure
struct CLUSTER
{
	PCA				mPCA;	//!< PCA object
	alg::FVECTOR	mCore;	//!< reference vector of each cluster
	alg::Matrix		mPi,	//!< a matrix variable
					mData;	//!< data of each cluster
	alg::LINDEX		mIndex; //!< data index of each cluster
};	

//!\brief Local PCA, partion data into clusters 
//!\warning the data should be stored in column matrix
class LocalPCAO
{
protected:
	bool			mbReset		;	//!< whether to reset reference vector
	unsigned int	mDataSize	,	//!< data size
					mDataDim	,	//!< data dimension
					mLatentDim	,	//!< latent dimension
					mClusterNo	;	//!< cluster size
	double **pData;					//!< data pointer
	std::vector< CLUSTER* >	mClusters;	//!< clusters
	alg::FVECTOR	mMinDis,		//!< distance to the reference vector
					mUpp,			//!< data range
					mLow;			//!< data range
public:
	//!\brief construct function
	LocalPCAO();

	//!\brief destruct function
	~LocalPCAO();

	//!\brief	initialize 
	//!\param	datadim data dimension
	//!\param	latentdim latent variable dimension
	//!\param	maxcluster maximum cluster number
	//!\return	void
	void	Initialize(unsigned int datadim, unsigned int latentdim, unsigned int maxcluster);
	
	//!\brief	train process
	//!\param	steps train iteration steps
	//!\param	size data number
	//!\param	pdata data pointer
	//!\return	void
	void	Train(unsigned int steps, unsigned int size, double **pdata);

	//!\brief	get cluster number
	//!\return	cluster number
	inline unsigned int ClusterSize() {return mClusterNo;}

	//!\brief	get the latent variable dimension
	//!\return	latent variable dimension
	inline unsigned int LatentDim() {return mLatentDim;}
	
	//!\brief	get data index of index-th cluster
	//!\brief	data index of index-th cluster
	inline alg::LINDEX& DataIndex(unsigned int index) {return (*mClusters[index]).mIndex;}

	//!\brief	get the distance to the reference vector of index-th point
	//!\return	the distance to the reference vector of index-th point
	inline double	DisToCore(unsigned int index) {return mMinDis[index];}

	//!\brief	get the mean of index-th cluster
	//!\return	mean of index-th cluster
	inline alg::FVECTOR& Mean(unsigned int index) {return (*mClusters[index]).mPCA.Mean();}
	
	//!\brief	get the eigenvector of index-th cluster
	//!\return	eigenvector of index-th cluster
	inline alg::Matrix& Eigenvector(unsigned int index) {return (*mClusters[index]).mPCA.Eigenvector();}

	//!\brief	get the eigenvalue of index-th cluster
	//!\return	eigenvalue of index-th cluster
	inline alg::FVECTOR& Eigenvalue(unsigned int index) {return (*mClusters[index]).mPCA.Eigenvalue();}
	
	//!\brief functions needed by Qusi-Newton method
	double f(double* r);
	void   g(double* px, double* pg);
protected:
	//!\brief	partition data according to current reference vectors
	//!\return	void
	void  Partition();

	//!\brief	calculate the distance from a point to a cluster
	//!\param	data data index
	//!\param	clu cluster index
	//!\return	the distance
	double Distance(unsigned int data, unsigned int clu);

	//!\brief	calculate the distance from a point to a cluster using K-means
	//!\param	data data index
	//!\param	clu cluster index
	//!\return	the distance
	double Distance_Eucdis(unsigned int data, unsigned int clu);

	//!\brief	calculate the distance from a point to a cluster using Local PCA
	//!\param	data data index
	//!\param	clu cluster index
	//!\return	the distance
	double Distance_Recdis(unsigned int data, unsigned int clu);

	//!\brief	update reference vector
	//!\param	clu cluster index
	//!\return	successful or not
	bool Update(unsigned int clu);

	//!\brief	update reference vector using Local PCA
	//!\param	clu cluster index
	//!\return	successful or not
	bool OptimizeRefe(unsigned int clu);

	//!\brief	update reference vector using K-means
	//!\param	clu cluster index
	//!\return	successful or not
	bool RebuildCenter(unsigned int clu);
	
	//!\brief	recalculate the distance from points to their reference vectors
	void  CalDistance();
};


} //namespace alg

} //namespace az

#endif //AZ_LOCALPCAO_H
