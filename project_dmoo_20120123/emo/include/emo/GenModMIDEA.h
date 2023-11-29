/*! \file	GenModMIDEA.h
	
	\brief	Model based generators for MOEA (MIDEA)
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Nov.30 2005 create
*/


#ifndef AZ_GENERATOR_MODEL_MIDEA_H
#define AZ_GENERATOR_MODEL_MIDEA_H

#include <vector>
#include "PopulationMO.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{

//!\brief	gen namespace, offspring generate strategies
namespace gen
{

//!\brief model based generator
namespace mod
{

//!\brief MIDEA
class MIDEA
{
protected:
	short	haveNextNextGaussian;

	int		**clusters,
			*clustersSize,
			popsize, 
			stringlength, 
			numberOfObjectives,
			maximumNumberOfClusters,
			RandomSeed,
			selsize, 
			useClustersCodingSpace,
			useClustersObjectiveSpace,
			clustersUsed, 
			clustersUsedOverall;

	double	**selected, 
			**objectiveValuesSelected,
			**mu, 
			**sigma,
			leaderThreshold,
			nextNextGaussian;
	double	mExtension;
public:
	MIDEA();

	//!\brief	set parameters
	//!\param	maxclu	maximum cluster
	//!\param	code	coding type
	//!\param	exten	extension
	//!\return	void
	void Set(unsigned int maxclu, unsigned int code, double exten, double threshold);

	//!\brief	Generator
	//!\param	sizenew	size of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	parent population
	//!\return	offspring population
	mea::CPopulationMO& Generate(unsigned int sizenew, mea::CPopulationMO& popnew, mea::CPopulationMO& popref);
protected:
	void FreeSpace();
	/*
	* Model selection in this implementation of MIDEA
	* involves using univariate factorizations in a
	* multiple of clusters. The samples are clustered
	* first. When no clustering should be used, a single
	* cluster is created containing all of the samples.
	*/
	void modelSelection( int useClustersCodingSpace, int useClustersObjectiveSpace,
						double leaderThreshold, int maximumNumberOfClusters );

	/*
	* Clustering. Clusters the samples into either
	* a single cluster or into possibly multiple clusters
	* using the BEND Leader algorithm.
	*/
	void clusterSelection( int useClustersCodingSpace, int useClustersObjectiveSpace,
						  double leaderThreshold, int maximumNumberOfClusters );

	/*
	* Allocates memory for the clusters.
	*/
	void allocateMemoryForGeneralStructures( int maximumNumberOfClusters );

	/*
	* Generates a single cluster containing all of
	* the samples in the sample vector.
	*/
	void singleCluster( void );

	/*
	* The BEND Leader algorithm.
	* This clustering algorithm is essentially the
	* Leader algorithm, using the BEND distance
	* (Bounding Box Euclidean Normalized Distance).
	* This distance is defined as Euclidean distance
	* normalized to the maximum and minimum distance
	* in each dimension, divided by an approximation
	* of the maximum distance in the sample set.
	* This approximation can be computed in linear
	* time with respect to the sample set by taking
	* the maximum distance over all points
	* to the minimal and maximal points in the
	* bounding box respectively. The leader algorithm
	* goes over the sample vector exactly once. For each
	* sample, it finds the cluster with the closest leader.
	* If this leader is closer than a given threshold, the
	* sample is added to that cluster. Otherwise, a new
	* cluster is created containing only this single
	* sample. More details are given here:
	*
	* Peter A.N. Bosman and Dirk Thierens. Advancing Continuous
	* IDEAs with Mixture Distributions and Factorization Selection
	* Metrics. In M. Pelikan and K. Sastry, organisers. Proceedings
	* of the Optimization by Building and Using Probabilistic Models
	* OBUPM Workshop at the Genetic and Evolutionary Computation
	* Conference - GECCO-2001, pages 208-212. 2001.
	*
	* This function clusters in the coding space.
	*/
	void BENDLeaderCodingSpace( double leaderThreshold, int maximumNumberOfClusters );

	/*
	* The BEND Leader algorithm in the objective space.
	*/
	void BENDLeaderObjectiveSpace( double leaderThreshold, int maximumNumberOfClusters );

	/*
	* The BEND Leader algorithm.
	*/
	void BENDLeader( double leaderThreshold, int maximumNumberOfClusters, double **samples, int numberOfSamples, int dimensions );

	/*
	* Returns the normalized Euclidean Distance between two samples.
	* Normalization is done with respect tothe range of the
	* samples in the solution vector in each dimension.
	*/
	double normalizedEuclideanDistance( int index1, int index2, double *sampleMin, double *sampleMax, double **samples, int dimensions );

	/*
	* Estimates the final parameters of the conditional
	* factorizations, including the ones required for sampling
	* new solutions.
	*/
	void modelFitting();

	/*
	* Estimates the final parameters of a single conditional
	* factorization in a single cluster.
	*/
	void modelFittingSingleCluster( int cluster );

	/*
	* Generates new solutions by sampling the
	* univariate binary factorization mixture.
	* The MIDEA framework is elitist in that it preserves
	* the selected individuals in the population.
	* The remainder of the population (n - tau*n) is
	* filled with new solutions.
	*/
	//void MIDEA::generateNewSolutions();

	/*
	* Generates a single new solution. First, a cluster is
	* selected randomly. Second, a new solution is actually
	* generated by sampling the univariate factorization
	* in that cluster.
	*/
	//void MIDEA::sampleSolutionFromMixture( int which  );

	/*
	* Generates a new solution from a univariate
	* factorization in a specified cluster.
	*/
	//void MIDEA::sampleSolutionFromUnivariateFactorization( int which, int cluster );

	/*
	* Returns a random, Gaussian ("normally") distributed
	* value with mean mu and standard deviation sigma.
	*/
	double sampleGaussian( double mu, double sigma );

	/*
	* Returns a random, Gaussian ("normally") distributed
	* value with mean 0 and standard deviation 1.
	*/
	double nextGaussian();

};//class MIDEA

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az

#endif //AZ_GENERATOR_MODEL_MIDEA_H
