//	GenModMIDEA.cpp

#include <ctime>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>
#include <float.h>
#include "emo/GenModMIDEA.h"

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

#define MASK			2147483647
#define PRIME			65539
#define SCALE			0.4656612875e-9

//check!
#define SAFESQRT(VALUE)						(VALUE <= 0 ? 0.0 : sqrt( VALUE ))
#define GETREALSELECTED(INDEX,MEMBER)		(selected[(MEMBER)][(INDEX)])

#define RANDOM01()			(( RandomSeed = ( (RandomSeed * PRIME) & MASK) ) * SCALE )
#define NEXTRANDOMNUMBER()	((int) (RANDOM01()*(RAND_MAX+1)))
#define RANDOM_MAX_VALUE	(RAND_MAX+1)
#define RANDOMNUMBER(MAX)	((long) (RANDOM01()*((double) MAX)))
#define CHANCEEVENT(CHANCE) (((CHANCE) == 0.0)?0:(RANDOM01()<CHANCE))
#define SETSEED(SEED)		{RandomSeed = (SEED); srand( (unsigned int) (SEED) );}
#define INITSEED()			{SETSEED((unsigned) time(NULL));}


MIDEA::MIDEA()
{
	clusters				= 0;
	clustersSize			= 0;
	selected				= 0; 
	objectiveValuesSelected	= 0;
	mu						= 0; 
	sigma					= 0;
	maximumNumberOfClusters	= 0;
	selsize					= 0;

	INITSEED();
}

//initialize the LPCA
void MIDEA::Set(unsigned int maxclu, unsigned int code, double exten, double threshold)
{
	mExtension = exten;
	if(code > 0) 
	{
		useClustersCodingSpace    = 1;	
		useClustersObjectiveSpace = 0;	
	}
	else
	{
		useClustersCodingSpace    = 0;	
		useClustersObjectiveSpace = 1;	
	}
	leaderThreshold			= threshold;
	maximumNumberOfClusters = maxclu;
}

//model-based generator
CPopulationMO& MIDEA::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	int i,j,cluster;	
	
	FreeSpace();

	popsize = selsize	= (int)popref.Size();
	stringlength		= popref.P().XSize();
	numberOfObjectives	= popref.P().FSize();

	objectiveValuesSelected = (double **) malloc( popsize*sizeof( double* ) );
	selected				= (double **) malloc( popsize*sizeof( double* ) );
	for(i=0; i<popsize; i++)
	{
		objectiveValuesSelected[i]	= (double *) malloc( numberOfObjectives*sizeof( double ) );
		selected[i]					= (double *) malloc( stringlength*sizeof( double ) );

		for(j=0; j<numberOfObjectives; j++) objectiveValuesSelected[i][j] = popref[i].F(j);
		for(j=0; j<stringlength; j++)		selected[i][j] = popref[i][j];
	}

	modelSelection( useClustersCodingSpace, useClustersObjectiveSpace, leaderThreshold, maximumNumberOfClusters );

	modelFitting();

#ifdef AZ_MODEL_OUT
	std::ofstream fhand("model.set");
	fhand<<"UMG"<<std::endl;
	fhand<<clustersUsed<<std::endl;
	for(int i = 0; i < clustersUsed; i++ )
	{
		fhand<<mu[i][0]<<"\t"<<mu[i][1]<<"\t"<<sigma[i][0]<<"\t"<<sigma[i][1]<<std::endl;
	}
	fhand.close();
#endif
	//sample
	popnew.Resize(sizenew);

	for(i=0; i<(int)sizenew; i++)
	{
		/* Select a factorization based on the mixing parameters */
		cluster = RANDOMNUMBER( clustersUsed );

		//== with extension ==//
		/* Sample from that factorization */
		for( j = stringlength-1; j >= 0; j-- )
			popnew[i][j] = sampleGaussian( mu[cluster][j], (1.0+mExtension)*sigma[cluster][j] );
	}

	FreeSpace();
	
	return popnew;
}

void MIDEA::FreeSpace()
{
	int i;

	if( objectiveValuesSelected )
		for( i = 0; i < selsize; i++ )
		{
			if( objectiveValuesSelected[i] )		free( objectiveValuesSelected[i] );
			if( selected[i] )						free( selected[i] );
		}
	if( objectiveValuesSelected )		free( objectiveValuesSelected );
    if( selected )						free( selected );

	if( mu )
		for( i = 0; i < maximumNumberOfClusters; i++ )
		{
			if( mu[i] )			free( mu[i] );
			if( sigma[i] )		free( sigma[i] );
			if( clusters[i] )	free( clusters[i] );
		}
	if( mu )			free( mu );
	if( sigma )			free( sigma );
	if( clustersSize )	free( clustersSize );
	if( clusters )		free( clusters );

	clusters				= 0;
	clustersSize			= 0;
	selected				= 0; 
	objectiveValuesSelected	= 0;
	mu						= 0; 
	sigma					= 0;
}

/*
* Model selection in this implementation of MIDEA
* involves using univariate factorizations in a
* multiple of clusters. The samples are clustered
* first. When no clustering should be used, a single
* cluster is created containing all of the samples.
*/
void MIDEA::modelSelection( int useClustersCodingSpace, int useClustersObjectiveSpace,
					double leaderThreshold, int maximumNumberOfClusters )
{
	clusterSelection( useClustersCodingSpace, useClustersObjectiveSpace, leaderThreshold, maximumNumberOfClusters );
}

/*
* Clustering. Clusters the samples into either
* a single cluster or into possibly multiple clusters
* using the BEND Leader algorithm.
*/
void MIDEA::clusterSelection( int useClustersCodingSpace, int useClustersObjectiveSpace,
					  double leaderThreshold, int maximumNumberOfClusters )
{
	/* If this is the first time, allocate memory for clusters and such */
	if( !clusters )
		allocateMemoryForGeneralStructures( maximumNumberOfClusters );

	if( useClustersCodingSpace == 0 && useClustersObjectiveSpace == 0 )
		singleCluster();
	else if( useClustersCodingSpace )
		BENDLeaderCodingSpace( leaderThreshold, maximumNumberOfClusters );
	else
		BENDLeaderObjectiveSpace( leaderThreshold, maximumNumberOfClusters );
}

/*
* Allocates memory for the clusters.
*/
void MIDEA::allocateMemoryForGeneralStructures( int maximumNumberOfClusters )
{
	int i;

	clustersUsedOverall = 0;
	clusters            = (int **) malloc( maximumNumberOfClusters*sizeof( int * ) );
	clustersSize        = (int *) malloc( maximumNumberOfClusters*sizeof( int ) );

	for( i = 0; i < maximumNumberOfClusters; i++ )
		clusters[i] = 0;
}

/*
* Generates a single cluster containing all of
* the samples in the sample vector.
*/
void MIDEA::singleCluster( void )
{
	int i;

	if( clustersUsedOverall == 0 )
	{
		clusters[0] = (int *) malloc( (2*popsize)*sizeof( int ) );
		clustersUsedOverall++;
	}

	clustersUsed = 1;
	for( i = 0; i < selsize; i++ )
		clusters[0][i] = i;
	clustersSize[0] = selsize;
}

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
void MIDEA::BENDLeaderCodingSpace( double leaderThreshold, int maximumNumberOfClusters )
{
	BENDLeader( leaderThreshold, maximumNumberOfClusters, selected, selsize, stringlength );
}

/*
* The BEND Leader algorithm in the objective space.
*/
void MIDEA::BENDLeaderObjectiveSpace( double leaderThreshold, int maximumNumberOfClusters )
{
	BENDLeader( leaderThreshold, maximumNumberOfClusters, objectiveValuesSelected,
		selsize, numberOfObjectives );
}

/*
* The BEND Leader algorithm.
*/
void MIDEA::BENDLeader( double leaderThreshold, int maximumNumberOfClusters, double **samples, int numberOfSamples, int dimensions )
{
	int    i, j, *Os, *Oc, qs, qc, cmin=0, *L;
	double value, dmin, distance, *sampleMin, *sampleMax, dist, dist2, maximumDistance;

	sampleMin = (double *) malloc( dimensions*sizeof( double ) );
	sampleMax = (double *) malloc( dimensions*sizeof( double ) );
	L         = (int *) malloc( numberOfSamples*sizeof( int ) );

	clustersUsed = 0;

	Os = (int *) malloc( numberOfSamples*sizeof( int ) );
	Oc = (int *) malloc( numberOfSamples*sizeof( int ) );

	/* Compute minimum and maximum among the samples in each dimension */
	for( i = 0; i < dimensions; i++ )
	{
		sampleMin[i] = samples[0][i];
		sampleMax[i] = samples[0][i];
		for( j = 1; j < numberOfSamples; j++ )
		{
			value = samples[j][i];
			if( value < sampleMin[i] )
				sampleMin[i] = value;
			if( value > sampleMax[i] )
				sampleMax[i] = value;
		}
	}

	/* Compute approximate maximum distance */
	maximumDistance = 0;
	for( j = 0; j < numberOfSamples; j++ )
	{
		dist2 = 0;
		for( i = 0; i < dimensions; i++ )
		{
			dist   = sampleMax[i] - sampleMin[i];
			value  = samples[j][i] - sampleMin[i];
			dist2 += dist == 0 ? 0 : (value*value)/(dist*dist);
		}
		dist2 = SAFESQRT( dist2 );
		if( dist2 > maximumDistance )
			maximumDistance = dist2;

		dist2 = 0;
		for( i = 0; i < dimensions; i++ )
		{
			dist   = sampleMax[i] - sampleMin[i];
			value  = samples[j][i] - sampleMax[i];
			dist2 += dist == 0 ? 0 : (value*value)/(dist*dist);
		}
		dist2 = SAFESQRT( dist2 );
		if( dist2 > maximumDistance )
			maximumDistance = dist2;
	}

	for( i = 0; i < numberOfSamples; i++ )
		Os[i] = i;

	/* Perform clustering by going over the samples in
	* a random order. */
	for( i = 0; i < numberOfSamples; i++ )
	{
		qs = RANDOMNUMBER( numberOfSamples - i );

		for( j = 0; j < clustersUsed; j++ )
			Oc[j] = j;

		dmin = leaderThreshold;

		for( j = 0; j < clustersUsed; j++ )
		{
			qc = RANDOMNUMBER( clustersUsed - j );

			distance = normalizedEuclideanDistance( Os[qs], L[Oc[qc]], sampleMin, sampleMax, samples, dimensions ) / maximumDistance;

			if( distance < dmin )
			{
				dmin = distance;
				cmin = Oc[qc];
			}

			Oc[qc] = Oc[clustersUsed-j-1];
		}

		if( clustersUsed == maximumNumberOfClusters )
		{
			dmin = 0.0;
			qc   = RANDOMNUMBER( clustersUsed - j );
			cmin = Oc[qc];
		}

		if( dmin < leaderThreshold )
		{
			clusters[cmin][clustersSize[cmin]] = Os[qs];
			clustersSize[cmin]++;
		}
		else
		{
			if( clustersUsed == clustersUsedOverall )
			{
				clusters[clustersUsed] = (int *) malloc( (2*popsize)*sizeof( int ) );
				clustersUsedOverall++;
			}

			L[clustersUsed]            = Os[qs];
			clusters[clustersUsed][0]  = Os[qs];
			clustersSize[clustersUsed] = 1;
			clustersUsed++;
		}

		Os[qs] = Os[numberOfSamples-i-1];
	}

	free( L );
	free( Os );
	free( Oc );
	free( sampleMin );
	free( sampleMax );
}

/*
* Returns the normalized Euclidean Distance between two samples.
* Normalization is done with respect tothe range of the
* samples in the solution vector in each dimension.
*/
double MIDEA::normalizedEuclideanDistance( int index1, int index2, double *sampleMin, double *sampleMax, double **samples, int dimensions )
{
	int    i;
	double value, dist, result;

	result = 0.0;
	for( i = 0; i < dimensions; i++ )
	{
		value   = samples[index1][i] - samples[index2][i];
		dist    = sampleMax[i] - sampleMin[i];
		result += dist == 0 ? 0 : (value*value)/(dist*dist);
	}
	result = SAFESQRT( result );

	return( result );
}

/*
* Estimates the final parameters of the conditional
* factorizations, including the ones required for sampling
* new solutions.
*/
void MIDEA::modelFitting()
{
	int i;

	for( i = 0; i < clustersUsed; i++ )
		if( clustersSize[i] > 0 )
			modelFittingSingleCluster( i );
}

/*
* Estimates the final parameters of a single conditional
* factorization in a single cluster.
*/
void MIDEA::modelFittingSingleCluster( int cluster )
{
	int    i, j;
	double value, value2;

	/* Allocate memory for new solution sampling parameters */
	if( !mu )
	{
		mu    = (double **) malloc( maximumNumberOfClusters*sizeof( double * ) );
		sigma = (double **) malloc( maximumNumberOfClusters*sizeof( double * ) );

		for( i = 0; i < maximumNumberOfClusters; i++ )
		{
			mu[i]    = NULL;
			sigma[i] = NULL;
		}
	}

	if( !mu[cluster] )
	{
		mu[cluster]    = (double *) malloc( stringlength*sizeof( double ) );
		sigma[cluster] = (double *) malloc( stringlength*sizeof( double ) );
	}

	for( i = 0; i < stringlength; i++ )
	{
		value = 1.0/((double) clustersSize[cluster]);

		mu[cluster][i] = 0.0;
		for( j = 0; j < clustersSize[cluster]; j++ )
			mu[cluster][i] += value*GETREALSELECTED(i,clusters[cluster][j]);

		sigma[cluster][i] = 0.0;
		for( j = 0; j < clustersSize[cluster]; j++ )
		{
			value2             = GETREALSELECTED(i,clusters[cluster][j]) - mu[cluster][i];
			sigma[cluster][i] += value*value2*value2;
		}

		sigma[cluster][i] = SAFESQRT( sigma[cluster][i] );
	}
}

/*
* Generates new solutions by sampling the
* univariate binary factorization mixture.
* The MIDEA framework is elitist in that it preserves
* the selected individuals in the population.
* The remainder of the population (n - tau*n) is
* filled with new solutions.
//*/
//void MIDEA::generateNewSolutions()
//{
//	int i, j;
//
//	/* Copy the top tau*n best samples to the population */
//	for( i = 0; i < selsize; i++ )
//		for( j = 0; j < stringlength; j++ )
//			SETREAL( GETREALSELECTED(j,i), j, i );
//	for( i = 0; i < selsize; i++ )
//	{
//		for( j = 0; j < numberOfObjectives; j++ )		objectiveValues[i][j]		= objectiveValuesSelected[i][j];
//		if( numberOfConstraints > 0 )
//			for( j = 0; j < numberOfConstraints; j++ )	constraintViolations[i][j]	= constraintViolationsSelected[i][j];
//	}
//
//	/* Fill the remainder of the n - tau*n positions in the
//	* offspring array. */
//	for( i = 0; i < offsize; i++ )
//		sampleSolutionFromMixture( i );
//}

/*
* Generates a single new solution. First, a cluster is
* selected randomly. Second, a new solution is actually
* generated by sampling the univariate factorization
* in that cluster.
*/
//void MIDEA::sampleSolutionFromMixture( int which  )
//{
//	int cluster;
//
//	/* Select a factorization based on the mixing parameters */
//	cluster = RANDOMNUMBER( clustersUsed );
//
//	/* Sample from that factorization */
//	sampleSolutionFromUnivariateFactorization( which, cluster );
//}

/*
* Generates a new solution from a univariate
* factorization in a specified cluster.
*/
//void MIDEA::sampleSolutionFromUnivariateFactorization( int which, int cluster )
//{
//	int i;
//
//	for( i = stringlength-1; i >= 0; i-- )
//		SETREALOFFSPRING( sampleGaussian( mu[cluster][i], sigma[cluster][i] ), i, which );
//}

/*
* Returns a random, Gaussian ("normally") distributed
* value with mean mu and standard deviation sigma.
*/
double MIDEA::sampleGaussian( double mu, double sigma )
{
	return( mu + sigma*nextGaussian() );
}

/*
* Returns a random, Gaussian ("normally") distributed
* value with mean 0 and standard deviation 1.
*/
double MIDEA::nextGaussian()
{
	double v1, v2, s, multiplier;

	if( haveNextNextGaussian )
	{
		haveNextNextGaussian = 0;

		return( nextNextGaussian );
	}
	else
	{
		do
		{
			v1 = 2 * (RANDOM01()) - 1;
			v2 = 2 * (RANDOM01()) - 1;
			s = v1 * v1 + v2 * v2;
		} while (s >= 1);

		multiplier           = SAFESQRT(-2 * log(s)/s);
		nextNextGaussian     = v2 * multiplier;
		haveNextNextGaussian = 1;

		return( v1 * multiplier );
	}
}

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az
