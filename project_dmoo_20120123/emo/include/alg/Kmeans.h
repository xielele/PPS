/*! \file	Kmeans.h
	
	\brief	Kmeans model
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Nov.18 2006 create
*/

#ifndef	AZ_KMEANS_H
#define	AZ_KMEANS_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "alg/Model.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

//!\brief Local PCA, partion data into clusters 
class Kmeans:public Model
{
public:
	//!\brief	train process
	//!\return	void
	void Train();
protected:
	//!\brief	calculate the distance between data m to cluster c
	//!\param	m	datat index
	//!\param	c	cluster index
	//!\return	distance
	double Distance(unsigned int m, unsigned int c);
};//class Kmeans


} //namespace alg

} //namespace az

#endif //AZ_KMEANS_H
