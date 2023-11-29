/*! \file	LocalPCA.h
	
	\brief	Local Principal Component Analysis (Local PCA) model
	
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
	\date	Nov.12 2006 rewrite
	\date	Nov.18 2006 reorganize
						combine Local PCA and Kmeans together
*/

#ifndef	AZ_LOCALPCA_H
#define	AZ_LOCALPCA_H

#include <vector>
#include "alg/Model.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{
//!\brief Local PCA, partion data into clusters 
class LocalPCA:public Model
{
protected:
	std::vector< std::vector< std::vector<double> > >	mvPI;	//!< matrix PI to each cluster
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
};//class LocalPCA

} //namespace alg

} //namespace az

#endif //AZ_LOCALPCA_H
