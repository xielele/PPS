/*! \file	HCSampler.h
	
	\brief	samplers in hyper cube
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	12.11 2006 create
*/

#ifndef	AZ_HCSAMPLER_H
#define	AZ_HCSAMPLER_H

#include <vector>

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{
	//!\breif	Latin Hyper Cube design
	//!\param	rand	rand matrix within [low, upp]
	//!\param	low		lower range
	//!\param	upp		upper range
	//!\return	void
	void LHC(std::vector< std::vector<double> >& rand, std::vector<double>& low, std::vector<double>& upp);

	//!\breif	uniform design
	//!\param	rand	rand matrix within [low, upp]
	//!\param	low		lower range
	//!\param	upp		upper range
	//!\return	void
	void Uniform(std::vector< std::vector<double> >& rand, std::vector<double>& low, std::vector<double>& upp);

} //namespace alg

} //namespace az

#endif //AZ_HCSAMPLER_H
