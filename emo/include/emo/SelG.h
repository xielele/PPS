/*! \file	SelG.h
	
	\brief	Guided selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Mar.25 2008 create
*/
#ifndef AZ_GUIDED_SELECTION_H
#define AZ_GUIDED_SELECTION_H

#include <vector>
#include <list>
#include "IndividualMO.h"
#include "PopulationMO.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief select strategies
namespace gsel
{
//!\brief 0 order selection
class Order0
{
protected:
	unsigned int  mPopSize;						//!< population size
	CPopulationMO *pPop;						//!< populaton
	std::vector< double >	mIdea;				//!< idea point	
	std::vector< std::vector<double> > mWeight;	//!< weight vectors
public:
	//!\brief	constructor
	//!\return	void
	Order0();

	//!\brief	destructor
	//!\return	void
	~Order0();

	//!\brief	initialize the selector
	//!\param	pop initial population
	//!\return	void
	void Init(CPopulationMO& pop);

	//!\brief	select from current population and offspring population
	//!\param	off offspring population
	//!\return	population
	CPopulationMO& Select(CPopulationMO& off);
protected:
	//!\brief	calculate the fitness according to ith sub-problem
	//!\param	ind individual
	//!\param	ith	sub-problem index
	//!\return	double
	double SopFit(CIndividualMO& ind, unsigned int ith);

	//!\brief	update the idea point
	//!\param	off offspring population
	//!\return	void
	void UpdateIdea(CPopulationMO& off);

	//!\brief	update the the population
	//!\param	off		offspring population
	//!\param	state	initialize the population (false) or update population (true)
	//!\return	void
	void UpdatePop(CPopulationMO& off, bool state);
};//class SCrowd

}//namespace gsel

} //namespace mea

} //namespace az


#endif //AZ_GUIDED_SELECTION_H
