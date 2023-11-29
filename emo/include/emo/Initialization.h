/*! \file	Initialization.h
	
	\brief	initialization strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Sep.27 2005 create
	\date	Aug.03 2006 add RandomSet
	\date   Mar.07 2008 redesign
*/
#ifndef AZ_INITIALIZATION_H
#define AZ_INITIALIZATION_H

#include "PopulationMO.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief initialization strategies
namespace ini
{
	//!\brief hybrid initialization strategy
	class Hybrid
	{
	protected:
		unsigned int		mEvas;		//!< evaluations
	public:
		Hybrid();

		~Hybrid();

		//!\brief	calculate the scalar obj
		//!\brief	obj value
		double obj(double* xy, int Dimension);

		//!\brief	get evaluations
		//!\return	evaluations
		inline unsigned int EvaTimes() {return mEvas;}

		//!\brief	initialize a population
		//!\param	pop population
		//!\param	size population size
		//!\return	population
		CPopulationMO& Initialize(CPopulationMO& pop, unsigned int size);	
	};//class Hybrid

	//!\brief random initialization strategy
	class Uniform
	{
	protected:
		unsigned int mEvas;	//!< evaluations
	public:
		//!\brief	get evaluations
		//!\return	evaluations
		inline unsigned int EvaTimes() {return mEvas;}

		//!\brief	initialize a population
		//!\param	pop population
		//!\param	size population size
		//!\return	population
		CPopulationMO& Initialize(CPopulationMO& pop, unsigned int size);	
	};//class Random

	//!\brief initialization with experiment design
	class LHC
	{
	protected:
		unsigned int mEvas;	//!< evaluations
	public:
		//!\brief	get evaluations
		//!\return	evaluations
		inline unsigned int EvaTimes() {return mEvas;}

		//!\brief	initialize a population
		//!\param	pop population
		//!\param	size population size
		//!\return	population
		CPopulationMO& Initialize(CPopulationMO& pop, unsigned int size);	
	};//class RandomSet

}//namespace ini

} //namespace mea

} //namespace az


#endif //AZ_INITIALIZATION_H
