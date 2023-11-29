/*! \file	Generator_LS.h

\brief	Local Search strategies

\author Aimin ZHOU
\author Department of Computer Science,
\author University of Essex, 
\author Colchester, CO4 3SQ, U.K
\author azhou@essex.ac.uk

\date	Apr.24 2006 create
*/

#ifndef AZ_GENERATOR_LS_H
#define AZ_GENERATOR_LS_H

#include "IndividualMO.h"
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

//!\brief local search strategies
namespace ls
{

	//!\brief base class
	class LSBase
	{
	protected:
		unsigned int mNoL,	//!< number of points to do local search
					 mNoE;	//!< maximum evaluations on each search point
	public:
		//!\brief	constractor
		LSBase();

		//!\brief	set parameters
		//!\param	noL	number of points to do local search
		//!\param	noE maximum evaluations on each search point	
		//!\return	void
		void Set( unsigned int noL, unsigned int noE);
	protected:
		//!\brief	select some points for local search
		//!\param	pop		reference population
		//!\param	start	start location in the reference population
		//!\param	end		end location in the reference population
		//!\param	size	the number of points to be selected
		//!\param	index	the index of the selected points
		//!\return	void	
		void SelectStartPoint(CPopulationMO& pop, unsigned int start, unsigned int end, unsigned int size, std::vector<unsigned int>& index);
	}; //class LSBase

	/////////////////////////////////////////////////////////////////////////////////////
	//!\brief random hill climb based on crossover & mutation
	class XHill:public LSBase
	{
	public:
		//!\brief	Generator
		//!\param	popnew	offspring population
		//!\param	pop		parent population
		//!\retun	offspring population
		CPopulationMO& Generate(CPopulationMO& popnew, CPopulationMO& pop);
	}; //class XHill

	/////////////////////////////////////////////////////////////////////////////////////
	//!\brief valley strategy
	class Valley:public LSBase
	{
	protected:
		unsigned int mNoG;	//!< number of points to build response surface model
	public:
		//!\brief	constractor
		Valley();

		//!\brief	Generator
		//!\param	popnew	offspring population
		//!\param	pop		parent population
		//!\retun	offspring population
		CPopulationMO& Generate(CPopulationMO& popnew, CPopulationMO& pop);
	}; //class XHill

} //namespace ls

} //namespace gen

} //namespace mea

} //namespace az

#endif //AZ_GENERATOR_LS_H
