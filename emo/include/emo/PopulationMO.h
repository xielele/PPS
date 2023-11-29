/*! \file	PopulationMO.h
	
	\brief	population class for MOEA
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Dec.03 2004 create
	\date	Dec.14 2004 redesign to a template class
	\date	Feb.06 2005 change
	\date	Apr.07 2005 change
	\date	Sep.27 2005 rewrite & reorganize structure
	\date	Apr.27 2006 add CPopulationMO(), P(CParameter)
	\date	May.09 2007 add Union operator
	\date	Oct.15 2007 add Divide operator
	\date	Mar.21 2008 change rank sort
	\date	Apr.30 2008 add DominateSort
*/

#ifndef	AZ_POPULATIONMO_H
#define	AZ_POPULATIONMO_H

#include <iostream>
#include <iomanip>
#include <vector>
#include "Parameter.h"
#include "IndividualMO.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//! \brief clear a pointer sequence
//! \param s reference to a swquence
template<class SEQ> void CLEAR(SEQ& s)
{
	typename SEQ::iterator i;
	for(i = s.begin(); i != s.end(); ++i) if(*i) delete *i;
	s.clear();
}

//!\brief population class for MOEA
class CPopulationMO
{
public:
	typedef CIndividualMO IND_TYPE;
protected:
	bool mbSort;					//!< whether the population has been sorted
	CParameter* pPar;				//!< pointer to the parameter object
	std::vector< IND_TYPE* > mvPop;	//!< the population vector stored individual pointers
public:
	//!\brief	constructor
	//!\param	par parameter objective
	//!\return	void
	CPopulationMO(CParameter& par);	

	//!\brief	constructor
	//!\param	pop another population
	//!\return	void
	CPopulationMO(CPopulationMO& pop);

	//!\brief	deconstrutor
	//!\return	void
	~CPopulationMO();	

	//!\brief	set parameters
	//!\param	par parameter object
	//!\return	void
	void P(CParameter& par);

	//!\brief	get the parameter object reference
	//!\return	the parameter object reference
	inline CParameter& P() {return *pPar;}
	
	//!\brief	see whether the population has been sorted
	//!\return	success if it has been sorted
	inline bool IsSort() {return mbSort;}

	//!\brief	set the sort state
	//!\param	sort new state
	//!\return	success if it has been sorted
	inline bool IsSort(bool sort ) {mbSort=sort; return mbSort;}
	
	//!\brief	get the population size
	//!\return	population size
	inline unsigned int	Size() {return (unsigned int)mvPop.size();}

	//!\brief	get the index-th individual reference
	//!\param	i individual index
	//!\return	individual reference
	inline IND_TYPE& operator[](unsigned int i) {return *mvPop[i];}

	//!\brief	get the index-th individual reference
	//!\param	i individual index
	//!\return	individual reference
	inline IND_TYPE& In(unsigned int i) {return *mvPop[i];}

	//!\brief	get the index-th individual pointer
	//!\param	i individual index
	//!\return	individual pointer
	inline IND_TYPE*& At(unsigned int i) {return mvPop[i];}

	//!\brief	see whether the individual is in the population
	//!\param	ind individual
	//!\return	success if the individual is in the population
	bool IsContain(IND_TYPE& ind);

	//!\brief	clear the population
	//!\return	void
	void Clear();

	//!\brief	evaluate the population
	//!\return	void
	void Evaluate();

	//!\brief	shuffle the population
	//!\return	void
	void Shuffle();	

	//!\brief	resize the population
	//!\param	s new size
	//!\return	void
	void Resize(unsigned int s);

	//!\brief	swap two individuals
	//!\param	i individual index
	//!\param	j individual index
	//!\return	void
	void Swap(unsigned int i, unsigned int j);
	
	//!\brief	erase individuals
	//!\param	i start individual index
	//!\return	void
	void Erase(unsigned int i);

	//!\brief	copy an individual to population
	//!\param	pind pointer to an individual
	//!\return	population
    CPopulationMO& Copy(IND_TYPE*& pind);

	//!\brief	copy an individual to population
	//!\param	ind reference to an individual
	//!\return	population
	CPopulationMO& Copy(IND_TYPE& ind);

	//!\brief	copy a population to population
	//!\param	pop reference to a population
	//!\return	population
	CPopulationMO& Copy(CPopulationMO& pop);

	//!\brief	combine an individual to population
	//!\param	pind pointer to an individual
	//!\return	population
	CPopulationMO& Combine(IND_TYPE*& pind);

	//!\brief	combine an individual to population
	//!\param	ind peference to an individual
	//!\return	population
	CPopulationMO& Combine(IND_TYPE& ind);

	//!\brief	combine a population to population
	//!\param	pop reference to a population
	//!\return	population
	CPopulationMO& Combine(CPopulationMO& pop);

	//!\brief	find the union of two sets (populations)
	//!\param	popA reference to a population
	//!\param	popB reference to a population
	//!\return	union of the two sets
	CPopulationMO& Unite(CPopulationMO& popA, CPopulationMO& popB);

	//!\brief	find the difference of two sets (populations)
	//!\param	popA reference to a population
	//!\param	popB reference to a population
	//!\return	difference of the two sets
	CPopulationMO& Minus(CPopulationMO& popA, CPopulationMO& popB);

	//!\brief	assign a population 
	//!\param	pop reference to a population
	//!\return	population
	CPopulationMO& operator=(CPopulationMO& pop);

	//!\brief	get a sub-population with Rank = r
	//!\param	pop sub-population
	//!\param	r rank
	//!\return	sub-population
	CPopulationMO& RankSub(CPopulationMO& pop, unsigned int r);

	//!\brief	get the maximum rank value
	//!\return	maximum rank value
	unsigned int RankMax();

	//!\brief	get the sub-population size with Rank=r
	//!\param	r rank value
	//!\param	start start index
	//!\param	end end index
	//!\return	sub-population size
	unsigned int RankSize(unsigned int r, unsigned int& start, unsigned int& end);

	//!\brief	assign rank value and sort population by general dominance
	//!\return	void
	void RankSort();

	//!\brief	assign rank value and sort population by e-dominance
	//!\return	void
	void ERankSort();

	//!\brief	sort the population to nondominated and dominated part
	//!\return	unsigned int the number of nondominated solutions
	unsigned int DominateSort();

	//!\brief	assign rank value and sort population by strategies
	//!\param	esort e-dominance or not
	//!\return	void
	void RankSort(bool esort);

	//!\brief	write the population to I/O stream
	//!\param	os output stream
	//!\return	output stream
	std::ostream& Write(std::ostream& os);

	//!\brief	write the population to a file
	//!\param	name file name
	//!\return	void
	void Write(std::string name);

	//!\brief	write the population to a file
	//!\param	name file name
	//!\return	void
	void Write(const char *name);

	//!\brief	write the population to I/O stream
	//!\param	os output stream
	//!\param	pop population
	//!\return	output stream
	friend std::ostream& operator<<(std::ostream& os, CPopulationMO& pop);

	//!\brief	read population from I/O stream
	//!\param	is input stream
	//!\return	input stream
	std::istream& Read(std::istream& is);

	//!\brief	read population from a file
	//!\param	name file name
	//!\return	void
	void Read(std::string name);

	//!\brief	read population from a file
	//!\param	name file name
	//!\return	void
	void Read(const char *name);

	//!\brief	read population from I/O stream
	//!\param	is input stream
	//!\param	pop population
	//!\return	void
	friend std::istream& operator>>(std::istream& is, CPopulationMO& pop);
};//class CPopulationMO

} //namespace mea

} //namespace az

#endif //AZ_POPULATIONMO_H
