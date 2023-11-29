/*! \file	Gen.h
	
	\brief	Evolutionary Aglorithm Generators
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Dec.08 2004 create
	\date	Sep.26 2005	redesign and rewrite
	\date	Oct.14 2005 add DE operator
	\date	Nov.28 2005 add XSimplex
	\date	Nov.08 2005 add ESimplex
*/

#ifndef AZ_GENERATOR_EA_H
#define AZ_GENERATOR_EA_H

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

//!\brief	Polynomial mutation operator (PM)
//\param	ind individual to be mutated
//\return	indivdual
CIndividualMO& PM(CIndividualMO& ind);

//!\brief	Simulated Binary crossover (SBX)
//\param	son1 offspring
//\param	son2 offspring
//\param	parent1 parent
//\param	parent2 parent
//\return	void
void SBX(CIndividualMO& son1, CIndividualMO& son2, CIndividualMO& parent1, CIndividualMO& parent2);

//!\brief	Parent-Centric Recombination (PCX)
//!\param	sizenew size of offpsring population
//!\param	popnew	offspring population
//!\param	pop		parent population
//\return	bool
bool PCX(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);

/////////////////////////////////////////////////////////////////////////////////////
//!\brief	DE operator
class XDE
{
protected:
	double	mF,				//!< setp length
			mCR;			//!< crossover probability
public:
	//!\brief	constructor
	//!\brief	no
	XDE();

	//!\brief	set parameters
	//!\param	f		 step length
	//!\param	cr		 crossover probability
	//!\return	void
	void Set(double f, double cr);

	//!\brief	Generator
	//!\param	sizenew size of offpsring population
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
};//class DE

/////////////////////////////////////////////////////////////////////////////////////
//!\brief	Fuzzy operator
class XFuzzy
{
public:
	//!\brief	constructor
	//!\brief	no
	XFuzzy();

	//!\brief	Generator
	//!\param	sizenew size of offpsring population
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
};//class Fuzzy


/////////////////////////////////////////////////////////////////////////////////////
//!\brief	GTX operator
class XGTX
{
protected:
	std::vector<unsigned int>	mvMatePool;			//!< mate pool
	std::vector<double>			mvCoefficient;		//!< coefficient
public:
	//!\brief	constructor
	//!\brief	no
	XGTX();

	//!\brief	Generator
	//!\param	sizenew size of offpsring population
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
};//class GTX

/////////////////////////////////////////////////////////////////////////////////////
//!\brief	Shape based crossover and mutation
class XSimplex
{
protected:
	double	mBeta,			//!< algorithm parameter
			mExten,			//!< extension
			mAlpha;			//!< neighborhood control parameter
	unsigned int mMaxGen,	//!< maximum generation
				 mCurGen,	//!< current generation
				 mType;		//!< neighborhood strategy
	std::vector< std::vector<double> > mMatDis;	//!< distance matrix
public:
	//!\brief	constructor
	//!\return	no
	XSimplex();

	//!\brief	set algorithm parameter
	//!\param	maxgen	maximum generation
	//!\param	beta	parameter of Gaussian noise
	//!\param	alpha	neighborhood control parameter
	//!\param	type	neighborhood strategy
	//!\return	no
	void Set(unsigned int maxgen, double beta, double alpha, unsigned int type);

	//!\brief	Generator
	//!\param	sizenew size of offpsring population
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
private:
	//!\brief	Generator for 1-D structure with linear models
	//!\param	popnew offspring population
	//!\param	pop	parent population
	//!\return	offspring population
	CPopulationMO& Generate1D(CPopulationMO& popnew, CPopulationMO& pop);

	//!\brief	Generator for 2-D structure with linear models
	//!\param	popnew offspring population
	//!\param	pop	parent population
	//!\return	offspring population
	CPopulationMO& Generate2D(CPopulationMO& popnew, CPopulationMO& pop);

	//!\brief	calculate the distance between two individuals
	//\param	ind1 individual 1
	//\param	ind2 individual 2
	//\return	distance(double)
	double Distance(CIndividualMO& ind1, CIndividualMO& ind2);

	//!\brief	MaxiMin sort 
	//!\param	pop		reference population
	//!\param	index	core point index
	//!\param	size	the number of points to be chosen out
	//!\return	void
	void MaxiMin(CPopulationMO& pop, std::vector<unsigned int>& index, unsigned int size);
};//class XSimplex

/////////////////////////////////////////////////////////////////////////////////////////////
//!\brief	Parent-Centric Recombination (PCX)
class XPCX
{
public:
	//!\brief	Generator
	//!\param	sizenew size of offpsring population
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
};//class PCX

/////////////////////////////////////////////////////////////////////////////////////////////
//!\brief	Parent-Centric Recombination (PCX)
class XSBX
{
public:
	//!\brief	Generator
	//!\param	sizenew size of offpsring population
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
};//class SBX

/////////////////////////////////////////////////////////////////////////////////////
//!\brief	NSDE operator
class XNSDE
{
protected:
	double	mF,				//!< setp length
			mCR;			//!< crossover probability
public:
	//!\brief	constructor
	//!\brief	no
	XNSDE();

	//!\brief	set parameters
	//!\param	f		 step length
	//!\param	cr		 crossover probability
	//!\return	void
	void Set(double f, double cr);

	//!\brief	Generator
	//!\param	sizenew size of offpsring population
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
};//class XNSDE

/////////////////////////////////////////////////////////////////////////////////////
//!\brief	Extended Simplex operator
class ESimplex
{
public:
	//!\brief	constructor
	//!\brief	no
	ESimplex();

	//!\brief	Generator
	//!\param	sizenew size of offpsring population
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop);
};//class ESimplex

} //namespace gen

} //namespace mea

} //namespace az

#endif //AZ_GENERATOR_EA_H
