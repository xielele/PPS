/*! \file	IndividualMO.h
	
	\brief	individual class for MOEA
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Dec.03 2004 2004  create
	\date	Sep.27 2005 rewrite & reorganize structure
	\date	Oct.14 2005 add CDominate()
	\date	Nov.15 2005 add ID
	\date	Mar.21 2008 add e-dominate
*/

#ifndef	AZ_INDIVIDUALMO_H
#define	AZ_INDIVIDUALMO_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include "Parameter.h"
#include "alg/Random.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{

//!\brief individual class for MOEA
class CIndividualMO
{
protected:
	unsigned int mRank,			//!< rank value
				 mID,			//!< ID
				 mOpt;			//!< ID of Optimizer, used for hybrid method
	double	mConstraint;		//!< constraint = sum(|mvEq|) + sum(mvIneq > 0)
	std::vector<double>	mvX	,	//!< variables
						mvF,	//!< objective values
						mvEq,	//!< equality values
						mvIneq;	//!< inequality values
	CParameter* pPar;			//!< the pointer to the parameter
public:
	//!\brief	Constractor
	//!\param	par parameters
	//!\return	void
	CIndividualMO(CParameter& par);

	//!\brief	Constractor
	//!\param	ind another individual
	//!\return	void
	CIndividualMO(CIndividualMO& ind);

	//!\brief	get the parameter object reference
	//!\return	the parameter object reference
	inline CParameter& P() {return *pPar;}

	//!\brief	check to see if it's feasible
	//!\return	whether the indiviudal is feasible
	inline bool IsFeasible() {return mConstraint <= P().TolC();}

	//!\brief	get the ith objective
	//!\param	i objective index
	//!\return	the ith objective reference
	inline double& F(unsigned int i) {return mvF[i];}

	//!\brief	get the ith variable
	//!\param	i variable index
	//!\return	the ith variable reference
	inline double& X(unsigned int i) {return mvX[i];}

	//!\brief	get the ith variable
	//!\param	i variable index
	//!\return	the ith variable reference
	inline double& operator[](unsigned int i) {return mvX[i];}

	//!\brief	get the constraint reference
	//!\return	the constraint reference
	inline double& C() {return mConstraint;}				

	//!\brief	get the objective vector reference
	//!\return	the objective vector reference
	inline std::vector<double>&	F() {return mvF;}

	//!\brief	get the variable vector reference
	//!\return	the variable vector reference
	inline std::vector<double>&	X() {return mvX;}

	//!\brief	get the variable vector reference
	//!\return	the variable vector reference
	inline std::vector<double>&	operator()() {return mvX;}

	//!\brief	set the rank value
	//!\param	r new rank value
	//!\return	rank value
	inline unsigned int Rank(unsigned int r) {mRank=r; return r;}

	//!\brief	get the rank value
	//!\return	rank value
	inline unsigned int Rank() {return mRank;}

	//!\brief	set ID
	//!\param	id new ID
	//!\return	new ID
	inline unsigned int ID(unsigned int id) {mID=id; return id;}

	//!\brief	get ID
	//!\return	ID
	inline unsigned int ID() {return mID;}

	//!\brief	set Optizer
	//!\param	opt new Optimizer
	//!\return	new Optimizer
	inline unsigned int OPT(unsigned int opt) {mOpt=opt; return opt;}

	//!\brief	get Optimizer
	//!\return	Optimizer
	inline unsigned int OPT() {return mOpt;}

	//!\brief	dominance check
	//!\param	ind another individual
	//!\return	dominate result: 1 dominate; -1 dominated; 0 non-dominated
	int Dominate(CIndividualMO& ind);

	//!\brief	constratint-dominance check
	//!\param	ind another individual
	//!\return	dominate result: 1 dominate; -1 dominated; 0 non-dominated
	int CDominate(CIndividualMO& ind);

	//!\brief	e-constratint-dominance check
	//!\param	ind another individual
	//!\param	e parameter vector
	//!\return	dominate result: 1 dominate; -1 dominated; 0 non-dominated
	int CEDominate(CIndividualMO& ind, std::vector<double>& e);

	//!\brief	evaluate the individual
	//!\return	void
	void Evaluate();

	//!\brief	make indiviudal to be a feasible one
	//!\brief	no
	void Check();

	//!\brief	set variable to be another individual
	//!\param	ind another individual
	//!\return	reference to the individual
	CIndividualMO& operator =(CIndividualMO& ind);

	//!\brief	check to see if two individuals are equal
	//!\param	ind another individual
	//!\return	whether they are equal
	bool operator ==(CIndividualMO& ind);

	//!\brief	check to see which is better
	//!\param	ind another individual
	//!\return	true if this individual is better than ind
	bool operator<(CIndividualMO& ind);

	//!\brief	write to a stream
	//!\param	os output stream
	//!\param	ind individual
	//!\return	output stream
	friend std::ostream& operator<<(std::ostream& os, CIndividualMO& ind);

	//!\brief	read from a stream
	//!\param	is input stream
	//!\param	ind individual
	//!\return	input stream
	friend std::istream& operator>>(std::istream& is, CIndividualMO& ind);
};//class CIndividualMO

} //namespace mea

} //namespace az

#endif //AZ_INDIVIDUALMO_H
