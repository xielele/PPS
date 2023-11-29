/*! \file	Sel.h
	
	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Sep.27 2005 create
	\date	Nov.16 2005 add SelectSort
	\date	Nov.22 2005 add CDistance2 selection which only uses the 2nd nearest distance
	\date	Nov.23 2005 add CStrength
	\date   Mar.11 2006 add CSelectionGTM
	\date	Mar.12 2006 add CSelectionMS
	\date	Mar.13 2006 add CSelectionMM
	\date	Feb.20 2007	add SCrowdXF
	\date	Feb.21 2007 add SMaxMinX
	\data	Apr.30 2007 add SRegF2
	\data	Mar.19 2008 add AMS
	\data	Mar.21 2008 rewrite SCrowdXF as Omni-Optimizer
	\data	Jun.22 2008 add SelectModel
	\date	Jul.06 2008 add SEmpty
*/
#ifndef AZ_SELECTION_H
#define AZ_SELECTION_H

#include <vector>
#include <list>
#include "IndividualMO.h"
#include "PopulationMO.h"
#include "alg/amoeba.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief select strategies
namespace sel
{

//!\brief an empty selection (noting to do)
class SEmpty
{
public:
	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size) {return pop;}
};//class SEmpty

//!\brief crowded selection(NSGA-II) strategy
class SCrowd
{
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
};//class SCrowd

//////////////////////////////////////////////////////////////////////////////////////////////////////
//!\brief crowded selection(NSGA-II) strategy
class SCrowd2
{
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
};//class SCrowd2

//!\brief K-nearest neighbor selection(SPEA2) strategy
class SNeighbor
{
protected:
	//!\brief structure descirbes each point
	struct NODE
	{
		bool bKeep;					//!< been kept
		std::vector<int> vMin;		//!< minimum distance point index to this point
		std::list<int> vInfluence;	//!< the points influnced by this point
	};

	int mSize,									//!< point number
		mRange;									//!< how many nearest points been considered
	std::vector< NODE > vNode;					//!< points vector
	std::vector< std::vector<double> > mDis;	//!< distance matrix
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);	
protected:
	//!\brief	sort the population [start, end)
	//!\param	pop population
	//!\param	start start point
	//!\param	end end point
	//!\param	remainsize remain size
	//!\return	population
	CPopulationMO& Sort(CPopulationMO& pop, int start, int end, int remainsize);

	//!\brief	get the distnace between two points
	//!\param	i point index
	//!\param	j point index
	//!\return	distance between two points
	double& DisMat(int i, int j);

	//!\brief	initialize some parameters
	//!\return	void
	void SetPar();

	//!\brief	erase one point
	//!\return	void
	void EraseOne();
};//class SNeighbor

//!\brief distance selection strategy like CDistance but only uses the 2nd nearest distance to estimate density
class SNeighbor2:public SNeighbor
{
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);
protected:
	//!\brief	erase one point
	//!\return	void
	void EraseOne();
};//class CDistance2

//!\brief strength selection(SPEA2) strategy
class SStrength: public SNeighbor
{
protected:
	std::vector< double > mFitness;
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);
protected:
	//!\brief	assign rank
	//!\param	pop population
	//!\return	population
	CPopulationMO& RankAssignment(CPopulationMO& pop);

	//!\brief	sort the population [start, end) according to strength
	//!\param	pop population
	//!\param	start start point
	//!\param	end end point
	//!\return	population
	CPopulationMO& StrengthSort(CPopulationMO& pop, unsigned int start, unsigned int end);
};//class SNeighbor2

//!\brief NSGA-II + SPEA2 strategy
class SCrowdNeighbor
{
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);	

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);	
};//class SCrowdNeighbor

//!\brief NSGA-II + SPEA2 strategy
class SCrowdNeighbor2
{
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);	

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);	
};//class SCrowdNeighbor2

//////////////////////////////////////////////////////////////////////////////////////////////////////
//!\brief Max-Min Selection (MMS) strategy
class SMaxMin
{
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);	

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);	
private:
};//class SMaxMin

//////////////////////////////////////////////////////////////////////////////////////////////////////
//!\brief Max-Min Selection (MMS) strategy in X-space
class SMaxMinX
{
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);	

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);	
private:
};//class SMaxMin

////////////////////////////////////////////////////////////////////////////////////////////////////////
////!\brief Selection based on Regularity in F-space for 2-objective problems (SRegF2)
//class SRegF2:public alg::amoeba
//{
//protected:
//	bool	mIsModel;									//!< the model is built or nont
//	std::vector< std::vector<double> > mvRef, mvRef0;	//!< reference points
//	std::vector<double> mvLen, mvLen0;					//!< length of each segment
//	double mLen, mLen0;									//!< the total length
//	unsigned int	mNoRef,								//!< number of reference points
//					mDim,								//!< dimension of decision pace (default is 2)
//					mCount;								//!< generations a model used
//	// for local search
//	std::vector<double> LDAB, LDOAB, LC;
//	CPopulationMO* pLPOP;
//public:
//	SRegF2() {mIsModel = false; mLen0 = 0.0; mCount = 0;}
//
//	//!\brief	select some reference central points from the population
//	//!\param	pop		population
//	//!\param	cenp	central points
//	//!\return	void
//	void Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp);
//
//	//!\brief	select some best individuals into next generation by using model 
//	//!\param	pop		population
//	//!\param	size	size of best ones
//	//!\param	flag	indicating selected or not
//	//!\param	fix		fix reference point or not
//	//!\return	population
//	CPopulationMO& Select(CPopulationMO& pop, unsigned int size, std::vector<bool>& flag);
//
//	double objective(std::vector<double>& x);
//	CPopulationMO& LocalSearch(CPopulationMO& popc, CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select some reference central points from the population
//	//!\param	pop		population
//	//!\param	cenp	central points
//	//!\param	inmodel	central points are from model
//	//!\return	void
//	void SelectCentre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp, bool inmodel);
//
//	//!\brief	select some best individuals into next generation by using (line) model 
//	//!\param	pop		population
//	//!\param	size	size of best ones
//	//!\param	flag	indicating selected or not
//	//!\param	fix		fix reference point or not
//	//!\return	population
//	CPopulationMO& SelectLine(CPopulationMO& pop, unsigned int size, std::vector<bool>& flag, bool fix);
//
//	////!\brief	select some best individuals into next generation by using (line-segment) model 
//	////!\param	pop		population
//	////!\param	size	size of best ones
//	////!\param	flag	indicating selected or not
//	////!\param	fix		fix reference point or not
//	////!\return	population
//	//CPopulationMO& SelectLineSegment(CPopulationMO& pop, unsigned int size, std::vector<bool>& flag, bool fix);
//
//	//!\brief	sort population, the first 'size' are best ones
//	//!\param	pop		population
//	//!\param	size	size of best ones
//	//!\return	population
//	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select from current population and offspring population
//	//!\param	pop		combined population
//	//!\param	size	population size
//	//!\return	population
//	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
//private:
//	//!\brief	build reference points
//	//!\param	pop		population
//	//!\return	void
//	void ChooseRef(CPopulationMO& pop);
//
//	////!\brief	set reference points
//	////!\param	iC current reference point
//	////!\param	iL		left reference point
//	////!\param	iR		right reference point
//	////!\param	vRef	set of reference points
//	////!\param	pop		population
//	////!\return	void
//	//void ChooseRef(unsigned int iC, unsigned int iL, unsigned int iR, std::vector< std::vector<double> >& vRef, CPopulationMO& pop);
//	////!\brief	set reference points
//	////!\param	vRef	set of reference points
//	////!\param	pop		population
//	////!\return	void
//	//void ChooseRef(std::vector< std::vector<double> >& vRef, CPopulationMO& pop);
//
//	////!\brief	select points to reference points
//	////!\param	pop		population
//	////!\param	vExist	indicates points are selected or not
//	////!\param	vRef	reference points
//	////!\param	size	number of points to select
//	////!\param	fix		fix reference point or not
//	////!\return	void
//	//void Select(CPopulationMO& pop, std::vector<bool>& vExist, std::vector< std::vector<double> >& vRef, unsigned int size, bool fix);
//};//class SRegF2
//
////////////////////////////////////////////////////////////////////////////////////////////////////////
////!\brief Approximation Model based Selection in Two Spaces
//class DualMM_OLD
//{
//protected:
//	unsigned int	mDimF,mDimX;	//!< dimensions of objective and decision spaces
//	std::vector<bool> mvSel;						//!< indicate a point is selected or not
//	std::vector< std::vector<double> > mvDisX;		//!< distance matrix in the search space
//	std::vector<double>	mvDisToSet;					//!< distance to a selected set
//	//parameters for 2-obj problems
//	std::vector<double>	mvTarLen,					//!< target point in the line segment
//						mvProL,						//!< projection on the line segment
//						mvProOL;					//!< projection orthogonal to the line segment
//	//parameters for 3-obj problems
//	std::vector< std::vector<double> > mvWeight;	//!< weight vectors
//	double Fa, Fb, Fc, Fd;							//!< plane function parameters Fa*x + Fb*y + Fc*z + Fd = 0
//	double Trans;									//!< translation along orthogonal direction
//	std::vector<double> A,B,C;						//!< vertex
//public:
//	//!\brief	sort population, the first 'size' are best ones
//	//!\param	pop		population
//	//!\param	size	size of best ones
//	//!\return	population
//	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select from current population and offspring population
//	//!\param	pop		combined population
//	//!\param	size	population size
//	//!\return	population
//	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select some reference central points from the population
//	//!\param	pop		population
//	//!\param	cenp	central points
//	//!\return	void
//	void Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp);
//private:
//	//!\brief	build approximation model and choose target points
//	//!\param	pop		population
//	//!\param	size	number of target points
//	//!\return	void
//	void BuildModel2(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	build approximation model and choose target points
//	//!\param	pop		population
//	//!\param	size	number of target points
//	//!\return	void
//	void BuildModel3(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	init distance matrix
//	//!\param	pop		population
//	//!\return	void
//	void InitMatrix(CPopulationMO& pop);
//
//	unsigned int SelX(CPopulationMO& pop);
//	unsigned int SelF(unsigned int k, CPopulationMO& pop);
//	void UpdateDistance(unsigned int index);
//};//class AMS

//!\brief crowded selection(NSGA-II) strategy in both the decision and the objective spaces (Omni-Optimizer)
class DualCrowd
{
public:
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
};//class SCrowd

////////////////////////////////////////////////////////////////////////////////////////////////////////
////!\brief MaxiMin Selection in Two Spaces
//class DualMaxiMin
//{
//protected:
//	unsigned int	mDimF,mDimX;	//!< dimensions of objective and decision spaces
//	std::vector<bool> mvSel;						//!< indicate a point is selected or not
//	std::vector< std::vector<double> >	mvDisX,		//!< distance matrix in the search space
//										mvDisF;		//!< distance matrix in the objective space
//	std::vector<double>	mvDisToSetX,				//!< distance to a selected set
//						mvDisToSetF;				
//public:
//	//!\brief	sort population, the first 'size' are best ones
//	//!\param	pop		population
//	//!\param	size	size of best ones
//	//!\return	population
//	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size, unsigned int method);
//
//	//!\brief	select from current population and offspring population
//	//!\param	pop		combined population
//	//!\param	size	population size
//	//!\return	population
//	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select some reference central points from the population
//	//!\param	pop		population
//	//!\param	cenp	central points
//	//!\return	void
//	void Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp);
//private:
//	//!\brief	init distance matrix
//	//!\param	pop		population
//	//!\return	void
//	void InitMatrix(CPopulationMO& pop);
//
//	unsigned int SelX(CPopulationMO& pop);
//	unsigned int SelF(CPopulationMO& pop);
//	void UpdateDistance(unsigned int index);
//};//class AMS

////////////////////////////////////////////////////////////////////////////////////////////////////////
////!\brief Approximation Model based Selection in Two Spaces
//class DualMM
//{
//protected:
//	unsigned int	mDimF,mDimX;	//!< dimensions of objective and decision spaces
//	std::vector<bool> mvSel;						//!< indicate a point is selected or not
//	std::vector< std::vector<double> > mvDisX;		//!< distance matrix in the search space
//	std::vector<double>	mvDisToSet;					//!< distance to a selected set
//	//parameters for 2-obj problems
//	std::vector< std::vector<double> >	mvSimplex,	//!< the points to form a simplex
//										mvWeight;	//!< weight vectors
//	std::vector<double>					mNormal;	//!< normal direction
//	double Trans;									//!< translation along normal direction
//public:
//	//!\brief	sort population, the first 'size' are best ones
//	//!\param	pop		population
//	//!\param	size	size of best ones
//	//!\return	population
//	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select from current population and offspring population
//	//!\param	pop		combined population
//	//!\param	size	population size
//	//!\return	population
//	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select some reference central points from the population
//	//!\param	pop		population
//	//!\param	cenp	central points
//	//!\return	void
//	void Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp);
//private:
//	//!\brief	build approximation model and choose target points
//	//!\param	pop		population
//	//!\param	size	number of target points
//	//!\return	void
//	void BuildSimplex(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	init distance matrix
//	//!\param	pop		population
//	//!\return	void
//	void InitMatrix(CPopulationMO& pop);
//
//	unsigned int SelX(CPopulationMO& pop);
//	unsigned int SelF(unsigned int k, CPopulationMO& pop);
//	void UpdateDistance(unsigned int index);
//};//class AMS
//
////!\brief Approximation Model based Selection in Two Spaces
//class DualM
//{
//protected:
//	unsigned int	mDimF,mDimX;	//!< dimensions of objective and decision spaces
//	std::vector<bool> mvSel;						//!< indicate a point is selected or not
//	std::vector< std::vector<double> > mvDisX;		//!< distance matrix in the search space
//	std::vector<double>	mvDisToSet;					//!< distance to a selected set
//	//parameters for 2-obj problems
//	std::vector< std::vector<double> >	mvSimplex,	//!< the points to form a simplex
//										mvWeight;	//!< weight vectors
//	std::vector<double>					mNormal;	//!< normal direction
//	double Trans;									//!< translation along normal direction
//public:
//	//!\brief	sort population, the first 'size' are best ones
//	//!\param	pop		population
//	//!\param	size	size of best ones
//	//!\return	population
//	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select from current population and offspring population
//	//!\param	pop		combined population
//	//!\param	size	population size
//	//!\return	population
//	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	select some reference central points from the population
//	//!\param	pop		population
//	//!\param	cenp	central points
//	//!\return	void
//	void Centre(CPopulationMO& pop, std::vector< std::vector<double> >& cenp);
//private:
//	//!\brief	build approximation model and choose target points
//	//!\param	pop		population
//	//!\param	size	number of target points
//	//!\return	void
//	void BuildSimplex(CPopulationMO& pop, unsigned int size);
//
//	//!\brief	init distance matrix
//	//!\param	pop		population
//	//!\return	void
//	void InitMatrix(CPopulationMO& pop);
//
//	unsigned int SelF(unsigned int k, CPopulationMO& pop);
//	void UpdateDistance(unsigned int index);
//};

//!\brief Approximation Model based Selection
class SelectModel
{
protected:
	unsigned int	mDimF,mDimX;	//!< dimensions of objective and decision spaces
	std::vector<bool> mvSel;						//!< indicate a point is selected or not
	std::vector< std::vector<double> > mvDisX;		//!< distance matrix in the search space
	std::vector<double>	mvDisToSet;					//!< distance to a selected set
	//parameters for 2-obj problems
	std::vector< std::vector<double> >	mvSimplex,	//!< the points to form a simplex
										mvWeight;	//!< weight vectors
	std::vector<double>					mNormal;	//!< normal direction
	double Trans;									//!< translation along normal direction

	// the following parameters are used to adaptively tune the extension rate
	std::vector< double >	mvExtreme,	//!< the location of the extreme points
							mvExtension;//!< the extension rate for each extreme point
	std::vector< unsigned int > mvAge;	//!< the age of each extreme point
public:
	//!\brief	constructor
	//!\return	none
	SelectModel();

	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop		population
	//!\param	size	size of best ones
	//!\return	population
	CPopulationMO& SelectSort(CPopulationMO& pop, unsigned int size);

	//!\brief	select from current population and offspring population
	//!\param	pop		combined population
	//!\param	size	population size
	//!\return	population
	CPopulationMO& Select(CPopulationMO& pop, unsigned int size);
private:
	//!\brief	build approximation model and choose target points
	//!\param	pop		population
	//!\param	size	number of target points
	//!\return	void
	void BuildSimplex(CPopulationMO& pop, unsigned int size);

	//!\brief	init distance matrix
	//!\param	pop		population
	//!\return	void
	void InitMatrix(CPopulationMO& pop);

	unsigned int SelF(unsigned int k, CPopulationMO& pop);
	unsigned int SelX(CPopulationMO& pop);
	void UpdateDistance(unsigned int index);
};

}//namespace sel

} //namespace mea

} //namespace az


#endif //AZ_SELECTION_H
