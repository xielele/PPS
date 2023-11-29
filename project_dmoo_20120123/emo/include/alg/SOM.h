/*! \file	SOM.h
	
	\brief	Self-Ogranizing Maps

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Mar.23 2005 create
	\date	Oct.04 2005 modify
*/

#ifndef AZ_SOM_H
#define	AZ_SOM_H

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

//\brief Batch version of SOM
class CBatchSOM
{
protected:
	int 
		mDim,		//!< the dimension of the data
		mRow,		//!< the row number of the toplogy
		mCol,		//!< the column number of the toplogy
		mWeightNo,	//!< the number of weights
		mDataNo,	//!< the number of data
		*pNear;		//!< the index of the nearest weight for all data point		
	double 
		mTheta,		//!< the current variance for neighbourhood function
		mThetaMin,	//!< the minimal variance for neighbourhood function
		mThetaMax,	//!< the maximal variance for neighbourhood function
		**pWeight,	//!< the weight of the toplogy
		**pData,	//!< the poINTer to the training data
		*pDis;		//!< the distance from a data to its core		
	bool bInit;		//!< whether it need to be initialized	
public:
	//!\brief	constructor
	//!\return	no
	CBatchSOM();

	//!\brief	destructor
	//!\return	no
	~CBatchSOM();

	//!\brief	initizalie parameters
	//!\param	row row number of neuro grid
	//!\param	col column number of neuro grid
	//!\param	dim data dimension
	//!\return	no
	void Initialize(int row, int col, int dim);

	//!\brief	main train step
	//!\param	steps train steps
	//!\param	size data size
	//!\param	pdata data pointer
	//!\return	no
	void Train(int steps, int size, double **pdata);

	//!\brief	get data size
	//!\return	data size
	inline int DataSize() {return mDataNo;}
	
	//!\brief	get weight matrix
	//!\return	weight matrix
	inline double** Weight() {return pWeight;}

	//!\brief	get a weight
	//!\param	row row index of the weitht matrix
	//!\param	col column index of the weight matrix
	//!\return	the weight	
	inline double* Weight(int row, int col) {return pWeight[row*mCol + col];}

	//!\brief	distance to the nearest neuron
	//!\param	k data index
	//!\return	the distance
	inline double DisToCore(int k) {return pDis[k];}

	//!\brief	get the nearest neruon index
	//!\param	k data index
	//!\return	the neuron index
	inline int NearCore(int k) {return pNear[k];}
protected:
	//!\brief	clear memory
	//!\return	no
	void	Clear();

	//!\brief	find the nearest neuron to data k
	//!\param	k data index
	//!\return	no
	void	NearWeight(int k);

	//!\brief	udate k-th weight
	//!\param	k weight inidex
	//!\return	no
	void	UpdateWeight(int k);

	//!\brief	calculate the distance between two neurons
	//!\param	k neuro index
	//!\param	n neuro index
	//!\return	the distance
	double	Neighbourhood(int k, int n);
}; //class CBatchSOM

} //namespace alg

} //namespace az

#endif	//AZ_SOM_H
