/*! \file	GenMod.h
	
	\brief	Model based generators for MOEA
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Aug.16 2005 create
	\date	Sep.25 2005 rewrite & reorganize structure
	\date	Nov.25 2005 add ROUND() and PERTUBATION() to ModelBase
	\date	Apr.07 2006 add class ModelPCA
	\date	Nov.12 2006 add class ModelLocalPCAU
	\date	Nov.27 2006 add class MixGauss
	\date	Jan.28 2007	add class ModelAL
	\date	May.10 2007 rewrite class ModelAL
	\date	Jun.25 2007 add class HybridRM
	\date	Jun.26 2007 rename ModelLocalPCAU -> RM
						rename ModelAL -> RM2
*/


#ifndef AZ_GENERATOR_MODEL_H
#define AZ_GENERATOR_MODEL_H

#include <vector>
#include "alg/PCA.h"
#include "alg/LocalPCAO.h"
#include "alg/SOM.h"
#include "alg/GTM.h"
#include "alg/Matrix.h"
#include "PopulationMO.h"
#include "Gen.h"
#include "Sel.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{

//!\brief	gen namespace, offspring generate strategies
namespace gen
{

//!\brief model based generator
namespace mod
{

//=========================================================================================
//!\brief Guided Crossover
class GuidedXOver
{
public:
	//!\brief	crossover
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& XOver(CPopulationMO& popnew, CPopulationMO& popref);
protected:
	//!\brief	find the nearest neighbor in reference population
	//!\param	ind individual
	//!\param	start the begin search index in reference population
	//!\param	end the last search index in reference population
	//!\param	popref reference population
	//!\return	nearest neighbor index
	unsigned int Neighbor(CPopulationMO::IND_TYPE& ind, unsigned int start, unsigned int end, CPopulationMO& popref);
};//class GuidedXOver

//=========================================================================================
//!\brief the basic class for model-based generators
class ModelBase
{
protected:
	double			**pData,	//!< pointer to data
					mExtension;	//!< Principal Curve(Surface) extension ratio
	unsigned int	mDataSize,	//!< data number
					mDataDim,	//!< data dimension
					mTrainSteps;//!< inner train steps
public:
	//!\brief	constructor
	//!\return	void
	ModelBase();

	//!\brief	destructor
	//!\return	void
	virtual ~ModelBase();
protected:
	//!\brief	clear data pool
	//!\return	void
	void Clear();
	
	//!\brief	convert a real number to a nearest positive integer number
	//!\param	X real number
	//!\return	integer number
	unsigned int ROUND(double X);

	//!\brief	pertubation around 1.5
	//!\return	pertubation number
	double PERTUBATION();
};//class ModelBase

//=========================================================================================
//!\brief EDA+GA hybrid generator(A), EDA and GA runs iteratively
class ModelHybridA:public ModelBase, public gen::XSBX, public GuidedXOver
{
protected:
	alg::LocalPCAO mLocalPCA;	//!< Local PCA object
	double			mRho;		//!< convergence condition
	unsigned int	mLatentDim,	//!< latent dimension
					mMaxCluster;//!< maximum cluster number
	bool			mbReset,	//!< flag to reset Local PCA
					mbGA;		//!< flag to use GA	
public:
	//!\brief	constructor
	//!\return	void
	ModelHybridA();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	initialize the Local PCA
	//!\param	latent latent dimension
	//!\param	cluster maximum cluster number
	//!\param	trainsteps Local PCA train steps
	//!\param	extension extension ratio
	//!\param	rho convergence condition
	//!\return	void
	void Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension, double rho);
		
	//!\brief	reset Local PCA
	//!\return	void
	void Reset();
protected:
	//!\brief	modle-based generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);	

	//!\brief	modle-based generator for 1D structure
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref);	

	//!\brief	modle-based generator for 2D structure
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref);	

	//!\brief	project point into a dimension
	//!\param	dim dimension index
	//!\param	index point index
	//!\param	clu cluster index
	//!\return	projection value
	double Project(unsigned int dim, unsigned int index, unsigned int clu);
}; //class ModelHybridA

//=========================================================================================
//!\brief EDA+GA hybrid generator(B), EDA and GA runs in a single step
class ModelHybridB:public ModelBase, public gen::XSBX, public GuidedXOver
{
protected:
	alg::LocalPCAO mLocalPCA;	//!< Local PCA object
	double			mRho;		//!< convergence shreshhold
	unsigned int	mLatentDim,	//!< latent dimension
					mMaxCluster;//!< maximum cluster number
	bool			mbReset;	//!< flag to reset Local PCA
public:
	//!\brief	constructor
	//!\return	void
	ModelHybridB();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	initialize the LPCA
	//!\param	latent latent dimension
	//!\param	cluster maximum cluster number
	//!\param	trainsteps Local PCA train steps
	//!\param	extension extension ratio
	//!\param	rho convergence shreshhold
	//!\return	void
	void Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension, double rho);
		
	//!\brief	reset Local PCA
	//!\return	void
	void Reset();
protected:
	//!\brief	model-based generator for 1D structure
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref);	

	//!\brief	model-based generator for 2D structure
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref);	
	
	//!\brief	project point into a dimension
	//!\param	dim dimension index
	//!\param	index point index
	//!\param	clu cluster index
	//!\return	projection value
	double Project(unsigned int dim, unsigned int index, unsigned int clu);
}; //class ModelHybridB

//=========================================================================================
//!\brief Local PCA based EDA generator
class ModelLocalPCA:public ModelBase, public GuidedXOver
{
protected:
	alg::LocalPCAO mLocalPCA;	//!< Local PCA object
	unsigned int	mLatentDim,	//!< latent dimension
					mMaxCluster;//!< maximum cluster number
	bool			mbReset;	//!< flag to reset Local PCA
public:
	//!\brief	constructor
	//!\return	void
	ModelLocalPCA();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	initialize the LPCA
	//!\param	latent latent dimension
	//!\param	cluster maximum cluster number
	//!\param	trainsteps Local PCA train steps
	//!\param	extension extension ratio
	//!\return	void
	void Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension);
		
	//!\brief	reset Local PCA
	//!\return	void
	void Reset();
protected:
	//!\brief	modle-based generator for 1D structure
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref);	

	//!\brief	modle-based generator for 2D structure
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref);	

	//!\brief	project point into a dimension
	//!\param	dim dimension index
	//!\param	index point index
	//!\param	clu cluster index
	//!\return	projection value
	double Project(unsigned int dim, unsigned int index, unsigned int clu);
}; //class ModelLocalPCA

//=========================================================================================
//!\brief SOM based model generator
class ModelSOM:public ModelBase, public GuidedXOver
{
protected:
	alg::CBatchSOM	mSom;		//!< SOM object
	unsigned int		mRowGrid,	//!< row number of SOM grid
						mColGrid;	//!< column number of SOM grid
	bool				mbWeight;	//!< need to reset the weights
public:
	//!\brief	constractor
	//!\return	void
	ModelSOM();

	//!\brief	destractor
	//!\return	void
	~ModelSOM();
	
	//!\brief	initialize the SOM
	//!\param	row row number of the neuro grid
	//!\param	col column numbe of the neuro grid
	//!\param	trainstpes train steps
	//!\param	extension extension ratio when sampling
	//!\return	void
	void Set(unsigned int row, unsigned int col, unsigned int trainsteps, double extension);
	
	//!\brief	reset the initial weights of SOM
	//!\return	void
	void Reset();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);
protected:
	//!\brief	Principal Curve Generator
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate1D(CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	Principal Surface Generator
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate2D(CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	the distance between two center points
	//!\param	p1 point 1
	//!\param	p2 point 2
	//!\return	the distance
	double DisCC(unsigned int cp1, unsigned int cp2);

	//!\brief	the distance between a center point and a variable
	//!\param	cp center point index
	//!\param	vp variable point index
	//!\return	distance
	double DisCV(unsigned int cp, unsigned int vp);
};//ModelSOM

//=========================================================================================
//!\brief GTM based model generator
class ModelGTM1:public ModelBase, public GuidedXOver
{
protected:
	alg::GTM	mGTM;			//!< GTM object
	alg::Matrix		mW,			//!< weight matrix
					mFI,		//!< matrix of the basis functions
					mT;			//!< data matrix
	double			mBeta;		//!< inversence variance
	unsigned int	mLatentDim,	//!< latent variable dimension
					mNoLatent,	//!< number of latent variables
					mNoBaseFun;	//!< number of base functions
	bool			mbWeight;	//!< need to reset the weights
public:
	//!\brief	constractor
	//!\return	void
	ModelGTM1();

	//!\brief	destractor
	//!\return	void
	~ModelGTM1();
	
	//!\brief	initialize the GTM
	//!\param	noLatent number of latent variables
	//!\param	noBaseFun number of base functions
	//!\param	dimLatent latent structure dimension(1 or 2)
	//!\param	trainstpes train steps
	//!\param	extension extension ratio when sampling
	//!\return	void
	void Set(unsigned int noLatent, unsigned int noBaseFun, unsigned int dimLatent, unsigned int trainsteps, double extension);
	
	//!\brief	reset the initial weights of SOM
	//!\return	void
	void Reset();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);
protected:
	//!\brief	Principal Curve Generator
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	Principal Surface Generator
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	the distance between two center points
	//!\param	p1 point 1
	//!\param	p2 point 2
	//!\return	the distance
	double DisCC(unsigned int cp1, unsigned int cp2);

	//!\brief	the distance between a center point and a variable
	//!\param	cp center point index
	//!\param	vp variable point index
	//!\return	distance
	double DisCV(unsigned int cp, unsigned int vp);
};//ModelGTM

//=========================================================================================
//!\brief GTM based model generator
class ModelGTM2:public ModelBase, public GuidedXOver
{
protected:
	alg::GTM	mGTM;			//!< GTM object
	alg::Matrix	mW,			//!< weight matrix
					mFI,		//!< matrix of the basis functions
					mT;			//!< data matrix
	double			mBeta,		//!< inversence variance
					mR;			//!< variance
	unsigned int	mLatentDim,	//!< latent variable dimension
					mNoLatent,	//!< number of latent variables
					mNoBaseFun;	//!< number of base functions
	bool			mbWeight;	//!< need to reset the weights
public:
	//!\brief	constractor
	//!\return	void
	ModelGTM2();

	//!\brief	destractor
	//!\return	void
	~ModelGTM2();
	
	//!\brief	initialize the GTM
	//!\param	noLatent number of latent variables
	//!\param	noBaseFun number of base functions
	//!\param	dimLatent latent structure dimension(1 or 2)
	//!\param	trainstpes train steps
	//!\param	extension extension ratio when sampling
	//!\return	void
	void Set(unsigned int noLatent, unsigned int noBaseFun, unsigned int dimLatent, unsigned int trainsteps, double extension);
	
	//!\brief	reset the initial weights of SOM
	//!\return	void
	void Reset();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);
protected:
	//!\brief	Principal Curve Generator
	//!\param	popnew offspring population
	//!\param	popref	parent population
	//!\return	offspring population
	CPopulationMO& ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	Principal Surface Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	new	offspring population
	CPopulationMO& ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref);
};//ModelGTM2

//=========================================================================================
//!\brief inner clustering without learning, PCA is used to build model
class ModelPCA:public ModelBase
{
protected:
	bool			mbConverged;		//!< converged or not
	unsigned int	mLatentDim,			//!< latent dimension
					mModelSize,			//!< the number of model points
					mCandidateSize,		//!< size of candidate solutions
					mLastCandidateSize,	//!< last candidate size used while not converged
					mMaxGen,			//!< the maximum generation
					mCurGen,			//!< current generation
					mMatingStrategy,	//!< neighborhood strategy
					mHisIndex,			//!< history memory index
					mDeltaT;			//!< test period
	double			mSigmoidBeta,		//!< parameter of sigmoid function
					mConThreshold;		//!< convergence threshold
	alg::Matrix	mMat;					//!< store data in current cluster
	std::vector< CPopulationMO* >	pvHis;	//!< history pool
	std::vector< unsigned int > mvIndex,	//index points 
								mvOffSize;	//number of points to generate around the seed
	std::vector< std::vector< unsigned int > > mvModelIndex;	//!< index of model points
public:
	//!\brief	constructor
	//!\return	void
	ModelPCA();

	//!\brief	deconstructor
	//!\return	void
	~ModelPCA();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	set parameters
	//!\param	extension	extension ratio
	//!\param	maxgen		maximum generation
	//!\param	modsize		the number of model points
	//!\param	mate		mating restriction strategy
	//!\param	threshold	convergence threshold
	//!\return	void
	void Set(double extension, double threshold, unsigned int maxgen, unsigned int modsize, unsigned int mate);
		
	//!\brief	reset Local PCA
	//!\return	void
	void Reset();
protected:
	//!\brief	model-based generator for 1D structure
	//!\param	popnew	offspring population
	//!\param	popref	parent population
	//!\param	pca		PCA object
	//!\param	index	the location of popnew to create new individuals
	//!\param	size	the number of new points to create
	//!\return	offspring population
	CPopulationMO& ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref, alg::PCA& pca, unsigned int& index, unsigned int size);	

	//!\brief	model-based generator for 2D structure
	//!\param	popnew	offspring population
	//!\param	popref	parent population
	//!\param	pca		PCA object
	//!\param	index	the location of popnew to create new individuals
	//!\param	size	the number of new points to create
	//!\return	offspring population
	CPopulationMO& ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref, alg::PCA& pca, unsigned int& index, unsigned int size);	
	
	//!\brief	project point into a dimension
	//!\param	dim dimension index
	//!\param	index point index
	//!\param	pca PCA object
	//!\return	projection value
	double Project(unsigned int dim, unsigned int index, alg::PCA& pca);

	//!\brief	assign offspring to each group with equal probability
	//!\param	pop		reference population
	//!\param	tolsize	total size of offspring
	//!\param	number	the number of points to generate around the seed
	//!\return	void
	void Assign_Naive(CPopulationMO& pop, unsigned int tolsize, std::vector<unsigned int>& number);

	//!\brief	assign offspring to each group according to density
	//!\param	pop		reference population
	//!\param	tolsize	total size of offspring
	//!\param	number	the number of points to generate around the seed
	//!\return	void
	void Assign_Density(CPopulationMO& pop, unsigned int tolsize, std::vector<unsigned int>& number);

	//!\brief	MaxiMin sort 
	//!\param	pop		reference population
	//!\param	index	core point index
	//!\param	size	the number of points to be chosen out
	//!\return	void
	void MaxiMin(CPopulationMO& pop, std::vector<unsigned int>& index, unsigned int size);

	//!\brief	Convergence Test
	//!\param	pop		current population
	//!\return	convergence measurement
	double ConvergenceTest(CPopulationMO& pop);
}; //class ModelPCA

//!\brief MixGauss, mixture Gaussian model-based generator
class MixGauss
{
protected:
	double			mExtension,		//!< extension
					mErrTol;		//!< error tolerance
	unsigned int	mNClu,			//!< number of clusters
					mNX,			//!< number of trainning data
					mDX,			//!< dimension of trainning data
					mMaxIter,		//!< maximal trainning steps
					mDLat,			//!< dimension of latent space
					mType;			//!< type: UMG, MMG, PPCA
public:
	//!\brief	set parameters
	//!<param	type	model type
	//!<param	dlat	dimension of latent space
	//!\param	nclu	number of cluster
	//!\param	exten	extension
	//!\param	maxiter	maximal trainning steps
	//!\param	errtol	error tolerance
	//!\return	void
	void Set(unsigned int type, unsigned int dlat, unsigned int nclu, double exten, unsigned int maxiter=100, double errtol=1.0E-10);

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);
};//class MixGaussian

//=========================================================================================
//!\brief Local PCA based EDA generator
class RM:public ModelBase
{
protected:
	unsigned int	mLatentDim,	//!< latent dimension
					mMaxCluster;//!< maximum cluster number
public:
	//!\brief	constructor
	//!\return	void
	RM();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	initialize the LPCA
	//!\param	latent latent dimension
	//!\param	cluster maximum cluster number
	//!\param	trainsteps Local PCA train steps
	//!\param	extension extension ratio
	//!\return	void
	void Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension);
}; //class RM

//=========================================================================================
//!\brief Local PCA based EDA generator
class RM_OLD:public ModelBase
{
protected:
	alg::LocalPCAO mLocalPCA;	//!< Local PCA object
	unsigned int	mLatentDim,	//!< latent dimension
					mMaxCluster;//!< maximum cluster number
public:
	//!\brief	constructor
	//!\return	void
	RM_OLD();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	initialize the LPCA
	//!\param	latent latent dimension
	//!\param	cluster maximum cluster number
	//!\param	trainsteps Local PCA train steps
	//!\param	extension extension ratio
	//!\return	void
	void Set(unsigned int latent, unsigned int cluster, unsigned int trainsteps, double extension);
protected:
	//!\brief	project point into a dimension
	//!\param	dim dimension index
	//!\param	index point index
	//!\param	clu cluster index
	//!\return	projection value
	double Project(unsigned int dim, unsigned int index, unsigned int clu);
}; //class RM_OLD

//=========================================================================================
//!\brief EDA offspring generator with Adaptive Local Model
class RM2:public ModelBase
{
protected:
	unsigned int	mMinM,				//!< minimum number of models
					mMaxM,				//!< maximum number of models
					mNoM,				//!< number of adaptive local models (reference points)
					mGen,				//!< current generation
					mNoNew;				//!< number of new solutions
	double			mUPFExtRatio,		//!< Utopian PF extension rate
					mModExtRatio,		//!< model extension rate
					mPCAThreshold;		//!< threshold to determine the number of principal components		
	std::vector< unsigned int > mvNoChi,	//!< number of child points to generate around the seed
								mvDim;		//!< dimension of each subset
	std::vector< double >		mvStd;		//!< standard devision of noise in each subset
	std::vector< std::vector< double > >	mvCenM;		//!< locations of central point of models	
	std::vector< alg::Matrix >	mvDatMat;				//!< matrix to store data for each 'cluster'
	std::vector< alg::PCA > mvPCA;						//!< PCA
public:
	//!\brief	constructor
	//!\return	void
	RM2();

	//!\brief	deconstructor
	//!\return	void
	~RM2();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	set parameters
	//!\param	uext		utopian model extension ratio
	//!\param	mext		model extension ratio
	//!\param	threshold	principal components threshold
	//!\param	minmodel	minimum number of models
	//!\param	maxmodel	maximum number of models
	//!\return	void
	void Set(double uext, double mext, double threshold, unsigned int minmodel, unsigned int maxmodel);		
protected:
	//!\brief	adaptive linear model-based generator
	//!\param	popnew	offspring population
	//!\param	index	index location in offspring population
	//!\param	popref	parent population
	//!\param	pca		PCA object
	//!\param	std		std. of noise
	//!\param	dim		dimension of manifold
	//!\param	size	the number of new points to create
	//!\return	offspring population
	CPopulationMO& ALModelGen(CPopulationMO& popnew, unsigned int& index, CPopulationMO& popref, alg::PCA& pca, double std, unsigned int dim, unsigned int size);	

	//!\brief	project point into a dimension
	//!\param	dim dimension index
	//!\param	index point index
	//!\param	pca PCA object
	//!\return	projection value
	double Project(unsigned int dim, unsigned int index, alg::PCA& pca);

	//!\brief	set model parameters
	//!\param	pop		reference population
	//!\param	sizenew size of new solutions
	//!\return	void
	void BuildModel(CPopulationMO& pop, unsigned int sizenew);

	//!\brief	partition the population
	//!\param	pop		reference population
	//!\return	void
	void Partition(CPopulationMO& pop);

	//!\brief	generate reference points
	//!\param	pop		 population
	//!\param	cenp	 reference points
	//!\return	void
	void Reference(CPopulationMO& pop, std::vector< std::vector<double> >& cenp);

}; //class ModelPCA

class HybridRM
{
protected:
	RM				gRM;				//!< RM generator
	RM2		gERM;				//!< enhanced RM generator
	unsigned int	mMinM,				//!< minimum number of models
					mMaxM,				//!< maximum number of models
					mNoM,				//!< number of adaptive local models (reference points)
					mMaxGen,			//!< maximum number of generations
					mGen,				//!< current generation
					mNoNew,				//!< number of new solutions
					mHybStr;			//!< hybrid strategy 0: iterated, 1: F->X, 2: adaptively with probability 
	double			mExtension;			//!< Principal Curve(Surface) extension ratio
public:
	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	set parameters
	//!\param	extension	extension ratio
	//!\param	threshold	principal components threshold
	//!\param	maxgen		maximum number of generations
	//!\param	minmodel	minimum number of models
	//!\param	maxmodel	maximum number of modelss
	//!\param	type		hybrid strategy
	//!\return	void
	void Set(double extension, double threshold, unsigned int maxgen, unsigned int minmodel, unsigned int maxmodel, unsigned int type);
}; //class HybridRM

//=========================================================================================
//!\brief EDA offspring generator with Regularity Models in both the objective and decision spaces
class RM3:public ModelBase
{
protected:
	unsigned int	mMaxM,				//!< maximum number of models
					mMaxG,				//!< maximum number of generations
					mNoM,				//!< number of adaptive local models (reference points)
					mGen,				//!< current generation
					mNoNew;				//!< number of new solutions
	double			mUPFExtension,		//!< Utopian PF extension
					mPCAThreshold;		//!< threshold to determine the number of principal components
	//std::vector< unsigned int > mvNoChi,	//!< number of child points to generate around the seed
	//							mvDim;		//!< dimension of each subset
	//std::vector< METHOD >		mvMethod;	//!< method to generate solutions
	//std::vector< double >		mvNoDo,		//!< number of dominate points in each cluster
	//							mvStd;		//!< standard devision of noise in each subset
	//std::vector< alg::Matrix >	mvDatMat;				//!< matrix to store data for each 'cluster'
	//std::vector< alg::PCA > mvPCA;						//!< PCA
	enum			METHOD{MODEL, NOMETHOD};//!< methods used in this algorithm
	struct NODE
	{
		unsigned int no;			//!< number of child to generate in the subpopulation
		unsigned int dim;			//!< dimension of the subpopulation
		METHOD		 method;		//!< method for this subpopulation to generate new trial solutions
		double		 ratio;			//!< ratio of dominated soltuions in this subpopulation
		double		 std;			//!< standard devasion of noise in this subpopulation
		alg::Matrix  data;			//!< data of this subpopulation
		alg::PCA	 pca;			//!< PCA
	};
	std::vector< NODE > mvPops;		//!< subpopulations 
public:
	//!\brief	constructor
	//!\return	void
	RM3();

	//!\brief	deconstructor
	//!\return	void
	~RM3();

	//!\brief	Generator
	//!\param	sizenew	number of new trial solutions
	//!\param	popnew	offspring population
	//!\param	popref	reference population(current population)
	//!\return	offspring population
	CPopulationMO& Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref);

	//!\brief	set parameters
	//!\param	extension	extension ratio
	//!\param	threshold	principal components threshold
	//!\param	maxmodel	maximum number of models
	//!\param	maxgen		maximum number of generations
	//!\return	void
	void Set(double extension, double threshold, unsigned int maxmodel, unsigned int maxgen);		
protected:

	//!\brief	adaptive linear model-based generator
	//!\param	popnew	offspring population
	//!\param	index	index location in offspring population
	//!\param	popref	parent population
	//!\param	pca		PCA object
	//!\param	std		std. of noise
	//!\param	dim		dimension of manifold
	//!\param	size	the number of new points to create
	//!\return	offspring population
	CPopulationMO& ALModelGen(CPopulationMO& popnew, unsigned int& index, CPopulationMO& popref, alg::PCA& pca, double std, unsigned int dim, unsigned int size);	

	//!\brief	project point into a dimension
	//!\param	dim dimension index
	//!\param	index point index
	//!\param	pca PCA object
	//!\return	projection value
	double Project(unsigned int dim, unsigned int index, alg::PCA& pca);

	//!\brief	set model parameters
	//!\param	pop		reference population
	//!\param	sizenew size of new solutions
	//!\return	void
	void BuildModel(CPopulationMO& pop, unsigned int sizenew);

	//!\brief	partition the population
	//!\param	pop		reference population
	//!\return	void
	void Partition(CPopulationMO& pop);
}; //class RM3

} //namespace mod

} //namespace gen

} //namespace mea

} //namespace az

#endif //AZ_GENERATOR_MODEL_H
