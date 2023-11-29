//	AlgD.cpp
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "AlgD.h"
#include "emo/Sel.h"
#include "alg/AR.h"
#include "alg/Matrix.h"

//#define SAVE_CEN 1
//#define PRE_POINT 1
#define PRE_MODEL 1

//!\brief	az namespace, the top namespace
namespace az
{
//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{
//!\brief namespace of dynamic evolutionary algoirhtm
namespace dea
{

const double PI = 3.141592653589793;
double		 T	= 0.0;
double		 TS = 0.0;

//centre moving functions
void mc(double t, unsigned int index, double& x1, double& x2)
{
	double r;
	switch(index)
	{
	case 1:	// line
		x1 = t*PI;
		x2 = t*PI;
		break;
	case 2: // sin
		x1 = t*PI;
		x2 = sin(t*PI);
		break;
	case 3: // sin with jump
		x1 = t*PI;
		x2 = sin(t*PI);
		x2+= (x2>=0) ? 0.5:-0.5;
		break;
	case 4: // heart curve (with circle)
        r  = 1;
        x1 = r*cos(2*t*PI)+2*r*cos(t*PI);
        x2 = r*sin(2*t*PI)+2*r*sin(t*PI);
		break;
	default:
		x1 = x2 = 0.0;
		break;
	}
}

//shape moving functions
double ms(double x1, double x2, double t, unsigned int index)
{
	double h;
	switch(index)
	{
	case 1:
        return x2 - 0.5*(sin(t*PI)+1);
	case 2:
        h  = 1.5+sin(t*PI);
        return x2 - pow(x1,h);
	default:
		return 0.0;
	}
}

//DMO problems
void DF1(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double x1, x2;
	mc(T, 1, x1, x2);
	F[0] = X[0]-x1;
	double gx= 1.0;
	for(unsigned int i=1; i<X.size(); i++)
		gx += pow(ms(F[0],X[i]-x2, T, 1),2.0);
	F[1] = gx*(1.0 -  pow(F[0]/gx, 1.5+sin(T*PI)));
}
void DF2(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double x1, x2;
	mc(T, 1, x1, x2);
	F[0] = X[0]-x1;
	double gx= 1.0;
	for(unsigned int i=1; i<X.size(); i++)
		gx += pow(ms(F[0],X[i]-x2, T, 2),2.0);
	F[1] = gx*(1.0 -  pow(F[0]/gx, 1.5+sin(T*PI)));
}
void DF3(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double x1, x2;
	mc(T, 2, x1, x2);
	F[0] = X[0]-x1;
	double gx= 1.0;
	for(unsigned int i=1; i<X.size(); i++)
		gx += pow(ms(F[0],X[i]-x2, TS, 1),2.0);
	F[1] = gx*(1.0 -  pow(F[0]/gx, 1.5+sin(TS*PI)));
}
void DF4(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double x1, x2;
	mc(T, 2, x1, x2);
	F[0] = X[0]-x1;
	double gx= 1.0;
	for(unsigned int i=1; i<X.size(); i++)
		gx += pow(ms(F[0],X[i]-x2, TS, 2),2.0);
	F[1] = gx*(1.0 -  pow(F[0]/gx, 1.5+sin(TS*PI)));
}
void DF5(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double x1, x2;
	mc(T, 3, x1, x2);
	F[0] = X[0]-x1;
	double gx= 1.0;
	for(unsigned int i=1; i<X.size(); i++)
		gx += pow(ms(F[0],X[i]-x2, T, 1),2.0);
	F[1] = gx*(1.0 -  pow(F[0]/gx, 1.5+sin(T*PI)));
}
void DF6(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double x1, x2;
	mc(T, 3, x1, x2);
	F[0] = X[0]-x1;
	double gx= 1.0;
	for(unsigned int i=1; i<X.size(); i++)
		gx += pow(ms(F[0],X[i]-x2, T, 2),2.0);
	F[1] = gx*(1.0 -  pow(F[0]/gx, 1.5+sin(T*PI)));
}
void DF7(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double x1, x2;
	mc(T, 4, x1, x2);
	F[0] = X[0]-x1;
	double gx= 1.0;
	for(unsigned int i=1; i<X.size(); i++)
		gx += pow(ms(F[0],X[i]-x2, T, 1),2.0);
	F[1] = gx*(1.0 -  pow(F[0]/gx, 1.5+sin(T*PI)));
}
void DF8(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double x1, x2;
	mc(T, 4, x1, x2);
	F[0] = X[0]-x1;
	double gx= 1.0;
	for(unsigned int i=1; i<X.size(); i++)
		gx += pow(ms(F[0],X[i]-x2, T, 2),2.0);
	F[1] = gx*(1.0 -  pow(F[0]/gx, 1.5+sin(T*PI)));
}

void DFRange(std::vector<double>& low, std::vector<double>& upp, std::string& name)
{
	double x1, x2;
	if(name == std::string("DF1") || name == std::string("DF2"))
		mc(T, 1, x1, x2);
	else if(name == std::string("DF3") || name == std::string("DF4"))
		mc(T, 2, x1, x2);
	else if(name == std::string("DF5") || name == std::string("DF6"))
		mc(T, 3, x1, x2);
	else if(name == std::string("DF7") || name == std::string("DF8"))
		mc(T, 4, x1, x2);

	low[0] = x1; upp[0] = x1+1.0;
	for(unsigned int i=1; i<(unsigned int)(low.size()); i++)
	{
		low[i] = x2;
		upp[i] = x2+1.0;
	}
}

//// modified FDA problem
//void ZZJ(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{
//	F[0] = X[0];
//	double Gt = sin(PI*(T+1.5));
//	double Ht = 1.5 + sin(PI*(T+1.5));
//	double gx = 1.0;
//	for(unsigned int i=1; i<X.size(); i++)
//		gx +=  pow(3*X[i]-1.0 + Gt - pow(X[0],Ht),2.0);
//	F[1] = gx*(1.0 -  pow(F[0]/gx, Ht));
//}
////	FDA1 test problem
//void FDA1(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{
//	double Gt= sin(0.5*PI*(T+1.0));
//	F[0] = X[0];
//	double gx= 1.0;
//	for(unsigned int i=1; i<X.size(); i++)
//		gx += pow(X[i]-Gt,2.0);
//	F[1] = gx*(1.0 -  pow(F[0]/gx, 0.5));
//}

DMOO::DMOO(
		unsigned int	popsize	 ,
		unsigned int	stepmax	 ,
		unsigned int	tinterval,
		double			t0		 ,
		double			tstep	 ,
		unsigned int	torder	 ,
		unsigned int	noLatent ,
		unsigned int	noBaseFun,
		unsigned int	dimLatent,
		unsigned int	trainsteps,
		double			extension,
		CParameter&		par		 )
		:mPop(par)
{
	mPopSize = popsize;
	mMaxStep = stepmax;
	mTaoT	 = tinterval;
	mT0		 = t0;
	mDeltT	 = tstep;
	mMaxOrder= torder;
	pPar	 = &par;
	mExtension	= extension;
	mLatentDim	= dimLatent;
	mNoLatent	= noLatent;
	mNoBaseFun	= noBaseFun;
	mTrainSteps = trainsteps;

	if(pPar->Problem() == std::string("DF1"))
		P().Evaluator( DF1 );
	if(pPar->Problem() == std::string("DF2"))
		P().Evaluator( DF2 );
	if(pPar->Problem() == std::string("DF3"))
		P().Evaluator( DF3 );
	if(pPar->Problem() == std::string("DF4"))
		P().Evaluator( DF4 );
	if(pPar->Problem() == std::string("DF5"))
		P().Evaluator( DF5 );
	if(pPar->Problem() == std::string("DF6"))
		P().Evaluator( DF6 );
	if(pPar->Problem() == std::string("DF7"))
		P().Evaluator( DF7 );
	if(pPar->Problem() == std::string("DF8"))
		P().Evaluator( DF8 );
	//if(pPar->Problem() == std::string("ZZJ"))
	//	P().Evaluator( ZZJ );
	//else
	//	P().Evaluator( FDA1 );
	P().FSize( 2 );
	P().ESize( 0 );
	P().ISize( 0 );

	Reset();
}

void DMOO::Reset()
{
	mStep	 = 0;
	mEvas	 = 0;
	T		 = mT0;
	TS		 = 0.0;
	hC.clear();
	DFRange(P().XLow(), P().XUpp(), P().Problem());	// set the original search space

	az::rnd::seed((long) time(NULL));

//============================================================================
#ifdef SAVE_CEN
	std::ofstream o1("pc.dat");
	o1.close();
	std::ofstream o2("rc.dat");
	o2.close();
#endif
//============================================================================
}

// main evolution step
unsigned int DMOO::Step()
{
	unsigned int i,j;
#ifdef PRE_POINT
	int p1,p2,p3;
#endif
	CPopulationMO popnew(P()),popold(P());

	mIsPre = false;

	switch(State())
	{
	// 0 history node (begining of the run, could be seen as a big jump)
	// => random initialize population
	case 0:
		Population().Resize(mPopSize);
		for(i=0; i<mPopSize; i++)
			for(j=0; j<P().XSize(); j++)
				Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
		Population().Evaluate();
		mEvas += Population().Size();
		break;
	// 1 history node (a big jump just happened)
	// => hypermutation + random initialization
	case 1:
		popold = Population();
		for(i=0; i<mPopSize; i++)
		{
			// random
			if(rnd::rand()<0.5)
			{
				for(j=0; j<P().XSize(); j++)
					Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
			}
			// hypermutation
			else
			{
				for(j=0; j<P().XSize(); j++)
				{
					Population()[i][j] = popold[i][j] + rnd::gaussian()*(P().XUpp(j)-P().XLow(j))*0.25;
					if(Population()[i][j] < P().XLow(j))		Population()[i][j] = P().XLow(j) + 0.5*rnd::rand(0.0,P().XUpp(j)-P().XLow(j));
					else if(Population()[i][j] > P().XUpp(j))	Population()[i][j] = P().XUpp(j) - 0.5*rnd::rand(0.0,P().XUpp(j)-P().XLow(j));
				}
			}
		}
		Population().Evaluate();
		mEvas += Population().Size();
		break;
	// 2 or more history nodes
	// => prediction with history nodes
	case 2:
		// Predict new centre point & variance
		Predict();	// mC & mPV
		mIsPre = true;
#ifdef PRE_POINT
		// 70% current + 24% random + 6% predict
		for(i=(unsigned int)0.7*Population().Size(); i<Population().Size()-6; i++)
		{
			for(j=0; j<P().XSize(); j++)
				Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
		}
		
		i = Population().Size()-6;
		for(j=0; j<P().XSize(); j++)
			Population()[i][j] = pC[j]+mPV*rnd::gaussian();
		i = Population().Size()-5;
		for(j=0; j<P().XSize(); j++)
			Population()[i][j] = pC[j]+mPV*rnd::gaussian();
		i = Population().Size()-4;
		for(j=0; j<P().XSize(); j++)
			Population()[i][j] = pC[j+P().XSize()]+mPV*rnd::gaussian();
		i = Population().Size()-3;
		for(j=0; j<P().XSize(); j++)
			Population()[i][j] = pC[j+P().XSize()]+mPV*rnd::gaussian();
		i = Population().Size()-2;
		for(j=0; j<P().XSize(); j++)
			Population()[i][j] = pC[j+2*P().XSize()]+mPV*rnd::gaussian();
		i = Population().Size()-1;
		for(j=0; j<P().XSize(); j++)
			Population()[i][j] = pC[j+2*P().XSize()]+mPV*rnd::gaussian();

		Check(Population());
#elif PRE_MODEL
		// set variance
		if(mPV<0.01)	mBeta = 1.0/0.01;
		else			mBeta = 1.0/mPV;

		// generate new trial solutions based on the centre point and variance
		ModelGen1D(Population(), mPopSize);
#endif

		Population().Evaluate();
		mEvas += Population().Size();
		break;
	// no change
	default:
		// calculate centre point
#ifdef PRE_POINT
		mC.resize(3*P().XSize());
		p1 = p2 = p3 = -1;
		for(i=0; i<Population().Size(); i++) 
		{
			if(p1<0 || Population()[i].F(0) < Population()[p1].F(0)) p1 = i;
			if(p2<0 || Population()[i].F(1) < Population()[p2].F(1)) p2 = i;
		}
		for(i=0; i<Population().Size(); i++) 
		{
			if(p3<0 || 
				(Population()[i].F(0) - Population()[p1].F(0))*(Population()[i].F(0) - Population()[p1].F(0)) + 
				(Population()[i].F(1) - Population()[p2].F(1))*(Population()[i].F(1) - Population()[p2].F(1)) <
				(Population()[p3].F(0) - Population()[p1].F(0))*(Population()[p3].F(0) - Population()[p1].F(0)) + 
				(Population()[p3].F(1) - Population()[p2].F(1))*(Population()[p3].F(1) - Population()[p2].F(1)))
				p3 = i;
		}
		for(i=0; i<P().XSize(); i++) mC[i]				 = Population()[p1][i];
		for(i=0; i<P().XSize(); i++) mC[i+P().XSize()]	 = Population()[p2][i];
		for(i=0; i<P().XSize(); i++) mC[i+2*P().XSize()] = Population()[p3][i];

		// assign new data to GTM
		mT.Resize(Population().Size(),P().XSize());
		for( i=0; i<Population().Size(); i++ ) for( j=0; j<P().XSize(); j++ ) mT(i,j) = Population()[i][j];
#elif PRE_MODEL
		mC.resize(P().XSize());
		for(i=0; i<P().XSize(); i++) mC[i] = 0.0;
		for(i=0; i<Population().Size(); i++) for(j=0; j<P().XSize(); j++) mC[j] += Population()[i][j];
		for(i=0; i<P().XSize(); i++) mC[i] /= double(Population().Size());

		// assign new data to GTM
		mT.Resize(Population().Size(),P().XSize());
		for( i=0; i<Population().Size(); i++ ) for( j=0; j<P().XSize(); j++ ) mT(i,j) = Population()[i][j]-mC[j];
#endif

		// reset the weights of GTM
		//if(mLatentDim==1)
		mGTM.Initialize1(mFI,mW,mBeta,mT,mNoLatent,mNoBaseFun,2.0);
		// train GTM model
		mGTM.Train(mW,mBeta,mT,mFI,mTrainSteps);
		// sample new solutions
		ModelGen1D(popnew, mPopSize);

		//check the range of new solutions
		Check(popnew);
		//remove these repeated solutions
		Check(popnew, Population());

		//evaluate new solutions
		popnew.Evaluate();
		mEvas += popnew.Size();

		//environmental select
		Population().Combine(popnew);
#define CES az::mea::sel::SCrowd2
		CES sel;
		sel.Select(Population(), mPopSize);
		break;
	}

	mStep++;

	//!!!! CHANGE STATE IF POSSIBLE
	// the change happens in (mTaoT*k + 1)th generation, k=1,2,...
	if(mStep > 1 && (mStep % mTaoT == 0))
	{
		TS+= 0.1;
		T += mDeltT;									// update time T
		DFRange(P().XLow(), P().XUpp(), P().Problem());	// update feasilbe search space
	}

	return mStep;
}

// check the range of each individual
void DMOO::Check(CPopulationMO& pop)
{
	for(unsigned int i=0; i<pop.Size(); i++) pop[i].Check();
}

// delete duplicate individuals
void DMOO::Check(CPopulationMO& popnew, CPopulationMO& pop)
{
	CPopulationMO tmp(P());
	for(unsigned int i=0; i<popnew.Size(); i++)
		if(!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
	popnew.Clear();
	popnew.Combine(tmp);
}

// check current state
unsigned int DMOO::State()
{
	unsigned int i;
	// state 0: start of the running
	if(mStep == 0) return 0;

	// check wheter there is a change: search space & fitness value
	unsigned int in = rnd::rand((unsigned int)0, Population().Size());
	bool change = false;
	for(i=0; i<P().XSize(); i++)
		if(Population()[in][i]<P().XLow(i) || Population()[in][i]>P().XUpp(i))
		{
			change=true;
			break;
		}
	// state 3: no change happens
	if(!change) // still in the range
	{
		CIndividualMO ind = Population()[in];
		ind.Evaluate();
		mEvas++;
		// no change detected
		if(fabs(ind.F(0)-Population()[in].F(0))+fabs(ind.F(1)-Population()[in].F(1))<1.0E-5) return 3;
	}

	// otherwise there is some change
	hC.insert(hC.begin(),mC);
//============================================================================
#ifdef SAVE_CEN
	std::ofstream o1("rc.dat",std::ios::app);
	for(unsigned int i=0; i<mC.size(); i++) o1<<mC[i]<<"\t"; o1<<std::endl;
	o1.close();
#endif
//============================================================================

	// state 1: a big jump happens
	if(Jump())
	{
		std::list< std::vector<double> >::iterator it = hC.begin();
		hC.erase(++it, hC.end());	// history nodes are useless when a big jump happens
		return 1;
	}

	// state 2: a normal change happens
	return 2;
}

// detect whether the last change is a big jump
bool DMOO::Jump()
{
	// only one history point
	if(hC.size()<2)  return true;

	// only two history points
	if(hC.size()==2) return false;

	//return false;

	double dp = 0.0, dr = 0.0;
#ifdef PRE_POINT
	unsigned int dim = 3*P().XSize();
#elif PRE_MODEL
	unsigned int dim = P().XSize();
#endif
	std::list< std::vector<double> >::iterator it, it0;
	it0 = it = hC.begin(); it0++;
	for(unsigned int i=0; i<dim; i++)
	{
		dp += ((*it0)[i]-pC[i])*((*it0)[i]-pC[i]);			// predicted moving variance
		dr += ((*it)[i]-(*it0)[i])*((*it)[i]-(*it0)[i]);	// real moving variance
	}

	return dp>0.0001 && dr>25*dp;
}

// predict the location of the next centre point
void DMOO::Predict()
{
	unsigned int i,k,d,dim,order;
	std::list< std::vector<double> >::iterator it, it0;

#ifdef PRE_POINT
	dim = 3*P().XSize();
#elif PRE_MODEL
	dim = P().XSize();
#endif

	// only two history nodes
	if(hC.size()==2)
	{
		it = it0 = hC.begin();
		it0++;
		mPV = 0.0;
		for(i=0; i<dim; i++)
		{
			mC[i] = 2.0*(*it)[i] - (*it0)[i];
			if(mC[i]>P().XUpp(i%P().XSize()))
				mC[i] = rnd::rand((*it)[i], P().XUpp(i%P().XSize()));
			else if(mC[i]<P().XLow(i%P().XSize()))
				mC[i] = rnd::rand(P().XLow(i%P().XSize()),(*it)[i]);
			mPV  += (mC[i]-(*it)[i])*(mC[i]-(*it)[i]);
		}
	}
	// more than two history nodes
	else
	{
		order = std::min((unsigned int)hC.size()-1, (unsigned int)mMaxOrder);
		double **px, **pa, *pv;
		px = new double*[dim]; for(i=0; i<dim; i++) px[i] = new double[hC.size()];
		pa = new double*[dim]; for(i=0; i<dim; i++) pa[i] = new double[order];
		pv = new double[dim];
		k  = (unsigned int)hC.size()-1;
		it = hC.begin();
		d  = 0;
		while(it!=hC.end())// && d<= 2*order+10) // maximum (order+3) history nodes are used
		{
			for(i=0; i<dim; i++) px[i][k] = (*it)[i];
			k--;
			it++;
			d++;
		}
		alg::aruv(px, dim, (unsigned int)hC.size(), order, pa, pv);

		std::vector<double> tc(dim);
		for(i=0; i<dim; i++) tc[i] = 0.0;
		it = hC.begin();
		for(k=0; k<order; k++)
		{
			for(i=0; i<dim; i++) tc[i] += (*it)[i]*pa[i][k];
			it++;
		}

		mPV = 0.0;
		it = it0 = hC.begin(); it0++;
		for(i=0; i<dim; i++)
		{
			mC[i] = tc[i];
			// not allowd too big prediction
			if(fabs(tc[i]-(*it)[i])>10*fabs((*it)[i]-(*it0)[i]))
				mC[i] = (*it)[i] + rnd::rand(0.0, 10*((*it)[i]-(*it0)[i]));
			// keep it in the boundary
			if(mC[i]>P().XUpp(i%P().XSize()))
				mC[i] = rnd::rand((*it)[i], P().XUpp(i%P().XSize()));
			else if(mC[i]<P().XLow(i%P().XSize()))
				mC[i] = rnd::rand(P().XLow(i%P().XSize()),(*it)[i]);

			mPV  += (mC[i]-(*it)[i])*(mC[i]-(*it)[i]);
		}

		for(i=0; i<dim; i++) delete []px[i]; delete []px;
		for(i=0; i<dim; i++) delete []pa[i]; delete []pa;
		delete []pv;
	}
	pC = mC;	// predicted centre

//============================================================================
#ifdef SAVE_CEN
	std::ofstream o2("pc.dat",std::ios::app);
	for(i=0; i<dim; i++) o2<<mC[i]<<"\t"; o2<<std::endl;
	o2.close();
#endif
//============================================================================
}

// generate new trial solution with GTM based method
CPopulationMO& DMOO::ModelGen1D(CPopulationMO& popnew, unsigned int size)
{
	unsigned int n,k,d;
	double dist;

	popnew.Resize(size);

	alg::Matrix X(size, 1), MU(mNoBaseFun, 1);

	//assign latent variables
	for(n=0; n<size; n++)
		X(n,0) = rnd::rand(double(n),double(n+1))*(2.0+2*mExtension)/double(size)-1.0-mExtension;
	//centre points of base functions
	for(n=0; n<mNoBaseFun; n++)
		MU(n, 0)= (-1.0 + double(n)*2.0/double(mNoBaseFun-1.0))*double(mNoBaseFun)/double(mNoBaseFun-1.0);
	//variance of base functions
	double sigma = 2.0 * (MU(1,0) - MU(0,0));

	alg::Matrix MDIS(X.RowSize(), MU.RowSize());

  	for (n = 0; n < MU.RowSize(); n++)
		for (k = 0; k < X.RowSize(); k++)
		{
			MDIS(k, n) = 0.0;
			for (d = 0; d<X.ColSize(); d++)
      		{
				dist = MU(n, d)-X(k, d);
				MDIS(k,n) += dist*dist;
			}
		}

	alg::Matrix FI(MDIS.RowSize(), MDIS.ColSize()+1);
    double tmp = -1.0/(2*sigma*sigma);
	for(n=0; n<MDIS.RowSize(); n++)
		for(k=0; k<MDIS.ColSize(); k++)
			FI(n,k) = exp(MDIS(n,k)*tmp);

	//Add bias basis function
	for(n=0; n<FI.RowSize(); n++) FI(n,FI.ColSize()-1) = 1.0;

	//create new solutions
	//Gaussian noise variance
	double radius  = 1.0 / sqrt(mBeta);
	for(n=0; n<popnew.Size(); n++)
	{
		for(d=0; d<popnew.P().XSize(); d++)
		{
#ifdef PRE_POINT
			popnew[n][d] = radius*rnd::gaussian();
#elif PRE_MODEL
			popnew[n][d] = mC[d] + radius*rnd::gaussian();
#endif
			for(k=0; k<FI.ColSize(); k++) popnew[n][d] += FI(n,k)*mW(k,d);
			if(popnew[n][d] < popnew.P().XLow(d))		popnew[n][d] = P().XLow(d) + 0.5*rnd::rand(0.0, P().XUpp(d)-P().XLow(d));
			else if(popnew[n][d] > popnew.P().XUpp(d))	popnew[n][d] = P().XUpp(d) - 0.5*rnd::rand(0.0, P().XUpp(d)-P().XLow(d));
		}
	}

	return popnew;
}

} //namespace dea

} //namespace mea

} //namespace az
