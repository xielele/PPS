//	AlgD.cpp
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "AlgD.h"
#include "emo/Sel.h"
#include "alg/AR.h"
#include "alg/Matrix.h"
#include "emo/GenMod.h"

//#define SAVE_CEN 1

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
double		 T	= 0.0, T0;

// FDAs are from M. Farina, K. Deb and P. Amato's paper
//FDA1
void FDA1(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, G = sin(0.5*PI*T);
	for(unsigned int i=1; i<X.size(); i++)
		gx += (X[i]-G)*(X[i]-G);
	F[0] = X[0];
	F[1] = gx*(1-sqrt(F[0]/gx));
}
//FDA2
void FDA2(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, hx = 0.0, H = 0.75+0.7*sin(0.5*PI*T);
	unsigned int i,xii = (unsigned int)(X.size()/2);
	for(i=1; i<xii; i++) gx += X[i]*X[i];
	for(i=xii; i<X.size(); i++) hx += (X[i]-H)*(X[i]-H); hx += H;
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx, 1.0/hx));
}
//FDA3
void FDA3(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, Ft = pow(10.0,2.0*sin(0.5*PI*T)), Gt = fabs(sin(0.5*PI*T));
	unsigned int i;
	for(i=1; i<X.size(); i++) gx += (X[i]-Gt)*(X[i]-Gt); gx += Gt;
	F[0] = pow(X[0],Ft);
	F[1] = gx*(1-pow(F[0]/gx, 0.5));
}
//FDA4
void FDA4(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 0.0, G = fabs(sin(0.5*PI*T));
	unsigned int i;
	for(i=2; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
	F[0] = (1.0+gx)*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
	F[1] = (1.0+gx)*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
	F[2] = (1.0+gx)*sin(0.5*PI*X[0]);
}
// == DMOPs are from C.K Goh and K.C Tan 's paper
//dMOP1
void DMOP1(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 0.0, H = 0.75*sin(0.5*PI*T)+1.25;
	unsigned int i;
	for(i=1; i<X.size(); i++) gx += X[i]*X[i]; gx = 1.0 + 9.0*gx;
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx, H));
}
//dMOP2
void DMOP2(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, H = 0.75*sin(0.5*PI*T)+1.25, G = sin(0.5*PI*T);
	unsigned int i;
	for(i=1; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx, H));
}
//dMOP3
void DMOP3(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	double gx = 1.0, H = 0.75*sin(0.5*PI*T)+1.25, G = sin(0.5*PI*T);
	unsigned int i;
	for(i=1; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
	F[0] = X[0];
	F[1] = gx*(1-pow(F[0]/gx, 0.5));
}

// newly designed DMOPs problems
//centre moving functions
void LT(double t, unsigned int index, double& x1, double& x2)
{
	t *= 0.5;
	switch(index)
	{	
	case 0:	// Lissajous curve
		x1 = (cos(t*PI)+1.0)*2.0;
		x2 = (sin(2.0*t*PI)+1.0)*2.0;
		break;
	case 1:	// Rose curve
		x1 = (cos(1.5*t*PI)*sin(0.5*t*PI)+1.0)*2.0;
		x2 = (cos(1.5*t*PI)*cos(0.5*t*PI)+1.0)*2.0;
		break;
	case 2: // Heart curve
		x1 = ((1.0-sin(t*PI))*sin(t*PI)+2.0)*1.7;
		x2 = ((1.0-sin(t*PI))*cos(t*PI)+1.5)*1.4;
		break;
	case 3: // discontinus Lissajous curve
		t  = t-floor(t);
		x1 = (cos(t*PI)+1.0)*2.0;
		x2 = (sin(2.0*t*PI)+1.0)*2.0;
		break;
	default:
		x1 = x2 = 0.0;
		break;
	}
}
//shape moving functions
double HT(double t)
{
	return 1.25+0.75*sin(t*PI);
}
// problem framework
void DMOP(std::vector< double >& F, std::vector< double >& X, unsigned int index, bool turn)
{
	double a, b, Gi, gx1=0.0, gx2=0.0, ht=HT(T);
	
	bool old = (((unsigned int)(T*10.0+0.001)) % 2 == 1);
	
	LT(T, index, a, b);
	
	for(unsigned int i=1; i<X.size(); i++)
	{
		//if(turn && old)
		//	Gi = pow(fabs(2.0*X[0]-2.0*a-1.0), ht+double(i-1.0)/double(X.size()-2.0));
		//else	
		//	Gi = 1.0 - pow(fabs(2.0*X[0]-2.0*a-1.0), ht+double(i-1.0)/double(X.size()-2.0));

		if(turn && old)
			Gi = pow(fabs(X[0]-a), ht+double(i+1.0)/double(X.size()));
		else	
			Gi = 1.0 - pow(fabs(X[0]-a), ht+double(i+1.0)/double(X.size()));


		if(i % 2 == 1)
			gx1 += pow(X[i] - b - Gi ,2.0);
		else 
			gx2 += pow(X[i] - b - Gi ,2.0);
	}
	F[0] = pow(fabs(X[0]-a), ht)     + 0.5*gx1;
	F[1] = pow(fabs(X[0]-a-1.0), ht) + 0.5*gx2;
}
void DMOPA(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	DMOP(F,X,0,false);
}
void DMOPB(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	DMOP(F,X,1,false);
}
void DMOPC(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	DMOP(F,X,2,false);
}
void DMOPD(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	DMOP(F,X,3,false);
}
void DMOPE(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{
	DMOP(F,X,0,true);
}
void DMOPF(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
{	
	double gx = 0.0, G = sin(0.5*PI*T);
	for(unsigned int i=2; i<X.size(); i++) gx += pow(X[i]-pow(0.5*(X[0]+X[1]), HT(T)+double(i+1.0)/double(2.0*X.size()))-G,2.0);
	F[0] = (1.0+gx)*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
	F[1] = (1.0+gx)*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
	F[2] = (1.0+gx)*sin(0.5*PI*X[0]);
}

void DFRange(std::vector<double>& low, std::vector<double>& upp, std::string& name)
{
	if(	name == std::string("FDA1")  || name == std::string("FDA2")  || name == std::string("FDA3")||
		name == std::string("DMOP1") || name == std::string("DMOP2") || name == std::string("DMOP3"))
	{
		low[0] = 0;upp[0] =  1;
		for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] = -1; upp[i] = 1;}
	}
	else if(name == std::string("FDA4"))
	{
		low[0] = 0;upp[0] =  1;
		for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] =  0; upp[i] = 1;}
	}
	else if(name == std::string("DMOPF"))
	{
		low[0] = 0; upp[0] =  1;
		low[1] = 0; upp[1] =  1;
		for(unsigned int i=2; i<(unsigned int)(low.size()); i++){low[i] = -1.0; upp[i] = 2.0;}
	}
	else
	{	//DMOPA - DMOPE
		for(unsigned int i=0; i<(unsigned int)(low.size()); i++){low[i] = 0.0; upp[i] = 5.0;}
	}
}

DMOO::DMOO(
        unsigned int	strategy,
		std::string&	optimizer,
		unsigned int	popsize	,
		unsigned int	stepmax	,
		unsigned int	taot	,
		unsigned int	nt		,
		unsigned int	torder	,
		double			t0		,	
		CParameter&		par		)
		:mPop(par),mBest(par)
{
	switch(strategy)
	{
	case 1:
		mStrategy = PRANDI;
		break;
	case 2:
		mStrategy = PPOINT;
		break;
	case 3:
		mStrategy = PLINES;
		break;
	case 4:
		mStrategy = PCENTR;
		break;
	default:
		mStrategy = PRANDI;
		break;
	}
	mOptimizer = optimizer;

	mPopSize = popsize;
	mMaxStep = stepmax;
	mTaoT	 = taot;
	mDelT	 = 1.0/nt;
	T0		 = mT0 = t0;
	mMaxOrder= torder;
	pPar	 = &par;

	if(pPar->Problem() == std::string("FDA1"))
	{
		P().Evaluator( FDA1 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("FDA2"))
	{
		P().Evaluator( FDA2 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("FDA3"))
	{
		P().Evaluator( FDA3 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("FDA4"))
	{
		P().Evaluator( FDA4 );
		P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOP1"))
	{
		P().Evaluator( DMOP1 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOP2"))
	{
		P().Evaluator( DMOP2 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOP3"))
	{
		P().Evaluator( DMOP3 );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPA"))
	{
		P().Evaluator( DMOPA );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPB"))
	{
		P().Evaluator( DMOPB );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPC"))
	{
		P().Evaluator( DMOPC );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPD"))
	{
		P().Evaluator( DMOPD );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPE"))
	{
		P().Evaluator( DMOPE );
		P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
	}
	else if(pPar->Problem() == std::string("DMOPF"))
	{
		P().Evaluator( DMOPF );
		P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
	}

	Reset();
}

void DMOO::Reset()
{
	mStep		= 0;
	mEvas		= 0;
	T			= mT0;
	mbToChange	= false;
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

// get the current and predicted centers
bool DMOO::Centers(std::vector<double>& cc, std::vector<double>& pc) 
{
	if(hC.size()<1) return false;

	cc = *(hC.begin());

	if(IsPrediction()) pc=pC;
	else
	{
		pc.resize(cc.size()); for(unsigned int i=0; i<cc.size(); i++) pc[i] = -1.0;
	}
	return true;
}

// main evolution step
unsigned int DMOO::Step()
{
	unsigned int i,j;
	CPopulationMO popnew(P()),popold(P());

	mIsPre	= false;
	
	mbChange= true;

	switch(State())
	{
	// 0 history node (begining of the run, could be seen as a big jump)
	// => random initialize population
	case 0:
	case 1:
		Population().Resize(mPopSize);
		for(i=0; i<mPopSize; i++)
			for(j=0; j<P().XSize(); j++)
				Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
		Population().Evaluate();
		mEvas += Population().Size();
		break;
	//// 1 history node (a big jump just happened)
	//// => hypermutation + random initialization
	//case 1:
	//	popold = Population();
	//	for(i=0; i<mPopSize; i++)
	//	{
	//		// half population by random initialization
	//		if(i<mPopSize/2)
	//		{
	//			for(j=0; j<P().XSize(); j++)
	//				Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
	//		}
	//		// half by hypermutation
	//		else
	//		{
	//			for(j=0; j<P().XSize(); j++)
	//			{
	//				Population()[i][j] = popold[i][j] + rnd::gaussian()*(P().XUpp(j)-P().XLow(j))*0.25;
	//				if(Population()[i][j] < P().XLow(j))		Population()[i][j] = P().XLow(j) + 0.5*rnd::rand(0.0,P().XUpp(j)-P().XLow(j));
	//				else if(Population()[i][j] > P().XUpp(j))	Population()[i][j] = P().XUpp(j) - 0.5*rnd::rand(0.0,P().XUpp(j)-P().XLow(j));
	//			}
	//		}
	//	}
	//	Population().Evaluate();
	//	mEvas += Population().Size();
	//	break;
	// 2 or more history nodes
	// => prediction with history nodes
	case 2:
		// Predict new centre point & variance
		Predict();	// mC & mPV
		mIsPre = true;

		switch(mStrategy)
		{
			case PRANDI:
				{
					Population().Resize(mPopSize);
					for(i=0; i<mPopSize; i++) for(j=0; j<P().XSize(); j++) Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
				}
				break;
			case PPOINT:
				{
					// permutate the population
					Population().Shuffle();
					// predict
					for(i=0; i<P().FSize()+1; i++)
					{
						for(j=0; j<P().XSize(); j++) Population()[3*i+0][j] = pC[j]+mSD*(rnd::rand()>0.5?1.0:-1.0);
						for(j=0; j<P().XSize(); j++) Population()[3*i+1][j] = pC[j]+mSD*(rnd::rand()>0.5?1.0:-1.0);
						for(j=0; j<P().XSize(); j++) Population()[3*i+2][j] = pC[j];
					}
					unsigned int rr = (unsigned int)(0.3*(Population().Size()-P().FSize()-1.0));
					// 70% current + 30% random
					for(i=P().FSize()+1; i<P().FSize()+1+rr; i++)
					{
						for(j=0; j<P().XSize(); j++) Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
					}
				}
				break;
			case PLINES:
				{
					double len1,len2; unsigned int left;
					len1 = len2 = 0;
					for(j=0; j<P().XSize(); j++)
					{
						len1 += (pC[j]               - pC[j+P().XSize()])*(pC[j]               - pC[j+P().XSize()]);
						len2 += (pC[j+2*P().XSize()] - pC[j+P().XSize()])*(pC[j+2*P().XSize()] - pC[j+P().XSize()]);
					}
					left = (unsigned int)(len1*Population().Size()/(len1+len2));
					if(left<2) left = 2; if(left>Population().Size()-2) left = Population().Size()-2;
					for(i=0; i<left; i++)
					{
						for(j=0; j<P().XSize(); j++)
						{
							Population()[i][j] = pC[j] + (i+0.0)/(left+0.0)*(pC[j+P().XSize()]-pC[j]) + rnd::gaussian()*mSD;
						}
					}
					for(i=left; i<Population().Size(); i++)
					{
						for(j=0; j<P().XSize(); j++)
						{
							Population()[i][j] = pC[j+P().XSize()] + (i-left+0.0)/(Population().Size()-left-1.0)*(pC[j+2*P().XSize()]-pC[j+P().XSize()]) + rnd::gaussian()*mSD;
						}
					}
				}
				break;
			case PCENTR:
				{
					// generate new trial solutions based on the centre point and variance
					for(i=0; i<Population().Size(); i++)
						for(j=0; j<P().XSize(); j++)
						{
							Population()[i][j] = pC[j] + mBest[i][j] + rnd::gaussian()*mSD;
						}
				}
				break;
		}
		Check(Population());
		Population().Evaluate();
		mEvas += Population().Size();
		break;
	// no change
	default:
		// no change
		mbChange	= false;
		// calculate centre point
		// sample new solutions
		Generate(popnew, mPopSize);

		//check the range of new solutions
		Check(popnew);
		//remove these repeated solutions
		Check(popnew, Population());

		//evaluate new solutions
		popnew.Evaluate();
		mEvas += popnew.Size();

		//environmental select
		Population().Combine(popnew);
		az::mea::sel::SCrowd2 sel;
		sel.Select(Population(), mPopSize);
		break;
	}

	mStep++;

	//char chr[100];
	//sprintf(chr,"data\\%d.pop",mStep);
	//std::string sss(chr);
	//Population().Write(sss.c_str());

	//!!!! CHANGE STATE IF POSSIBLE
	// the change happens in (mTaoT*k + 1)th generation, k=1,2,...
	if(mStep > 1 && (mStep % mTaoT == 0))
	{
		T += mDelT;										// update time T
		DFRange(P().XLow(), P().XUpp(), P().Problem());	// update feasilbe search space
		mbToChange = true;

		//char chr[100];
		//sprintf(chr,"data\\%d.pop",mStep);
		//std::string sss(chr);
		//Population().Write(sss.c_str());
	}
	else
	{
		mbToChange = false;
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
	for(unsigned int i=0; i<popnew.Size(); i++) if(!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
	popnew.Clear();
	popnew.Combine(tmp);
}

// check current state
unsigned int DMOO::State()
{
	unsigned int i, j;

	// state 0: start of the running
	if(mStep == 0) return 0;

	// check wheter there is a change: search space & fitness value
	bool change = false;
	unsigned int in = rnd::rand((unsigned int)0, Population().Size());
	for(i=0; i<P().XSize(); i++)
		if(Population()[in][i]<P().XLow(i) || Population()[in][i]>P().XUpp(i))
		{
			change=true;
			break;
		}
	// state 3: no change happens
	if(!change) // still in the range
	{
		unsigned int count = 0; unsigned int cnum = (unsigned int)(Population().Size()*0.05); if(cnum<1) cnum = 1;
		for(i=0; i<cnum; i++)
		{
			CIndividualMO ind = Population()[in];
			ind.Evaluate();
			mEvas++;
			// no change detected
			if(fabs(ind.F(0)-Population()[in].F(0))+fabs(ind.F(1)-Population()[in].F(1))<1.0E-30) count++;
			else break;
		}
		if(count == cnum) return 3;
	}

	// otherwise there is some change
	switch(mStrategy)
	{
	case PRANDI:
		return 0;
		break;
	case PPOINT:
	case PLINES:
		{
			double dismin, dis;
			std::vector<unsigned int> pp(P().FSize()+1); for(i=0; i<pp.size(); i++) pp[i] = 0;
			mC.resize((P().FSize()+1)*P().XSize());
			for(i=1; i<Population().Size(); i++) 
			{
				for(j=0; j<pp.size()-1; j++) if(Population()[i].F(j) < Population()[pp[j]].F(j)) pp[j] = i;	// obj j 
			}
			dismin = 1.0E100;
			for(i=0; i<Population().Size(); i++) 
			{
				dis = 0;
				for(j=0; j<pp.size()-1; j++) dis += (Population()[i].F(j) - Population()[pp[j]].F(j))*(Population()[i].F(j) - Population()[pp[j]].F(j));
				if(dis<dismin) {pp[pp.size()-1] = i; dismin = dis;} // CTI (close-to-idea) point
			}
			for(i=0; i<pp.size(); i++) for(j=0; j<P().XSize(); j++) mC[i*P().XSize()+j] = Population()[pp[i]][j];
		}
		break;
	case PCENTR:
		{
			unsigned int j;
			mC.resize(P().XSize());
			for(i=0; i<P().XSize(); i++) mC[i] = 0.0;
			for(i=0; i<Population().Size(); i++) for(j=0; j<P().XSize(); j++) mC[j] += Population()[i][j];
			for(i=0; i<P().XSize(); i++) mC[i] /= double(Population().Size());
			
			mBest = Population();
			for(i=0; i<Population().Size(); i++) for(j=0; j<P().XSize(); j++) mBest[i][j] -= mC[j];
		}
		break;
	}
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
	// 
	return false;

	// only two history points
	if(hC.size()==2) return false;

	unsigned int dim;
	switch(mStrategy)
	{
		case PRANDI:
			return false;
			break;
		case PPOINT:
		case PLINES:
			dim = (P().FSize()+1)*P().XSize();
			break;
		case PCENTR:
			dim = P().XSize();
			break;
	}
	std::list< std::vector<double> >::iterator it2, it1, it0;
	it0 = it1 = it2 = hC.begin(); it1 ++; it2 ++; it2++;
	double cm = 0.0, lm = 0.0;
	for(unsigned int i=0; i<dim; i++)
	{
		cm += ((*it0)[i]-(*it1)[i])*((*it0)[i]-(*it1)[i]);		// current movement
		lm += ((*it2)[i]-(*it1)[i])*((*it2)[i]-(*it1)[i]);		// last movement
	}
	return cm>0.0001 && cm>9*lm;
}

// predict the location of the next centre point
void DMOO::Predict()
{
	unsigned int i,k,d,dim,order;
	std::list< std::vector<double> >::iterator it, it0;
	switch(mStrategy)
	{
		case PRANDI:
			dim = 0;
			break;
		case PPOINT:
		case PLINES:
			dim = (P().FSize()+1)*P().XSize();
			break;
		case PCENTR:
			dim = P().XSize();
			break;
	}

	// only two history nodes
	if(hC.size()==2)
	{
		it = it0 = hC.begin();
		it0++;
		mSD = 0.0;
		for(i=0; i<dim; i++)
		{
			mC[i] = 2.0*(*it)[i] - (*it0)[i];
			if(mC[i]>P().XUpp(i%P().XSize()))
				mC[i] = rnd::rand((*it)[i], P().XUpp(i%P().XSize()));
			else if(mC[i]<P().XLow(i%P().XSize()))
				mC[i] = rnd::rand(P().XLow(i%P().XSize()),(*it)[i]);
			mSD  += (mC[i]-(*it)[i])*(mC[i]-(*it)[i]);
		}
		mSD *= 2.0;
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

		mSD = 0.0;
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

			mSD  += (mC[i]-(*it)[i])*(mC[i]-(*it)[i]);
		}

		for(i=0; i<dim; i++) delete []px[i]; delete []px;
		for(i=0; i<dim; i++) delete []pa[i]; delete []pa;
		delete []pv;
	}
	mSD= 0.25*sqrt(mSD/double(dim)); //sqrt(mSD)/double(dim);//
	pC = mC;	// predicted centre

//============================================================================
#ifdef SAVE_CEN
	std::ofstream o2("pc.dat",std::ios::app);
	for(i=0; i<dim; i++) o2<<mC[i]<<"\t"; o2<<std::endl;
	o2.close();
#endif
//============================================================================
}

// generate new trial solutions by any offspring generator
CPopulationMO& DMOO::Generate(CPopulationMO& popnew, unsigned int size)
{
	if(mOptimizer == std::string("GTM"))
	{
		az::mea::gen::mod::ModelGTM2 gen;
		gen.Set(25,2,popnew.P().FSize()-1,30,0.25);
		gen.Generate(size, popnew, Population());
	}
	else if(mOptimizer == std::string("PCA"))
	{
		az::mea::gen::mod::RM gen;
		gen.Set(popnew.P().FSize()-1, 5, 30, 0.25);
		gen.Generate(size, popnew, Population());
	}
	else if(mOptimizer == std::string("NSDE"))
	{
		popnew.P().ETA_PM() = 20.0;
		popnew.P().ETA_SBX()= 20.0;
		popnew.P().Pc()		= 0.8;
		popnew.P().Pm()		= 0.05;
		az::mea::gen::XNSDE	gen;
		gen.Set(0.5,1.0);
		gen.Generate(size, popnew, Population());
	}
	else if(mOptimizer == std::string("NSGA"))
	{
		popnew.P().ETA_PM() = 20.0;
		popnew.P().ETA_SBX()= 20.0;
		popnew.P().Pc()		= 0.8;
		popnew.P().Pm()		= 0.05;
		az::mea::gen::XSBX gen;
		gen.Generate(size, popnew, Population());
	}

	return popnew;
}

} //namespace dea
} //namespace mea
} //namespace az

//==============================================================================================
// removed in Jul.30 2008
//==============================================================================================
//// newly designed DMOPs problems
////centre moving functions
//void LT(double t, unsigned int index, double& x1, double& x2)
//{
//	double r; unsigned int kk;
//	switch(index)
//	{	
//	case 0:	// line
//		x1 = 0;
//		x2 = sin(0.5*t*PI);
//		break;
//	case 1:	// rotated line
//		x1 = 0.5*t*PI;
//		x2 = 0.5*t*PI;
//		break;
//	case 2: // sin
//		x1 = 0.5*t*PI;
//		x2 = sin(0.5*t*PI);
//		break;
//	case 3: // heart curve (with circle)
//        r  = 0.5;
//        x1 = r*cos(t*PI)+2*r*cos(0.5*t*PI);
//        x2 = r*sin(t*PI)+2*r*sin(0.5*t*PI);
//		break;
//	case 4: // sin with jump
//		x1 = 0.5*t*PI;
//		x2 = sin(0.5*t*PI);
//		x2+= x2>=0 ? -1.0:1.0;
//		//kk = (unsigned int)(0.5*(t-T0));
//		//x2+= ((kk / 2) % 2 == 0) ? 1:-1; 
//		break;
//	default:
//		x1 = x2 = 0.0;
//		break;
//	}
//}
////shape moving functions
//double HT(double t)
//{
//	return 1.25+0.75*sin(0.5*t*PI);
//}
//// problem framework
//void DMOP(std::vector< double >& F, std::vector< double >& X, unsigned int index)
//{
//	double x1, x2;
//	LT(T, index, x1, x2);	
//	double gx = 0.0;
//	F[0] = X[0]-x1;
//	for(unsigned int i=1; i<X.size(); i++)
//	{
//		gx += pow(X[i]+pow(F[0], HT(T)+double(i-1.0)/double(X.size()-2)) - x2 - 1.0 ,2.0);
//	}
//	gx = 1.0 + gx/double(5);
//	F[1] = gx*(1.0 -  pow(F[0]/gx, HT(T)));
//}
//void DMOPA(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{
//	DMOP(F,X,0);
//}
//void DMOPB(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{
//	DMOP(F,X,1);
//}
//void DMOPC(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{
//	DMOP(F,X,2);
//}
//void DMOPD(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{
//	DMOP(F,X,3);
//}
//void DMOPE(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{
//	DMOP(F,X,4);
//}
//void DMOPF(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{	
//	double x1, x2, gx = 0.0;
//	LT(T, 2, x1, x2);	
//	unsigned int i;
//	for(i=2; i<X.size(); i++) gx += pow(X[i]+pow(X[0]-x1, HT(T)+double(i-2.0)/double(X.size()-3)) - x2 - 1.0 ,2.0);
//	F[0] = (1.0+gx)*cos(0.5*PI*(X[0]-x1))*cos(0.5*PI*(X[1]-x2));
//	F[1] = (1.0+gx)*cos(0.5*PI*(X[0]-x1))*sin(0.5*PI*(X[1]-x2));
//	F[2] = (1.0+gx)*sin(0.5*PI*(X[0]-x1));
//}
//void DMOPG(std::vector< double >& F, std::vector< double >& E, std::vector< double >& I, std::vector< double >& X)
//{	
//	double x1, x2, Gi;
//	bool old = (((unsigned int)(T*10.0)) % 2 == 1); 
//	LT(T, 2, x1, x2);	
//	double gx = 1.0;
//	F[0] = X[0]-x1;
//	for(unsigned int i=1; i<X.size(); i++)
//	{
//		if(old)
//			Gi = pow(F[0], HT(T)+double(i-1.0)/double(X.size()-2));
//		else
//			Gi = 1.0 - pow(F[0], HT(T)+double(i-1.0)/double(X.size()-2));
//		gx += pow(X[i]-x2-Gi,2.0);
//	}
//	F[1] = gx*(1.0 -  pow(F[0]/gx, HT(T)));
//}
//
//void DFRange(std::vector<double>& low, std::vector<double>& upp, std::string& name)
//{
//	if(	name == std::string("FDA1")  || name == std::string("FDA2")  || name == std::string("FDA3")||
//		name == std::string("DMOP1") || name == std::string("DMOP2") || name == std::string("DMOP3"))
//	{
//		low[0] = 0;upp[0] =  1;
//		for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] = -1; upp[i] = 1;}
//	}
//	else if(name == std::string("FDA4"))
//	{
//		low[0] = 0;upp[0] =  1;
//		for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] =  0; upp[i] = 1;}
//	}
//	else
//	{	//DMOPA - DMOPE
//		unsigned int index/* = name[4]-char('A')*/; double x1, x2;
//		if(name == std::string("DMOPA"))		index = 0;
//		else if(name == std::string("DMOPB"))	index = 1;
//		else if(name == std::string("DMOPC"))	index = 2;
//		else if(name == std::string("DMOPD"))	index = 3;
//		else if(name == std::string("DMOPE"))	index = 4;
//		else if(name == std::string("DMOPF"))	index = 2;
//		else if(name == std::string("DMOPG"))	index = 2;
//		LT(T,index,x1,x2);
//		low[0] = x1; upp[0] = x1+1.0;
//		for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] = x2; upp[i] = x2+1.0;}
//	}
//}
//==============================================================================================
//==============================================================================================