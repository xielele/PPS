//	AlgD.cpp
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "AlgD.h"
#include "emo/Sel.h"
#include "alg/AR.h"
#include "alg/BPTSeries.h"
#include "alg/Fitting.h"
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

//std::ofstream tmp;//

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
		double			alpha   ,
		CParameter&		par		)
		:mPop(par),mBest0(par),mBest1(par)
{
	switch(strategy)
	{
	case 1:
		mStrategy = INIRIS;
		break;
	case 2:
		mStrategy = INIFPS;
		break;
	case 3:
		mStrategy = INIPPS;
		break;
	case 4:
		mStrategy = INIPPSBP;
		break;
	case 5:
		mStrategy = INIPPSNR;
		break;
	default:
		mStrategy = INIRIS;
		break;
	}
	mOptimizer = optimizer;

	mPopSize = popsize;
	mMaxStep = stepmax;
	mTaoT	 = taot;
	mDelT	 = 1.0/nt;
	T0		 = mT0 = t0;
	mMaxOrder= torder;
	mAlpha   = alpha;
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

	//tmp.open("tmp.dat");
}

DMOO::~DMOO()
{
//	tmp.close();
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
}

// main evolution step
unsigned int DMOO::Step()
{
	// environment change
	if(EnvironmentChange())
	{
		switch(mStrategy)
		{
			case INIRIS:
				InitRIS(true);
				break;
			case INIFPS:
				InitFPS();
				break;
			case INIPPS:
			case INIPPSBP:
			case INIPPSNR:
				InitPPS();
				break;
		}
		Population().Evaluate(); mEvas += Population().Size();
	}
	// no environment change
	else
	{
		CPopulationMO popnew(P());

		// generate new solutions
		Generate(popnew, mPopSize);

		//check the range of new solutions
		Check(popnew);
		//remove these repeated solutions
		Check(popnew, Population());

		//evaluate new solutions
		popnew.Evaluate(); mEvas += popnew.Size();

		//environmental select
		Population().Combine(popnew);
		az::mea::sel::SCrowd2 sel;
		sel.Select(Population(), mPopSize);
	}
	mStep++;

	//!!!! CHANGE STATE IF POSSIBLE
	// the change happens in (mTaoT*k + 1)th generation, k=1,2,...
	if(mStep > 1 && (mStep % mTaoT == 0))
	{
		T += mDelT;										// update time T
		mbToChange = true;
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

bool DMOO::EnvironmentChange()
{
	if(mStep == 0) return true;

	CIndividualMO ind(P());
	unsigned int i, j, in, cnum = (unsigned int)(Population().Size()*0.05); if(cnum<1) cnum = 1;
	for(i=0; i<cnum; i++)
	{
		in  = rnd::rand((unsigned int)0, Population().Size());
		ind = Population()[in];
		ind.Evaluate(); mEvas++;
		for(j=0; j<Population().P().FSize(); j++) if(fabs(ind.F(j)-Population()[in].F(j))>1.0E-10) return true;
	}
	return false;
}

// predict the location of the next centre point
void DMOO::Predict()
{
	switch(mStrategy)
	{
		case INIPPS:
		{
			unsigned int i, k, d, order, dim=(unsigned int)(mC.size());
			std::list< std::vector<double> >::iterator it, it0;

			while(hC.size()>20+mMaxOrder) hC.pop_back();

			order = std::min((unsigned int)hC.size()-1, (unsigned int)mMaxOrder);
			double **px, **pa, *pv;
			px = new double*[dim]; for(i=0; i<dim; i++) px[i] = new double[hC.size()];
			pa = new double*[dim]; for(i=0; i<dim; i++) pa[i] = new double[order];
			pv = new double[dim];
			k  = (unsigned int)hC.size()-1;
			it = hC.begin();
			d  = 0;
			while(it!=hC.end())
			{
				for(i=0; i<dim; i++) px[i][k] = (*it)[i];
				k--;
				it++;
				d++;
			}
			alg::aruv(px, dim, (unsigned int)hC.size(), order, pa, pv);
			mStdC.resize(dim); for(i=0; i<dim; i++) mStdC[i] = pv[i];

			pC.resize(mC.size()); for(i=0; i<dim; i++) pC[i] = 0.0;
			it = hC.begin();
			for(k=0; k<order; k++)
			{
				for(i=0; i<dim; i++) pC[i] += (*it)[i]*pa[i][k];
				it++;
			}
			it = hC.begin();
			for(i=0; i<dim; i++)
			{
				// keep it in the boundary
				if(     pC[i]>P().XUpp(i%P().XSize())) pC[i] = (*it)[i]; //P().XUpp(i%P().XSize());
				else if(pC[i]<P().XLow(i%P().XSize())) pC[i] = (*it)[i]; //P().XLow(i%P().XSize());
				if(     pC[i]>P().XUpp(i%P().XSize())) pC[i] = P().XUpp(i%P().XSize());
				else if(pC[i]<P().XLow(i%P().XSize())) pC[i] = P().XLow(i%P().XSize());
			}

			for(i=0; i<dim; i++) delete []px[i]; delete []px;
			for(i=0; i<dim; i++) delete []pa[i]; delete []pa;
			delete []pv;
		}
		break;
		case INIPPSBP:
		{
			unsigned int i, d, dim=(unsigned int)(mC.size());
			std::list< std::vector<double> >::iterator it, it0;
			while(hC.size()>20+mMaxOrder) hC.pop_back();
			mStdC.resize(dim);
			pC.resize(dim);
			double *pd = new double[hC.size()];
			for(d=0; d<dim; d++)
			{
				// copy data
				it = hC.begin(); i = 0;
				while(it != hC.end())
				{
					pd[i++] = (*it)[d];
					it++;
				}
				// train
				pC[d] = bpts::BPPredict(pd, hC.size(), mStdC[d]);
				// boundary checking
				it = hC.begin();
				if(     pC[d]>P().XUpp(d%P().XSize())) pC[d] = (*it)[d];
				else if(pC[d]<P().XLow(d%P().XSize())) pC[d] = (*it)[d];
				if(     pC[d]>P().XUpp(d%P().XSize())) pC[d] = P().XUpp(d%P().XSize());
				else if(pC[d]<P().XLow(d%P().XSize())) pC[d] = P().XLow(d%P().XSize());
			}
			delete []pd;
		}
		break;
		case INIPPSNR:
		{
			unsigned int i, k, d, order, num, dim=(unsigned int)(mC.size());
			std::list< std::vector<double> >::iterator it, it1, it2;
			while(hC.size()>20+mMaxOrder) hC.pop_back();
			mStdC.resize(dim);
			pC.resize(dim);
			order	= 2;
			num		= hC.size();
			double *pd	= new double[num-2];
			double *pt1 = new double[num-2];
			double *pt2	= new double[num-2];
			double *pc	= new double[6];
			for(d=0; d<dim; d++)
			{
				// copy data
				it = it1 = it2 = hC.begin(); i = 0;
				it2++; it++; it++;
				while(it != hC.end())
				{

					pt1[i]	= (*it1)[d];
					pt2[i]	= (*it2)[d];
					pd[i]	= (*it)[d];
					it++; it1++; it2++; i++;
				}
				// train
				az::alg::poly_fit(pt1, pt2, pd, num-2, order, pc);
				// prediction
				mStdC[d] = 0.0;
				for(i=0; i<num-2; i++)
				{
					double val = pc[0] + pc[1] * pt1[i] + pc[2] * pt2[i] + pc[3] * pt1[i] * pt1[i] + pc[4] * pt1[i] * pt2[i] + pc[5] * pt2[i] * pt2[i];
					mStdC[d] += (val-pd[i])*(val-pd[i]);
				}
				mStdC[d] = sqrt(mStdC[d]/num);
				pC[d] = 0.0;
				double t1 = pt2[num-2], t2 = pd[num-2];
				double val = pc[0] + pc[1] * t1 + pc[2] * t2 + pc[3] * t1 * t1 + pc[4] * t1 * t2 + pc[5] * t2 * t2;

				for(k=0; k<=order; k++) pC[d] += pow(num+1.0, k+0.0)*pc[k];
				// boundary checking
				it = hC.begin();
				if(     pC[d]>P().XUpp(d%P().XSize())) pC[d] = (*it)[d];
				else if(pC[d]<P().XLow(d%P().XSize())) pC[d] = (*it)[d];
				if(     pC[d]>P().XUpp(d%P().XSize())) pC[d] = P().XUpp(d%P().XSize());
				else if(pC[d]<P().XLow(d%P().XSize())) pC[d] = P().XLow(d%P().XSize());
			}
			delete []pd;
			delete []pt1;
			delete []pt2;
			delete []pc;
		}
	}
}

void DMOO::InitPPS()
{	
	//Step 0: initialize the population
	if(mStep == 0) {InitRIS(true); return;}

	unsigned int i, j, k, xdim=P().XSize(), psize = Population().Size();

	//Step 1: save center point and shape points
	mC.resize(xdim); 
	for(k=0; k<xdim; k++) 
	{
		mC[k]  = 0.0;
		for(i=0; i<psize; i++) mC[k] += Population()[i][k]; 
		mC[k] /= double(psize);
		//tmp<<mC[k]<<" ";
	}
	//tmp<<std::endl;

	hC.insert(hC.begin(), mC);
	
	mBest1 = mBest0;
	mBest0 = Population();
	for(i=0; i<psize; i++) for(j=0; j<xdim; j++) mBest0[i][j] -= mC[j];

	//Step 2: if the history info. is not enough for prediction, then use random initialization
	if(hC.size()<2*mMaxOrder) {InitRIS(false); return;}

	//Step 3: find the parent point for each point in mBest0
	std::vector<unsigned int> pindex(psize);
	double dismin, dis;
	for(i=0; i<psize; i++)
	{
		pindex[i] = 0; dismin = 1.0E100;
		for(j=0; j<psize; j++)
		{
			dis = 0.0; for(k=0; k<xdim; k++) dis += (mBest0[i][k] - mBest1[j][k])*(mBest0[i][k] - mBest1[j][k]);
			if(dis<dismin) {dismin = dis; pindex[i] = j;}
		}
	}
	//std::vector<double> mu(xdim); 
	//for(k=0; k<xdim; k++) 
	//{
	//	mu[k]  = 0.0; 
	//	for(i=0; i<psize; i++) mu[k] += mBest0[i][k] - mBest1[pindex[i]][k];
	//	mu[k] /= double(psize);
	//}
	double mstd = 0.0;
	for(k=0; k<xdim; k++)
	{
		//for(i=0; i<psize; i++) mstd += (mBest0[i][k] - mBest1[pindex[i]][k] - mu[k])*(mBest0[i][k] - mBest1[pindex[i]][k] - mu[k]);
		for(i=0; i<psize; i++) mstd += (mBest0[i][k] - mBest1[pindex[i]][k])*(mBest0[i][k] - mBest1[pindex[i]][k]);
	}
	//mstd = sqrt(mstd/double(psize*xdim));
	mstd = mstd/double(psize*xdim);

	//Step 4: predict the center point
	Predict();	
	
	CPopulationMO pop(Population());
	//Step 5: sample new trial solutions
	for(k=0; k<xdim; k++)
	{
		double std = sqrt(mStdC[k]*mStdC[k] + mstd);
		for(i=0; i<psize; i++) 
		{
			Population()[i][k] = pC[k] + mBest0[i][k] + rnd::gaussian() * std;
			//Population()[i][k] = pC[k] + mBest0[i][k] + rnd::gaussian() * mAlpha*mStdC[k];// mAlpha*(mStdC[k] + mstd);
			//Population()[i][k] = pC[k] + 2*mBest0[i][k] - mBest1[pindex[i]][k] + rnd::gaussian() * mAlpha*mStdC[k];// mAlpha*(mStdC[k] + mstd);
			if(     Population()[i][k]>P().XUpp(k)) Population()[i][k] = 0.5*(pop[pindex[i]][k] + P().XUpp(k));
			else if(Population()[i][k]<P().XLow(k)) Population()[i][k] = 0.5*(pop[pindex[i]][k] + P().XLow(k));
		}
	}
}

void DMOO::InitFPS()
{
	//Step 0: initialize the population
	if(mStep == 0) {InitRIS(true); return;}

	unsigned int i, j, k, xdim=P().XSize(), cdim = (P().FSize()+1)*P().XSize(), psize = Population().Size();

	//Step 1: save extreme points
	mC.resize(cdim);
	double dismin, dis;
	std::vector<unsigned int> pp(P().FSize()+1); for(i=0; i<pp.size(); i++) pp[i] = 0;
	for(i=1; i<Population().Size(); i++) 
	{
		for(j=0; j<pp.size()-1; j++) if(Population()[i].F(j) < Population()[pp[j]].F(j)) pp[j] = i;
	}
	dismin = 1.0E100;
	for(i=0; i<Population().Size(); i++) 
	{
		dis = 0;
		for(j=0; j<P().FSize(); j++) dis += (Population()[i].F(j) - Population()[pp[j]].F(j))*(Population()[i].F(j) - Population()[pp[j]].F(j));
		if(dis<dismin) {pp[pp.size()-1] = i; dismin = dis;} // CTI (close-to-idea) point
	}
	for(i=0; i<pp.size(); i++) for(j=0; j<xdim; j++) mC[i*xdim+j] = Population()[pp[i]][j];
	hC.insert(hC.begin(), mC);

	//Step 2: if the history info. is not enough for prediction, then use random initialization
	if(hC.size()<2*mMaxOrder) {InitRIS(false); return;}

	//Step 3: predict the extreme points
	Predict();	

	//Step 4: generate new points
	// permutate the population
	Population().Shuffle();
	CPopulationMO pop(Population());
	// predict
	for(i=0; i<pp.size(); i++)
	{
		for(j=0; j<xdim; j++) Population()[3*i+0][j] = pC[i*xdim+j] + 1.6449*mStdC[j]*(rnd::rand()>0.5?1.0:-1.0);
		for(j=0; j<xdim; j++) Population()[3*i+1][j] = pC[i*xdim+j] + 1.6449*mStdC[j]*(rnd::rand()>0.5?1.0:-1.0);
		for(j=0; j<xdim; j++) Population()[3*i+2][j] = pC[i*xdim+j];
		for(j=0; j<xdim; j++) for(k=0; k<3; k++)
		{
			if(     Population()[3*i+k][j]>P().XUpp(j)) Population()[3*i+k][j] = 0.5*(pop[3*i+k][j] + P().XUpp(j));
			else if(Population()[3*i+k][j]<P().XLow(j)) Population()[3*i+k][j] = 0.5*(pop[3*i+k][j] + P().XLow(j));
		}
	}
	// 70% current + 30% random
	for(i=0; i<(unsigned int)(0.3*(Population().Size()-pp.size()*3)); i++)
	{
		for(j=0; j<xdim; j++) Population()[(unsigned int)(pp.size()*3+1+i)][j] = rnd::rand(P().XLow(j),P().XUpp(j));
	}
}

void DMOO::InitRIS(bool all)
{
	unsigned int i,j;
	if(all)
	{
		Population().Resize(mPopSize);
		for(i=0; i<mPopSize; i++) for(j=0; j<P().XSize(); j++) Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
	}
	else
	{
		for(i=0; i<mPopSize/2; i++) for(j=0; j<P().XSize(); j++) Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
	}
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
	else if(mOptimizer == std::string("RM2"))
	{
		az::mea::gen::mod::RM2 gen;
		gen.Set(1.0, 1.0, 0.8, 1, 10);
		gen.Generate(size, popnew, Population());
	}

	return popnew;
}

} //namespace dea
} //namespace mea
} //namespace az
