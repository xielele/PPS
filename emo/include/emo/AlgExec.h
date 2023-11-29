/*! \file	AlgExec.h
	
	\brief	implements of multi-objective evoluationary algorithms
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Dec.02 2006 create
	\date	May.22 2007 redesign and rewrite
	\date	Jun.30 2008 redesign and rewrite
*/

#ifndef	AZ_ALG_EXEC_H
#define	AZ_ALG_EXEC_H

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

// include files
#include "alg/Random.h"
#include "Problem.h"

#include "emo/Parameter.h"
#include "emo/Config.h"
#include "emo/Alg.h"
#include "emo/Initialization.h"
#include "emo/IndividualMO.h"
#include "emo/PopulationMO.h"

// define the basic class
class EAExec
{
public:
	az::mea::CParameter	mPar;		//parameter object
	unsigned int		mSta,		//flag to indicate weather temporary results
						mGen,		//current generation
						mMaxGen,	//max generation
						*pIndex,	// index of storage space
						*pNo;		// number of nondominated solutions in each generation
	double				mTa,		//run time
						**pPF,		// Pareto Front
						*pT;		// run time
	az::mea::EABasic	*pEA;		// moea object
public:
	EAExec()
	{
		pEA		= 0;
	}

	virtual ~EAExec()
	{
		Clear();
	}

	void Clear()
	{
		if( pEA != 0) {delete pEA; pEA = 0;}
	}

	void SetArray( double* pt, double** ppf, unsigned int* pindex, unsigned int* pno )
	{
		pIndex= pindex;
		pNo   = pno;
		pPF   = ppf;
		pT    = pt;
	}

	void Save(std::string& path)
	{
		pEA->Save(path);
	}
	
	int Step( int s )
	{
		pEA->Step();

		mGen++;

		Record();

		if(pEA->IsTerminate())
		{
   			// record run time
			if(mSta > 0)
			{
				double tm,ts;
				pEA->Time(tm,ts);
				mTa = tm+ts;            
				*pT = mTa;
			}
			return -1;
		}
		else
			return 1;
	}

	// initialize algorithm parameters
	int Init( az::mea::Config& config )
	{
		std::string methodstr;
		double bound;
		unsigned int i,dimension,popsize,xcoding;
		std::string problem;
		
		// random seed
		az::rnd::seed((long) time(NULL));

		config.Get(std::string("COMMONP"),std::string("TOLERANCEF"),0,mPar.TolF());
		config.Get(std::string("COMMONP"),std::string("TOLERANCEX"),0,mPar.TolX());
		config.Get(std::string("COMMONP"),std::string("TOLERANCEC"),0,mPar.TolC());
		config.Get(std::string("SBX"),std::string("PC"),0,mPar.Pc());
		config.Get(std::string("SBX"),std::string("PM"),0,mPar.Pm());
		config.Get(std::string("SBX"),std::string("ETAPC"),0,mPar.ETA_SBX());
		config.Get(std::string("SBX"),std::string("ETAPM"),0,mPar.ETA_PM());
		config.Get(std::string("PCX"),std::string("ETA"),0,mPar.ETA_PCX());

		config.Get(std::string("COMMONA"),std::string("STA"),0,mSta);
		config.Get(std::string("COMMONA"),std::string("XCODING"),0,xcoding);
		config.Get(std::string("COMMONA"),std::string("INI"),0,methodstr);
		
		config.Get(std::string("COMMONA"),std::string("PROBLEM"),0,problem);
		PROBLEM::BuildProblems();
		mPar.FSize( PROBLEM::FSize(problem.c_str()) ); 
		mPar.ESize( PROBLEM::ESize(problem.c_str()) ); 
		mPar.ISize( PROBLEM::ISize(problem.c_str()) ); 
		mPar.Evaluator( PROBLEM::Evaluator(problem.c_str()) );

		config.Get(std::string("COMMONP"),std::string("DIMENSION"),0,dimension);
		mPar.XSize(dimension);
		for( i=0; i<config.GetSize(std::string("COMMONP"),std::string("BOUNDUPPER")); i++ )
		{
			config.Get(std::string("COMMONP"),std::string("BOUNDUPPER"),i,bound);
			mPar.XRealUpp(i) = bound;
		}
		for( i=config.GetSize(std::string("COMMONP"),std::string("BOUNDUPPER")); i<mPar.XSize(); i++ )
			mPar.XRealUpp(i) = mPar.XRealUpp(i-1);
		for( i=0; i<config.GetSize(std::string("COMMONP"),std::string("BOUNDLOWER")); i++ )
		{
			config.Get(std::string("COMMONP"),std::string("BOUNDLOWER"),i,bound);
			mPar.XRealLow(i) = bound;
		}
		for( i=config.GetSize(std::string("COMMONP"),std::string("BOUNDLOWER")); i<mPar.XSize(); i++ )
			mPar.XRealLow(i) = mPar.XRealLow(i-1);
		if(xcoding>0)
		{
			mPar.XCoding() = true;
			for(i=0; i<mPar.XSize(); i++) { mPar.XLow(i) = 0.0; mPar.XUpp(i) = 1.0; }
		}
		else
		{
			mPar.XCoding() = false;
			for(i=0; i<mPar.XSize(); i++) { mPar.XLow(i) = mPar.XRealLow(i); mPar.XUpp(i) = mPar.XRealUpp(i); }
		}

		config.Get(std::string("COMMONA"),std::string("POPSIZE"),0,popsize);
		config.Get(std::string("COMMONA"),std::string("ITERATIONS"),0,mMaxGen);

		// create initial population
		az::mea::CPopulationMO pop(mPar);
		unsigned int inieva;
		if(methodstr==std::string("LHC"))
		{
			az::mea::ini::LHC init;
			init.Initialize(pop, popsize);
			inieva = init.EvaTimes();
		}
		else if(methodstr==std::string("HYBRID"))
		{
			az::mea::ini::Hybrid init;
			init.Initialize(pop, popsize);
			inieva = init.EvaTimes();
		}
		else
		{
			az::mea::ini::Uniform init;
			init.Initialize(pop, popsize);
			inieva = init.EvaTimes();
		}

		/////////////////////////////////////	
		//create MA object
		Clear();
		mGen	= 0;
		pEA		= NewEA(mMaxGen-1, mPar, pop);
		pEA->EvaTimes() += inieva;

		SetGenerator(config);

		//record
		Record();

		return 1;
	}

protected:
	// set generator parameters
	virtual void SetGenerator(az::mea::Config& config) = 0;

	virtual az::mea::EABasic* NewEA( unsigned int	stepmax, az::mea::CParameter& par, az::mea::CPopulationMO& pop) = 0;

	// record population 
	void Record()
	{
		unsigned int nondom = 0, dom = 0, k = 0, s = 0, redim = pEA->P().XSize() < 4 ? pEA->P().XSize() : 4;
		// staticsic run, save data into an array		
		if(mSta>0)
		{
			//// save every 10 generations
			//if(mGen % 10 != 0) 
			//{
			//	pNo[mGen] = 0;
			//	return;
			//}
			// sort the population: nondominate solutions -> dominate solutions
			pEA->Population().RankSort();
			pEA->Population().RankSize(1,dom, nondom);
			dom = pEA->Population().Size() - nondom;

			for(k=0; k<nondom; k++)
			{
				for(s=0; s<pEA->Population().P().FSize(); s++) pPF[(*pIndex)][s] = pEA->Population()[k].F(s);
				for(s=0; s<redim; s++) pPF[(*pIndex)][pEA->Population().P().FSize()+s] = pEA->Population()[k][s];
				(*pIndex)++;
			}
			pNo[mGen]   = nondom;
		}
	}
};

#endif //AZ_ALG_EXEC_H
