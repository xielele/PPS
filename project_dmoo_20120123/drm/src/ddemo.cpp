#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "emo/Parameter.h"
#include "emo/Config.h"
#include "emo/Sel.h"
#include "AlgD.h"

az::mea::CParameter mPar;
az::mea::dea::DMOO* pEA;

#ifdef WIN32
BOOL APIENTRY DllMain( HANDLE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch(ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
		break;
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
    return TRUE;
}
#endif

#ifdef	WIN32
extern "C" __declspec(dllexport) int  Init( const char *parfile )
#elif	_LINUX
extern "C"{ int Init( const char *parfile )
#endif
{
		double bound;
		unsigned int i,dimension,popsize,xcoding, mSta, mMaxGen, strategy;
		std::string problem, method;
		
		// random seed
		az::rnd::seed((long) time(NULL));

		az::mea::Config config;
		config.Load(std::string(parfile), '#');

		config.Get(std::string("COMMONP"),std::string("TOLERANCEF"),0, mPar.TolF());
		config.Get(std::string("COMMONP"),std::string("TOLERANCEX"),0, mPar.TolX());
		config.Get(std::string("COMMONP"),std::string("TOLERANCEC"),0, mPar.TolC());
		config.Get(std::string("COMMONA"),std::string("STA"),		0, mSta);
		config.Get(std::string("COMMONA"),std::string("XCODING"),	0, xcoding);
		config.Get(std::string("COMMONA"),std::string("STRATEGY"),	0, strategy);
		config.Get(std::string("COMMONA"), std::string("METHOD"),	0, method);		
		config.Get(std::string("COMMONA"),std::string("PROBLEM"),	0, mPar.Problem());
		config.Get(std::string("COMMONP"),std::string("DIMENSION"),	0, dimension);
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

		unsigned int torder, taot, nt; double t0, alpha;
		config.Get(std::string("DGTM"),std::string("TAOT"),			0, taot);
		config.Get(std::string("DGTM"),std::string("TORDER"),		0, torder);
		config.Get(std::string("DGTM"),std::string("TINI"),			0, t0);
		config.Get(std::string("DGTM"),std::string("NT"),			0, nt);
		config.Get(std::string("DGTM"),	  std::string("ALPHA"),		0, alpha);

		pEA = new az::mea::dea::DMOO(strategy, method, popsize, mMaxGen, taot, nt, torder, t0, alpha, mPar); 

		return 0;
}
#ifdef _LINUX
}
#endif

#ifdef WIN32
extern "C" __declspec(dllexport) int Step( int s )
#elif  _LINUX
extern "C"{ int Step( int s )
#endif
{
	if(pEA == 0)
	{
		std::cout<<"Error"<<std::endl;
		return -1;
	}
	else
	{
		int ss = (*pEA).Step();
		std::string str = std::string("dmoo.pop");
		pEA->Write(str);
		if(pEA->IsTerminate()) {delete pEA; pEA=0; return -1;}
		return ss;
	}
}
#ifdef _LINUX
}
#endif
