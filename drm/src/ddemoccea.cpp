#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "problem.h"
#include "COEA.h"
#include "emo/Config.h"

COEA *pEA;

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
	int dimension, popsize, maxgen, tt, nt, t0;
	std::string name;

	az::mea::Config config;
	config.Load(std::string(parfile), '#');
	
	config.Get(std::string("COMMONA"),std::string("PROBLEM"),	0, name);
	config.Get(std::string("COMMONP"),std::string("DIMENSION"),	0, dimension);
	config.Get(std::string("COMMONA"),std::string("POPSIZE"),	0, popsize);
	config.Get(std::string("COMMONA"),std::string("ITERATIONS"),0, maxgen);
	config.Get(std::string("DGTM"),std::string("TAOT"),			0, tt);
	config.Get(std::string("DGTM"),std::string("TINI"),			0, t0);
	config.Get(std::string("DGTM"),std::string("NT"),			0, nt);

	IND::evalFp = OBJ_FUN;

	RND.seed((long) time(NULL));

	pEA = new COEA(tt, nt, dimension, maxgen, popsize, t0, name);

	pEA->Init();

	std::string str = std::string("dmoo.pop");
	pEA->Write(str);

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