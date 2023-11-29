#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include "emo/Config.h"

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
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}

    return TRUE;
}
#endif

#ifdef	WIN32
extern "C" __declspec(dllexport) int Problem( )
#elif	_LINUX
extern "C"{int  Problem( )
#endif
{
	std::string str;
	az::mea::Config mConfig;
	mConfig.Load("dmoo.set", '#');
	mConfig.Get("COMMONA", "PROBLEM", 0, str);
	if(str==std::string("FDA1"))		return 1;
	else if(str==std::string("FDA2"))	return 2;
	else if(str==std::string("FDA3"))	return 3;
	else if(str==std::string("FDA4"))	return 4;
	else if(str==std::string("DMOP1"))	return 5;
	else if(str==std::string("DMOP2"))	return 6;
	else if(str==std::string("DMOP3"))	return 7;
	else if(str==std::string("DMOPA"))	return 8;
	else if(str==std::string("DMOPB"))	return 9;
	else if(str==std::string("DMOPC"))	return 10;
	else if(str==std::string("DMOPD"))	return 11;
	else if(str==std::string("DMOPE"))	return 12;
	else if(str==std::string("DMOPF"))	return 13;
	else return 0;
}
#ifdef _LINUX 
}
#endif

#ifdef	WIN32
extern "C" __declspec(dllexport) int TInterval( )
#elif	_LINUX
extern "C"{ int  TInterval( )
#endif
{
	int inte;
	az::mea::Config mConfig;
	mConfig.Load("dmoo.set", '#');
	mConfig.Get("DGTM", "TAOT", 0, inte);
	return inte;
}
#ifdef _LINUX 
}
#endif

#ifdef	WIN32
extern "C" __declspec(dllexport) double TInit( )
#elif	_LINUX
extern "C"{ double TInit( )
#endif
{
	double tini;
	az::mea::Config mConfig;
	mConfig.Load("dmoo.set", '#');
	mConfig.Get("DGTM", "TINI", 0, tini);
	return tini;
}
#ifdef _LINUX 
}
#endif

#ifdef	WIN32
extern "C" __declspec(dllexport) double TStep( )
#elif	_LINUX
extern "C"{ double TStep( )
#endif
{
	unsigned int nt;
	az::mea::Config mConfig;
	mConfig.Load("dmoo.set", '#');
	mConfig.Get("DGTM", "NT", 0, nt);
	return 1.0/nt;
}
#ifdef _LINUX 
}
#endif
