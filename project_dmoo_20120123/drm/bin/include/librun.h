#ifndef AZ_MOO_DLL_H
#define AZ_MOO_DLL_H

extern "C" __declspec(dllexport) int Init( const char *problem );
extern "C" __declspec(dllexport) int Step( int s );

#endif	//AZ_MOO_DLL_H
