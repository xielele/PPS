#ifndef AZ_DMOO_DLL_H
#define AZ_DMOO_DLL_H

extern "C" __declspec(dllexport) int Problem();
extern "C" __declspec(dllexport) int TInterval();
extern "C" __declspec(dllexport) double TInit();
extern "C" __declspec(dllexport) double TStep();

#endif	//AZ_DMOO_DLL_H
