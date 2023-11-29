// =====================================================================
//
// Purpose:			EDA with local search
//
// Author:			A. Zhou
//
// Last Modified:   2008/03/07
//
// Description:     implementation of EDA/L based on J. Sun's work
//
// =====================================================================
#ifndef AZ_EDAL_H
#define AZ_EDAL_H

namespace edal
{
	typedef double (*PObj)(double*, int);	// objective function pointer

	// population initialization
	void UDM_INITIALIZE(PObj calobj, int method, double a, double b, int dimension, int Popsize, double**x, double* f);

	// select population to build model
	void selection(double** population, double* popcost, int dimension, int Popsize, double** selectedpop, int selsize, int* selected);

	// umda offspring generator
	int umda(PObj calobj, double** y, double*, double& cost, int M, int L, int K, double sigma, double a, double b, int dimension);
	int umda(PObj calobj, double** y, double* ycost, double*, double& cost, int M, int L, int K, double sigma, double a, double b, int dimension);

	// simplex local search
	int	amoeba(PObj calobj, double* start, double lamda, int fcall, double& cost, int L, double a, double b, int dimension);

	// EDA/L optimizer
	int optimizer(PObj CalObj, int MaxEva, int MaxLS, int PopSize, int Dimension, double BoundLow, double BoundUpp, double* BestX);
}

#endif
