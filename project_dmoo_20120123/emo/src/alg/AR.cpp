// AR.cpp

#include <cmath>
#include <vector>
#include "alg/AR.h"
#include "alg/Matrix.h"

namespace az
{

namespace alg
{
	// input (output) matrix and vector are column-wise
	// pX: col_t(1), col_t(2), ..., col_t(size)
	bool aruv(double** pX, unsigned int dim, unsigned int size, unsigned int order, double** pA, double* pV)
	{
		if(size <= order) return false;

		unsigned int i,j,g;
		double pre;
		alg::Matrix X(order,order);
		std::vector<double> Y(order);

		for(unsigned int d=0; d<dim; d++)
		{
			// change model to Y=A*X
			for(i=0; i<order; i++)
			{
				for(j=i; j<order; j++)
				{
					X(i,j) = 0.0;
					for(g=size-2; g>=order-1; g--) X(i,j) += pX[d][g-i]*pX[d][g-j];
					X(j,i) = X(i,j);
				}

				Y[i] = 0.0;
				for(g=size-1; g>=order; g--) Y[i] += pX[d][g]*pX[d][g-i-1];
			}

			//
			X.Inv();

			// coefficient
			for(i=0; i<order; i++)
			{
				pA[d][i] = 0.0;
				for(j=0; j<order; j++) pA[d][i] += Y[j]*X(j,i);
			}

			// variance
			pV[d] = 0.0;
			for(g=size-1; g>=order; g--)
			{
				pre = 0.0;
				for(i=0; i<order; i++) pre += pA[d][i]*pX[d][g-1-i];
				pV[d] += (pX[d][g]-pre)*(pX[d][g]-pre);
			}
			pV[d] = sqrt(pV[d]/double(size-order));

			//// variance
			//double mm = 0.0; std::vector<double> err(size-order);
			//for(g=size-1; g>=order; g--)
			//{
			//	pre = 0.0;
			//	for(i=0; i<order; i++) pre += pA[d][i]*pX[d][g-1-i];
			//	err[g-order] = pX[d][g]-pre;
			//	mm += pX[d][g]-pre;
			//}
			//mm /= double(size-order);
			//pV[d] = 0.0; for(g=size-1; g>=order; g--) pV[d] += (err[g-order] - mm)*(err[g-order] - mm);
			//pV[d] = sqrt(pV[d]/double(size-order));
		}

		return true;
	}

	// input (output) matrix and vector are column-wise
	// pX: col_t(1), col_t(2), ..., col_t(size)
	bool aruv(double** pX, unsigned int dim, unsigned int size, unsigned int order, double** pA, double* pC, double* pV)
	{
		if(size <= order) return false;

		unsigned int i,j,g;
		double pre;
		alg::Matrix X(order,order);
		std::vector<double> Y(order);

		for(unsigned int d=0; d<dim; d++)
		{
			// mean
			pC[d] = 0.0;
			for(i=0; i<size; i++) pC[d] += pX[d][i]; pC[d] /= double(size);

			// change model to Y=A*X
			for(i=0; i<order; i++)
			{
				for(j=i; j<order; j++)
				{
					X(i,j) = 0.0;
					for(g=size-2; g>=order-1; g--) X(i,j) += (pX[d][g-i]-pC[d])*(pX[d][g-j]-pC[d]);
					X(j,i) = X(i,j);
				}

				Y[i] = 0.0;
				for(g=size-1; g>=order; g--) Y[i] += (pX[d][g]-pC[d])*(pX[d][g-i-1]-pC[d]);
			}

			//
			X.Inv();

			// coefficient
			for(i=0; i<order; i++)
			{
				pA[d][i] = 0.0;
				for(j=0; j<order; j++) pA[d][i] += Y[j]*X(j,i);
			}

			// variance
			pV[d] = 0.0;
			for(g=size-1; g>=order; g--)
			{
				pre = pC[d];
				for(i=0; i<order; i++) pre += pA[d][i]*(pX[d][g-1-i]-pC[d]);
				pV[d] += (pX[d][g]-pre)*(pX[d][g]-pre);
			}
			pV[d] = sqrt(pV[d]/double(size-order));
		}

		return true;
	}

	// input (output) matrix and vector are column-wise
	// pX: col_t(1), col_t(2), ..., col_t(size)
	bool armv(double** pX, unsigned int dim, unsigned int size, unsigned int order, double** pA, double* pC, double* pV)
	{
	    unsigned int d,o,s,ir,ic,ord,dr,oc,dc;
		double pre;
		if(size <= order) return false;

		// mean
		for(d=0; d<dim; d++)
		{
			pC[d]  = 0.0;
			for(s=0; s<size; s++) pC[d] += pX[d][s];
			pC[d] /= double(size);
		}

		// change model to Y = AX
		alg::Matrix X(order*dim, order*dim);
		alg::Matrix Y(      dim, order*dim);

		for(ord=0; ord<order; ord++)
		{
			for(dr=0; dr<dim; dr++)
			{
				ir = ord*dim+dr;
				for(oc=0; oc<order; oc++)
					for(dc=0; dc<dim; dc++)
					{
						ic = oc*dim+dc;
						X(ir,ic) = 0.0;
						for(s=size-2; s>=order-1; s--) X(ir,ic) += (pX[dr][s-ord]-pC[dr])*(pX[dc][s-oc]-pC[dc]);
					}
			}
		}
		for(dr=0; dr<dim; dr++)
		{
			ir = dr;
			for(oc=0; oc<order; oc++)
				for(dc=0; dc<dim; dc++)
				{
					ic = oc*dim+dc;
					Y(ir,ic) = 0.0;
					for(s=size-2; s>=order-1; s--) Y(ir,ic) += (pX[dr][s+1]-pC[dr])*(pX[dc][s-oc]-pC[dc]);
				}
		}

		X.Inv();

		for(d=0; d<dim; d++)
		{
			for(s=0; s<dim*order; s++)
			{
				pA[d][s] = 0.0;
				for(o=0; o<dim*order; o++) pA[d][s] += Y(d,o)*X(o,s);
			}
		}

		//
		for(d=0; d<dim; d++) pV[d] = 0.0;
		for(s=size-1; s>=order; s--)
		{
			for(d=0; d<dim; d++)
			{
				pre = pC[d];
				for(ord=0; ord<order; ord++)
					for(dr=0; dr<dim; dr++)
						pre += pA[d][ord*dim+dr]*(pX[dr][s-ord-1]-pC[dr]);
				pV[d] += (pre-pX[d][s])*(pre-pX[d][s]);
			}
		}
		for(d=0; d<dim; d++) pV[d] = sqrt(pV[d]/double(size-order));

		return true;
	}

} //namespace alg

} //namespace az
