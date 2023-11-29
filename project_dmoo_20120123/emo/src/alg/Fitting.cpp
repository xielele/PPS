// Fitting.cpp

#include <cmath>
#include "alg/Matrix.h"
#include "alg/Fitting.h"

namespace az
{

namespace alg
{
	//C.T=X
	//C=(c0,c1,c2,...,cn), X=(x1, x2, ..., xm)
	void poly_fit(double* pT, double* pX, unsigned int size, unsigned int order, double* pC)
	{
		unsigned int i,j;
		alg::Matrix T(order+1, size);
		for(i=0; i<size; i++)
			for(j=0; j<=order; j++)
				T(j, i) = pow(pT[i], double(j));
		alg::Matrix TT(T);
		TT.Trans();
		alg::Matrix TTT;
		T.Multiply(TT, TTT);
		TTT.Inv();
		alg::Matrix tmp;
		TT.Multiply(TTT , tmp);
		for(i=0; i<=order; i++)
		{
			pC[i] = 0.0;
			for(j=0; j<size; j++) pC[i] += pX[j] * tmp(j, i);
		}
	}

	// C.T=X
	// C=(c0,c1,c2,...,cn), X=(x1, x2, ..., xm)
	void poly_fit(std::vector<double>& T, std::vector<double>& X, unsigned int order, std::vector<double>& C)
	{
		poly_fit(&(*(T.begin())), &(*(X.begin())), (unsigned int)(T.size()), order, &(*(C.begin())));
	}

	// C.T=X
	// C=(c0,c1,c2,...,cn), X=(x1, x2, ..., xm)
	void poly_fit(double* pT1, double* pT2, double* pX, unsigned int size, unsigned int order, double* pC)
	{
		unsigned int i,j,k,row=0, p1, p2, r;
		for(i=1; i<=order+1; i++) row += i;
		alg::Matrix T(row, size);

		r=0;
		for(i=0; i<=order; i++)
		{
			p1=i; p2=0;
			for(j=0; j<i+1; j++)
			{
				for(k=0; k<size; k++) T(r,k) = pow(pT1[k], double(p1))*pow(pT2[k], double(p2));
				p1--; p2++;r++;
			}
		}

		alg::Matrix TT(T);
		TT.Trans();
		alg::Matrix TTT;
		T.Multiply(TT, TTT);
		TTT.Inv();
		alg::Matrix tmp;
		TT.Multiply(TTT , tmp);
		for(i=0; i<row; i++)
		{
			pC[i] = 0.0;
			for(j=0; j<size; j++) pC[i] += pX[j] * tmp(j, i);
		}
	}

	// C.T=X
	// C=(c0,c1,c2,...,cn), X=(x1, x2, ..., xm)
	void poly_fit(std::vector< double >& T1, std::vector< double >& T2, std::vector< double >& X, unsigned int order, std::vector< double >& C)
	{
		poly_fit(&(*(T1.begin())), &(*(T2.begin())), &(*(X.begin())), (unsigned int)(T1.size()), order, &(*(C.begin())));
	}

} //namespace alg

} //namespace az
