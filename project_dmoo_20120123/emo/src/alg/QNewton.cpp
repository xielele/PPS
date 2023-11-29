// QNewton.cpp

#include <math.h>
#include <algorithm>
#include "alg/QNewton.h"

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

//constractor
QNewton::QNewton(					
	int		dim,
	FTYPE	fun,
	GTYPE	gun,
	double	*x0,
	double	*low,
	double	*up,
	double	tolx,
	double	tolg)
{
	int i,j; 
	double sum;
	mDim	= dim;
	pFun	= fun;
	pGun	= gun;
	mTolG	= tolg;
	mTolX	= tolx;
	mStepMax	= 100;
	mIter		= 0;
	mEval		= 0;
	pBoundUpp	= new double[dim];
	pBoundLow	= new double[dim];
	pG			= new double[dim];
	pDG			= new double[dim];
	pHDG		= new double[dim];
	pP			= new double[dim];
	pScild		= new double[dim];
	pX			= new double[dim];
	pXBest		= new double[dim];
	pHessin		= new double*[dim];

	for(i=0; i<dim; i++)
	{
		pBoundUpp[i]	= up[i];
		pBoundLow[i]	= low[i];
		pX[i]			= pXBest[i] 
						= pScild[i]
						= x0[i];
	}

	mFBest = mF = FEva(pScild);
	GEva(pScild, pG);

	sum = 0.0;
	for(i=0; i<dim; i++)
	{
		pHessin[i] = new double[dim];
		for(j=0; j<dim; j++) pHessin[i][j] = 0.0;
		pHessin[i][i] = 1.0;
		pP[i]	= -pG[i];
		sum += pScild[i]*pScild[i];
	}
	sum		= sqrt(sum);
	mStepMax= std::max<double>(sum,dim);
}

//destractor
QNewton::~QNewton()
{
	int i;
	delete []pBoundUpp;
	delete []pBoundLow;
	delete []pG;
	delete []pDG;
	delete []pHDG;
	delete []pP;
	delete []pScild;
	delete []pX;
	delete []pXBest;
	for(i=0; i<mDim; i++) delete []pHessin[i];
	delete []pHessin;
}

//one step
bool QNewton::Step()
{
	int i,j;
	double test, temp, den;

	mIter ++;
	
	lnsrch(pX, mF, pG, pP);

	for(i=0; i<mDim; i++) 
	{
		pP[i]		= pX[i] - pScild[i];
		pScild[i]	= pX[i];
	}

	test = 0.0;
	for(i=0; i<mDim; i++)
	{
		temp=fabs(pP[i]) / std::max<double>(fabs(pScild[i]), 1.0);
		if(temp>test) test = temp;
	}

	if(test < mTolX) 
	{
		return false;
	}

	for(i=0; i<mDim; i++) pDG[i] = pG[i];

	GEva(pScild, pG);
	
	test=0.0;
	den	= std::max(mF, 1.0);
	for(i=0; i<mDim; i++) 
	{
		temp=fabs(pG[i])*std::max<double>(fabs(pScild[i]), 1.0)/den;
		if(temp > test) test=temp;
	}
	if (test < mTolG) 
	{
		return false;
	}

	for(i=0; i<mDim; i++) pDG[i] = pG[i] - pDG[i];
	for(i=0; i<mDim; i++) 
	{
		pHDG[i]=0.0;
		for(j=0; j<mDim; j++) pHDG[i] += pHessin[i][j]*pDG[j];
	}

	double fac = 0.0, fae = 0.0, sumdg = 0.0, sump = 0.0, fad;
	for( i=0; i<mDim; i++ ) 
	{
		fac		+= pDG[i]*pP[i];
		fae		+= pDG[i]*pHDG[i];
		sumdg	+= pDG[i]*pDG[i];
		sump	+= pP[i]*pP[i];
	}
	if(fac*fac > 3.0e-8*sumdg*sump) 
	{
		fac=1.0/fac;
		fad=1.0/fae;
		for( i=0; i<mDim; i++) pDG[i] = fac*pP[i] - fad*pHDG[i];
		for( i=0; i<mDim; i++) 
			for( j=0; j<mDim; j++ ) 
				pHessin[i][j] += fac*pP[i]*pP[j] - fad*pHDG[i]*pHDG[j] + fae*pDG[i]*pDG[j];
	}
	for(i=0; i<mDim; i++) 
	{
		pP[i]=0.0;
		for(j=0; j<mDim; j++) pP[i] -= pHessin[i][j]*pG[j];
	}

	return true;
}

//evaluate a vector
double QNewton::FEva(double* px)
{
	mEval++;
	double f = (*pFun)(px);
	if(f < mFBest) 
	{ 
		mFBest = f; 
		for(int i=0; i<mDim; i++) pXBest[i] = px[i]; 
	}
	return f;
}

//calculate the gradient of a vector
void QNewton::GEva(double* px, double* pg)
{
	(*pGun)(px,pg);
}

//linear search
bool QNewton::lnsrch(
	double*& x, 
	double & f, 
	double*& g, 
	double*& p)
{
	int i;
	double a,alam,alam2=0,alamin,b,disc,f2=0,fold2=0,rhs1,rhs2,slope,sum,temp,test,tmplam,fold;
	
	double* scild = new double[mDim];
	
	for(i=0; i<mDim; i++) scild[i] = x[i];
	fold = f;

	for(sum=0.0,i=0; i<mDim; i++) sum += p[i]*p[i];
	sum=sqrt( sum );
	if(sum > mStepMax) for (i=0; i<mDim; i++) p[i] *= mStepMax/sum;
	
	for (slope=0.0,i=0; i<mDim; i++) slope += g[i]*p[i];
	test=0.0;
	for( i=0 ;i<mDim; i++ ) 
	{
		temp=fabs(p[i])/std::max<double>(fabs(scild[i]), 1.0);
		if(temp > test) test=temp;
	}
	alamin = mTolX/test;
	alam=1.0;
	for (;;) 
	{
		//assign new value
		for(i=0; i<mDim; i++)
		{
			x[i] = scild[i] + alam*p[i];
			if(x[i] < pBoundLow[i])			x[i] = pBoundLow[i] + mTolX;
			else if( x[i] > pBoundUpp[i] )	x[i] = pBoundUpp[i] - mTolX;
		}
		f = FEva(x);

		// unsuccess
		if(alam < alamin) 
		{
			for(i=0; i<mDim; i++) x[i] = scild[i];
			delete []scild;
			return false;
		} 
		// success
		else if(f <= fold+1.0e-3*alam*slope) 
		{
			delete []scild;
			return true;
		}
		else 
		{
			if(alam == 1.0) tmplam = -slope/( 2.0*( f-fold-slope ) );
			else 
			{
				rhs1	= f-fold-alam*slope;
				rhs2	= f2-fold2-alam2*slope;
				a		= (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b		= (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else
				{
					disc=b*b-3.0*a*slope;
					if (disc<0.0)	
					{
						delete []scild;
						return false;//throw error( "Roundoff problem in lnsrch." );
					}
					else tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam) tmplam=0.5*alam;
			}
		}
		alam2	= alam;
		f2		= f;
		fold2	= fold;
		alam	= std::max(tmplam, 0.1*alam);
	}
	delete []scild;
}

} //namespace alg

} //namespace az

