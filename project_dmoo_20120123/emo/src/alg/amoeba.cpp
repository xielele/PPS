// amoeba.cpp

#include <algorithm>
#include <cmath>
#include "alg/amoeba.h"

namespace az
{
namespace alg
{
	// get the row sum
	void amoeba::get_psum(MAT &p, VEC &psum)
	{
		unsigned int i,j;
		for (j=0;j<p[0].size();j++) {
			for (psum[j]=0.0,i=0;i<p.size();i++)
				psum[j] += p[i][j];
		}
	}
	
	// main function of simplex method
	void amoeba::opt(MAT &p, VEC &y, const double ftol, unsigned int &nfunk, const unsigned int NMAX)
	{
		const double TINY=1.0e-10;
		unsigned int i,j,
					ihi,					// index of vertice with highest value
					ilo,					// index of vertice with lowest value
					inhi,					// index of vertice with the 2nd highest value
					mpts = (unsigned int)p.size(),			// number of vertices
					ndim = (unsigned int)p[0].size();		// dimension
		double rtol,ysave,ytry;

		VEC psum(ndim);	get_psum(p,psum);	// sum of all vertices
		nfunk = 0;
		for (;;) 
		{
			// find the lowest, 1st and 2nd highest values
			ilo = 0;
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			for(i=0;i<mpts;i++) 
			{
				if (y[i] <= y[ilo])  ilo = i;
				if (y[i] >  y[ihi]) {inhi= ihi; ihi=i;} 
				else if (y[i] > y[inhi] && i != ihi) inhi = i;
			}

			// stop conditions
			rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
			if(rtol < ftol || nfunk >= NMAX) 
			{
				std::swap(y[0],y[ilo]);
				std::swap(p[0],p[ilo]);
				break;					
			}
			
			// try simplex
			nfunk += 2;
			// reflection
			ytry = amotry(p,y,psum,ihi,-1.0);
			// reflection and expansion
			if(ytry <= y[ilo]) ytry = amotry(p,y,psum,ihi,2.0);
			// contraction
			else if(ytry >= y[inhi]) 
			{
				ysave = y[ihi];
				ytry  = amotry(p,y,psum,ihi,0.5);
				// multiple contraction
				if(ytry >= ysave) 
				{
					for(i=0;i<mpts;i++) 
					{
						if(i != ilo) 
						{
							for(j=0;j<ndim;j++)
								p[i][j] = psum[j] = 0.5*(p[i][j]+p[ilo][j]);
							y[i] = objective(psum);
						}
					}
					nfunk += ndim;
					get_psum(p,psum);
				}
			} 
			else --nfunk;
		}	
	}

	// one try of simplex method
	double amoeba::amotry(MAT &p, VEC &y, VEC &psum, const unsigned int ihi, const double fac)
	{
		unsigned int j,ndim = (unsigned int)p[0].size();
		double fac1,fac2,ytry;
		VEC ptry(ndim);

		fac1 = (1.0-fac)/ndim;
		fac2 = fac1-fac;
		for(j=0;j<ndim;j++) ptry[j] = psum[j]*fac1-p[ihi][j]*fac2;
		ytry = objective(ptry);
		if(ytry < y[ihi]) 
		{
			y[ihi] = ytry;
			for(j=0;j<ndim;j++) 
			{
				psum[j]  += ptry[j]-p[ihi][j];
				p[ihi][j] = ptry[j];
			}
		}
		return ytry;
	}
} //namespace alg
} //namespace az
