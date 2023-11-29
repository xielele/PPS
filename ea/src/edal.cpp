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

#include <cmath>
#include <ctime>
#include <iostream>
#include "edal.h"
#include "Random.h"

namespace edal
{
#define n		100	// maximal problem dimension
#define tol		1e-10		// used for bfgs

//-----------------------------------------------------------------------
int swap(double& tempx, double& tempy)
{
	double auxf;
	
	auxf = tempx;
	tempx = tempy;
	tempy = auxf;

	return 0;

}
void Quicksort(double *array, int length)
{
	const int M1=length, NSTACK = 50;

	int i,ir,j,k,jstack=-1,l=0;

	double astar;

	int *istack;

	istack = new int[NSTACK];

	ir = length-1;

	for (;;){
		if (ir-l < M1)
		{
			for (j=l+1; j<=ir; j++)
			{
				astar = array[j];
				for (i=j-1; i>=l; i--)
				{
					if (array[i] <= astar) break;
					array[i+1] = array[i];
				}
				array[i+1] = astar;
			}
			if (jstack < 0) break;
			ir = istack[jstack--];
			l  = istack[jstack--];
		}
		else
		{
			k = (l+ir) >> 1;
			swap(array[k],array[l+1]);
			if (array[l] > array[ir])
			{
				swap(array[l],array[ir]);
			}
			if (array[l+1] > array[ir])
			{
				swap(array[l+1],array[ir]);
			}
			if (array[l] > array[l+1])
			{
				swap(array[l],array[l+1]);
			}
			i=l+1;
			j=ir;
			astar = array[l+1];
			for (;;)
			{
				do i++; while(array[i] < astar);
				do j--; while(array[j] > astar);
				if (j < i) break;
				swap(array[i],array[j]);
			}
			array[l+1] = array[j];
			array[j] = astar;
			jstack += 2;
			if (jstack >= NSTACK ) std::cout<<"NSTACK too small in sort"<<std::endl;
			if (ir-i+1 >= j-1)
			{
				istack[jstack] =ir;
				istack[jstack-1] = i;
				ir = j-1;
			}
			else
			{
				istack[jstack] = j-1;
				istack[jstack-1]= l;
				l=i;
			}
		}
	}

	delete []istack;
}
int intswap(int& tempx, int& tempy)
{
	int auxf;
	     
	auxf = tempx;
	tempx = tempy;
	tempy = auxf;
	       
	return 0;      
}
void indexx(double *arr, int *indx, int length)
// Indexs an array arr[0,..length-1], i.e., outputs the array indx[0..length-1]
// such that arr[indx[j]] is in ascending order for j=0,1,...length-1.
// The input array arr[length] is not changed
{
   const int M1=length, NSTACK = 50;
   int i, indxt,ir,j,k,jstack=-1,l=0;
   double astar;
   int *istack;

   istack = new int[NSTACK];

   ir = length -1;

   for (j=0; j<length; j++) indx[j] =j;
   for (;;)
     {
         if (ir -l < M1)
         {
            for (j=l+1; j<=ir; j++)
            {
               indxt = indx[j];

               astar = arr[indxt];

               for (i=j-1; i>=l; i--)
               {
                  if (arr[indx[i]] <= astar) break;
                  indx[i+1] = indx[i];
               }
               indx[i+1]=indxt;
            }
            if (jstack < 0) break;
            ir = istack[jstack--];
            l  = istack[jstack--];
         }
         else
         {
            k = (l+ir)>>1;
            intswap(indx[k],indx[l+1]);
            if (arr[indx[l]] > arr[indx[ir]])
            {
               intswap(indx[l],indx[ir]);
            }
            if (arr[indx[l+1]] > arr[indx[ir]])
            {
               intswap(indx[l+1],indx[ir]);
            }
            if (arr[indx[l]] > arr[indx[l+1]])
            {
               intswap(indx[l],indx[l+1]);
            }  
            
            i = l+1;
            j = ir;
            indxt = indx[l+1];
            astar = arr[indxt];
            for (;;)
            {
               do i++; while (arr[indx[i]] < astar);
               do j--; while (arr[indx[j]] > astar);
               if (j<i) break;
               intswap(indx[i],indx[j]); 
            }
            indx[l+1] = indx[j];
            indx[j]   = indxt;
               
            jstack += 2;
			if (jstack >= NSTACK)std::cout<<"NSTACK too small in indexx"<<std::endl;
            if (ir-i+1 >= j-1)
            {
               istack[jstack] = ir;   
               istack[jstack-1] = i;
               ir = j-1;
            }
            else
            {
               istack[jstack] = j-1;
               istack[jstack-1] = l;
               l = i;
            }
        }
    }
    delete []istack;
}

int SWAP(double& tempx, double& tempy)
{
	double auxf;
	
	auxf = tempx;
	tempx = tempy;
	tempy = auxf;

	return 0;

}
double amotry(double p[n+1][n], double yp[n+1], double psum[n], PObj calobj, const int ihi, const double fac, int dimension, double a, double b)
{
	int j,ndim;
	double fac1, fac2, ytry;
	double ptry[n];

	ndim = dimension;

	fac1 = (1.0-fac)/ndim;
	fac2 = fac1 - fac;

	for (j=0; j<ndim; j++)
	{
		ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;

		if (ptry[j] <= a) ptry[j] = a;

		if (ptry[j] >= b) ptry[j] = b;
	}

	ytry = (*calobj)(ptry, dimension);

	if (ytry < yp[ihi])
	{
		yp[ihi] = ytry;
		for (j=0; j<ndim; j++)
		{
			psum[j] += ptry[j] - p[ihi][j];
			p[ihi][j] = ptry[j];
		}
	}

	return ytry;
}
//-----------------------------------------------------------------------

// population initialization
void UDM_INITIALIZE(PObj calobj, int method, double a, double b, int dimension, int Popsize, double**x, double* f)
{
	int l,k;
	double r1;
	
	for (l=0; l<Popsize; l++)
	{
		for (k=0; k<dimension; k++)
		{
			do r1 = az::rnd::rand(); while(r1 < 1.0e-30);
			x[l][k] = r1*(b-a) + a;
		}
	}

	for (l=0; l<Popsize; l++) f[l] = (*calobj)(x[l],dimension);	
}

// select population to build model
void selection(double** population, double* popcost, int dimension, int Popsize, double** selectedpop, int selsize, int* selected)
{
    int i,j,k;

    int* index;
    index = new int[Popsize];
    
    indexx(popcost, index, Popsize);
    
    // select the best half solutions
    for (i = 0; i<selsize; i++)
    {
            k = index[i];
            for (j=0; j<dimension; j++)
                    selectedpop[i][j] = population[k][j];
            selected[k] = 1;
            //selcost[i] = popcost[k];
    }       
    delete[] index;
}

int umda(PObj calobj, double **y, double* point, double& cost, int M, int L, int K, double sigma, double a, double b, int dimension)
{
	int	i,j,k,l,fin,tj;
	double	*tempx=new double[M+2], *prob=new double[M+1], r2, sum; //r3, 
	int *PM;
	double seglength,tempa,tempb,tempmin,tempfit;
	double tpopulation[n][n];

	PM = new int[M+1];

	for (l=0; l<K; l++)
	{
		for (i=0; i<dimension; i++)
		{
			tempa = y[0][i];
			tempb = y[0][i];

			for (j=1; j<= M; j++)
			{
				tempx[j] = y[j-1][i];
				if (tempx[j] < tempa) tempa = tempx[j];
				if (tempx[j] > tempb) tempb = tempx[j];
			}
			
			seglength = (tempb - tempa)/L;

			tempx[0]   = (tempa-seglength >= a)?tempa-seglength:a;
			tempx[M+1] = (tempb+seglength <= b)?tempb+seglength:b;

			tempa = tempx[0];
			tempb = tempx[M+1];
			
			// sort tempx, the temporary of the i variable
			Quicksort(tempx, M+2);

			// generate new point
			// first assign the prob. in each interval, there are M+1 intervals
			
			for (j=0; j<M+1; j++)
				prob[j] = (double) 1/(M+1);

			sum = 0.0;
			for (j=0; j<M+1; j++)
				sum += prob[j];

			for (j=0; j<M+1; j++)
				prob[j] = (double) prob[j]/sum;

			for (j=1; j<M+1; j++)
				prob[j] = prob[j] + prob[j-1];

			k=0; 
			fin = 0;
			r2 = az::rnd::rand();
					
			do
			{
				if (r2 < prob[k])
				{
					//printf(" k = %d, tempx[k] = %f\n",k,tempx[k]);
					if (tempx[k] == tempx[k+1])
					{
						//printf("k = %d, tempa = %f, tempb = %f\n", k, tempa, tempb);

						//tpopulation[l][i] = tempx[k] + sigma*noise(&rnd_uni_init);
						//double tp;
						//do tp = tempx[k] + sigma*noise(&rnd_uni_init); while(tp <= tempa || tp >= tempb);

						tpopulation[l][i] = tempx[k];
						
						if (tpopulation[l][i] <= tempa) {tpopulation[l][i] = tempa; }
						if (tpopulation[l][i] >= tempb) {tpopulation[l][i] = tempb; }
					}
					else
					{
						//r3 = noise(&rnd_uni_init);

						//x[xj][i] = (tempx[k]+ tempx[k+1])/2.0 + sigma*r3*(tempx[k+1]-tempx[k]);
						//tpopulation[l][i] = (tempx[k]+ tempx[k+1])/2.0 + sigma*r3*(tempx[k+1]-tempx[k]);
						//tpopulation[l][i] = tempx[k] + sigma*r3*(tempx[k+1]-tempx[k]);
						tpopulation[l][i] = (tempx[k]+ tempx[k+1])/2.0;
					
						if (tpopulation[l][i] <= tempx[k]) {tpopulation[l][i] = tempx[k];}
						if (tpopulation[l][i] >= tempx[k+1]) {tpopulation[l][i] = tempx[k+1];}
					}

					fin = 1;
				}
				else
					k++;
			}while(fin==0);
			//printf(" genrated tpopulation[%d][%d] = %f\n",l,i,tpopulation[l][i]);
			//getchar();
		}
	}

	tj = 0;
	tempmin = (*calobj)(tpopulation[0], dimension);

	for (j=1; j<K; j++)
	{
		tempfit = (*calobj)(tpopulation[j],dimension);
		if (tempfit < tempmin)
		{
			tempmin = tempfit;
			tj = j;
		}
	}
	for (j=0; j<dimension; j++)
		point[j] = tpopulation[tj][j];

	cost = tempmin;

	delete []PM;
	delete []tempx;
	delete []prob;

	return 0;
}

int umda(PObj calobj, double **y, double *ycost, double* point, double& cost, int M, int L, int K, double sigma, double a, double b, int dimension)
{
	int	i,j,k,l,fin,tj;
	double	*tempx=new double[M+2], *prob=new double[M+1], r2, sum; //r3, 
	int *PM;
	double seglength,tempa,tempb,tempmin,tempfit;
	double tpopulation[n][n];

	//double tmprob[M+1];
	//double suml = 0.0;
	//for(i=0; i<M; i++)
	//	suml += ycost[i];

	PM = new int[M+1];

	for (l=0; l<K; l++)
	{
	for (i=0; i<dimension; i++)
	{
		tempa = y[0][i];
		tempb = y[0][i];

		//tmprob[0] = (double) (1.0/(M+1));
		//printf("temprob[9] = %f\n", tmprob[0]);
		//getchar();

		for (j=1; j<= M; j++)
		{
			tempx[j] = y[j-1][i];
			//tmprob[j] = ycost[j-1]/suml;

			if (tempx[j] < tempa)
				tempa = tempx[j];

			if (tempx[j] > tempb)
				tempb = tempx[j];
		}
		
		seglength = (tempb - tempa)/L;

		tempx[0]   = (tempa-seglength >= a)?tempa-seglength:a;
		tempx[M+1] = (tempb+seglength <= b)?tempb+seglength:b;

		tempa = tempx[0];
		tempb = tempx[M+1];
		
		// sort tempx, the temporary of the i variable
		Quicksort(tempx, M+2);

		//for (j=0; j<M+1; j++)
		//	printf("%1.3f ", tempx[j]);
		//printf("\n");
		//getchar();

		// generate new point
		// first assign the prob. in each interval, there are M+1 intervals
		
		for (j=0; j<M+1; j++)
			prob[j] = (double) 1/(M+1);

		sum = 0.0;
		for (j=0; j<M+1; j++)
			sum += prob[j];

		for (j=0; j<M+1; j++)
			prob[j] = (double) prob[j]/sum;

		for (j=1; j<M+1; j++)
			prob[j] = prob[j] + prob[j-1];

		k=0; 
		fin = 0;
		r2 = az::rnd::rand();
				
		do
		{
			if (r2 < prob[k])
			{
				//printf(" k = %d, tempx[k] = %f\n",k,tempx[k]);
				if (tempx[k] == tempx[k+1])
				{
					//printf("k = %d, tempa = %f, tempb = %f\n", k, tempa, tempb);

					//tpopulation[l][i] = tempx[k] + sigma*noise(&rnd_uni_init);
					//double tp;
					//do tp = tempx[k] + sigma*noise(&rnd_uni_init); while(tp <= tempa || tp >= tempb);

					tpopulation[l][i] = tempx[k];
					
					if (tpopulation[l][i] <= tempa) {tpopulation[l][i] = tempa; }
					if (tpopulation[l][i] >= tempb) {tpopulation[l][i] = tempb; }
				}
				else
				{
					//r3 = noise(&rnd_uni_init);

					//x[xj][i] = (tempx[k]+ tempx[k+1])/2.0 + sigma*r3*(tempx[k+1]-tempx[k]);
					//tpopulation[l][i] = (tempx[k]+ tempx[k+1])/2.0 + sigma*r3*(tempx[k+1]-tempx[k]);
					//tpopulation[l][i] = tempx[k] + sigma*r3*(tempx[k+1]-tempx[k]);
					tpopulation[l][i] = (tempx[k]+ tempx[k+1])/2.0;
				
					if (tpopulation[l][i] <= tempx[k]) {tpopulation[l][i] = tempx[k];}
					if (tpopulation[l][i] >= tempx[k+1]) {tpopulation[l][i] = tempx[k+1];}
				}

				fin = 1;
			}
			else
				k++;
		}while(fin==0);
		//printf(" genrated tpopulation[%d][%d] = %f\n",l,i,tpopulation[l][i]);
		//getchar();
	}
	}

	tj = 0;
	tempmin = (*calobj)(tpopulation[0], dimension);

	for (j=1; j<K; j++)
	{
		tempfit = (*calobj)(tpopulation[j],dimension);
		if (tempfit < tempmin)
		{
			tempmin = tempfit;
			tj = j;
		}
	}
	for (j=0; j<dimension; j++)
		point[j] = tpopulation[tj][j];

	cost = tempmin;

	delete []PM;
	delete []tempx;
	delete []prob;

	return 0;
}

// simplex search
int amoeba(PObj calobj, double* start, double lamda, int fcall, double& cost, int L, int dimension, double a, double b)
// start[n] - the starting point, lamda - p+lamda*e_i
// fcall - the fitness evaluation calls,
{
	const double TINY = 1e-10;	
	int i,ihi,ilo,inhi,j;
	double rtol,ysave,ytry;
	double p[n+1][n], yp[n+1], sum, tmpmin;
	int mpts,ndim,indicator,stop;
	int nfunk;

	double psum[n];

	mpts = dimension+1;
	ndim = dimension;
	nfunk = 0;


	// first generate the additional points from the starting point start[n]
	// and assign the fitness value to this points
	for (i=0; i<dimension; i++)	p[dimension][i] = start[i];

	yp[dimension] = cost;

	tmpmin		= yp[dimension];
	indicator	= 0;
	stop		= 0;

//loop:

	for (j=0; j<dimension; j++)
	{
		for (i=0; i<dimension; i++)
		{
			if (j == i)
			{
				if (lamda <= 1e-10)	lamda = (double) (b-a)/((double)L);
				
				p[j][i] = start[i] + lamda;

				//if (p[j][i] <= a) p[j][i] = a+2*lamda;
				//if (p[j][i] >= b) p[j][i] = b-2*lamda;
			}
			else p[j][i] = start[i];
		}
		yp[j] = (*calobj)(p[j], dimension);
	}

	// then do downhill simplex search from the points p[n+1][n]
	// calculate psum[n]
	for (j=0; j<ndim; j++)
	{
		sum = 0.0;
		for (i=0; i<mpts; i++)	sum += p[i][j];

		psum[j] = sum;
	}

	nfunk += ndim;

	for(;;)
	{
		ilo = 0;

		// first we must determine which point is th highest(worst)
		// next-highest, and lowest(best), by looping over the points in the simplex

		ihi = yp[0]>yp[1] ? (inhi=1,0):(inhi=0,1);

		for (i=0; i<mpts; i++)
		{
			if (yp[i] <= yp[ilo]) ilo = i;
			if (yp[i] > yp[ihi])
			{
				inhi = ihi;
				ihi  = i;
			}
			else if (yp[i] > yp[inhi] && i!= ihi) inhi = i;
		}

		rtol = 2.0*fabs(yp[ihi]-yp[ilo])/(fabs(yp[ihi])+fabs(yp[ilo])+TINY);

		// compute the fractional range from highest to lowest and return if satisfactory

		if (rtol < tol)			// if returning, put best point and value in slot 1
		{
			double temp;
			temp = yp[0];
			yp[0] = yp[ilo];
			yp[ilo] = temp;
			//SWAP(yp[0],yp[ilo]);

			for (i=0; i<ndim; i++)
			{
				double tmp;

				tmp = p[0][i];
				p[0][i] = p[ilo][i];
				p[ilo][i] = tmp;
				//SWAP(p[0][i],p[ilo][i]);
			}
			
			double miny		= yp[0];
			int minyindex	= 0;

			for (i=0; i<mpts; i++)
			{
				if (yp[i] < miny) 
				{
					miny		= yp[i];			
					minyindex	= i;
				}
			}

			cost = yp[minyindex];

			for (i=0; i<ndim; i++) start[i] = p[minyindex][i];

			break;
		}


		if (nfunk >= fcall || stop == 1)  
		{
			double miny = yp[0];
			int minyindex = 0;

			for (i=0; i<mpts; i++)
			{
				if (yp[i] < miny) 
				{
					miny		= yp[i];			
					minyindex	= i;
				}
			}

			cost = yp[minyindex];

			for (i=0; i<ndim; i++)	start[i] = p[minyindex][i];

			return nfunk;

			break;
		}

		nfunk += 2;

		// begin a new iteration. First extrapolate by a factor-1 through the face
		// of the simplex across from the high point, i.e., reflect the simpelx from the high point

		ytry = amotry(p,yp,psum,calobj, ihi,-1.0,dimension,a,b);

		if(ytry <= yp[ilo])
			// Givers a result better thant the best point, so try an additional extrapolation 
			// by a factor 2
		{
			ytry=amotry(p,yp,psum,calobj,ihi,2.0,dimension,a,b);
		}
		else if (ytry >= yp[inhi])
		{
			// the reflected point is worse than the second highest, so look for an intermediate
			// lower point,i.e., do a one-dimensional contraction

			ysave	= yp[ihi];
			ytry	= amotry(p,yp,psum,calobj,ihi,0.5,dimension,a,b);
			if(ytry >= ysave)
			{
				for(i=0; i<mpts; i++)
				{
					if (i != ilo)
					{
						for (j=0; j<ndim; j++)
						{
							p[i][j] = psum[j]=0.5*(p[i][j]+p[ilo][j]);

							if (p[i][j] <= a) p[i][j] = a;

							if (p[i][j] >= b) p[i][j] = b;
						}

						yp[i] = (*calobj)(psum,dimension);
					}
				}

				nfunk += ndim;
				// recompute psum
				for (j=0; j<ndim; j++)
				{
					sum = 0.0;
					for (i=0; i<mpts; i++)
						sum += p[i][j];
					psum[j] = sum;
				}
			}
		}
		else
			--nfunk;
	}

	return nfunk;
}

int optimizer(PObj CalObj, int MaxEva, int MaxLS, int PopSize, int Dimension, double BoundLow, double BoundUpp, double* BestX)
{
	int		xk = 3;
	int		i,j,k,s;
	int		gen = 0;
	double	fbest = 1.0E100,fold;
	int		L = PopSize;
	int		suc = 0;
	int		count = 0;

	double lambda = (BoundUpp-BoundLow)/((double) PopSize);

	az::rnd::seed((long) time(NULL));

	double** Population	= new double*[PopSize]; for (i=0; i<PopSize; i++)	Population[i] = new double[Dimension];
	double* PopCost		= new double[PopSize];			
	double** SelPop		= new double*[PopSize/2]; for(i = 0; i<PopSize/2; i++) SelPop[i] = new double[Dimension];
	int* SelIndex		= new int[PopSize];
	int* index			= new int[PopSize];	
	int* NTIX			= new int[PopSize];	
	for (j=0; j<PopSize; j++)	NTIX[j] = 0;
	
	//initialize population
	UDM_INITIALIZE(CalObj, 3, BoundLow, BoundUpp, Dimension, PopSize, Population, PopCost);

	for (j=0; j<PopSize; j++)
	{
		count ++;
		if (PopCost[j] < fbest)
		{
			fbest = PopCost[j];
			for(s=0; s<Dimension; s++) BestX[s] = Population[j][s];
		}
	}

	do
	{
		fold = fbest;
		
		for (i=0; i<PopSize; i++) SelIndex[i] = 0;

		// selection
		selection(Population, PopCost, Dimension, PopSize, SelPop, PopSize/2, SelIndex);

		// sort population
		indexx(PopCost, index, PopSize);

		// update the solution's status
		for (i=0; i<PopSize; i++)
		{
			if (SelIndex[i] == 0)	// the i-th solution is discard, no matter its previous status
				NTIX[i] = 0;
			else
				NTIX[i] = NTIX[i] + 1;
		}
		
		// do expensive local search
		if (suc >= 4)
		{
			// only the best two solutions maybe carry out local search
			for (k=0; k<1; k++)	
			{
				// if the i-th solution survives for th times, then we think it is a good starting point
				// for expensive local search
				i = index[k];

				// There are two options: 1. downhill simplex (amoeba) 2. UOBYQB 
				// The programs may have some bugs, will modify later
				// but you can use the program
				count += amoeba(CalObj, Population[i], lambda, MaxLS, PopCost[i], L, Dimension, BoundLow, BoundUpp);
				count += Dimension;

				// local search 
				//count += UOBYQB(fun_type, X,1.0e-1,MAXFUN,PopCost[i],a,b,Dimension);

				NTIX[i] = 0;

				if (PopCost[i] < fbest)
				{
					fbest = PopCost[i];
					for(s=0; s<Dimension; s++) BestX[s] = Population[i][s];
				}
			}
		}

		// eda search
		for (j=0; j<PopSize; j++)
		{
			if (SelIndex[j] == 0 || NTIX[j] ==0) // means that the solution is not SelPop, need to be replaced
			//if (SelIndex[j] == 0) // means that the solution is not SelPop, need to be replaced
			{
				double* tmpx = new double[Dimension];
				double tmpcost;

				umda(CalObj, SelPop, tmpx, tmpcost, PopSize/2, L, xk, 0.01, BoundLow, BoundUpp, Dimension);
				count += xk;

				// replace the current one
				for (i=0; i<Dimension; i++)	Population[j][i] = tmpx[i];

				PopCost[j] = tmpcost;

				if (PopCost[j] < fbest)
				{
					fbest = PopCost[j];
					for(s=0; s<Dimension; s++) BestX[s] = Population[j][s];
				}
				delete[] tmpx;
			}
			// else we reserve the solution
		}

		gen++;

		if (fbest < fold)
			suc = 0;
		else
			suc ++;	
	}while(count<MaxEva);// End Generation Loop

	delete[] SelIndex;
	delete[] index;	
	delete[] PopCost;
	for(i=0; i<PopSize; i++)  delete[] Population[i];
	for(i=0; i<PopSize/2; i++)delete[] SelPop[i];
	delete[] SelPop;
	delete[] Population;
	delete[] NTIX;

	return count;
}
}//namespace edal
