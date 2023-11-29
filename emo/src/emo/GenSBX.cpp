/*! \file	Generator_SBX.cpp
	
	\brief	Evolutionary Aglorithm Generators
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Mar.30 2006 redesign
*/

#include <algorithm>
#include <float.h>
#include "emo/Gen.h"

#if defined(WIN32)
    #define wxFinite(n) _finite(n)
#elif defined(_LINUX)
    #define wxFinite(n) finite(n)
#else
    #define wxFinite(n) ((n)==(n))
#endif

namespace az
{
namespace mea
{
namespace gen
{
// Polynomial mutation operator (PM)
CIndividualMO& PM(CIndividualMO& ind)
{
	double rnd,delta,deltaq,lown,uppn;
	double val,xy;

	for(unsigned int n = 0; n < ind.P().XSize(); n++)
	{
		lown = ind.P().XLow(n); uppn = ind.P().XUpp(n);
		if(rnd::rand() <= ind.P().Pm())
		{
			if(ind[n] > lown)
			{			
				rnd = rnd::rand();
				
				if(rnd <= 0.5)
				{
					delta	= (ind[n] - lown)/(uppn - lown);
					xy		= 1.0 - delta;
					val		= 2*rnd + (1.0 - 2.0*rnd)*(pow(xy, (ind.P().ETA_PM()+1.0)));
					deltaq	= pow(val,1.0/(ind.P().ETA_PM() +1.0)) - 1.0;
				}
				else
				{
					delta	= (uppn - ind[n])/(uppn - lown);
					xy		= 1.0 - delta;
					val		= 2.0*(1.0 - rnd)+2.0*(rnd - 0.5)*(pow(xy, (ind.P().ETA_PM()+1.0)));
					deltaq	= 1.0 - (pow(val, 1.0/(ind.P().ETA_PM() +1.0)));
				}

				ind[n] = ind[n] + deltaq * (uppn - lown);

				if(ind[n] < lown) ind[n] = lown;
				else if(ind[n] > uppn) ind[n] = uppn;
			}
			else // ind[n]  == LowBounds[n]
			{
				ind[n] = rnd::rand(lown, uppn) ;
			}
			//if( _finite(ind[n]) == 0 )
			//{
			//	int ssss=0;
			//}
		}//end of if
	} // end of loop
	return ind;
}

// Simulated Binary crossover (SBX)
void SBX(CIndividualMO& son1, CIndividualMO& son2, CIndividualMO& parent1, CIndividualMO& parent2)
{
	double betaq,beta,alpha,pmax,pmin;
	
	if(rnd::rand() < parent1.P().Pc())
	{
		for(unsigned int n = 0; n < parent1.P().XSize(); n++)
		{
			double lown = parent1.P().XLow(n);
			double uppn = parent1.P().XUpp(n);
			//Crossover
			if(rnd::rand() < 0.5)
			{
				if(fabs(parent1[n] - parent2[n]) > 1.0E-6) 
				{
					pmax	= std::max<double>(parent1[n] , parent2[n]);
					pmin	= std::min<double>(parent1[n] , parent2[n]);

					//Find beta value*/
					if((pmin - lown) > (uppn - pmax))
						beta = 1 + (2.0 * (uppn - pmax) / (pmax-pmin));
					else
						beta = 1 + (2.0 * (pmin - lown) / (pmax-pmin));

					//*Find alpha*/
					double expp  = parent1.P().ETA_SBX() + 1.0;
					beta  = 1.0/beta;
					alpha = 2.0 - pow(beta, expp);

					double rnd3 = rnd::rand();
					if(rnd3 <= 1.0/alpha)
					{
						alpha = alpha*rnd3;
						expp  = 1.0/(parent1.P().ETA_SBX()+1.0);
						betaq = pow(alpha,expp);
					}
					else
					{
						alpha = alpha*rnd3;
						alpha = 1.0/(2.0-alpha);
						expp  = 1.0/(parent1.P().ETA_SBX()+1.0);
						betaq = pow(alpha,expp);
					}

					if( wxFinite(betaq) == 0 )
					{
						betaq = 0.0;
					}

				}
				else
				{
					betaq= 1.0;
					pmin = parent1[n]; 
					pmax = parent2[n];
				}
				//Generating two children
				if(rnd::rand()<0.5)
				{
					son1[n]	= 0.5*((pmin+pmax) - betaq*(pmax-pmin));
					son2[n]	= 0.5*((pmin+pmax) + betaq*(pmax-pmin));
				}
				else
				{
					son1[n]	= 0.5*((pmin+pmax) + betaq*(pmax-pmin));
					son2[n]	= 0.5*((pmin+pmax) - betaq*(pmax-pmin));
				}

				if(son1[n] < lown) son1[n] = lown;
				else if(son1[n] > uppn) son1[n] = uppn;
				if(son2[n] < lown) son2[n] = lown;
				else if(son2[n] > uppn) son2[n] = uppn;

				//if(parent1[n] < parent2[n])
				//{
				//	son1[n]	= 0.5*((pmin+pmax) - betaq*(pmax-pmin));
				//	son2[n]	= 0.5*((pmin+pmax) + betaq*(pmax-pmin));
				//}
				//else
				//{
				//	son1[n]	= 0.5*((pmin+pmax) + betaq*(pmax-pmin));
				//	son2[n]	= 0.5*((pmin+pmax) - betaq*(pmax-pmin));
				//}
			}
			//Copy
			else
			{
				son1[n]	= parent1[n];
				son2[n]	= parent2[n];
			}
		}
	}
	else
	{
		for(unsigned int n = 0; n < parent1.P().XSize(); n++)
		{
			son1[n]	= parent1[n];
			son2[n]	= parent2[n];
		}
	}
}

// SBX generator
CPopulationMO& XSBX::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{
	unsigned int i,par;

	popnew.Resize(sizenew);
	
	std::vector<unsigned int> a(sizenew),b(sizenew);
	for(i=0; i<sizenew; i++) a[i] = b[i] = i % pop.Size();
	std::random_shuffle(a.begin(),a.end());
	std::random_shuffle(b.begin(),b.end());

	for(par=0; par+1<sizenew; par+=2)
	{
		//Crossover
		mea::gen::SBX( popnew[par], popnew[par+1], pop[ pop[a[par]] < pop[a[par+1]] ? a[par] : a[par+1] ], pop[ pop[b[par]] < pop[b[par+1]] ? b[par] : b[par+1] ] );
		
		//Mutation
		mea::gen::PM( popnew[par] );
		mea::gen::PM( popnew[par+1] );
	}
	if(sizenew % 2 != 0)
	{
		//Copy
		popnew[sizenew-1] = pop[ a[sizenew-1] ];
		
		//Mutation
		mea::gen::PM(popnew[sizenew-1]);
	}
	a.clear();b.clear();
	return popnew;
}

} //namespace gen
} //namespace mea
} //namespace az
