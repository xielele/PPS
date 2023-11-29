/*! \file	Generator_LS_Valley.cpp
	
	\brief	Local Search strategy with valley modeling
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Apr.24 2006 create
*/

#include "algorithm/Fitting.h"
#include "mea/Generator_LS.h"

//!\brief EA Evolutionary Algorithms
namespace EA
{
//!\brief offspring generate strategies
namespace GEN
{
//!\brief local search strategies
namespace LS
{
	//!\brief	constractor
	Valley::Valley()
	{
		mNoG = 5;
	}

	//!\brief	Generator
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\retun	offspring population
	CPopulationMO& Valley::Generate(CPopulationMO& popnew, CPopulationMO& pop)
	{
		unsigned int i,k,s,r0,r1,ibest,iworst,iL,iG,iN;
		double dis,dis1,fit;

		if(mNoG < pop.Size()) mNoG = std::min((unsigned int)10, pop.Size()/2);

		std::vector<unsigned int>	index;			//the point index to be selected to do LS

		std::vector< double >	fita(pop.Size()),	//fitness along search direction
								fitat(mNoG),		//fitness along search direction
								fito(pop.Size()),	//fitness orthogonal to search direction
								fitot(mNoG),		//fitness orthogonal to search direction
								direction(pop.P().FSize()),	//normalized search direction
								nadir(pop.P().FSize()),		//a nearest point which can dominate all current best points
								Y(mNoG),			//fitness
								X(mNoG),			//decision variable
								C(3);				//coefficient
		
		CPopulationMO popt(pop.P()); popt.Resize(mNoG);

		//Step 1: select the pionts to do LS
		pop.RankSort();
		pop.RankSize(1, r0, r1);
		LSBase::SelectStartPoint(pop, r0, r1, std::min(r1-r0, mNoL), index);
		
		//Step 2: find the nadir point
		for(k=0; k<nadir.size(); k++)
		{
			nadir[k] = -1.0E200;
			for(i=r0; i<r1; i++) if(pop[i].F(k) > nadir[k]) nadir[k] = pop[i].F(k);
		}

		//Step 3: LS
		popnew.Resize((unsigned int)(index.size())*mNoE);
		// LS one by one
		iN = 0;
		for(iL=0; iL<index.size(); iL++)
		{
			//Step 3.1: define the search direction and normalize it
			dis = 0.0;
			for(k=0; k<nadir.size(); k++) 
			{
				direction[k] = pop[index[iL]].F(k) - nadir[k];
				dis			+= direction[k] * direction[k];
			}
			dis = sqrt(dis);
			for(k=0; k<direction.size(); k++) direction[k] /= dis;

			//Step 3.2: recalculate the fitness along and orthogonal to the search direction
			for(i=0; i<pop.Size(); i++)
			{
				fita[i] = 0.0; dis = 0.0;
				for(k=0; k<direction.size(); k++) 
				{
					fita[i] += (pop[i].F(k) - pop[index[iL]].F(k))*direction[k];
					dis		+= (pop[i].F(k) - pop[index[iL]].F(k))*(pop[i].F(k) - pop[index[iL]].F(k));
				}
				fito[i] = sqrt(dis - fita[i]*fita[i]);
			}

			//Step 3.3: select the nearest points in orthogonal directions
			// is it necessary to add some panelty when selecting??????
			for(k=0; k<mNoG; k++)
			{
				dis = (k==0) ? -1.0E200:fitot[k-1];
				dis1= 1.0E200; s = 0;
				for(i=0; i<pop.Size(); i++) if(fito[i]>dis && fito[i]<dis1) {s=i; dis1=fito[i];}
				popt[k]  = pop[s];
				fitat[k] = fita[s];
				fitot[k] = fito[s];
			}

			//Step 3.4: build response surface model
			ibest = iworst = 0;
			for(k=0; k<mNoG; k++) 
			{
				Y[k] = fitat[k] - fitot[k]*fitot[k]/sqrt(fitot[k]*fitot[k]+fitat[k]*fitat[k]);

				if(Y[k]>Y[ibest])  ibest  = k;
				if(Y[k]<Y[iworst]) iworst = k;
			}
			iG=0;
			while(iG++<mNoG)
			{
				for(s=0; s<pop.P().XSize(); s++)
				{
					for(i=0; i<mNoG; i++) X[i] = popt[i][s];
					LEARN::poly_fit(X,Y,2,C);
					if(C[2]<0) popnew[iN][s] = -0.5*C[1]/C[2];
					else popnew[iN][s] = popt[ibest][s];
					if(popnew[iN][s]>popt.P().XUpp(k))		popnew[iN][s] = 0.5*(popt[rnd::rand((unsigned int)0,mNoG)][s] + popt.P().XUpp(k));
					else if(popnew[iN][s]<popt.P().XLow(k))	popnew[iN][s] = 0.5*(popt[rnd::rand((unsigned int)0,mNoG)][s] + popt.P().XLow(k));
				}
				popnew[iN].Evaluate();
				fit = 0.0; dis = 0.0;
				for(k=0; k<direction.size(); k++) 
				{
					fit += (popnew[iN].F(k) - pop[index[iL]].F(k))*direction[k];
					dis += (popnew[iN].F(k) - pop[index[iL]].F(k))*(popnew[iN].F(k) - pop[index[iL]].F(k));
				}
				fit -= (dis - fit*fit)/sqrt(dis);
				iN++;
				if(fit>Y[iworst])
				{
					Y[iworst]	= fit;
					popt[iworst]= popnew[iN-1];
					for(i=0; i<mNoG; i++) 
					{
						if(i==0) { ibest = iworst = 0; }
						else { if(Y[i]>Y[ibest])  ibest  = i; if(Y[i]<Y[iworst]) iworst = i; }
					}
				}
				else
				{
					break;
				}
			}
		}
		popnew.Erase(iN);
		return popnew;
	}
} //namespace LS

} //namespace GEN

} //namespace EA
