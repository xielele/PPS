/*! \file	Generator_LS_XHill.cpp
	
	\brief	Local Search strategy with XHill
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Apr.24 2006 create
*/

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
	//!\brief	Generator
	//!\param	popnew	offspring population
	//!\param	pop		parent population
	//!\retun	offspring population
	CPopulationMO& XHill::Generate(CPopulationMO& popnew, CPopulationMO& pop)
	{
		//the point index to be selected to do LS
		std::vector<unsigned int> index;
	
		//Step 1: select the pionts to do LS
		unsigned int r0,r1,r2;
		pop.RankSort();
		pop.RankSize(1, r0, r1);
		LSBase::SelectStartPoint(pop, r0, r1, std::min(r1-r0, mNoL), index);
		
		//Step 2: LS
		popnew.Resize((unsigned int)(index.size())*mNoE);
		// define the range of reference points: [r1, r2)
		r2 = std::min(r1 + 2*index.size(), pop.Size());
		while(r2<pop.Size() && pop[r2].Rank() == pop[r2-1].Rank()) r2++;
		if(r2-r1 < index.size()) r1 = r0;
		// LS one by one
		unsigned int i,j,k,p; double dis,disp,t;
		for(i=0; i<index.size(); i++)
		{
			//find the nearest reference point
			disp = 1.0E200; p=r1;
			for(j=r1; j<r2; j++)
			{
				dis = 0.0;
				for(k=0; k<pop.P().XSize(); k++) dis += (pop[index[i]][k] - pop[j][k])*(pop[index[i]][k] - pop[j][k]);
				if(dis>0 && dis<disp) {disp=dis; p=j;}
			}

			//a local model built by point index[i] and point p
			disp = sqrt(disp)*0.2/sqrt(double(pop.P().XSize())); 
			for(j=0; j<mNoE; j++)
			{
				t = -1.0 + 2.0*rnd::rand(double(j), double(j+1.0))/double(mNoE);
				for(k=0; k<pop.P().XSize(); k++)
				{
					popnew[i*mNoE+j][k] = pop[index[i]][k] + (pop[index[i]][k]- pop[p][k])*t + disp*rnd::gaussian();
					//the border handeling is very important
					if(popnew[i*mNoE+j][k]>pop.P().XUpp(k))		 popnew[i*mNoE+j][k] = 0.5*(pop[index[i]][k] + pop.P().XUpp(k));
					else if(popnew[i*mNoE+j][k]<pop.P().XLow(k)) popnew[i*mNoE+j][k] = 0.5*(pop[index[i]][k] + pop.P().XLow(k));
				}
			}
		}
		
		popnew.Evaluate();

		return popnew;
	}
} //namespace LS

} //namespace GEN

} //namespace EA
