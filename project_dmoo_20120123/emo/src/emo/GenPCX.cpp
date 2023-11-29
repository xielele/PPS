/*! \file	Generator_PCX.cpp
	
	\brief	Evolutionary Aglorithm Generators
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Mar.30 2006 redesign
*/

#include <cmath>
#include <algorithm>
#include <vector>
#include <float.h>
#include "emo/Parameter.h"
#include "alg//Random.h"
#include "emo/IndividualMO.h"
#include "emo/PopulationMO.h"
#include "emo/Gen.h"

namespace az
{
namespace mea
{
namespace gen
{

// Parent-Centric Recombination (PCX)
bool PCX(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{
	unsigned int i,k,s;
	double tmp1,tmp2,meanD;
	std::vector< double >	centroid(pop.P().XSize()),
							dist(pop.Size()),
							D(pop.Size()),
							tmpV1(pop.P().XSize()),
							tmpV2(pop.P().XSize());
	std::vector< std::vector<double> > diff(pop.Size());

	//calculate the centroid
	for(k=0; k<pop.P().XSize(); k++)
	{
		centroid[k] = 0.0;
		for(i=0; i<pop.Size(); i++) centroid[k] += pop[i][k];
		centroid[k] /= double(pop.Size());
	}

	// calculate the distace (d) from centroid to the index parent arr1[0]  
	// also distance (diff) between index and other parents are computed
	for(i=0; i<diff.size(); i++)
	{
		dist[i] = 0.0;
		diff[i].resize(pop.P().XSize());
		for(k=0; k<pop.P().XSize(); k++)
		{
			if(i==0) diff[i][k] = centroid[k] - pop[0][k];
			else	 diff[i][k] = pop[i][k] - pop[0][k];
			dist[i] += diff[i][k]*diff[i][k];
		}
		//dist[i] = sqrt(dist[i]);
		if(dist[i]<1.0E-40) return false;
	}

	// orthogonal directions are computed (see the paper)
	meanD = 0.0;
	for(i=1; i<D.size(); i++)
	{
		tmp1 = 0; for(k=0; k<pop.P().XSize(); k++) tmp1 += diff[i][k]*diff[0][k];
		tmp2 = tmp1*tmp1/dist[0];
		if(tmp2>=dist[i]) return false;
		D[i] = sqrt(dist[i]-tmp2);
		meanD += D[i];
	}
	meanD /= double(D.size()-1.0);

	//int sss;
	popnew.Resize(sizenew);

	for(i=0; i<popnew.Size(); i++)
	{
		// Next few steps compute the child, by starting with a random vector
		for(k=0;k<pop.P().XSize();k++)
			tmpV1[k] = rnd::gaussian()*meanD*pop.P().ETA_PCX();

		tmp1 = 0.0; for(s=0; s<pop.P().XSize(); s++) tmp1 += tmpV1[s]*diff[0][s];
		for(k=0;k<pop.P().XSize();k++)
		{
			tmpV2[k] = tmpV1[k] - tmp1*diff[0][k] / dist[0];//(dist[0]*dist[0]);
		}

		tmp1 = rnd::gaussian()*pop.P().ETA_PCX();

		for(k=0;k<pop.P().XSize();k++)
		{
			popnew[i][k] = pop[i][k] + tmpV2[k] + tmp1*diff[0][k];
			if(popnew[i][k] < pop.P().XLow(k))		popnew[i][k] = 0.5*(pop.P().XLow(k)+pop[i][k]);
			else if(popnew[i][k] > pop.P().XUpp(k)) popnew[i][k] = 0.5*(pop.P().XUpp(k)+pop[i][k]);
		}
	}
	return true;
}

//=============================================================================================
// PCX generator
CPopulationMO&  XPCX::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{
	unsigned int i,index,i1,i2,nson;
	
	popnew.Resize(sizenew);

	CPopulationMO parent(pop.P()),son(pop.P());
	parent.Resize(3); son.Resize(3);
	
	std::vector<unsigned int> a(pop.Size()+2),b(pop.Size()+2);
	for( i=0; i<pop.Size(); i++ ) a[i]=b[i]=i;
	a[pop.Size()]   = b[pop.Size()]   = rnd::rand((unsigned int)0,pop.Size());
	a[pop.Size()+1] = b[pop.Size()+1] = rnd::rand((unsigned int)0,pop.Size());
	std::random_shuffle(a.begin(),a.end());
	std::random_shuffle(b.begin(),b.end());

	index = 0; nson = 3;
	while(index<pop.Size())
	{
		for(i=0; i<3; i++) 
			if(pop[a[index+i]] < pop[b[index+i]]) 
				parent[i] = pop[a[index+i]]; 
			else 
				parent[i] = pop[b[index+i]];

		if(index > sizenew - 3) 
		{
			nson = sizenew-index;  
			son.Resize(nson);
		}

		//Crossover
		while(!mea::gen::PCX(nson, son, parent))
		{
			for(i=0; i<3; i++) 
			{
				i1 = rnd::rand((unsigned int)0,pop.Size());
				i2 = i1; while(i2==i1) i2 = rnd::rand((unsigned int)0,pop.Size());
				if(pop[i1] < pop[i2]) 
					parent[i] = pop[i1]; 
				else 
					parent[i] = pop[i2];
			}
		}
		
		//Mutation
		for(i=0; i<son.Size(); i++)
		{
			mea::gen::PM(son[i]);
			popnew[index+i] = son[i];
		}
		index += son.Size();
	}
	return popnew;
}

} //namespace gen
} //namespace mea
} //namespace az
