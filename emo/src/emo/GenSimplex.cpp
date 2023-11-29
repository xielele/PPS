/*! \file	Generator_Simplex.cpp
	
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
#include "alg//Random.h"
#include "alg//Matrix.h"
#include "emo/Parameter.h"
#include "emo/Gen.h"

namespace az
{
namespace mea
{
namespace gen
{

// constructor
XSimplex::XSimplex()
{
	mBeta  = 0.1;
	mAlpha = 5.0;
	mExten = 0.0;
	mType  = 0;
	mMaxGen= 1;
}

// set algorithm parameter
void XSimplex::Set(unsigned int maxgen, double beta, double alpha, unsigned int type)
{
	mMaxGen= maxgen;
	mAlpha = alpha;
	mBeta  = beta;
	mAlpha = alpha;
	mType  = type;
	mMaxGen= maxgen;
}

// Generator
CPopulationMO& XSimplex::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& pop)
{	
	unsigned int i,j;

	popnew.Resize(sizenew);

	mMatDis.resize(pop.Size());
	for(i=0; i<pop.Size(); i++) mMatDis[i].resize(pop.Size());
	for(i=0; i<pop.Size(); i++)
	{
		mMatDis[i][i] = 1.0E100;
		for(j=i+1; j<pop.Size(); j++)
			mMatDis[j][i] = mMatDis[i][j] = Distance(pop[i], pop[j]);
	}

	switch(pop.P().FSize())
	{
	case 2:
		mExten = 0.5;
		Generate1D(popnew,pop);
		break;
	case 3:
		mExten = 1.414;
		Generate2D(popnew,pop);
		break;
	default:
		break;
	}

	return popnew;
}

// generator for 1-D structure
CPopulationMO& XSimplex::Generate1D(CPopulationMO& popnew, CPopulationMO& pop)
{
	unsigned int i,j,index, seedindex, neiNo=0;
	double delta;
	std::vector<unsigned int> p(2); 
	std::vector<double> t;
	std::vector<unsigned int> seeds(pop.Size()),candidate(pop.Size());

	//initialize the neighborhood parameters
	switch(mType)
	{
	case 0:
		break;
	case 1:
	case 2:
		//initialize the distance matrix
		mMatDis.resize(pop.Size());
		for(i=0; i<pop.Size(); i++) mMatDis[i].resize(pop.Size());
		for(i=0; i<pop.Size(); i++)
		{
			mMatDis[i][i] = 1.0E100;
			for(j=i+1; j<pop.Size(); j++)
				mMatDis[j][i] = mMatDis[i][j] = Distance(pop[i], pop[j]);
		}
		
		//find the index points (seed points)
		for(i=0; i<pop.Size(); i++) candidate[i] = i;
		MaxiMin(pop, seeds, pop.Size());
		break;
	default:
		break;
	}
	switch(mType)
	{
	case 0:
		break;
	case 1:
		neiNo	= 10 + (unsigned int)((pop.Size()-20)*pow(mCurGen++/double(mMaxGen), mAlpha));
		break;
	case 2:
		neiNo	= 10 + (unsigned int)((pop.Size()-20)*(1.0-pow(mCurGen++/double(mMaxGen), mAlpha)));
		break;
	default:
		break;
	}

	index = seedindex = 0;
	while(index < popnew.Size())
	{
		t.clear();
		switch(mType)
		{
		case 0:
		default:
			do{p[0]=rnd::rand((unsigned int)(0), pop.Size());}while(false);
			do{p[1]=rnd::rand((unsigned int)(0), pop.Size());}while(p[1]==p[0]);
			break;
		case 1:
			p[0] = seeds[seedindex++];
			std::random_shuffle(candidate.begin(), candidate.end());
			for(i=0; i<10; i++)
				for(j=i+1; j<neiNo; j++)
					if(mMatDis[p[0]][candidate[i]]>mMatDis[p[0]][candidate[j]])
						std::swap(candidate[i],candidate[j]);
			std::random_shuffle(candidate.begin(), candidate.begin()+10);
			p[0] = candidate[0];
			p[1] = candidate[1];
			break;
		case 2:
			p[0] = seeds[seedindex++];
			for(i=0; i<neiNo; i++)
				for(j=i+1; j<pop.Size(); j++)
					if(mMatDis[p[0]][candidate[i]]>mMatDis[p[0]][candidate[j]])
						std::swap(candidate[i],candidate[j]);
			std::random_shuffle(candidate.begin(), candidate.begin()+neiNo);
			p[0] = candidate[0];
			p[1] = candidate[1];
			break;
		}

		// sort two parents
		if(pop[p[0]].Rank()>pop[p[1]].Rank()) std::swap(p[0],p[1]);

		//strategy: how many offspring points to create
		//  ----o----|----o----
		if( pop[p[0]].Rank()<pop[p[1]].Rank() ||	//one non-dominated point
			index == popnew.Size()-1 )				//create one offspring point
		{
			t.resize(1);
			t[0] = rnd::rand(-mExten,0.5);
		}
		else										//two non-dominated points
		{
			t.resize(2);
			t[0] = rnd::rand(-mExten,0.5);
			t[1] = rnd::rand(    0.5,1.0+mExten);
		}

		//mutation variance
		delta  = Distance(pop[p[0]],pop[p[1]]);
		delta *= mBeta/(sqrt(double(pop.P().XSize())));

		//create offspring points
		for(j=0; j<t.size(); j++)
		{
			for(i=0; i<pop.P().XSize(); i++) 
			{
				popnew[index][i] = pop[p[0]][i] + (pop[p[1]][i]- pop[p[0]][i])*t[j] + delta*rnd::gaussian();

				//the border handeling is very important
				if(popnew[index][i]>pop.P().XUpp(i))		popnew[index][i] = 0.5*(pop[p[j]][i] + pop.P().XUpp(i));
				else if(popnew[index][i]<pop.P().XLow(i))	popnew[index][i] = 0.5*(pop[p[j]][i] + pop.P().XLow(i));
			}
			index++;
		}
	}
	return popnew;
}

// generator for 2-D structure
CPopulationMO& XSimplex::Generate2D(CPopulationMO& popnew, CPopulationMO& pop)
{
	unsigned int i,j,j1,size,index, seedindex, neiNo=0;
	double delta,t1,t2;
	std::vector<double> O(pop.P().XSize());
	std::vector<unsigned int> p(3);
	std::vector<unsigned int> seeds(pop.Size()),candidate(pop.Size());

	//initialize the neighborhood parameters
	switch(mType)
	{
	case 0:
		break;
	case 1:
	case 2:
		//initialize the distance matrix
		mMatDis.resize(pop.Size());
		for(i=0; i<pop.Size(); i++) mMatDis[i].resize(pop.Size());
		for(i=0; i<pop.Size(); i++)
		{
			mMatDis[i][i] = 1.0E100;
			for(j=i+1; j<pop.Size(); j++)
				mMatDis[j][i] = mMatDis[i][j] = Distance(pop[i], pop[j]);
		}
		
		//find the index points (seed points)
		for(i=0; i<pop.Size(); i++) candidate[i] = i;
		MaxiMin(pop, seeds, pop.Size());
		break;
	default:
		break;
	}
	switch(mType)
	{
	case 0:
		break;
	case 1:
		neiNo	= 10 + (unsigned int)((pop.Size()-20)*pow(mCurGen++/double(mMaxGen), mAlpha));
		break;
	case 2:
		neiNo	= 10 + (unsigned int)((pop.Size()-20)*(1.0-pow(mCurGen++/double(mMaxGen), mAlpha)));
		break;
	default:
		break;
	}

	index = seedindex = 0;
	while(index < popnew.Size())
	{
		switch(mType)
		{
		case 0:
		default:
			do{p[0]=rnd::rand((unsigned int)(0), pop.Size());}while(false);
			do{p[1]=rnd::rand((unsigned int)(0), pop.Size());}while(p[1]==p[0]);
			do{p[2]=rnd::rand((unsigned int)(0), pop.Size());}while(p[2]==p[0]||p[2]==p[1]);
			break;
		case 1:
			p[0] = seeds[seedindex++];
			std::random_shuffle(candidate.begin(), candidate.end());
			for(i=0; i<10; i++)
				for(j=i+1; j<neiNo; j++)
					if(mMatDis[p[0]][candidate[i]]>mMatDis[p[0]][candidate[j]])
						std::swap(candidate[i],candidate[j]);
			std::random_shuffle(candidate.begin(), candidate.begin()+10);
			p[0] = candidate[0];
			p[1] = candidate[1];
			p[2] = candidate[2];
			break;
		case 2:
			p[0] = seeds[seedindex++];
			for(i=0; i<neiNo; i++)
				for(j=i+1; j<pop.Size(); j++)
					if(mMatDis[p[0]][candidate[i]]>mMatDis[p[0]][candidate[j]])
						std::swap(candidate[i],candidate[j]);
			std::random_shuffle(candidate.begin(), candidate.begin()+neiNo);
			p[0] = candidate[0];
			p[1] = candidate[1];
			p[2] = candidate[2];
			break;
		}

		for(i=0; i<3; i++) for(j=i+1; j<3; j++) if(pop[p[j]].Rank()<pop[p[i]].Rank()) std::swap(p[i],p[j]);

		//one offspring point
		if(	pop[p[0]].Rank()< pop[p[1]].Rank() ||
			pop[p[0]].Rank()==pop[p[1]].Rank() && index == popnew.Size()-1)
		{
			size = 1;
		}
		//two offspring points
		else if(pop[p[0]].Rank()==pop[p[1]].Rank() && pop[p[1]].Rank()< pop[p[2]].Rank() ||
				pop[p[0]].Rank()==pop[p[1]].Rank() && pop[p[1]].Rank()==pop[p[2]].Rank() && index == popnew.Size()-2)
		{
			size = 2;
		}
		//three offspring points
		else
		{
			size = 3;
		}

		//mutation variance
		delta  = Distance(pop[p[0]],pop[p[1]])+Distance(pop[p[0]],pop[p[2]])+Distance(pop[p[1]],pop[p[2]]);
		delta *= mBeta/(3.0*sqrt(double(pop.P().XSize())));

		//center point
		for(i=0; i<O.size(); i++) O[i] = (pop[p[0]][i] + pop[p[1]][i] + pop[p[2]][i])/3.0;

		for(j=0; j<size; j++)
		{
			t1 = rnd::rand(0.0,mExten); t2 = rnd::rand(0.0,mExten);
			j1 = (j + (rnd::rand(0.0,1.0)<0.5 ? 1:2)) % 3;

			for(i=0; i<pop.P().XSize(); i++) 
			{
				popnew[index][i] =	O[i] 
									+ t1*(pop[p[j]][i] - O[i]) 
									+ t2*(mExten-t1)*((pop[p[j]][i] + pop[p[j1]][i])/2.0 - O[i])
									+ delta*rnd::gaussian();

				//border handeling
				if(popnew[index][i]>pop.P().XUpp(i))		popnew[index][i] = 0.5*(pop[p[j]][i] + pop.P().XUpp(i));
				else if(popnew[index][i]<pop.P().XLow(i))	popnew[index][i] = 0.5*(pop[p[j]][i] + pop.P().XLow(i));
			}
			index++;
		}
	}
	return popnew;
}

// calculate the distance between two individuals
double XSimplex::Distance(CIndividualMO& ind1, CIndividualMO& ind2)
{
	double dis=0.0; unsigned int i;
	for(i=0; i<ind1.P().XSize(); i++) dis += (ind1[i]-ind2[i])*(ind1[i]-ind2[i]);
	return sqrt(dis);
}

// MaxiMin sort 
void XSimplex::MaxiMin(CPopulationMO& pop, std::vector<unsigned int>& index, unsigned int size)
{
	unsigned int i, j, k, iindex;
	double fmin;
	std::vector< double >	dismin;						// minimum distance
	std::vector< bool >		exist;						// flag to show whether the point is selected

	if(size<pop.P().FSize()) size = pop.P().FSize();
	index.resize(size);

	//Step 1: calculate the distance matrix
	dismin.resize(pop.Size());
	exist.resize(pop.Size());
	for(i=0; i<pop.Size(); i++)
	{
		exist[i]	= false;
		dismin[i]	= 1.0E200;
	}
	
	//Step 3: select the solutions with minimum Fi(x)
	for(k=0; k<pop.P().FSize(); k++)
	{
		fmin = 1.0E200; iindex = 0;
		for(i=0; i<pop.Size(); i++) 
			if(!exist[i] && pop[i].F(k)<fmin)
			{
				fmin  = pop[i].F(k);
				iindex= i;
			}
		exist[iindex] = true; index[k] = iindex; 
		for(i=0; i<pop.Size(); i++) if(!exist[i] && mMatDis[i][iindex]<dismin[i]) dismin[i] = mMatDis[i][iindex];
	}

	//Step 4: select point one by one
	for(i=pop.P().FSize(); i<size; i++)
	{
		//find the one to be selected
		iindex = 0; fmin = -1.0E200;
		for(j=0; j<pop.Size(); j++) if(!exist[j] && dismin[j] > fmin) {fmin=dismin[j]; iindex=j;} 
		exist[iindex] = true; index[i] = iindex;

		//update the minimum distance of other points
		for(j=0; j<pop.Size(); j++) if(!exist[j] && dismin[j] > mMatDis[j][iindex]) dismin[j] = mMatDis[j][iindex];
	}
}

} //namespace gen
} //namespace mea
} //namespace az
