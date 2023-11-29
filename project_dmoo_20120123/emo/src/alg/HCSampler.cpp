// HCSampler.cpp

#include <algorithm>
#include "alg/HCSampler.h"
#include "alg/Random.h"

namespace az
{

namespace alg
{
	void LHC(std::vector< std::vector<double> >& rand, std::vector<double>& low, std::vector<double>& upp)
	{
		unsigned int i,j;
		for(i=0; i<rand.size(); i++) 
		{
			for(j=0; j<rand[i].size(); j++) rand[i][j] = low[i] + (upp[i]-low[i])*(j+rnd::rand(0.0,1.0))/double(rand[i].size());
			std::random_shuffle(rand[i].begin(), rand[i].end());
		}
	}

	void Uniform(std::vector< std::vector<double> >& rand, std::vector<double>& low, std::vector<double>& upp)
	{
		unsigned int i,j;
		for(i=0; i<rand.size(); i++) 
		{
			for(j=0; j<rand[i].size(); j++) rand[i][j] = rnd::rand(low[i],upp[i]);
		}
	}

} //namespace alg

} //namespace az
