#include <fstream>
#include "alg/Random.h"
#include "alg/HCSampler.h"
#include "emo/Parameter.h"
#include "emo/Initialization.h"

#ifdef INI_HYBRID
#include "edal.h"
#endif

namespace az
{
namespace mea
{
namespace ini
{

//-------------------------------------------------------------------------------------
//Random Initialization
CPopulationMO& LHC::Initialize(CPopulationMO& pop, unsigned int size)
{
	unsigned int i,j,dim = pop.P().XSize();
	std::vector< std::vector<double> > xx(dim);
	for(i=0; i<dim; i++) xx[i].resize(size);

	az::alg::LHC(xx, pop.P().XLow(), pop.P().XUpp());

	pop.Resize(size);
	for(i=0; i<size; i++)
		for(j=0; j<pop.P().XSize(); j++)
			pop[i][j] = xx[j][i];

	pop.Evaluate();
	pop.RankSort();
	mEvas = size;
	return pop;
}

//-------------------------------------------------------------------------------------
//Random Initialization
CPopulationMO& Uniform::Initialize(CPopulationMO& pop, unsigned int size)
{
	unsigned int i,j,dim = pop.P().XSize();
	std::vector< std::vector<double> > xx(dim);
	for(i=0; i<dim; i++) xx[i].resize(size);
	
	az::alg::Uniform(xx, pop.P().XLow(), pop.P().XUpp());

	pop.Resize(size);
	for(i=0; i<size; i++)
		for(j=0; j<dim; j++)
			pop[i][j] = xx[j][i];

	pop.Evaluate();
	pop.RankSort();
	mEvas = size;
	return pop;
}

//-------------------------------------------------------------------------------------
#ifdef INI_HYBRID
CIndividualMO		*pInd;
std::vector<double>	vWeight;
double objective(double* xy, int Dimension)
{
	unsigned int i;
	for(i=0; i<Dimension; i++) (*pInd)[i] = xy[i];
	(*pInd).Evaluate();
	double fit = 0.0;
	for(i=0; i<vWeight.size(); i++) fit += (*pInd).F(i)*vWeight[i];
	return fit;
}

Hybrid::Hybrid() {pInd = 0;}

Hybrid::~Hybrid(){if(pInd!=0) delete pInd;}

// initialize a population
CPopulationMO& Hybrid::Initialize(CPopulationMO& pop, unsigned int size)
{
	unsigned int i,j,dim = pop.P().XSize(), nobj = pop.P().FSize();
	std::vector< std::vector<double> > xx(dim);
	for(i=0; i<dim; i++) xx[i].resize(size-nobj);

	az::alg::Uniform(xx, pop.P().XLow(), pop.P().XUpp());

	// random initialization
	pop.Resize(size);
	for(i=0; i<size-nobj; i++)
		for(j=0; j<pop.P().XSize(); j++)
			pop[i+nobj][j] = xx[j][i];

	// biased initialization
	mEvas = 0;
	if(pInd!=0) delete pInd; pInd = new CIndividualMO(pop.P());
	vWeight.resize(nobj);
	double *indbest = new double[dim];
	for(i=0; i<nobj; i++)
	{
		vWeight[i] = 0.98; for(j=0; j<nobj; j++) if(j!=i) vWeight[j] = 0.02/double(nobj-1);
		mEvas += edal::optimizer(objective, 50000, 500, 100, dim, 0, 1, indbest);
		for(j=0; j<dim; j++) pop[i][j] = indbest[j];
	}

	pop.Evaluate();
	pop.RankSort();
	mEvas += size-nobj;
	delete []indbest;
	return pop;
}
#else
CPopulationMO& Hybrid::Initialize(CPopulationMO& pop, unsigned int size)
{
	unsigned int i,j,dim = pop.P().XSize();
	std::vector< std::vector<double> > xx(dim);
	for(i=0; i<dim; i++) xx[i].resize(size);
	
	az::alg::Uniform(xx, pop.P().XLow(), pop.P().XUpp());

	pop.Resize(size);
	for(i=0; i<size; i++)
		for(j=0; j<dim; j++)
			pop[i][j] = xx[j][i];

	pop.Evaluate();
	pop.RankSort();
	mEvas = size;
	return pop;
}
#endif

}//namespace ini
} //namespace mea
} //namespace az
