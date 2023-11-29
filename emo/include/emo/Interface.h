/*! \file	Interface.h
	
	\brief	friend interfaces to use generators and selectors in the package

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	Aug.22 2007 create
*/

#ifndef AZ_INTERFACE_H
#define AZ_INTERFACE_H

#include "emo/Gen.h"
#include "emo/GenMod.h"
#include "emo/Sel.h"
#include "emo/Parameter.h"

//!\brief	interface to call offspring generator
//!\param	index	the ID of offspring generator
//!				- 1 Local PCA based generator (see Q. Zhang, A. Zhou and Y. Jin, RM-MEDA: A Regularity Model Based Multiobjective Estimation of Distribution Algorithm, IEEE Trans. on Evol. Comp., 2007)
//!				- 2 GTM based generator (see A. Zhou, Q. Zhang, Y. Jin, B. Sendhoff and E. Tsang, Modelling the Population Distribution in Multi-objective Optimization by Generative Topographic Mapping, PPSN 2006)
//!				- 3 SBX (real coded) (see K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, A Fast and Elitist Multiobjective Genetic Algorithm:NSGA-II, IEEE Trans. on Evol. Comp., 2002)
//!\param	popf	pointer to fitness values of parent population, ndat x nobj
//!\param	popx	pointer to decision variables of parent population, ndat x ndec
//!\param	nobj	number of objectives
//!\param	ndec	number of decision variables
//!\param	ndat	number of data used (size of parent population)
//!\param	xupp	upper boundary of decision variables
//!\param	xlow	lower boundary of decision variables
//!\param	par		pointer ot parameters
//!\param	npar	number of parameters
//!\param	newx	offspring population (!!! IMPORTANT: IT MUST BE CREATED BEFORE THE CALLING), nnew x ndec
//!\param	nnew	number of new trial solutions to generate
//!\return	void (newx)
void Generator(unsigned int index,
			   double** popf, double** popx, unsigned int nobj, unsigned int ndec, unsigned int ndat, 
			   double*  xupp, double*  xlow,
			   double*  par,  unsigned int npar,
			   double** newx, unsigned int nnew)
{
	az::rnd::seed((long) time(NULL));

	unsigned int i,j;
	az::mea::CParameter mPar;
	mPar.FSize(nobj);
	mPar.XSize(ndec);
	for(i=0; i<ndec; i++)
	{
		mPar.XUpp(i) = xupp[i];
		mPar.XLow(i) = xlow[i];
	}
	az::mea::CPopulationMO pop(mPar), popnew(mPar);
	pop.Resize(ndat);
	for(i=0; i<ndat; i++)
	{
		for(j=0; j<nobj; j++) pop[i].F(j) = popf[i][j];
		for(j=0; j<ndec; j++) pop[i][j]   = popx[i][j];
	}

	// generate new trial solutions
	switch(index)
	{
	// RM-MEDA
	case 1:
		{
			//!!!IMPORTANT: PARAMETERS
			unsigned int latdim = nobj-1;					// dimension of latent variables
			unsigned int noclus = (unsigned int)par[0];		// number of clusters used in Local PCA algorithm
															// make sure there are at least 10-20 points in each cluster
			unsigned int trains = (unsigned int)par[1];		// Local PCA algorithm training steps
			double       extens = par[2];					// extension ratio in RM-MEDA sampling process
															// 0.1-0.25 for 2-obj problems
			az::mea::gen::mod::RM gen;
			gen.Set(latdim, noclus, trains, extens);
			gen.Generate(nnew, popnew, pop);
		}
		break;
	// GTM based method
	case 2:
		{
			//!!!IMPORTANT: PARAMETERS
			unsigned int latdim = nobj-1;					// dimension of latent variables
			unsigned int latvar	= (unsigned int)par[0];		// the number latent variables
			unsigned int basisf = (unsigned int)par[1];		// the number basis functions
			unsigned int trains = (unsigned int)par[2];		// GTM algorithm training steps
			double       extens = par[3];					// extension ratio in RM-MEDA sampling process
															// 0.1-0.25 for 2-obj problems
			az::mea::gen::mod::ModelGTM2 gen;
			gen.Set(latvar, basisf, latdim, trains, extens);
			gen.Generate(nnew, popnew, pop);
		}
		break;
	// SBX
	case 3:
		{
			//!!!IMPORTANT: PARAMETERS
			mPar.Pc()		= par[0];		// crossover probability
			mPar.Pm()		= par[1];		// mutation probability
			mPar.ETA_SBX()	= par[2];		// eta for SBX(simulated binary crossover)
			mPar.ETA_PM()	= par[3];		// eta for PM(polynomial mutation)
			
			az::mea::gen::XSBX gen;
			gen.Generate(nnew, popnew, pop);
		}
		break;
	default:
		break;
	}

	for(i=0; i<nnew; i++)
	{
		for(j=0; j<ndec; j++) newx[i][j] = popnew[i][j];
	}
}

//!\brief	interface to call selection operator
//!\param	index	the ID of selection operator
//!				- 1 modified nondominated crowding sortring (see Q. Zhang, A. Zhou and Y. Jin, RM-MEDA: A Regularity Model Based Multiobjective Estimation of Distribution Algorithm, IEEE Trans. on Evol. Comp., 2007
//!                   and K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, A Fast and Elitist Multiobjective Genetic Algorithm:NSGA-II, IEEE Trans. on Evol. Comp., 2002)
//!				- 2 MaxiMin selection (see E.J. Solteiro Pires, P.B. de Moura Oliveira and J.A. Tenreiro Machado, Multi-objective MaxiMin Sorting Scheme, EMO 2005)
//!\param	popf1	pointer to fitness values of parent population, ndat x nobj
//!\param	popx1	pointer to decision variables of parent population, ndat x ndec
//!\param	popf2	pointer to fitness values of parent population, ndat x nobj
//!\param	popx2	pointer to decision variables of parent population, ndat x ndec
//!\param	nobj	number of objectives
//!\param	ndec	number of decision variables
//!\param	ndat1	number of data used (size of parent population)
//!\param	ndat2	number of data used (size of parent population)
//!\param	newf	offspring population (!!! IMPORTANT: IT MUST BE CREATED BEFORE THE CALLING), nnew x nobj
//!\param	newx	offspring population (!!! IMPORTANT: IT MUST BE CREATED BEFORE THE CALLING), nnew x ndec
//!\param	nnew	number of new trial solutions to generate
//!\return	void (newx)
void Selector(unsigned int index,
			   double** popf1, double** popx1, double** popf2, double** popx2,
			   unsigned int nobj, unsigned int ndec, unsigned int ndat1, unsigned int ndat1, 
			   double** newf, double** newx, unsigned int nnew)
{
	az::rnd::seed((long) time(NULL));

	unsigned int i,j;
	az::mea::CParameter mPar;
	mPar.FSize(nobj);
	mPar.XSize(ndec);
	az::mea::CPopulationMO pop1(mPar), pop2(mPar), popnew(mPar);
	pop1.Resize(ndat1);pop2.Resize(ndat2);
	for(i=0; i<ndat1; i++)
	{
		for(j=0; j<nobj; j++) pop1[i].F(j) = popf1[i][j];
		for(j=0; j<ndec; j++) pop1[i][j]   = popx1[i][j];
	}
	for(i=0; i<ndat2; i++)
	{
		for(j=0; j<nobj; j++) pop2[i].F(j) = popf2[i][j];
		for(j=0; j<ndec; j++) pop2[i][j]   = popx2[i][j];
	}

	popnew.Combine(pop1);
	popnew.Combine(pop2);
	
	// select next generation
	switch(index)
	{
	// nondominated crowding distance sorting
	case 1:
		{
			az::mea::sel::SCrowd2 sel;
			sel.Select(popnew, nnew);
		}
		break;
	// MaxiMin selection
	case 2:
		{	
			az::mea::sel::SMaxMin sel;
			sel.Select(popnew, nnew);
		}
		break;
	default:
		break;
	}

	for(i=0; i<nnew; i++)
	{
		for(j=0; j<nobj; j++) newf[i][j] = popnew[i].F(j);
		for(j=0; j<ndec; j++) newx[i][j] = popnew[i][j];
	}
}

#endif //AZ_INTERFACE_H