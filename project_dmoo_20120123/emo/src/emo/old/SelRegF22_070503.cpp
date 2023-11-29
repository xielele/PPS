/*! \file	SelRegF22.cpp
	
	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	May.01 2007 create
*/

#include <fstream>
#include <vector>
#include <list>
#include "mea/Sel.h"

namespace az
{
namespace mea
{
namespace sel
{
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//!\brief	sort population, the first 'size' are best ones
	//!\param	pop population
	//!\param	size size of best ones
	//!\return	population
	CPopulationMO& SRegF22::SelectSort(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()<=size) return pop;

		//pop.Read("F.pop");
		pop.Write("F.pop");

		unsigned int i, k, s, iRef = 20;

		std::vector< bool > vExist(pop.Size());
		for(i=0; i<vExist.size(); i++) vExist[i] = true;

		std::vector< std::vector<double> > vRef(iRef);
		for(i=0; i<iRef; i++) vRef[i].resize(pop.P().FSize());

		// find two extreme points
		unsigned int iA, iB;
		iA = iB = 0;
		for(i=1; i<pop.Size(); i++)
		{ 
			if( (0.95*pop[i].F(0)+0.05*pop[i].F(1)) < (0.95*pop[iA].F(0)+0.05*pop[iA].F(1)) )		iA = i;
			else if((0.95*pop[i].F(1)+0.05*pop[i].F(0)) < (0.95*pop[iB].F(1)+0.05*pop[iB].F(0)))	iB = i;
		}
		// set two extreme reference points
		vRef[0][0]		= pop[iA].F(0) - 0.10*fabs(pop[iB].F(0)-pop[iA].F(0));
		vRef[0][1]		= pop[iA].F(1) + 0.05*fabs(pop[iB].F(1)-pop[iA].F(1));
		vRef[iRef-1][0] = pop[iB].F(0) + 0.05*fabs(pop[iB].F(0)-pop[iA].F(0));
		vRef[iRef-1][1] = pop[iB].F(1) - 0.10*fabs(pop[iB].F(1)-pop[iA].F(1));

		for(i=1; i<iRef-1; i++) ChooseRef(vRef, i, pop);
		
		{
			std::ofstream os("F.mod");	
			os<<std::scientific<<std::setprecision(10);
			for(i=0; i<iRef; i++)
			{
				for(k=0; k<vRef[0].size(); k++) os<<vRef[i][k]<<"\t";
				os<<std::endl;
			}
			os.close();
		}

		std::vector<double> vDis(iRef-1); double tD = 0.0;
		std::vector<unsigned int> vNo(iRef-1); unsigned int tN=0;
		for(i=0; i<vDis.size(); i++) 
		{
			vDis[i] = 0.0;
			for(k=0; k<vRef[0].size(); k++) vDis[i] += (vRef[i+1][k]-vRef[i][k])*(vRef[i+1][k]-vRef[i][k]);
			vDis[i] = sqrt(vDis[i]);
			tD += vDis[i];
		}
		for(i=0; i<vDis.size(); i++)
		{
			vNo[i] = (unsigned int)(size*vDis[i]/tD);
			tN += vNo[i];
		}
		for(; tN<size; tN++) vNo[rnd::rand((unsigned int)0, (unsigned int)vNo.size())]++;

		std::vector<double> C(vRef[0]);
		for(i=0; i<iRef-1; i++)
		{
			for(k=0; k<vNo[i]; k++)
			{
				for(s=0; s<vRef[0].size(); s++) C[s] = vRef[i][s] + double(k+1.0)/double(vNo[i])*(vRef[i+1][s]-vRef[i][s]);

				Select(pop, vExist, vRef[i], vRef[i+1], C);
			}
		}

		//{
		//	CPopulationMO popt(pop.P());
		//	popt.Copy(pop[iA]); popt.Copy(pop[iB]);
		//	popt.Write("popt.pop");
		//	popt.Clear();
		//}

		//// select
		//Select(pop, end, size, vExist, rA, rB);
	
		// sort
		int a=0,b=(int)pop.Size()-1;
		while(a<b)
		{
			while(a<(int)pop.Size()	&& !vExist[a]) a++;
			while(b>=0				&&  vExist[b]) b--;
			if(a<b)
			{
				pop.Swap(a,b);
				a++; b--;
			}
		}

		return pop;
	}

	//!\brief	select from current population and offspring population
	//!\param	pop combined population
	//!\param	size population size
	//!\return	population
	CPopulationMO& SRegF22::Select(CPopulationMO& pop, unsigned int size)
	{
		if(pop.Size()>size) SelectSort(pop, size).Erase(size);

		return pop;
	}

	//!\brief	choose a reference point
	//!\param	vRef	set of reference points
	//!\param	iR		index of reference point
	//!\param	pop		current population
	//!\return	void
	void SRegF22::ChooseRef(std::vector< std::vector<double> >& vRef, unsigned int iR, CPopulationMO& pop)
	{
		unsigned int i, index = 0; double tgt, tg = 1.0E100;
		for(i=0; i<pop.Size(); i++) if( InSubRange(vRef, iR-1, pop[i]) && pop[i].F(0) > vRef[iR-1][0])
		{
			tgt = (pop[i].F(1) - vRef[iR-1][1])/(pop[i].F(0) - vRef[iR-1][0]);
			if(tgt<tg) {tg=tgt; index=i;}
		}
		std::vector<double> L(vRef[0]), R(vRef[vRef.size()-1]);
		if(iR>1)
		{
			for(i=0; i<L.size(); i++) L[i] = vRef[0][i] + (double)(iR-1.0)/(double)(vRef.size()-1.0)*(vRef[vRef.size()-1][i] - vRef[0][i]);
		}
		if(iR<vRef.size()-1)
		{
			for(i=0; i<R.size(); i++) R[i] = vRef[0][i] + (double)(iR)/(double)(vRef.size()-1.0)*(vRef[vRef.size()-1][i] - vRef[0][i]);
		}

		double PRO=0.0, LR = 0.0;
		for(i=0; i<R.size(); i++) LR += (L[i]-R[i])*(L[i]-R[i]);
		LR	 = sqrt(LR);
		for(i=0; i<R.size(); i++) PRO += (pop[index].F(i)-vRef[iR-1][i])*(R[i]-L[i]);
		PRO /= LR;

		for(i=0; i<R.size(); i++) 
		{
			vRef[iR][i] = vRef[iR-1][i] + LR/PRO*(pop[index].F(i)-vRef[iR-1][i]);
		}
		PRO = 0.0;
		for(i=0; i<R.size(); i++) PRO += (pop[index].F(i)-vRef[iR-1][i])*(vRef[iR][i]-R[i]);
		for(i=0; i<R.size(); i++) 
		{
			vRef[iR][i] += (PRO>0 ? 1.0:-1.0)*0.1*(R[i] - vRef[iR][i]);
		}

	}

	//!\brief	determine whether a point is in a range
	//!\param	vRef	set of reference points
	//!\param	iR		index of reference point
	//!\param	ind		individual
	//!\return	bool
	bool SRegF22::InSubRange(std::vector< std::vector<double> >& vRef, unsigned int iR, CIndividualMO& ind)
	{
		unsigned int i;
		bool left, right;
		if(iR == 0) left = true;
		else
		{
			std::vector<double> C(vRef[0]);
			for(i=0; i<C.size(); i++) C[i] += (double)(iR)/(double)(vRef.size()-1.0)*(vRef[vRef.size()-1][i] - vRef[0][i]);
			double pro = 0.0;
			for(i=0; i<C.size(); i++) pro += (ind.F(i)-C[i])*(vRef[vRef.size()-1][i] - C[i]);
			left = (pro>=0);
		}
		if(iR == vRef.size()-1) right = true;
		else
		{
			std::vector<double> C(vRef[0]);
			for(i=0; i<C.size(); i++) C[i] += (double)(iR+1.0)/(double)(vRef.size()-1.0)*(vRef[vRef.size()-1][i] - vRef[0][i]);
			double pro = 0.0;
			for(i=0; i<C.size(); i++) pro += (ind.F(i)-C[i])*(vRef[0][i] - C[i]);
			right = (pro>=0);
		}

		return left & right;
	}

	//!\brief	select a point to a reference point
	//!\param	pop population
	//!\param	vExist indicates points are selected or not
	//!\param	L reference point
	//!\param	R reference point
	//!\param	C reference point
	//!\return	void
	void SRegF22::Select(CPopulationMO& pop, std::vector<bool>& vExist, std::vector<double>& L, std::vector<double>& R, std::vector<double>& C)
	{
		unsigned int i;

		unsigned int index = 0;
		if(L[1]>R[1])
		{
			//||rB-rA||
			double dAB = sqrt( (L[0]-R[0])*(L[0]-R[0]) + (L[1]-R[1])*(L[1]-R[1]) );
			//max distance to C and index
			double proLR, proO, dis, mDis = -1.0E100; 
			//find the optimal point
			for(i=0; i<pop.Size(); i++) if(vExist[i])
			{
				proLR = (pop[i].F(0)-C[0])*(R[0]-L[0]) + (pop[i].F(1)-C[1])*(R[1]-L[1]); // \ dAB
				proO  = (pop[i].F(0)-C[0])*(R[1]-L[1]) - (pop[i].F(1)-C[1])*(R[0]-L[0]); // \ dAB
				dis   = proO - proLR*proLR/fabs(proO);
				dis  /= dAB;
				if(dis>mDis) {mDis = dis; index=i;}
			}
		}
		else
		{
			double dis, mDis = 1.0E100;
			//find the optimal point
			for(i=0; i<pop.Size(); i++) if(vExist[i])
			{
				dis = (pop[i].F(0)-C[0])*(pop[i].F(0)-C[0]) + (pop[i].F(1)-C[1])*(pop[i].F(1)-C[1]); 
				if(dis<mDis) {mDis = dis; index=i;}
			}
		}
		// add the point to the next generation
		vExist[index] = false;
	}
}//namespace sel
} //namespace mea
} //namespace az
