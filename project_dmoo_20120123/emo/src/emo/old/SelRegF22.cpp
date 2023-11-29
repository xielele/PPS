/*! \file	SelRegF22.cpp
	
	\brief	MOEA selection strategies

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	May.01 2007 create
	\date	May.03 2007 redesign
	\date	May.04 2007 redesign
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

		unsigned int i, k, s;
		
		unsigned int iDim = pop.P().FSize();					// dimension of the space (2: default)
		unsigned iSeg = 8;										// iSeg segments
		std::vector< std::vector<double> >	vRef(iSeg+1),		// reference points
											vSegL(iSeg),		// left point of each segment
											vSegR(iSeg);		// right point of each segment
		for(i=0; i<iSeg; i++) { vRef[i].resize(iDim); vSegL[i].resize(iDim); vSegR[i].resize(iDim);}
		vRef[iSeg].resize(iDim);
		std::vector< double >	vLenS(iSeg),						// length of each segment
								vLR(iDim),							// direction of LR
								vOLR(iDim);							// orthogonal direction to LR
		double LenS = 0.0;											// length of PF
		double LenR = 0.0;											// length of reference segment		
		std::vector< unsigned int > vInS(pop.Size()),				// which segment of each point belongs to
									vNoS(iSeg);						// how many points in each segment
		std::vector< bool > vExist(pop.Size());						// select or not
		for(i=0; i<vExist.size(); i++) vExist[i] = true;

		// Step 1: set reference points
		// find two extreme points
		unsigned int iA, iB;
		iA = iB = 0;
		for(i=1; i<pop.Size(); i++)
		{ 
			if( (0.95*pop[i].F(0)+0.05*pop[i].F(1)) < (0.95*pop[iA].F(0)+0.05*pop[iA].F(1)) )		iA = i;
			else if((0.95*pop[i].F(1)+0.05*pop[i].F(0)) < (0.95*pop[iB].F(1)+0.05*pop[iB].F(0)))	iB = i;
		}
		// set two extreme reference points
		// !!! here we use the default dimension 2
		vRef[0][0]		= pop[iA].F(0) - 0.05*fabs(pop[iB].F(0)-pop[iA].F(0));
		vRef[0][1]		= pop[iA].F(1) + 0.05*fabs(pop[iB].F(1)-pop[iA].F(1));
		vRef[iSeg][0]	= pop[iB].F(0) + 0.05*fabs(pop[iB].F(0)-pop[iA].F(0));
		vRef[iSeg][1]	= pop[iB].F(1) - 0.05*fabs(pop[iB].F(1)-pop[iA].F(1));
		// set the other reference points
		for(i=1; i<iSeg; i++) {for(k=0; k<iDim; k++) vRef[i][k] = vRef[0][k] + double(i)/double(iSeg)*(vRef[iSeg][k]-vRef[0][k]);}

		// Step 2: set directions LR and OLR
		//    /^
		// L /
		//   \
		//	  \> R	
		for(k=0; k<iDim; k++) {vLR[k] = vRef[iSeg][k]-vRef[0][k]; LenR += vLR[k]*vLR[k];}
		LenR = sqrt(LenR);
		for(k=0; k<iDim; k++) {vLR[k] /= LenR;}
		LenR/= double(iSeg);
		// !!! here we use the default dimension 2
		vOLR[0] = -vLR[1]; vOLR[1] = vLR[0];

		//Step 3: assign each point to the corrosponding subset (segment)
		for(i=0; i<vNoS.size(); i++) vNoS[i] = 0;
		for(i=0; i<pop.Size(); i++)
		{
			for(k=0; k<iSeg; k++) if(InSegment(k, pop[i], vRef, vLR)) break;
			vInS[i] = k;
			vNoS[k]++;
		}

		//Step 4: segment
		std::vector<double> vProX(pop.Size()),vProY(pop.Size()); double mProY = 1.0E100;
		for(i=0; i<pop.Size(); i++)
		{
			vProX[i] = vProY[i] = 0.0;
			for(k=0; k<iDim; k++) {vProX[i] += (pop[i].F(k)-vRef[0][k])*vLR[k]; vProY[i] += (pop[i].F(k)-vRef[0][k])*vOLR[k];}
			if(vProY[i]<mProY) mProY = vProY[i];
		}
		for(k=0; k<iSeg; k++)
		{
			// no point in the segment
			//if(vNoS[k]<1) {vSegL[k]=vRef[k]; vSegR[k]=vRef[k+1]; vLenS[k]=0.0; continue;}

			//Strategy 1 & 2
			//unsigned int iA = pop.Size(), iB = pop.Size(); double dA=1.0E100, dB=1.0E100,d;
			//for(i=0; i<pop.Size(); i++) if(vInS[i]==k)
			//{
			//	//if(iA == pop.Size()) iA = i; 
			//	//else if(vProY[i]<vProY[iA]) {iB = iA; iA = i;}
			//	//else if(iB == pop.Size() || vProY[i]<vProY[iB]) iB = i;

			//	d = vProY[i] - (vProX[i] - LenR*k)*(vProX[i] - LenR*k)/fabs(vProY[i]);
			//	if(d<dA) {dA=d; iA=i;}
			//	d = vProY[i] - (vProX[i] - LenR*(k+1))*(vProX[i] - LenR*(k+1))/fabs(vProY[i]);
			//	if(d<dB) {dB=d; iB=i;}
			//}

			//double rL=0.0, rR=0.0, t;
			//// only one point in the segment
			//if(vNoS[k] == 1) 
			//{
			//	t=0.0;
			//}
			//// more than one point in the segment
			//else
			//{
			//	if(fabs(vProX[iA]-vProX[iB]) < (vProY[iB]-vProY[iA])) t = (vProX[iB]>vProX[iA]) ? 1.0:-1.0;
			//	else t = (vProY[iB]-vProY[iA])/(vProX[iB]-vProX[iA]);
			//}
			//rR = (LenR*(k+1)-vProX[iA])*t + vProY[iA];
			//rL = (LenR*(k+0)-vProX[iA])*t + vProY[iA];
			//for(i=0; i<iDim; i++) {vSegL[k][i]=vRef[k][i]+rL*vOLR[i]; vSegR[k][i]=vRef[k+1][i]+rR*vOLR[i];}

			//Strategy 3: using reference points
			for(i=0; i<iDim; i++){vSegL[k][i]=vRef[k][i] + (mProY-LenR)*vOLR[i]; vSegR[k][i]=vRef[k+1][i] + (mProY-LenR)*vOLR[i];}

			//Strategy 4:
			//if(k==0)		vSegL[k] = vRef[0];
			//else			vSegL[k] = vSegR[k-1];
			//if(k==iSeg-1)	vSegR[k] = vRef[iSeg];
			//else
			//{
			//	double proy = 0.0, t = 0.0; unsigned int index = pop.Size();
			//	for(i=0; i<iDim; i++) proy += (vSegL[k][i]-vRef[0][i])*vOLR[i];
			//	for(i=0; i<pop.Size(); i++) if( (vProX[i]> k*LenR) && (vProX[i]<= (k+1.0)*LenR) )
			//	{
			//		if(vProX[i]> k*LenR + fabs(vProY[i]-proy))
			//		{
			//			if(index==pop.Size() || (vProY[i]-proy)/(vProX[i]-k*LenR) < t ) { index = i; t = (vProY[index]-proy)/(vProX[index]-k*LenR); }
			//		}
			//		else
			//		{
			//			if(index==pop.Size()) t = vProY[i] > proy ? 1:-1; 
			//		}
			//	}
			//	if(fabs(t)>1) t = t>0 ? 1:-1;
			//	
			//	for(i=0; i<iDim; i++) vSegR[k][i] = vRef[k+1][i] + (LenR*t+proy)*vOLR[i];
			//}

			vLenS[k]=0.0;
			for(i=0; i<iDim; i++) vLenS[k] += (vSegL[k][i]-vSegR[k][i])*(vSegL[k][i]-vSegR[k][i]);
			vLenS[k] = sqrt(vLenS[k]);
			LenS += vLenS[k];
		}

		//Step 5: assign the size of reference points to each segment
		unsigned int tN=0;
		for(k=0; k<iSeg; k++)
		{
			vNoS[k] = (unsigned int)(size*vLenS[k]/LenS);
			tN     += vNoS[k];
		}
		while(tN<size){ i = rnd::rand((unsigned int)0, iSeg); if(vNoS[i]>0) vNoS[i]++; tN++;}	

		//Step 6: select point 
		std::vector<double> C(iDim);
		for(k=0; k<iSeg; k++)
		{
			for(i=0; i<vNoS[k]; i++)
			{
				if(vNoS[k]>1){for(s=0; s<iDim; s++) C[s] = vSegL[k][s] + double(i)/double(vNoS[k]-1.0)*(vSegR[k][s]-vSegL[k][s]);}
				else		 {for(s=0; s<iDim; s++) C[s] = vSegL[k][s] + 0.5*(vSegR[k][s]-vSegL[k][s]);}
				Select(pop, vExist, vSegL[k], vSegR[k], C);
			}
		}
		
		{
			std::ofstream os("F.mod");	
			os<<std::scientific<<std::setprecision(10);
			os<<vRef[0][0]<<"\t"<<vRef[0][1]<<"\t"<<vRef[iSeg][0]<<"\t"<<vRef[iSeg][1]<<std::endl;
			for(i=0; i<iSeg; i++)
			{
				for(k=0; k<iDim; k++) os<<vSegL[i][k]<<"\t";
				for(k=0; k<iDim; k++) os<<vSegR[i][k]<<"\t";
				os<<std::endl;
			}
			os.close();
		}
	
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

	//!\brief	determine whether a point is in a range
	//!\param	iR		index of segment
	//!\param	ind		individual
	//!\param	vRef	set of reference points
	//!\param	vLR		direction from left to right
	//!\return	bool
	bool SRegF22::InSegment(unsigned int iR, CIndividualMO& ind, std::vector< std::vector<double> >& vRef, std::vector<double>& vLR)
	{
		unsigned int i;
		bool left, right;
		if(iR == 0) left = true;
		else
		{
			double pro = 0.0;
			for(i=0; i<vLR.size(); i++) pro += (ind.F(i)-vRef[iR][i])*vLR[i];
			left = (pro>=0);
		}

		if(!left) return false;

		if(iR == vRef.size()-2) right = true;
		else
		{
			double pro = 0.0;
			for(i=0; i<vLR.size(); i++) pro += (ind.F(i)-vRef[iR+1][i])*vLR[i];
			right = (pro<0);
		}

		return right;
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

		// add the point to the next generation
		vExist[index] = false;
	}
}//namespace sel
} //namespace mea
} //namespace az
