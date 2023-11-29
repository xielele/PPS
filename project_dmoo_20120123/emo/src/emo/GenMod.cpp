// Generator_Model.cpp

#include "emo/GenMod.h"

namespace az
{
namespace mea
{
namespace gen
{
namespace mod
{

// crossover
CPopulationMO& GuidedXOver::XOver(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j,ne;
	unsigned int start,end;
	double beta;

	//find the best solutions
	popref.RankSize(1,start,end);

	//maybe we should check whether the Guided Crossover is needed
	if(end-start>0.2*popref.Size()) return popnew;

	CPopulationMO::IND_TYPE ind1(popref.P()),ind2(popref.P());

	for(i=0; i<(unsigned int)popnew.Size(); i++)
	{
		////make individual leagel(needed by SBX operator)
		//popnew[i].Check();	

		//find a near good point in neighborhood
		//ne = Neighbor(popnew[i], start, end, popref);
		ne		= rnd::rand(start,end);

		beta	= rnd::rand(0.0,1.0); 
		for(j=0; j<popnew.P().XSize(); j++) 
		{	
			popnew[i][j] = popnew[i][j] + beta*(popref[ne][j] - popnew[i][j]) ;
			// border check strategy 1: DE strategy
			if(popnew[i][j] < popnew.P().XLow(j))		popnew[i][j] = 0.5*(popnew.P().XLow(j)+popref[i][j]);
			else if(popnew[i][j] > popnew.P().XUpp(j))	popnew[i][j] = 0.5*(popnew.P().XUpp(j)+popref[i][j]);
		}

		PM(popnew[i]);


		////===============================================================================
		////SBX crossover
		////create two offsprings by SBX
		//mea::gen::SBX(ind1,ind2,popref[ne],popnew[i]);
		////randomly replace i-th offspring by ind1 or ind2
		//if( rnd::rand(0.0,1.0)<0.5 )
		//	popnew[i] = ind2;
		//else 
		//	popnew[i] = ind1;
		////===============================================================================
		
		////===============================================================================
		////uniform crossover
		//unsigned int j;
		//for(j=0; j<popnew.P().XSize(); j++) if(rnd::rand()<0.5) popnew[i][j] = popref[ne][j];
		////===============================================================================


		////===============================================================================
		////uniform crossover 2
		//unsigned int j;
		//for(j=0; j<popnew.P().XSize(); j++) popnew[i][j] = rnd::rand(std::min<double>(popnew[i][j],popref[ne][j]),std::max<double>(popnew[i][j],popref[ne][j]));
		////===============================================================================

		////===============================================================================
		////uniform crossover 3
		//if(rnd::rand()<0.5)
		//{
		//	//make individual leagel(needed by SBX operator)
		//	popnew[i].Check();	
		//	//find a near good point in neighborhood
		//	ne = Neighbor(popnew[i], start, end, popref);
		//	unsigned int j;
		//	for(j=0; j<popnew.P().XSize(); j++) popnew[i][j] += rnd::gaussian()*0.5*(popref[ne][j]-popnew[i][j]);
		//}
		////===============================================================================
	}
	return popnew;
}	

//find the nearest neighbor in reference population
unsigned int GuidedXOver::Neighbor(CPopulationMO::IND_TYPE& ind, unsigned int start, unsigned int end, CPopulationMO& popref)
{
	return rnd::rand(start,end);

	//unsigned int i,j,index1=0,index2=0;
	//double dis1=1.0E100,dis2=1.0E100,dist;

	//for( i=start; i<end; i++ )
	//{
	//	dist = 0.0;
	//	for(j=0; j<(unsigned int)popref.P().XSize(); j++) dist += (ind[j] - popref[i][j])*(ind[j] - popref[i][j]);
	//	if(dist < dis1) 
	//	{
	//		index2 = index1; dis2=dis1;
	//		index1 = i;		 dis1=dist;
	//	}
	//	else if(dist < dis2) 
	//	{
	//		index2 = i; dis2=dist;
	//	}
	//}
	//return (rnd::rand(0.0,1.0)<0.5) ? index1:index2;

	unsigned int i,j,index=0,r1,r2;
	double mindis=1.0E100,dist;

	if(end > start+2) 
	{
		r1=rnd::rand(start,end); 
		do{r2=rnd::rand(start,end);}while(r2==r1); 
		if(r1>r2) std::swap(r1,r2);
	}
	else 
	{
		r1=start; r2=end;
	}

	for(i=r1; i<r2; i++)
	{
		dist = 0.0;
		for(j=0; j<(unsigned int)popref.P().XSize(); j++) dist += (ind[j] - popref[i][j])*(ind[j] - popref[i][j]);
		if(dist < mindis) 
		{
			index = i; mindis = dist;
		}
	}
	return index;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//constructor
ModelBase::ModelBase()
{
	mDataSize	= mDataDim	= 0;
	pData		= 0;
	mExtension	= 0.2;
}

//destructor
ModelBase::~ModelBase()
{
	Clear();
}
//clear data pool
void ModelBase::Clear()
{
	unsigned int i;
	if( pData != 0 )
	{
		for( i=0; i<mDataSize; i++ ) delete []pData[i];
		delete []pData;
		pData	  = 0;
		mDataSize = 0;
	}
}

// convert a real number to a nearest positive integer number 
unsigned int ModelBase::ROUND(double X) 
{ 
	if(X-(unsigned int)(X)<0.5 && (unsigned int)(X)>0 ) 
		return (unsigned int)(X); 
	else 
		return (unsigned int)(X)+1; 
}

// pertubation
double ModelBase::PERTUBATION()	
{ 
	double r;
	do{ r = 1.5 + rnd::gaussian()*0.5; }while(r<0.8 || r>2.2);
	return r;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az
