/*! \file	Generator_Model_GTM1.cpp
	
	\brief	Evolutionary Aglorithm Generator with GTM (finding center points)
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Apr.10 2006 redesign
*/
#include <cmath>
#include <fstream>
#include "alg/Matrix.h"
#include "emo/GenMod.h"

namespace az
{
namespace mea
{
namespace gen
{
namespace mod
{

//constractor
ModelGTM1::ModelGTM1()
{
	mbWeight	= true;
}

//destractor
ModelGTM1::~ModelGTM1()
{
}

//initialize the GTM
void ModelGTM1::Set(unsigned int noLatent, unsigned int noBaseFun, unsigned int dimLatent, unsigned int trainsteps, double extension)
{
	mTrainSteps	= trainsteps;
	mExtension	= extension;
	mLatentDim	= dimLatent;
	mNoLatent	= noLatent;
	mNoBaseFun	= noBaseFun;
	if(mLatentDim==1 && mNoLatent % 2 == 0) mNoLatent++;
}
	
//reset the initial weights of GTM
void ModelGTM1::Reset()
{
	mbWeight = true;
}

//build GTM model and sample new solutions
CPopulationMO& ModelGTM1::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j;

	Reset();
	
	//Step 1: clear the return population
	popnew.Resize(sizenew);
	
	//Step 2: assign new data
	if( popref.Size() != mDataSize ) 
	{
		
		mDataSize	= popref.Size();
		mDataDim	= popref.P().XSize();
		mT.Resize(mDataSize,mDataDim);
	}
	for( i=0; i<mDataSize; i++ ) for( j=0; j<mDataDim; j++ ) mT(i,j) = popref[i][j];
	
	//Step 3: reset the weights of GTM if necessary
	if( mbWeight )
	{
		if(mLatentDim==1)
			mGTM.Initialize1(mFI,mW,mBeta,mT,mNoLatent,mNoBaseFun,2.0);
		else
			mGTM.Initialize2(mFI,mW,mBeta,mT,mNoLatent,mNoBaseFun,2.0);
		mbWeight = false;
	}

	//Step 4: train GTM model
	mGTM.Train(mW,mBeta,mT,mFI,mTrainSteps);

	//Step 5: build Principal Curve(Surface) model and sample new solutions
    //1D
	if(mLatentDim==1)
		ModelGen1D(popnew, popref);
	//2D
	else ModelGen2D(popnew, popref);

	return popnew;
}

// the hidden structure is 1D, use polynominal function to estimate new solutions
// each segment has new solutions according to its arc length
// Nov.10 2005
CPopulationMO& ModelGTM1::ModelGen1D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j,nodeNo,neuro,popCur,size;
	double t0,t1,t2,c0,c1,c2,start,end,step,totalLength,radius;

	nodeNo  = mGTM.Project().RowSize();	 //center node number = 2N+1

	std::vector< double >	t,							//latent variable when sampling
							arcLength(nodeNo-1);		//arc length
			
	// segements model
	//  o----o----o----o----o----o

	std::ofstream file("center.txt");
	file<<nodeNo<<std::endl;
	for(i=0; i<nodeNo; i++) file<<mGTM.Project()(i,0)<<"\t"<<mGTM.Project()(i,1)<<std::endl;
	file.close();

	//Step2: calculate the arc length
	totalLength=0.0;
	for(neuro=0; neuro<nodeNo-1; neuro++)
	{
		arcLength[neuro] = DisCC(neuro, neuro+1);
		totalLength += arcLength[neuro];
	}

	//Step3: build quadratic models with Gaussian noise
	//Gaussian noise variance
	radius  = 1.0 / sqrt(mBeta);// * PERTUBATION();
	//X=C0 + C1*t + C2*t*t + N(0,r)
	popCur = 0;
	for( neuro=0; neuro<nodeNo-1; neuro+=2 )
	{
		//reference independent variables
		//  o----o----o
		// t0    t1   t2
		t0 = 0.0;
		t1 = arcLength[neuro];
		t2 = arcLength[neuro+1] + t1;
		// at the begining(extension) 
		//  ---o----o----o
		//     t0    t1   t2
		//start = (neuro==0) ? -t2*mExtension : 0.0;			
		start = -t2*mExtension;

		// at the end(extension)
		//  o----o----o---
		// t0    t1   t2		
		//end   = ((neuro+3)==nodeNo) ? t2*(1.0+mExtension) : t2;
		end   = t2*(1.0+mExtension);

		//offspring size 
		//         arc length 
		// size = ------------ * popsize
		//        total length
		size  = ((neuro+3)==nodeNo) ? (popnew.Size() - popCur) : (unsigned int)(popref.Size()*t2/totalLength);
		
		if(size<1) continue;

		// offspring independent variables
		// uniformally selected variables in [start, end]
		t.resize( size );
		step  = (end-start)/double(size);
		for(i=0; i<size; i++) { t[i] = rnd::rand(start, start+step); start += step; }

		//sample: interplate
		for(i=0; i<mDataDim; i++)
		{
			//c0 + c1*t0 + c2*t0*t0 = x0 => c0 = x0
			//c0 + c1*t1 + c2*t1*t1 = x1 => c1 = (x1-x0-c2*t1*t1)/t1
			//c0 + c1*t2 + c2*t2*t2 = x2 => c2 = ((x1-x0)t2 - (x2-x0)t1)/(t1t1t2-t1t2t2)
			c0 = mGTM.Project()(neuro,i);
			c2 = (mGTM.Project()(neuro+1,i)-c0)*t2 - (mGTM.Project()(neuro+2,i)-c0)*t1;
			c2/= t1*t2*(t1-t2);
			c1 = (mGTM.Project()(neuro+1,i)-c0)/t1 - c2*t1;

			//x=c0 + c1*t + c2*t*t + N(0,r)
			for(j=0; j<size; j++)
				popnew[popCur+j][i] = c0 + c1*t[j] + c2*t[j]*t[j] + radius*rnd::gaussian();
		}

		popCur += size;
	}
	
	return popnew;
}

// the hidden structure is 2D, use triangles to estimate new solutions
// new solutions are created by area ratio
// Nov.10 2005
CPopulationMO& ModelGTM1::ModelGen2D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int nodeNo, gridWidth, triNo, triNo1, triNo2, i, j, k, s, nodeA, nodeB, nodeC, nodeD;
	double totalArea, area, t1, t2, radius;
	nodeNo		= mGTM.Project().RowSize();					// center points number
	gridWidth	= (unsigned int)(sqrt(double(nodeNo+0.01)));	// the grid is square 
	triNo		= (gridWidth-1)*(gridWidth-1)*2;			// triangle number in the mesh grid
	std::vector<double> triArea(triNo);						// all triangle areas(used to selection)
	std::vector< std::vector<double> > sideLength(nodeNo);	// all the side lengthes of triangles

	for(i=0; i<nodeNo; i++)
	{
		sideLength[i].resize(nodeNo);
		for(j=0; j<nodeNo; j++) sideLength[i][j] = 0.0;
	}
	
	// triangle mesh grid
	//  _ _ _
	// |/|/|/|
	//  - - -
	// |/|/|/|
	//  - - -
	// |/|/|/|
	//  - - -

	//Step1: calculate the side lengthes
	for(i=0; i<gridWidth-1; i++)
	{
		for(j=0; j<gridWidth-1; j++)
		{
			// A B
			//  _ 
			// |/|
			//  - 
			// C D
			nodeA=i*gridWidth+j;	 nodeB=nodeA+1;
			nodeC=(i+1)*gridWidth+j; nodeD=nodeC+1;
			//AB
			sideLength[nodeA][nodeB] = sideLength[nodeB][nodeA] = DisCC(nodeA, nodeB);
			//AC
			sideLength[nodeA][nodeC] = sideLength[nodeC][nodeA] = DisCC(nodeA, nodeC);
			//BC
			sideLength[nodeC][nodeB] = sideLength[nodeB][nodeC] = DisCC(nodeC, nodeB);
			//BD
			if(j == gridWidth-2)
				sideLength[nodeD][nodeB] = sideLength[nodeB][nodeD] = DisCC(nodeD, nodeB);
			//CD
			if(i == gridWidth-2)
				sideLength[nodeC][nodeD] = sideLength[nodeD][nodeC] = DisCC(nodeC, nodeD);
		}
	}

	//Step2: calculate triangle areas
	totalArea = 0.0;
	for(i=0; i<gridWidth-1; i++)
	{
		for(j=0; j<gridWidth-1; j++)
		{
			// A B
			//  _
			// |/|
			//  -
			// C D
			nodeA=i*gridWidth+j;	 nodeB=nodeA+1;
			nodeC=(i+1)*gridWidth+j; nodeD=nodeC+1;
			triNo1 = i*(gridWidth-1)*2 + j*2;
			triNo2 = triNo1+1;
			//ABC
			totalArea += 0.25*sqrt( ( sideLength[nodeA][nodeB]+sideLength[nodeA][nodeC]+sideLength[nodeB][nodeC])*
									( sideLength[nodeA][nodeB]+sideLength[nodeA][nodeC]-sideLength[nodeB][nodeC])*
									( sideLength[nodeA][nodeB]-sideLength[nodeA][nodeC]+sideLength[nodeB][nodeC])*
									(-sideLength[nodeA][nodeB]+sideLength[nodeA][nodeC]+sideLength[nodeB][nodeC]));
			triArea[triNo1] = totalArea;
			 //BCD
			totalArea += 0.25*sqrt( ( sideLength[nodeD][nodeB]+sideLength[nodeD][nodeC]+sideLength[nodeB][nodeC])*
									( sideLength[nodeD][nodeB]+sideLength[nodeD][nodeC]-sideLength[nodeB][nodeC])*
									( sideLength[nodeD][nodeB]-sideLength[nodeD][nodeC]+sideLength[nodeB][nodeC])*
									(-sideLength[nodeD][nodeB]+sideLength[nodeD][nodeC]+sideLength[nodeB][nodeC]));
			triArea[triNo2] = totalArea;
		}
	}

	// Step3: create new solutions
	// Gaussian noise variance
	radius  = 1.0 / sqrt(mBeta);// * PERTUBATION();
	for(s=0; s<popnew.Size(); s++)
	{
		// Step3.1: find a triangle(route wheel by triangle areas)
		area = rnd::rand(0.0, totalArea);
		for(k=0; k<triNo; k++) if(area <= triArea[k]) break;
		
		// Step3.2: calculate the four corners of this cell
		i = k/(2*gridWidth-2); j = (k % (2*gridWidth-2))/2;
		nodeA=i*gridWidth+j;	 nodeB=nodeA+1;
		nodeC=(i+1)*gridWidth+j; nodeD=nodeC+1;

		t1 = rnd::rand(-mExtension, 1.0+mExtension);
		t2 = rnd::rand(-mExtension, (t1<0)? 1.0+mExtension : 1.0+mExtension-t1);

		//if(k%2==0) nodeA = nodeA;
		// A B
		//  _
		// |/
		// C 
		if(k%2==1) nodeA = nodeD;
		//   B
		//  /|
		//  -
		// C D
		
		for(j=0; j<mDataDim; j++)
			popnew[s][j] = (mGTM.Project()(nodeB,j) - mGTM.Project()(nodeA,j))*t1 +
						   (mGTM.Project()(nodeC,j) - mGTM.Project()(nodeA,j))*t2 +
						    mGTM.Project()(nodeA,j) + 
						    radius*rnd::gaussian();
	}

	return popnew;
}

//the distance between two center points
double ModelGTM1::DisCC(unsigned int cp1, unsigned int cp2)
{
	unsigned int i;
	double dis = 0.0;
	for(i=0; i<mDataDim; i++) dis += (mGTM.Project()(cp1,i)-mGTM.Project()(cp2,i))*(mGTM.Project()(cp1,i)-mGTM.Project()(cp2,i));
	return sqrt(dis);
}

//the distance between a center point and a variable
double ModelGTM1::DisCV(unsigned int cp, unsigned int vp)
{
	unsigned int i;
	double dis=0.0;
	for(i=0; i<mDataDim; i++) dis += (mGTM.Project()(cp,i)-mT(vp,i))*(mGTM.Project()(cp,i)-mT(vp,i));
	return sqrt(dis);
}

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az

