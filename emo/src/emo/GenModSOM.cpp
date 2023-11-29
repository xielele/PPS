/*! \file	Generator_Model_SOM.cpp
	
	\brief	Evolutionary Aglorithm Generator with SOM
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Apr.10 2006 redesign
*/

#include <cmath>
#include <fstream>
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
ModelSOM::ModelSOM()
{
	mbWeight	= true;
}

//destractor
ModelSOM::~ModelSOM()
{
}

//initialize the SOM
void ModelSOM::Set( unsigned int row, unsigned int col, unsigned int trainsteps, double extension )
{
	mRowGrid	= row;
	mColGrid	= col;
	mTrainSteps	= trainsteps;
	mExtension	= extension;
	if(mRowGrid==1 && mColGrid%2==0) mColGrid++;
	if(mColGrid==1 && mRowGrid%2==0) mRowGrid++;
}
	
//reset the initial weights of SOM
void ModelSOM::Reset()
{
	mbWeight = true;
}

//build SOM model and sample new solutions
CPopulationMO& ModelSOM::Generate(unsigned int sizenew, CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j;

	Reset();
	
	//Step 1: clear the return population
	popnew.Resize(sizenew);
	
	//Step 2: check the reference population size
	if( (unsigned int)(popref.Size()) < mRowGrid*mColGrid ) return popnew;
	
	//Step 3: assign new data
	if( popref.Size() != mDataSize ) 
	{
		Clear();	
		mDataSize	= popref.Size();
		mDataDim	= popref.P().XSize();
		pData		= new double*[mDataSize];
		for( i=0; i<mDataSize; i++ ) pData[i] = new double[mDataDim];
	}
	for( i=0; i<mDataSize; i++ ) for( j=0; j<mDataDim; j++ ) pData[i][j] = popref[i][j];
	
	//Step 4: reset the weights of SOM if necessary
	if( mbWeight )
	{
		mSom.Initialize( mRowGrid, mColGrid, mDataDim );
		mbWeight = false;
	}

	//Step 5: train SOM model
	mSom.Train( mTrainSteps, mDataSize, pData );

	//Step 6: build Principal Curve(Surface) model and sample new solutions
    //1D
	if( mRowGrid == 1 || mColGrid == 1 )
		Generate1D(popnew, popref);
	//2D
	else Generate2D(popnew, popref);

	//Step 7: guided crossover if necessary
	//GuidedXOver::XOver(popnew, popref);

	return popnew;
}

// the hidden structure is 1D, use polynominal function to estimate new solutions
// each segment has new solutions according to its arc length
// Nov.11 2005
CPopulationMO& ModelSOM::Generate1D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j,nodeNo,neuro,popCur,size;
	double t0,t1,t2,c0,c1,c2,start,end,step,totalLength, radius;

	nodeNo  = mRowGrid*mColGrid;						//center node number = 2N+1

	std::vector< double >	t,							//latent variable when sampling
							arcLength(nodeNo-1);		//arc length
			
	// segements model
	//  o----o----o----o----o----o

	std::ofstream file("center.txt");
	file<<nodeNo<<std::endl;
	for(i=0; i<nodeNo; i++) file<<mSom.Weight()[i][0]<<"\t"<<mSom.Weight()[i][1]<<std::endl;
	file.close();

	//Step1: calculate the arc length
	totalLength=0.0;
	for(neuro=0; neuro<nodeNo-1; neuro++)
	{
		arcLength[neuro] = DisCC(neuro, neuro+1);
		totalLength += arcLength[neuro];
	}

	//Step2: gaussian noise variance
	radius = 0.0;
	for(i=0; i<mDataSize; i++)  radius += mSom.DisToCore(i);
	radius  = radius / double(mDataSize*sqrt(double(mDataDim)));// * PERTUBATION();

	//Step3: sample
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
		start = -t2*mExtension;			
		// at the end(extension)
		//  o----o----o---
		// t0    t1   t2		
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
			c0 = mSom.Weight()[neuro][i];
			c2 = (mSom.Weight()[neuro+1][i]-c0)*t2 - (mSom.Weight()[neuro+2][i]-c0)*t1;
			c2/= t1*t2*(t1-t2);
			c1 = (mSom.Weight()[neuro+1][i]-c0)/t1 - c2*t1;

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
CPopulationMO& ModelSOM::Generate2D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int nodeNo, triNo, triNo1, triNo2, i, j, k, s, nodeA, nodeB, nodeC, nodeD;
	double totalArea, area, t1, t2, radius;
	nodeNo		= mRowGrid*mColGrid;						// center points number
	triNo		= (mRowGrid-1)*(mColGrid-1)*2;				// triangle number in the mesh grid
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
	for(i=0; i<mRowGrid-1; i++)
	{
		for(j=0; j<mColGrid-1; j++)
		{
			// A B
			//  _ 
			// |/|
			//  - 
			// C D
			nodeA=i*mColGrid+j;	 nodeB=nodeA+1;
			nodeC=(i+1)*mColGrid+j; nodeD=nodeC+1;
			//AB
			sideLength[nodeA][nodeB] = sideLength[nodeB][nodeA] = DisCC(nodeA, nodeB);
			//AC
			sideLength[nodeA][nodeC] = sideLength[nodeC][nodeA] = DisCC(nodeA, nodeC);
			//BC
			sideLength[nodeC][nodeB] = sideLength[nodeB][nodeC] = DisCC(nodeC, nodeB);
			//BD
			if(j == mColGrid-2)
				sideLength[nodeD][nodeB] = sideLength[nodeB][nodeD] = DisCC(nodeD, nodeB);
			//CD
			if(i == mRowGrid-2)
				sideLength[nodeC][nodeD] = sideLength[nodeD][nodeC] = DisCC(nodeC, nodeD);
		}
	}

	//Step2: calculate triangle areas
	totalArea = 0.0;
	for(i=0; i<mRowGrid-1; i++)
	{
		for(j=0; j<mColGrid-1; j++)
		{
			// A B
			//  _
			// |/|
			//  -
			// C D
			nodeA=i*mColGrid+j;		nodeB=nodeA+1;
			nodeC=(i+1)*mColGrid+j; nodeD=nodeC+1;
			triNo1 = i*(mColGrid-1)*2 + j*2;
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

	//Step3: gaussian noise variance
	radius = 0.0;
	for(i=0; i<mDataSize; i++)  radius += mSom.DisToCore(i);
	radius  = radius / double(mDataSize*sqrt(double(mDataDim)));// * PERTUBATION();

	//Step4: sample
	for(s=0; s<popnew.Size(); s++)
	{
		// Step3.1: find a triangle(route wheel by triangle areas)
		area = rnd::rand(0.0, totalArea);
		for(k=0; k<triNo; k++) if(area <= triArea[k]) break;
		
		// Step3.2: calculate the four corners of this cell
		i = k/(2*mColGrid-2); j = (k % (2*mColGrid-2))/2;
		nodeA=i*mColGrid+j;		nodeB=nodeA+1;
		nodeC=(i+1)*mColGrid+j; nodeD=nodeC+1;

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
			popnew[s][j] =	(mSom.Weight()[nodeB][j] - mSom.Weight()[nodeA][j])*t1 +
							(mSom.Weight()[nodeC][j] - mSom.Weight()[nodeA][j])*t2 +
							 mSom.Weight()[nodeA][j] +
						     radius*rnd::gaussian();
	}

	return popnew;
}

//the distance between two center points
double ModelSOM::DisCC(unsigned int cp1, unsigned int cp2)
{
	unsigned int i;
	double dis = 0.0;
	for(i=0; i<mDataDim; i++) dis += (mSom.Weight()[cp1][i]-mSom.Weight()[cp2][i])*(mSom.Weight()[cp1][i]-mSom.Weight()[cp2][i]);
	return sqrt(dis);
}

//the distance between a center point and a variable
double ModelSOM::DisCV(unsigned int cp, unsigned int vp)
{
	unsigned int i;
	double dis=0.0;
	for(i=0; i<mDataDim; i++) dis += (mSom.Weight()[cp][i]-pData[vp][i])*(mSom.Weight()[cp][i]-pData[vp][i]);
	return sqrt(dis);
}

} //namespace mod
} //namespace gen
} //namespace mea
} //namespace az



/***
//the hidden structure is 1D, use polynominal function to estimate new solutions
CPopulationMO& ModelSOM::Generate1D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int i,j,neuroNo,curveNo,neuro,popCur,size;
	double t0,t1,t2,c0,c1,c2,r,start,end,step;
	
	neuroNo = mRowGrid*mColGrid;			//be sure neuroNo = 2*N+1
	curveNo = unsigned int( neuroNo/2 );	//curve number

	std::vector< unsigned int > notoNeuro( neuroNo ),	//number of points attached to each neruo
								sortNeuro( neuroNo ),	//sorted index of neuro
								nearNeuro( neuroNo );	//the nearest index to each neuro
	std::vector< double >		distoNeuro( neuroNo );	//distance to each neuro
	std::vector< double >		t;

	//Step 1: sort neuros according to their objective 
	for(i=0; i<neuroNo; i++)
	{	
		sortNeuro[i] = i; nearNeuro[i] = neuroNo+1;
	}
	for(i=0; i<mDataSize; i++)
		if(nearNeuro[mSom.NearCore(i)]>neuroNo || mSom.DisToCore( mSom.NearCore(i) ) < mSom.DisToCore( nearNeuro[mSom.NearCore(i)] ) )
			nearNeuro[mSom.NearCore(i)] = i;
	for(i=0; i<neuroNo-1; i++)
		for(j=i+1; j<neuroNo; j++)
			if(popref[nearNeuro[sortNeuro[i]]].F(0)>popref[nearNeuro[sortNeuro[j]]].F(0))
			{
				neuro		 = sortNeuro[j];
				sortNeuro[j] = sortNeuro[i];
				sortNeuro[i] = neuro;
			}

	//Step2: statistic the number of points and the distance attached to each neruo
	for( i=0; i<neuroNo; i++ ) 
	{ 
		notoNeuro[i]  = 0; 
		distoNeuro[i] = 0.0; 
	}
	for( i=0; i<mDataSize; i++ ) 
	{ 
		notoNeuro[mSom.NearCore(i)]++; 
		distoNeuro[mSom.NearCore(i)] += mSom.DisToCore(i); 
	}

	//Step3: build quadratic models with gaussian noise
	//X=C0 + C1*t + C2*t*t + N(0,r)
	popCur = 0;
	for( neuro=0; neuro<neuroNo-1; neuro+=2 )
	{
		//reference independent variables
		t0 = 0.0;
		t1 = NormL2( sortNeuro[neuro+0], sortNeuro[neuro+1] );
		t2 = NormL2( sortNeuro[neuro+1], sortNeuro[neuro+2] ) + t1;
	
		//offspring size
		size  = notoNeuro[sortNeuro[neuro+1]];
		size += (neuro==0) ? notoNeuro[sortNeuro[neuro]] :  notoNeuro[sortNeuro[neuro]]-notoNeuro[sortNeuro[neuro]]/2;
		size += ((neuro+3)==neuroNo) ? notoNeuro[sortNeuro[neuro+2]] :  notoNeuro[sortNeuro[neuro+2]]/2;
		
		//offspring independent variables
		t.resize( size );
		//start = ((neuro==0) ? -mExtension*t2 : 0.0);
		//end   = ((neuro+3)==neuroNo) ? (1.0+mExtension)*t2 : t2;
		start = -mExtension*t2;
		end   = (1.0+mExtension)*t2;
		step  = (end-start)/double(size);
		for(i=0; i<size; i++) { t[i] = rnd::rand(start, start+step); start += step; }

		//gaussian noise variance
		r  = distoNeuro[sortNeuro[neuro]  ]/double(notoNeuro[sortNeuro[neuro  ]]) + 
				distoNeuro[sortNeuro[neuro+1]]/double(notoNeuro[sortNeuro[neuro+1]]) +
				distoNeuro[sortNeuro[neuro+2]]/double(notoNeuro[sortNeuro[neuro+2]]);
		r /= 3.0*sqrt(double(mDataDim));
		r *= PERTUBATION();

		//sample
		for(i=0; i<mDataDim; i++)
		{
			//c0 + c1*t0 + c2*t0*t0 = x0 => c0 = x0
			//c0 + c1*t1 + c2*t1*t1 = x1 => c1 = (x1-x0-c2*t1*t1)/t1
			//c0 + c1*t2 + c2*t2*t2 = x2 => c2 = ((x1-x0)t2 - (x2-x0)t1)/(t1t1t2-t1t2t2)
			c0 = mSom.Weight()[sortNeuro[neuro]][i];
			c2 = (mSom.Weight()[sortNeuro[neuro+1]][i]-c0)*t2 - (mSom.Weight()[sortNeuro[neuro+2]][i]-c0)*t1;
			c2/= t1*t2*(t1-t2);
			c1 = (mSom.Weight()[sortNeuro[neuro+1]][i]-c0)/t1 - c2*t1;

			//x=c0 + c1*t + c2*t*t + N(0,r)
			for(j=0; j<size; j++)
				popnew[popCur+j][i] = c0 + c1*t[j] + c2*t[j]*t[j] + r*rnd::gaussian();
		}

		popCur += size;
	}
	
	return popnew;
}
*********/
//CPopulationMO& ModelSOM::Generate1D1(CPopulationMO& popnew, CPopulationMO& popref)
//{
//	unsigned int i,j;
//	double r = 0.0;
//	
//	//gaussian noise variance
//	for( i=0; i<mDataSize; i++ ) r += mSom.DisToCore(i);
//	r /= double( mDataSize )*sqrt( double( mDataDim ) );
//	r *= PERTUBATION();

//	unsigned int num = mRowGrid*mColGrid;

//	FVECTOR C( 3 ), T( num ), X( num ), t( popnew.Size() );
//	T[0] = 0.0;
//	for( i=1; i<num; i++ )
//		T[i] = NormL2( mSom.Weight()[i], mSom.Weight()[i-1] ) + T[i-1];
//	for( i=0; i<num; i++ ) T[i] /= T[num-1];	
//	
//	double step = (1.0+2*mExtension)/double(popnew.Size());
//	double start= -mExtension;
//	for( i=0; i<popnew.Size(); i++ ) 
//	{
//		t[i]   = rnd::rand<double>(start,start+step); 
//		start += step;
//	}

//	for( i=0; i<mDataDim; i++ )
//	{
//		for( j=0; j<num; j++ ) X[j] = mSom.Weight()[j][i];
//		fit::polynomial1< FLOAT >( T, X, 2, C );
//		for( j=0; j<popnew.Size(); j++ ) popnew[j][i] = C[0] + C[1]*t[j] + C[2]*t[j]*t[j] + r*rnd::gaussian();
//	}

//	return popnew;
//}

/*************
//the hidden structure is 2D, use triangles to estimate new solutions
CPopulationMO& ModelSOM::Generate2D(CPopulationMO& popnew, CPopulationMO& popref)
{
	unsigned int row,col,i,j;
	double t1,t2;
	
	//distance to each neuro and the number of points attached to it
	double **pdis = new double*[mRowGrid];
	unsigned int **pnum = new unsigned int*[mRowGrid];
	for(row=0; row<mRowGrid; row++)
	{
		pdis[row] = new double[mColGrid];
		pnum[row] = new unsigned int[mColGrid];
		for(col=0; col<mColGrid; col++)
		{
			pdis[row][col] = 0.0;
			pnum[row][col] = 0;
		}
	}
	for( i=0; i<mDataSize; i++ ) 
	{ 
		pdis[mSom.NearCore(i)/mColGrid][mSom.NearCore(i)%mColGrid] += mSom.DisToCore(i); 
		pnum[mSom.NearCore(i)/mColGrid][mSom.NearCore(i)%mColGrid]++; 
	}

	//sphere radius of each trangle and its density
	double **pradius = new double*[2*(mRowGrid-1)];
	unsigned int **pdensity = new unsigned int*[2*(mRowGrid-1)];
	unsigned int total = 0;
	for(row=0; row<2*(mRowGrid-1); row++)
	{
		pradius[row] = new double[mColGrid-1];
		pdensity[row]= new unsigned int[mColGrid-1];
		for(col=0; col<mColGrid-1; col++)
		{
			pradius[row][col] = pdis[row/2  ][col  ]/double(pnum[row/2  ][col  ]) +
								pdis[row/2+1][col+1]/double(pnum[row/2+1][col+1]) +
								pdis[row/2+row%2][col+1-row%2]/double(pnum[row/2+row%2][col+1-row%2]);
			pradius[row][col]/= 3.0*sqrt(double(mDataDim));
			pradius[row][col]*= PERTUBATION();
			total += pnum[row/2  ][col  ] +
				    	pnum[row/2+1][col+1] +
						pnum[row/2+row%2][col+1-row%2];
			pdensity[row][col] = total;
        }
	}
	
	//sample offsprings
	unsigned int density;
	for(i=0; i<unsigned int(popnew.Size()); i++)
	{
		//find a trangle
		density = rnd::rand(unsigned int(1), unsigned int(total+1));
		for(row=0; row<2*(mRowGrid-1); row++)
		{
			for(col=0; col<mColGrid-1; col++)
				if(density<=pdensity[row][col])	break;
			if( col<mColGrid-1 && density<=pdensity[row][col]) break;
		}
		
		t1 = rnd::rand(-mExtension, 1.0+mExtension);
		t2 = rnd::rand(-mExtension, (t1<0)? 1.0+mExtension : 1.0+mExtension-t1);
		//lower trangle /|
		//				-
		if(row%2==0)
		{			
			//if((row/2)>0 && col==(mColGrid-2))
			//{
			//	t1 = rnd::rand<double>(-mExtension, 1.0);	//_
			//	t2 = rnd::rand<double>(0.0, 1.0-t1);		// |
			//}
			//else
			//{
			//	t2 = rnd::rand<double>( (row/2)==0 ? -mExtension:0.0, 1.0 );
			//	t1 = rnd::rand<double>( col==(mColGrid-2) ? -mExtension:0.0, 1.0-t2);
			//}
			for(j=0; j<mDataDim; j++)
				popnew[i][j] = (mSom.Weight(row/2  ,col  )[j] - mSom.Weight(row/2,col+1)[j])*t1 +
							(mSom.Weight(row/2+1,col+1)[j] - mSom.Weight(row/2,col+1)[j])*t2 +
								mSom.Weight(row/2  ,col+1)[j] +
							pradius[row][col]*rnd::gaussian();
		}
		//upper trangle  _
		//				|/
		else
		{	
			//if(col>0 && (row/2)==(mRowGrid-2))
			//{
			//	t2 = rnd::rand<double>(-mExtension, 1.0);	// |
			//	t1 = rnd::rand<double>(0.0, 1.0-t2);		//_
			//}
			//else
			//{
			//	t1 = rnd::rand<double>(col==0 ? -mExtension:0.0, 1.0);
			//	t2 = rnd::rand<double>((row/2)==(mRowGrid-2) ? -mExtension:0.0, 1.0-t1);
			//}

			for(j=0; j<mDataDim; j++)
				popnew[i][j] = (mSom.Weight(row/2+1,col+1)[j] - mSom.Weight(row/2+1,col)[j])*t1 +
							(mSom.Weight(row/2  ,col  )[j] - mSom.Weight(row/2+1,col)[j])*t2 +
								mSom.Weight(row/2+1,col  )[j] +
							pradius[row][col]*rnd::gaussian();
		}
	}

	//free space
	for(row=0; row<mRowGrid  ; row++) delete[]pdis[row]; delete[]pdis;
	for(row=0; row<mRowGrid  ; row++) delete[]pnum[row]; delete[]pnum;
	for(row=0; row<2*(mRowGrid-1); row++) delete[]pradius[row];  delete[]pradius;
	for(row=0; row<2*(mRowGrid-1); row++) delete[]pdensity[row]; delete[]pdensity;

	return popnew;
}
***/
