// LocalPCA.cpp

#include <algorithm>
#include <cmath>
#include "alg/Random.h"
#include "alg/LocalPCA.h"

namespace az
{

namespace alg
{

void LocalPCA::Train()
{
	unsigned int	c,	// cluster index
					m,	// row index
					n;	// col index
	unsigned int i;

	// initialize the PI matrix
	mvPI.resize(mNClu);
	for(c=0; c<mNClu; c++)
	{
		mvPI[c].resize(mDX);
		for(m=0; m<mDX; m++) mvPI[c][m].resize(mDX);
		for(m=0; m<mDX; m++) for(n=m; n<mDX; n++) mvPI[c][m][n] = mvPI[c][m][n] = (m==n) ? 1.0 : 0.0;
	}

	// initialize the center points
	std::vector<unsigned int> index(mNX);
	for(i=0; i<mNX; i++) index[i] = i;
	for(c=0; c<mNClu; c++)
	{
		std::swap(index[c], index[rnd::rand(c, mNX)]);
		mvMean[c] = mvX[index[c]];
	}

	// trainning
	unsigned int failupdate = 0;
	double	dis, mindis, pro, err;
	std::vector<double> meanold(mDX);
	while( mIter++ < mMaxIter && failupdate<mNClu )
	{
		// partition trainning data
		for(c=0; c<mNClu; c++) mvNo[c] = 0;	// number point in each cluster
		for(m=0; m<mNX; m++)
		{
			mindis		= 1.0E100;
			mvIndex[m]	= 0;				// to which cluster it belongs
			for(c=0; c<mNClu; c++)
				if( (dis = Distance(m,c))< mindis)
				{
					mindis		= dis;
					mvIndex[m]	= c;
				}
			mvNo[mvIndex[m]]++;
		}

		// update parameters
		failupdate	= 0;
		for(c=0; c<mNClu; c++) 
		{
			// save old mean
			meanold = mvMean[c];

			// no data assigned to cluster c
			if(mvNo[c] < 1)
			{
				mvMean[c] = mvX[rnd::rand((unsigned int)0, (unsigned int)mNX)];
				for(m=0; m<mDX; m++) 
				{
					mvPI[c][m][m] = 1.0;
					for(n=m+1; n<mDX; n++) mvPI[c][m][n] = mvPI[c][n][m] = 0.0;
				}
			}
			// only one data assigned to cluster c
			else if(mvNo[c] == 1)
			{
				for(m=0; m<mNX; m++) if(mvIndex[m] == c) break;
				mvMean[c] = mvX[m];
				for(m=0; m<mDX; m++) 
				{
					mvPI[c][m][m] = 1.0;
					for(n=m+1; n<mDX; n++) mvPI[c][m][n] = mvPI[c][n][m] = 0.0;
				}
			}
			// more than one data assigned to cluster c
			else if(mvNo[c] > 1)
			{
				// find data index
				index.resize(mvNo[c]); n=0;
				for(m=0; m<mNX; m++) if(mvIndex[m] == c) index[n++] = m;
				
				// calculate the eigenvalues and eigenvectors
				Eigen(mvMean[c], mvEigenvalue[c], mvEigenvector[c], index);

				// calculate pi
				for(m=0; m<mDX; m++) for(n=0; n<mDX; n++)
				{
					mvPI[c][m][n] = 0.0;
					for(i=mDLat; i<mDX; i++) mvPI[c][m][n] += mvEigenvector[c][i][m]*mvEigenvector[c][i][n];
				}
			}// else if(mvNo[c] > 1)

			// calculate the error
			err = 0.0;
			for(m=0; m<mDX; m++) err += (meanold[m] - mvMean[c][m])*(meanold[m] - mvMean[c][m]);
			failupdate += (sqrt(err)<mErrTol) ? 1:0;
		}// for(c=0; c<mNClu; c++)
	}// while( mIter++ < mMaxIter && failupdate<mNClu )

	// calculate the projection
	for(c=0; c<mNClu; c++) if(mvNo[c] > 1)
	{
		// copy data
		index.resize(mvNo[c]); n=0;
		for(m=0; m<mNX;   m++) if(mvIndex[m] == c) index[n++] = m;
		for(m=0; m<mDLat; m++) {mvProMin[c][m] = 1.0E100; mvProMax[c][m] = -1.0E100;}
		for(m=0; m<index.size(); m++) for(n=0; n<mDLat; n++)
		{
			pro = 0.0;
			for(i=0; i<mDX; i++) pro += (mvX[index[m]][i]-mvMean[c][i])*mvEigenvector[c][n][i];
			if(pro>mvProMax[c][n]) mvProMax[c][n] = pro;
			if(pro<mvProMin[c][n]) mvProMin[c][n] = pro;
		}
	}
}

double LocalPCA::Distance(unsigned int m, unsigned int c)
{
	unsigned int i,j;
	double dis = 0.0;
	std::vector<double> tmpdis(mDX);
	for(i=0; i<mDX; i++)
	{
		tmpdis[i] = 0.0;
		for(j=0; j<mDX; j++) tmpdis[i] += (mvX[m][j]-mvMean[c][j])*mvPI[c][j][i];
	}
	for(i=0; i<mDX; i++) dis += tmpdis[i]*(mvX[m][i]-mvMean[c][i]);
	return dis;
}

} //namespace alg

} //namespace az
