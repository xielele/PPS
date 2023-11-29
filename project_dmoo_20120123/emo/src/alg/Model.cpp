// Model.cpp
#ifdef AZ_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#else
#include "alg/Matrix.h"
#endif
#include <iomanip>
#include "alg/Model.h"

namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{

void Model::Set(unsigned int nclu, unsigned int nx, unsigned int dx, unsigned int dlat, unsigned int maxiter, double errtol)
{
	unsigned int i,j;
	
	mNClu	= nclu;
	mNX		= nx;
	mDX		= dx; 
	mDLat	= dlat;
	mMaxIter= maxiter;
	mErrTol	= errtol;
	mIter	= 0;

	mvIndex.resize(mNX);
	mvNo.resize(mNClu);
	
	mvX.resize(mNX);
	for(i=0; i<mNX; i++) mvX[i].resize(mDX);

	mvMean.resize(mNClu);
	for(i=0; i<mNClu; i++) mvMean[i].resize(mDX);

	mvEigenvalue.resize(mNClu);
	for(i=0; i<mNClu; i++) mvEigenvalue[i].resize(mDX);

	mvProMin.resize(mNClu);
	mvProMax.resize(mNClu);
	for(i=0; i<mNClu; i++)
	{
		mvProMin[i].resize(mDLat);
		mvProMax[i].resize(mDLat);
	}

	mvEigenvector.resize(mNClu);
	for(i=0; i<mNClu; i++) 
	{
		mvEigenvector[i].resize(mDX);
		for(j=0; j<mDX; j++) mvEigenvector[i][j].resize(mDX);
	}
}

void Model::Write(std::string file)
{
	unsigned int i,j,k;
	std::ofstream out(file.c_str());
	out<<std::scientific<<std::setprecision(5);
	
	out<<"Train Steps "<<mIter<<std::endl;

	for(i=0; i<mNClu; i++)
	{
		out<<std::endl<<"===========cluster "<<i<<"==========="<<std::endl;
		out<<"data"<<std::endl;
		for(j=0;j<mNX; j++) if(mvIndex[j]==i)
		{
			for(k=0; k<mDX; k++) out<<mvX[j][k]<<"\t";
			out<<std::endl;
		}
		out<<"mean"<<std::endl;	for(j=0; j<mDX; j++) out<<mvMean[i][j]<<"\t";out<<std::endl;
		out<<"eigenvalue"<<std::endl;	for(j=0; j<mDX; j++) out<<mvEigenvalue[i][j]<<"\t";out<<std::endl;
		out<<"eigenvector"<<std::endl;
		for(j=0; j<mDX; j++)
		{
			for(k=0; k<mDX; k++) out<<mvEigenvector[i][j][k]<<"\t"; out<<std::endl;
		}
		//out<<"PI"<<std::endl;
		//for(j=0; j<mDX; j++)
		//{
		//	for(k=0; k<mDX; k++) out<<mvPI[i][j][k]<<"\t"; out<<std::endl;
		//}
	}

	out.close();
}

void Model::Eigen(std::vector<double>& mean, std::vector<double>& eva, std::vector< std::vector<double> >& eve, std::vector< unsigned int >& index)
{
	unsigned int i,j,k;
	//calculate the mean
	for(i=0; i<mDX; i++)
	{
		mean[i] = 0.0;
		for(j=0; j<index.size(); j++) mean[i] += mvX[index[j]][i];
		mean[i] /= double(index.size());
	}

	//calulate the covariance
	std::vector< std::vector<double> > cov(mDX); for(i=0; i<mDX; i++) cov[i].resize(mDX);
	for(i=0; i<mDX; i++)
	{
		for(j=i; j<mDX; j++)
		{
			cov[i][j]  = 0.0;
			for(k=0; k<index.size(); k++) cov[i][j] += (mvX[index[k]][i] - mean[i])*(mvX[index[k]][j] - mean[j]);
			cov[i][j] /= double(index.size()-0.0);
			cov[j][i]  = cov[i][j];
		}
	}
	alg::Eigen(eva, eve, mDX, cov);
	cov.clear();
}

void Eigen(std::vector<double>& eva, std::vector< std::vector<double> >& eve, unsigned int no, std::vector< std::vector<double> >& cov)
{
	unsigned int i,j,mDX = (unsigned int)cov.size();
#ifdef AZ_GSL
//============================================================================================================
	gsl_matrix* cov1 	= gsl_matrix_alloc(mDX,mDX);
	gsl_vector* eval 	= gsl_vector_alloc(mDX);
	gsl_matrix* evec	= gsl_matrix_alloc(mDX,mDX);
	
	for(i=0; i<mDX; i++)
		for(j=0; j<mDX; j++)
			gsl_matrix_set(cov1, i, j, cov[i][j]);

	gsl_set_error_handler_off();
	gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(mDX);
	gsl_eigen_symmv(cov1, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);
	
	for(i=0; i<no; i++)
	{
		eva[i] = fabs(gsl_vector_get(eval,i));
		for(j=0; j<mDX; j++) eve[i][j] = gsl_matrix_get(evec,j,i);
	}
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_matrix_free(cov1);
#else
//============================================================================================================
	alg::Matrix cov1(mDX,mDX), eve1(mDX,mDX);
	std::vector<double> eva1(mDX);
	for(i=0; i<mDX; i++)
		for(j=0; j<mDX; j++)
			cov1(i,j)  = cov[i][j];
	cov1.Eig(eva1, eve1);
	for(i=0; i<no; i++)
	{
		eva[i] = eva1[i];
		for(j=0; j<mDX; j++) eve[i][j] = eve1(j,i);
	}
	eva1.clear();
//============================================================================================================
#endif
}

} //namespace alg
} //namespace az
