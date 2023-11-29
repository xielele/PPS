// PCA.cpp

#include <cmath>
#include <fstream>
#include "alg/PCA.h"

namespace az
{

namespace alg
{

//constructor
PCA::PCA()
{
	pData = 0;
}

//constructor
PCA::PCA(alg::Matrix& data)
{
	pData = &data;
}

//destructor
PCA::~PCA()
{
	pData = 0;
}

//initialize colvariance matrix to be an identity
void PCA::Initialize(unsigned int dim)
{
	mCov.Identity(dim);
	mEigenvector.Identity(dim);
	mEigenvalue.resize(dim);
	for(unsigned int i=0; i<dim; i++) mEigenvalue[i] = 1.0;
}

//get the range of the projections in a dimension 
void PCA::PrimaryRange(unsigned int dim, double& min, double& max)
{
	unsigned int i, j;
	double tmp;
	min = 1.0e20; max = -1.0e20;
	for(i=0; i<(*pData).ColSize(); i++)
	{
		tmp = 0.0;
		for(j=0; j<mMean.size(); j++)
			tmp += ((*pData)(j, i) - mMean[j])*mEigenvector(j, dim);
		if(tmp < min) min = tmp;
		if(tmp > max) max = tmp;
	}
}	
//calculate the mean, eigenvalue and eigenvectors 
void PCA::Train()
{ 
	//Step 1: calculate the mean
	(*pData).ColMean(mMean);

	//Step 2: get the covariance matrix
	alg::Matrix datatmp = (*pData);
	datatmp.ColSub(mMean);
	alg::Matrix datatmp1 = datatmp;
	datatmp1.Trans();
	datatmp.Multiply(datatmp1, mCov);
	mCov.Divide(double( (*pData).ColSize()>1 ? (*pData).ColSize()-1.0 : (*pData).ColSize() ));

	//Step 3: calculate the eigenvalue and eigenvector
	mCov.Eig(mEigenvalue, mEigenvector);
}

//write results into a stream
std::ostream& operator<<(std::ostream& os, PCA& pca)
{
	pca.Write(os);
	return os;
}

//write results into a stream
void PCA::Write( std::ostream& os )
{
	unsigned int i;
	//datas
	os<<(*pData).RowSize()<<"\t"<<(*pData).ColSize()<<std::endl;
	os<<(*pData)<<std::endl<<std::endl;
	//cov
	os<<"covariance"<<std::endl;
	os<<mCov<<std::endl<<std::endl;
	//mean
	os<<"mean"<<std::endl;
	for(i=0; i<mMean.size(); i++)
		os<<mMean[i]<<"\t";
	os<<std::endl<<std::endl;
	//eigenvalue
	os<<"eigenvalue"<<std::endl;
	for(i=0; i<mEigenvalue.size(); i++)
		os<<mEigenvalue[i]<<"\t";
	os<<std::endl<<std::endl;
	//eigenvector
	os<<"eigenvector"<<std::endl;
	os<<mEigenvector<<std::endl;
	//range
	os<<"range in primary dimension"<<std::endl;
	double min, max;
	PrimaryRange(0, min, max);
	os<<min<<"\t"<<max<<std::endl;
}

//write results into a file
void PCA::Write(std::string& name)
{
	std::ofstream os;
	os.open(name.c_str());
	
	Write(os);

	os.close();
}

//write results into a file
void PCA::Write( const char* name )
{
	std::string str = std::string(name);
	Write(str);
}

} //namespace alg

} //namespace az

