// Matrix.cpp

#include <float.h>
#include <math.h>
#include "alg/Matrix.h"

#if defined(WIN32)
    #define wxFinite(n) _finite(n)
#elif defined(_LINUX)
    #define wxFinite(n) finite(n)
#else
    #define wxFinite(n) ((n)==(n))
#endif

namespace az
{

namespace alg
{

//constructor
Matrix::Matrix(unsigned int row, int unsigned col)
	:mRow(row), mCol(col)
{
	if(mRow* mCol > 0)
	{
		pData = new double[mRow * mCol];
		for(unsigned int i=0; i<mRow*mCol; i++) pData[i]=0.0;
	}
	else 
		pData = 0;
}

//constructor
Matrix::Matrix(const Matrix& mat)
	:mRow(mat.mRow), mCol(mat.mCol)
{
	if(mRow* mCol > 0)
	{
		pData = new double[mRow * mCol];
		memcpy(pData, mat.pData, mRow*mCol*sizeof(double));
	}
	else 
		pData = 0;
}

//destructor
Matrix::~Matrix()
{
	if(pData) delete []pData;
}

//reset the matrix size
Matrix& Matrix::Resize(unsigned int row, unsigned int col)
{
	if(row != mRow || col != mCol)
	{
		if(pData) delete []pData;
		mRow	= row;
		mCol	= col;
		pData = new double[mRow*mCol];
		for(unsigned int i=0; i<mRow*mCol; i++)	pData[i] = 0.0;
	}
	return *this;
}

//create an identity matrix
Matrix& Matrix::Identity(unsigned int size)
{
	(*this).Resize(size, size);
	for(unsigned int i=0; i<size; i++)	(*this)(i, i) = 1.0;
	return *this;
}

//get an element
double& Matrix::operator()(unsigned int row, unsigned int col)
{
	CHECK(row < mRow && col < mCol, "Matrix::(row,col)");
	return pData[col*mRow + row];
}

//reset to another matrix
Matrix& Matrix::operator= (const Matrix& mat)
{
	Resize(mat.mRow, mat.mCol);
	for(unsigned int i=0; i<mRow*mCol; i++) pData[i] = mat.pData[i];
	return *this;
}

//get a row
FVECTOR& Matrix::Row(unsigned int row, FVECTOR& value)
{
	CHECK(row < mRow, "Matrix::Row()");
	value.resize(mCol);
	for(unsigned int i=0; i<mCol; i++) value[i] = (*this)(row, i);
	return value;
}

//get a column
FVECTOR& Matrix::Column(unsigned int col, FVECTOR& value)
{
	CHECK(col < mCol, "Matrix::Column()");
	value.resize(mRow);
	for(unsigned int i=0; i<mRow; i++) value[i] = (*this)(i, col);
	return value;
}

//get a sub-matrix except a row and a column
Matrix& Matrix::Sub(unsigned int row, unsigned int col, Matrix& mat)
{
	CHECK(row < mRow && mRow > 1 && col < mCol && mCol > 1, "Matrix::Sub()");

	mat.Resize(mRow-1, mCol-1);

	unsigned int i,j;

	std::vector<unsigned int> ir(mRow-1), jc(mCol-1);
	for(i=0; i<mRow-1; i++) ir[i] = i<row ? i : i+1;
	for(j=0; j<mCol-1; j++) jc[j] = j<col ? j : j+1;

	for(i=0; i<mat.RowSize(); i++)
		for(j=0; j<mat.ColSize(); j++)
			mat(i,j) = (*this)(ir[i],jc[j]);
	return mat;
}

//calculate the determinant of a square matrix
double Matrix::Det()
{
	CHECK(mRow == mCol, "Matrix::Det()");

	// compute the determinant by definition
	//if(mRow == 1) return (*this)(0,0);

	//if(mRow == 2) return (*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0);
	//
	//if(mRow == 3) return  (*this)(0,0) * (*this)(1,1) * (*this)(2,2) 
	//					 +(*this)(0,1) * (*this)(1,2) * (*this)(2,0) 
	//					 +(*this)(0,2) * (*this)(1,0) * (*this)(2,1)
	//					 -(*this)(0,0) * (*this)(1,2) * (*this)(2,1) 
	//					 -(*this)(0,1) * (*this)(1,0) * (*this)(2,2) 
	//					 -(*this)(0,2) * (*this)(1,1) * (*this)(2,0);
	//
	//double tmp = 0.0;
	//Matrix submatrix;
	//for(unsigned int i=0; i<mCol; i++)
	//	tmp += (i % 2 == 0 ? 1:-1) * (*this)(0,i) * (*this).Sub(0, i, submatrix).Det();
	//return tmp;

	// modified Nov.19 2006
	Matrix mat = *this;
	std::vector<unsigned int> indx(mCol); 
	double d;
	LUdcmp(mat, indx, d);
	for(unsigned int i=0; i<mCol; i++) d *= mat(i,i);
	return d;
}

//translate the matrix
Matrix& Matrix::Trans()
{
	if(mRow<2 || mCol<2) 
	{
		std::swap((*this).mCol, (*this).mRow);
		return *this;
	}
	Matrix mat(*this);
	std::swap((*this).mCol, (*this).mRow);
	for(unsigned int i=0; i<mRow; i++)
		for(unsigned int j=0; j<mCol; j++)
			(*this)(i, j) = mat(j, i);				
	return *this;
}

//inverse the matrix
Matrix& Matrix::Inv()
{
	CHECK(mRow == mCol, "Matrix::Inv()");

	//Adjoint method
	//unsigned int i,j;
	//double det=0.0;
	//Matrix adjmat(mRow,mCol), submatrix;
	//for(i=0; i<mRow; i++)
	//	for(j=0; j<mCol; j++)
	//		adjmat(i, j) = ((i + j)%2 == 0 ? 1 : -1)*(*this).Sub(j, i, submatrix).Det();
	//for(i=0; i<mRow; i++) det += (*this)(0,i)*adjmat(i,0);
	//*this = adjmat; this->Divide(det);

	//LU decomposition version
	Matrix lumat(*this);
	std::vector<unsigned int> indx(mRow);
	std::vector<double> col(mRow);
	double d;
	unsigned i,j;
	LUdcmp(lumat,indx,d);
	for(j=0; j<mRow; j++)
	{
		for(i=0; i<mRow; i++) col[i]=0.0;
		col[j]=1.0;
		LUbksb(lumat,indx,col);
		for(i=0; i<mRow; i++) (*this)(i,j)=col[i];
	}
	
	return *this;
}

//calculate the eigenvalue and egienvectors
void Matrix::Eig(FVECTOR& eigvalue, Matrix& eigvector)
{
	CHECK(mRow == mCol, "Matrix::Eig()");
	eigvalue.resize(mRow);
	FVECTOR interm;
	interm.resize(mRow);
	eigvector = *this;
	
	tred2(eigvalue, interm, eigvector);
	tqli(eigvalue, interm, eigvector);
	//Eigenvector in columns
	Sort(eigvalue, eigvector);
}

//multiply a matrix
Matrix& Matrix::Multiply(Matrix& mat, Matrix& result)
{
	unsigned int i,j,k;

	CHECK(mCol == mat.RowSize(), "Matrix::Myltiply()");
	
	result.Resize((*this).RowSize(), mat.ColSize());
	
	for(i=0; i<result.RowSize(); i++)
		for(j=0; j<result.ColSize(); j++)
		{
			result(i, j) = 0;
			for(k=0; k<mat.RowSize(); k++)
				result(i, j) += (*this)(i, k) * mat(k, j);
		}
	return result;
}

//left multiply a vector
FVECTOR& Matrix::LeftMultiply(FVECTOR& vec, FVECTOR& result)
{
	CHECK(mRow == vec.size(), "Matrix::LeftMultiply()");
	result.resize((*this).ColSize());
	for(unsigned int i=0; i<(*this).ColSize(); i++)
	{
		result[i] = 0;
		for(unsigned int j=0; j<(*this).RowSize(); j++)
			result[i] += (*this)(j, i) * vec[j];
	}
	return result;
}

//right multiply a vector
FVECTOR& Matrix::RightMultiply(FVECTOR& vec, FVECTOR& result)
{
	CHECK(mCol == vec.size(), "Matrix::RightMultiply()");
	result.resize((*this).RowSize());
	for(unsigned int i=0; i<(*this).RowSize(); i++)
	{
		result[i] = 0;
		for(unsigned int j=0; j<(*this).ColSize(); j++)
			result[i] += (*this)(i, j) * vec[j];
	}
	return result;
}

//divide a scalar
Matrix& Matrix::Divide(double sca)
{
	unsigned int i,j;
	for(i=0; i<(*this).RowSize(); i++)
		for(j=0; j<(*this).ColSize(); j++)
			(*this)(i, j) /= sca;
	return *this;
}

//get the mean of all columns
FVECTOR& Matrix::ColMean(FVECTOR& mean)
{
	mean.resize((*this).RowSize());
	unsigned int i,j;
	for(i=0; i<mean.size(); i++)
	{
		mean[i] = 0.0;
		for(j=0; j<(*this).ColSize(); j++)
			mean[i] += (*this)(i, j);
		mean[i] /= (double) (*this).ColSize() ;
	}
	return mean;
}

//get the mean of all rows
FVECTOR& Matrix::RowMean(FVECTOR& mean)
{
	mean.resize((*this).ColSize());
	unsigned int i,j;
	for(i=0; i<mean.size(); i++)
	{
		mean[i] = 0;
		for(j=0; j<(*this).RowSize(); j++)
			mean[i] += (*this)(j, i);
		mean[i] /= (double) (*this).RowSize() ;
	}	
	return mean;
}

//get standard variation of all columns
FVECTOR& Matrix::ColStd(FVECTOR& std)
{
	std.resize((*this).RowSize());
	FVECTOR mean;
	ColMean(mean);
	unsigned int i,j;
	for(i=0; i<std.size(); i++)
	{
		std[i] = 0.0;
		for(j=0; j<(*this).ColSize(); j++)
			std[i] += ((*this)(i,j) - mean[i])*((*this)(i,j) - mean[i]);
		std[i] /= double((*this).ColSize() -1);
		std[i] = sqrt(std[i]);
	}
	return std;
}

//get standard variation of all rows
FVECTOR& Matrix::RowStd(FVECTOR& std)
{
	std.resize((*this).ColSize());
	FVECTOR mean;
	RowMean(mean);
	unsigned int i,j;
	for(i=0; i<std.size(); i++)
	{
		std[i] = 0.0;
		for(j=0; j<(*this).RowSize(); j++)
			std[i] += ((*this)(j,i) - mean[i])*((*this)(j,i) - mean[i]);
		std[i] /= double((*this).RowSize() -1);
		std[i] = sqrt(std[i]);
	}
	return std;
}

//subtract a row vector
Matrix& Matrix::RowSub(FVECTOR& value)
{
	unsigned int i,j;
	for(i=0; i<(*this).RowSize(); i++)
		for(j=0; j<(*this).ColSize(); j++)
			(*this)(i, j) -= value[j];
	return *this;
}

//subtract a column vector
Matrix& Matrix::ColSub(FVECTOR& value)
{
	unsigned int i,j;
	for(i=0; i<(*this).RowSize(); i++)
		for(j=0; j<(*this).ColSize(); j++)
			(*this)(i, j) -= value[i];
	return *this;
}

//read a matrix
std::istream& operator>>(std::istream& is, Matrix& mat)
{
	unsigned int row,col;
	is>>row>>col;
	mat.Resize(row,col);
	for(unsigned int i=0; i<mat.mRow; i++)
		for(unsigned int j=0; j<mat.mCol; j++)
			is>>mat(i, j);
	return is;
}

//write a matrix
std::ostream& operator<< (std::ostream& os, Matrix& mat)
{
	os<<mat.RowSize()<<"\t"<<mat.ColSize()<<std::endl;
	os.setf(std::ios::right|std::ios::scientific);
	os.precision(5);   
	for(unsigned int i=0; i<mat.mRow; i++)
	{
		for(unsigned int j=0; j<mat.mCol; j++)
			os<<mat(i, j)<<"\t";
		os<<std::endl;
	}
	return os;
}

//solve A X = b, mat if the decomposition of A
void LUbksb(Matrix& mat, std::vector<unsigned int>& indx, std::vector<double>& b)
{
	int i,ii=0,ip,j, n=mat.RowSize();
	double sum;

	for(i=0;i<n;i++) 
	{
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if(ii!=0)
			for(j=ii-1;j<=i;j++) sum -= mat(i,j)*b[j];
		else if(sum!=0.0) ii=i+1;
		b[i]=sum;
	}

	for(i=n-1;i>=0;i--) 
	{
		sum=b[i];
		for(j=i+1;j<n;j++) sum -= mat(i,j)*b[j];
		b[i]=sum/mat(i,i);
	}
}

//LU decomposition 
void LUdcmp(Matrix& mat, std::vector<unsigned int>& indx, double& d)
{
	const double TINY = 1.0e-20;
	unsigned int i,imax=0,j,k;
	double big,dum,sum,temp;
	unsigned int n = mat.RowSize();
	std::vector<double> vv(n);

	d=1.0;
	for(i=0;i<n;i++) 
	{
		big=0.0;
		for(j=0;j<n;j++)
			if((temp=fabs(mat(i,j)))>big) big=temp;
		//if(big == 0.0) std::cout<<"error"<<std::endl;
		vv[i]=1.0/big;
	}

	for(j=0;j<n;j++) 
	{
		for(i=0;i<j;i++) 
		{
			sum = mat(i,j);
			for(k=0;k<i;k++) sum -= mat(i,k)*mat(k,j);
			mat(i,j) = sum;
		}
		big=0.0;
		for (i=j;i<n;i++) 
		{
			sum=mat(i,j);
			for(k=0;k<j;k++)
				sum -= mat(i,k)*mat(k,j);
			mat(i,j) = sum;
			if((dum=vv[i]*fabs(sum)) >= big) 
			{
				big=dum;
				imax=i;
			}
		}
		if(j != imax) 
		{
			for(k=0;k<n;k++) 
				std::swap(mat(imax,k), mat(j,k));
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (mat(j,j) == 0.0) mat(j,j)=TINY;
		if (j != n-1) 
		{
			dum=1.0/mat(j,j);
			for (i=j+1;i<n;i++) mat(i,j) *= dum;
		}
	}
}

//Householder reduction of Matrix a to tridiagonal form.
//
// Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
// Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
// Springer-Verlag, 1976, pp. 489-494.
// W H Press et al., Numerical Recipes in C, Cambridge U P,
// 1988, pp. 373-374. 
void Matrix::tred2(FVECTOR& eigenvalue, FVECTOR& interm, Matrix& eigenvector)
{
	double scale, hh, h, g, f;
	int k, j, i,l;

	for(i = mRow-1; i >= 1; i--)
	{
		l = i - 1;
		h = scale = 0.0;
		if(l > 0)
		{
			for (k = 0; k <= l; k++)
				scale += fabs(eigenvector(i, k));
			if (scale == 0.0)
				interm[i] = eigenvector(i, l);
			else
			{
				for (k = 0; k <= l; k++)
				{
					eigenvector(i, k) /= scale;
					h += eigenvector(i, k) * eigenvector(i, k);
				}
				f = eigenvector(i, l);
				g = f>0 ? -sqrt(h) : sqrt(h);
				interm[i] = scale * g;
				h -= f * g;
				eigenvector(i, l) = f - g;
				f = 0.0;
				for (j = 0; j <= l; j++)
				{
					eigenvector(j, i) = eigenvector(i, j)/h;
					g = 0.0;
					for (k = 0; k <= j; k++)
						g += eigenvector(j, k) * eigenvector(i, k);
					for (k = j+1; k <= l; k++)
						g += eigenvector(k, j) * eigenvector(i, k);
					interm[j] = g / h;
					f += interm[j] * eigenvector(i, j);
				}
				hh = f / (h + h);
				for (j = 0; j <= l; j++)
				{
					f = eigenvector(i, j);
					interm[j] = g = interm[j] - hh * f;
					for (k = 0; k <= j; k++)
						eigenvector(j, k) -= (f * interm[k] + g * eigenvector(i, k));
				}
			}
		}
		else
			interm[i] = eigenvector(i, l);
		eigenvalue[i] = h;
	}
	eigenvalue[0] = 0.0	;	
	interm[0] = 0.0;
	for(i = 0; i < (int)mRow; i++)
	{
		l = i - 1;
		if(eigenvalue[i])
		{
			for (j = 0; j <= l; j++)
			{
				g = 0.0;
				for (k = 0; k <= l; k++)
					g += eigenvector(i, k) * eigenvector(k, j);
				for (k = 0; k <= l; k++)
					eigenvector(k, j) -= g * eigenvector(k, i);
			}
		}
		eigenvalue[i] = eigenvector(i, i);
		eigenvector(i, i) = 1.0;
		for (j = 0; j <= l; j++)
			eigenvector(j, i) = eigenvector(i, j) = 0.0;
	}
}

#define SIGN(a, b) ((b) < 0 ? -fabs(a) : fabs(a))
//Tridiagonal QL algorithm -- Implicit 
void Matrix::tqli(FVECTOR& eigenvalue, FVECTOR& interm, Matrix& eigenvector)
{
	int m, l, iter, i, k;
	double s, r, p, g, f, dd, c, b;
	
	for (i = 1; i < (int)mRow; i++)
		interm[i-1] = interm[i];
	interm[mRow-1] = 0.0;
	for (l = 0; l < (int)mRow; l++)
	{
		iter = 0;
		do
		{
			if (iter++ > 100) break;

			for (m = l; m < (int)mRow-1; m++)
			{
				dd = fabs(eigenvalue[m]) + fabs(eigenvalue[m+1]);
				if (fabs(interm[m]) + dd == dd || !wxFinite(dd)) break;
			}
			if (m != l)
			{
				//if (iter++ == 30) erhand("No convergence in TLQI.");
				g = (eigenvalue[l+1] - eigenvalue[l]) / (2.0 * interm[l]);
				r = sqrt((g * g) + 1.0);
				g = eigenvalue[m] - eigenvalue[l] + interm[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m-1; i >= l; i--)
				{
					f = s * interm[i];
					b = c * interm[i];
					if (fabs(f) >= fabs(g))
					{
						c = g / f;
						r = sqrt((c * c) + 1.0);
						interm[i+1] = f * r;
						c *= (s = 1.0/r);
					}
					else
					{
						s = f / g;
						r = sqrt((s * s) + 1.0);
						interm[i+1] = g * r;
						s *= (c = 1.0/r);
					}
					g = eigenvalue[i+1] - p;
					r = (eigenvalue[i] - g) * s + 2.0 * c * b;
					p = s * r;
					eigenvalue[i+1] = g + p;
					g = c * r - b;
					for (k = 0; k < (int)mRow; k++)
					{
						f = eigenvector(k, i+1);
						eigenvector(k, i+1)= s * eigenvector(k, i) + c * f;
						eigenvector(k, i) = c * eigenvector(k, i) - s * f;
					}
				}
				eigenvalue[l] = eigenvalue[l] - p;
				interm[l] = g;
				interm[m] = 0.0;
			}
		}  while (m  != l);
	}
}

//sort the eigenvalue by decreasing order
void Matrix::Sort(FVECTOR& eigenvalue, Matrix& eigenvector)
{
	unsigned int i,j,k;
	unsigned int m = (unsigned int)eigenvalue.size();

	//repair 
	for(i=0; i<m; i++) if(!wxFinite(eigenvalue[i])) eigenvalue[i] = 1.0E-100;
	for(i=0; i<m; i++) 
	{
		j=0;
		for(k=0; k<m; k++) if(!wxFinite(eigenvector(k, i))) j = 1;
		if(j>0) for(k=0; k<m; k++) eigenvector(k, i) = 1.0/sqrt(m+0.0);
	}

	for(i=0; i<m; i++) if(eigenvalue[i]<0)
	{
		eigenvalue[i] *= -1;
		for(j=0; j<m; j++) eigenvector(j,i) *= -1;
	}
	for(i=0; i<m - 1; i++)
		for(j=i+1; j<m; j++)
			if(eigenvalue[j] > eigenvalue[i])
			{
				std::swap(eigenvalue[j], eigenvalue[i]);
				for(k=0; k<m; k++)
					std::swap(eigenvector(k, j), eigenvector(k, i));
			}
}

} //namespace alg

} //namespace az

