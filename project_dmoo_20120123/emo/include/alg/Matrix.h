/*! \file	Matrix.h
	
	\brief	Matrix: contains 2-D data
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Dec.30 2004 create
	\date	Sep.25 2005 rewrite & reorganize structure
	\date	Oct.20 2005 remove bugs in Sub,Det,Inv
	\date	Oct.21 2005 add LU decomposition
	\date	Oct.23 2005 add linear algebra functions(Cholesky,CholeskySolve,SVD,pinv)
*/

#ifndef	AZ_MATRIX_H
#define	AZ_MATRIX_H

#include <list>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	alg namespace, contains algorithms
namespace alg
{
	//!\brief error process for DATA namespace
	class error : public std::exception 
	{
	public:
		//!\brief	constructor
		//!\param	msg error message
		//!\return	void
		error(std::string const& msg) throw()
			: msg_(std::string("DATA error: ") + msg)
		{}
		
		//!\brief	destructor
		virtual ~error() throw() {}

		//!\brief	look up the error reason
		//!\return	error message
		virtual const char* what() const throw() { return msg_.c_str(); }
	protected:	
		std::string msg_;	//!< error message
	};

	//! check to ensure the expression is right
	#ifdef _DEBUG
		#define CHECK(cond, str)	if(!(cond))	{ throw error(str); }
	#else
		#define CHECK(cond, str)	{}	//if(!(cond))	{ throw error(str); }
	#endif

	//!\brief real vector
	typedef std::vector<double>	FVECTOR;

	//! index structure based on list
	typedef std::list< unsigned int >	LINDEX;

	//! index structure based on vector
	typedef std::vector< unsigned int > VINDEX;

	//!\brief Matrix class
	class Matrix
	{
	protected:
		unsigned int mRow,	//!< row size 
					 mCol;	//!< column size
		double*		 pData;	//!< pointer to the data		
	public:
		//!\brief	constructor
		//!\param	row row size
		//!\param	col column size
		//!\return	void
		Matrix(unsigned int row = 0, int unsigned col = 0);

		//!\brief	constructor
		//!\brief	mat reference matrix
		//!\return	void
		Matrix(const Matrix& mat);

		//!\brief	destructor
		//!\return	void
		~Matrix();

		//!\brief	reset the matrix size
		//!\param	row row size
		//!\param	col column size
		//!\return	reference of the matrix
		Matrix& Resize(unsigned int row, unsigned int col);

		//!\brief	create an identity matrix
		//!\param	size row and column size
		//!\return	reference of the matrix
		Matrix& Identity(unsigned int size);

		//!\brief	get the row size
		//!\return	row size
		inline unsigned int RowSize() {return mRow;}

		//!\brief	get the column size
		//!\return	column size
		inline unsigned int ColSize() {return mCol;}

		//!\brief	get the pointer to the data
		//!\return	the pointer to the data
		inline double* operator()() {return pData;}
		
		//!\brief	get an element
		//!\param	row row number
		//!\param	col column number
		//!\return	reference to the element
		double& operator()(unsigned int row, unsigned int col);

		//!\brief	reset to another matrix
		//!\param	mat matrix reference
		//!\return	reference of the matrix
		Matrix& operator= (const Matrix& mat);
		
		//!\brief	get a row
		//!\param	row row number
		//!\param	value a vector to store the row
		//!\return	reference to value
		FVECTOR& Row(unsigned int row, FVECTOR& value);

		//!\brief	get a column
		//!\param	col column number
		//!\param	value a vector to store the column
		//!\return	reference to value
		FVECTOR& Column(unsigned int col, FVECTOR& value);

		//!\brief	get a sub-matrix except a row and a column
		//!\param	row row number
		//!\param	col column number
		//!\param	mat a matrix to store the sub-matrix
		//!\return	sub-matrix
		Matrix& Sub(unsigned int row, unsigned int col, Matrix& mat);

		//!\brief	calculate the determinant of a square matrix
		//!\return	the deterministic
		double Det();

		//!\brief	translate the matrix
		//!\return	reference to the matrix
		Matrix& Trans();

		//!\brief	inverse the matrix
		//!\return	reference to the matrix
		Matrix& Inv();
		
		//!\brief	calculate the eigenvalue and egienvectors
		//!\param	eigvalue eigenvalue
		//!\param	eigvector eigenvector
		//!\return	void
		void Eig(FVECTOR& eigvalue, Matrix& eigvector);

		//!\brief	multiply a matrix
		//!\param	mat another matrix
		//!\param	result reuslt matrix
		//!\return	reuslt matrix
		Matrix& Multiply(Matrix& mat, Matrix& result);

		//!\brief	left multiply a vector
		//!\param	vec vector
		//!\param	result reuslt vector
		//!\return	reuslt vector
		FVECTOR& LeftMultiply(FVECTOR& vec, FVECTOR& result);

		//!\brief	right multiply a vector
		//!\param	vec vector
		//!\param	result reuslt vector
		//!\return	reuslt vector
		FVECTOR& RightMultiply(FVECTOR& vec, FVECTOR& result);

		//!\brief	divide a scalar
		//!\param	sca scalar
		//!\return	reference to the matrix
		Matrix& Divide(double sca);

		//!\brief	get the mean of all columns
		//!\param	mean the mean vector
		//!\return	mean vector
		FVECTOR& ColMean(FVECTOR& mean);

		//!\brief	get the mean of all rows
		//!\param	mean the mean vector
		//!\return	mean vector
		FVECTOR& RowMean(FVECTOR& mean);

		//!\brief	get standard variation of all columns
		//!\param	std the std vector
		//!\return	std vector
		FVECTOR& ColStd(FVECTOR& std);

		//!\brief	get standard variation of all rows
		//!\param	std the std vector
		//!\return	std vector
		FVECTOR& RowStd(FVECTOR& std);

		//!\brief	subtract a row vector
		//!\param	value row vector
		//!\return	reference to the matrix
		Matrix& RowSub(FVECTOR& value);

		//!\brief	subtract a column vector
		//!\param	value row vector
		//!\return	reference to the matrix
		Matrix& ColSub(FVECTOR& value);

		//!\brief	read a matrix
		//!\param	is input stream
		//!\param	mat matrix
		//!\return	reference to input stream
		friend std::istream& operator>> (std::istream& is, Matrix& mat);

		//!\brief	write a matrix
		//!\param	os output stream
		//!\param	mat matrix
		//!\return	reference to output stream
		friend std::ostream& operator<< (std::ostream& os, Matrix& mat);

		//!\brief	solve A X = b(Numerical Recipes in C++ pp.50-51)
		//!\param	mat LU docomposition of A
		//!\param	indx input vector that records the row permutation by LUdcmp
		//!\param	b right hand of equation
		//!\return	void
		friend void LUbksb(Matrix& mat, std::vector<unsigned int>& indx, std::vector<double>& b);

		//!\brief	LU decompostion of a rowwise permutation(Numerical Recipes in C++ pp.49-50)
		//!\param	mat input and output matrix
		//!\param	indx output vector that records the row permutation
		//!\param	d +-1 depending on whether the number of row interchanges was even or odd
		//!\return	void
		friend void LUdcmp(Matrix& mat, std::vector<unsigned int>& indx, double& d);
	protected:
		//!\brief Householder reduction of Matrix a to tridiagonal form.
		//!
		//! Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
		//! Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
		//! Springer-Verlag, 1976, pp. 489-494.
		//! W H Press et al., Numerical Recipes in C, Cambridge U P,
		//! 1988, pp. 373-374. 
		//!\param	eigenvalue eigenvalue
		//!\param	interm temporal variable
		//!\param	eigenvector eigenvector
		//!\return	void
		void tred2(FVECTOR& eigenvalue, FVECTOR& interm, Matrix& eigenvector);

		//!\brief Tridiagonal QL algorithm -- Implicit 
		//!\param	eigenvalue eigenvalue
		//!\param	interm temporal variable
		//!\param	eigenvector eigenvector
		//!\return	void
		void tqli(FVECTOR& eigenvalue, FVECTOR& interm, Matrix& eigenvector);

		//!\brief	sort the eigenvalue by decreasing order
		//!\param	eigenvalue eigenvalue
		//!\param	eigenvector eigenvector
		//!\return	void
		void Sort(FVECTOR& eigenvalue, Matrix& eigenvector);
	};//class Matrix

	//!\brief linear algebra functions

	//!\brief	cholesky factorization of A: L*L'
	//!\param	L factorization matrix(output)
	//!\param	A a square matrix(input)
	//!\return	success or not
	bool Cholesky(Matrix&L, Matrix&A);

	//!\brief 	Solve a linear system A*X = B, using cholesky factorization of A: L*L'
	//!\param	X a matrix so that L*L'*X = B(output)
	//!\param	A a square matrix(input)
	//!\param	B righthand matrix(input)
	//!\return	success or not
	bool CholeskySolve(Matrix& X, Matrix& A, Matrix& B);

	//!\brief	For an m-by-n matrix A with m >= n, so that A = U*S*V'.
	//!\param	U m-by-n orthogonal matrix(output)
	//!\param	S n-by-n diagonal matrix(output)
	//!\param	V n-by-n orthogonal matrix V(output)
	//!\param	A m-by-n matrix(input)
	//!\param	no
	void SVD(Matrix& U, Matrix&S, Matrix&V, Matrix& A);

	//!\brief	find Pseudo inverse matrix by SVD
	//!\param	inA inverse A(output)
	//!\param	A m-by-n matrix(input)
	//!\param	no
	void pinv(Matrix& inA, Matrix& A);

} //namespace alg

} //namespace az

#endif //AZ_MATRIX_H
