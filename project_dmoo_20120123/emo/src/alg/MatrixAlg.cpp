// MatrixAlgebra.cpp

#include <cmath>
#include <vector>
#include "alg/Matrix.h"

namespace az
{

namespace alg
{

	//!\brief	cholesky factorization of A1: L*L'
	// this algorithm comes from the version of TNT(jama_cholesky.h)
	bool Cholesky(Matrix&L, Matrix&A1)
	{
		int i,j,k;
		double d,s;

		if(A1.RowSize() != A1.ColSize()) return false;

		bool isspd=true;

		int n	= A1.RowSize();
		L.Resize(n,n);

		//Main loop.
		for(j=0; j<n; j++) 
		{
			d=0.0;
			for(k=0; k<j; k++) 
			{
				s=0.0;
				for(i=0; i<k; i++) s += L(k,i)*L(j,i);
				L(j,k)=s =(A1(j,k)-s)/L(k,k);
				d=d + s*s;
				isspd=isspd &&(A1(k,j) == A1(j,k)); 
			}
			d=A1(j,j) - d;
			isspd=isspd &&(d > 0.0);
			L(j,j)=sqrt(d > 0.0 ? d : 0.0);
			for(k=j+1; k<n; k++) L(j,k)=0.0;
		}
		return isspd;
	}

	//!\brief 	Solve a linear system A1*X=B, using cholesky factorization of A1: L*L'
	// this algorithm comes from the version of TNT(jama_cholesky.h)
	bool CholeskySolve(Matrix& X, Matrix& A1, Matrix& B)
	{
		int i,j,k;
		Matrix L_;

		X.Resize(0,0);
		if(A1.RowSize() != A1.ColSize() || A1.RowSize() != B.RowSize() || !Cholesky(L_,A1)) return false;

		int n	= A1.RowSize();
		int nx	= B.ColSize();

		//Step 1: Solve L*Y=B;
		X=B;
		for(j=0; j<nx; j++)
		{
			for(k=0; k<n; k++) 
			{
				for(i=0; i<k; i++) 
					X(k,j) -= X(i,j)*L_(k,i);
				X(k,j) /= L_(k,k);
			}
		}

		//Step 2: Solve L'*X=Y;
		for(j=0; j<nx; j++)
		{
			for(k=n-1; k>=0; k--) 
			{
				for(i=k+1; i<n; i++) 
					X(k,j) -= X(i,j)*L_(i,k);
				X(k,j) /= L_(k,k);
			}
		}

		return true;
	}

	#define MAX(a,b) ((a)<(b) ? (b):(a))
	#define MIN(a,b) ((a)<(b) ? (a):(b))
	double hypot(const double &a, const double &b)
	{
		if(a== 0) return fabs(b);
		else
		{
			double c=b/a;
			return fabs(a) * sqrt(1 + c*c);
		}
	}

	// an m-by-n matrix A1 with m>=n => A1=U*S*V'.
	// this algorithm comes from the version of TNT(jama_svd.h)
	void SVD(Matrix& U, Matrix&S, Matrix&V, Matrix& A)
	{
		int m =A.RowSize();
		int	n =A.ColSize();
		int nu=MIN(m,n);
		U.Resize(m, nu);
		V.Resize(n, n);
		S.Resize(MIN(m+1,n),MIN(m+1,n));
		std::vector<double> s(MIN(m+1,n)), e(n), work(m);
		Matrix A1(A);
		int wantu=1;		/* boolean */
		int wantv=1;		/* boolean */
		int i=0, j=0, k=0;

		// Reduce A1 to bidiagonal form, storing the diagonal elements
		// in s and the super-diagonal elements in e.

		int nct=MIN(m-1,n);
		int nrt=MAX(0,MIN(n-2,m));
		for(k=0; k < MAX(nct,nrt); k++) 
		{
			if(k<nct) 
			{
				// Compute the transformation for the k-th column and
				// place the k-th diagonal in s[k].
				// Compute 2-norm of k-th column without under/overflow.
				s[k]=0;
				for(i=k; i<m; i++) 
				{
					s[k]=hypot(s[k],A1(i,k));
				}
				if(s[k] != 0.0)
				{
					if(A1(k,k)<0.0) 
					{
						s[k]=-s[k];
					}
					for(i=k; i<m; i++) 
					{
						A1(i,k) /= s[k];
					}
					A1(k,k) += 1.0;
				}
				s[k]=-s[k];
			}

			for(j=k+1; j<n; j++) 
			{
				if((k<nct) &&(s[k] != 0.0))  
				{
					// Apply the transformation.
					double t=0;
					for(i=k; i<m; i++) 
					{
						t += A1(i,k)*A1(i,j);
					}
					t=-t/A1(k,k);
					for(i=k; i<m; i++) 
					{
						A1(i,j) += t*A1(i,k);
					}
				}
				// Place the k-th row of A1 into e for the
				// subsequent calculation of the row transformation.
				e[j]=A1(k,j);
			}
			if(wantu &(k<nct)) 
			{
				// Place the transformation in U for subsequent back
				// multiplication.
				for(i=k; i<m; i++) 
				{
					U(i,k)=A1(i,k);
				}
			}
			if(k<nrt) 
			{
				// Compute the k-th row transformation and place the
				// k-th super-diagonal in e[k].
				// Compute 2-norm without under/overflow.
				e[k]=0;
				for(i=k+1; i<n; i++) 
				{
					e[k]=hypot(e[k],e[i]);
				}
				if(e[k] != 0.0) 
				{
					if(e[k+1]<0.0) 
					{
						e[k]=-e[k];
					}
					for(i=k+1; i<n; i++) 
					{
						e[i] /= e[k];
					}
					e[k+1] += 1.0;
				}
				e[k]=-e[k];
				if((k+1<m) &(e[k] != 0.0)) 
				{
					// Apply the transformation.
					for(i=k+1; i<m; i++) 
					{
						work[i]=0.0;
					}
					for(j=k+1; j<n; j++) 
					{
						for(i=k+1; i<m; i++) 
						{
							work[i] += e[j]*A1(i,j);
						}
					}
					for(j=k+1; j<n; j++) 
					{
						double t=-e[j]/e[k+1];
						for(i=k+1; i<m; i++) 
						{
							A1(i,j) += t*work[i];
						}
					}
				}
				if(wantv) 
				{
					// Place the transformation in V for subsequent
					// back multiplication.
					for(i=k+1; i<n; i++) 
					{
						V(i,k)=e[i];
					}
				}
			}
		}

		// Set up the final bidiagonal matrix or order p.
		int p=MIN(n,m+1);
		if(nct<n) 
		{
			s[nct]=A1(nct,nct);
		}
		if(m<p) 
		{
			s[p-1]=0.0;
		}
		if(nrt+1<p) 
		{
			e[nrt]=A1(nrt,p-1);
		}
		e[p-1]=0.0;

		// If required, generate U.
		if(wantu) 
		{
			for(j=nct; j<nu; j++) 
			{
				for(i=0; i<m; i++) 
				{
					U(i,j)=0.0;
				}
				U(j,j)=1.0;
			}
			for(k=nct-1; k>=0; k--) 
			{
				if(s[k] != 0.0) 
				{
					for(j=k+1; j<nu; j++) 
					{
						double t=0;
						for(i=k; i<m; i++) 
						{
							t += U(i,k)*U(i,j);
						}
						t=-t/U(k,k);
						for(i=k; i<m; i++) 
						{
							U(i,j) += t*U(i,k);
						}
					}
					for(i=k; i<m; i++ ) 
					{
						U(i,k)=-U(i,k);
					}
					U(k,k)=1.0 + U(k,k);
					for(i=0; i<k-1; i++) 
					{
						U(i,k)=0.0;
					}
				} 
				else 
				{
					for(i=0; i<m; i++) 
					{
						U(i,k)=0.0;
					}
					U(k,k)=1.0;
				}
			}
		}

		// If required, generate V.
		if(wantv) 
		{
			for(k=n-1; k>=0; k--) 
			{
				if((k<nrt) &(e[k] != 0.0)) 
				{
					for(j=k+1; j<nu; j++) 
					{
						double t=0;
						for(i=k+1; i<n; i++) 
						{
							t += V(i,k)*V(i,j);
						}
						t=-t/V(k+1,k);
						for(i=k+1; i<n; i++) 
						{
							V(i,j) += t*V(i,k);
						}
					}
				}
				for(i=0; i<n; i++) 
				{
					V(i,k)=0.0;
				}
				V(k,k)=1.0;
			}
		}

		// Main iteration loop for the singular values.
		int pp=p-1;
		int iter=0;
		double eps=pow(2.0,-52.0);
		while(p > 0) 
		{
			int k=0;
			int kase=0;
			// Here is where a test for too many iterations would go.

			// This section of the program inspects for
			// negligible elements in the s and e arrays.  On
			// completion the variables kase and k are set as follows.

			// kase=1     if s(p) and e[k-1] are negligible and k<p
			// kase=2     if s(k) is negligible and k<p
			// kase=3     if e[k-1] is negligible, k<p, and
			//              s(k), ..., s(p) are not negligible(qr step).
			// kase=4     if e(p-1) is negligible(convergence).

			for(k=p-2; k>=-1; k--) 
			{
				if(k == -1) 
				{
					break;
				}
				if(fabs(e[k]) <= eps*(fabs(s[k]) + fabs(s[k+1]))) 
				{
					e[k]=0.0;
					break;
				}
			}
			if(k == p-2) 
			{
				kase=4;
			}
			else 
			{
				int ks;
				for(ks=p-1; ks>=k; ks--) 
				{
					if(ks == k) 
					{
						break;
					}
					double t =(ks != p ? fabs(e[ks]) : 0.) + 
						(ks != k+1 ? fabs(e[ks-1]) : 0.);
					if(fabs(s[ks]) <= eps*t)  
					{
						s[ks]=0.0;
						break;
					}
				}
				if(ks == k) 
				{
					kase=3;
				} 
				else if(ks == p-1) 
				{
					kase=1;
				} 
				else 
				{
					kase=2;
					k=ks;
				}
			}
			k++;

			// Perform the task indicated by kase.

			switch(kase) 
			{
				// Deflate negligible s(p).
			case 1: 
				{
					double f=e[p-2];
					e[p-2]=0.0;
					for(j=p-2; j>=k; j--) 
					{
						double t=hypot(s[j],f);
						double cs=s[j]/t;
						double sn=f/t;
						s[j]=t;
						if(j != k) 
						{
							f=-sn*e[j-1];
							e[j-1]=cs*e[j-1];
						}
						if(wantv) 
						{
							for(i=0; i<n; i++) 
							{
								t=cs*V(i,j) + sn*V(i,p-1);
								V(i,p-1)=-sn*V(i,j) + cs*V(i,p-1);
								V(i,j)=t;
							}
						}
					}
				}
				break;

				// Split at negligible s(k).
			case 2: 
				{
					double f=e[k-1];
					e[k-1]=0.0;
					for(j=k; j<p; j++)
					{
						double t=hypot(s[j],f);
						double cs=s[j]/t;
						double sn=f/t;
						s[j]=t;
						f=-sn*e[j];
						e[j]=cs*e[j];
						if(wantu) 
						{
							for(i=0; i<m; i++) 
							{
								t=cs*U(i,j) + sn*U(i,k-1);
								U(i,k-1)=-sn*U(i,j) + cs*U(i,k-1);
								U(i,j)=t;
							}
						}
					}
				}
				break;
				
				// Perform one qr step.
			case 3: 
				{
					// Calculate the shift.
					double scale=MAX(MAX(MAX(MAX(
						fabs(s[p-1]),fabs(s[p-2])),fabs(e[p-2])), 
						fabs(s[k])),fabs(e[k]));
					double sp=s[p-1]/scale;
					double spm1=s[p-2]/scale;
					double epm1=e[p-2]/scale;
					double sk=s[k]/scale;
					double ek=e[k]/scale;
					double b =((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
					double c =(sp*epm1)*(sp*epm1);
					double shift=0.0;
					if((b != 0.0) ||(c != 0.0)) 
					{
						shift=sqrt(b*b + c);
						if(b<0.0) 
						{
							shift=-shift;
						}
						shift=c/(b + shift);
					}
					double f =(sk + sp)*(sk - sp) + shift;
					double g=sk*ek;

					// Chase zeros.

					for(j=k; j<p-1; j++) 
					{
						double t=hypot(f,g);
						double cs=f/t;
						double sn=g/t;
						if(j != k) 
						{
							e[j-1]=t;
						}
						f=cs*s[j] + sn*e[j];
						e[j]=cs*e[j] - sn*s[j];
						g=sn*s[j+1];
						s[j+1]=cs*s[j+1];
						if(wantv) 
						{
							for(i=0; i<n; i++) 
							{
								t=cs*V(i,j) + sn*V(i,j+1);
								V(i,j+1)=-sn*V(i,j) + cs*V(i,j+1);
								V(i,j)=t;
							}
						}
						t=hypot(f,g);
						cs=f/t;
						sn=g/t;
						s[j]=t;
						f=cs*e[j] + sn*s[j+1];
						s[j+1]=-sn*e[j] + cs*s[j+1];
						g=sn*e[j+1];
						e[j+1]=cs*e[j+1];
						if(wantu &&(j<m-1)) 
						{
							for(i=0; i<m; i++) 
							{
								t=cs*U(i,j) + sn*U(i,j+1);
								U(i,j+1)=-sn*U(i,j) + cs*U(i,j+1);
								U(i,j)=t;
							}
						}
					}
					e[p-2]=f;
					iter=iter + 1;
				}
				break;

				// Convergence.
			case 4: 
				{
					// Make the singular values positive.
					if(s[k] <= 0.0) 
					{
						s[k] =(s[k]<0.0 ? -s[k] : 0.0);
						if(wantv) 
						{
							for(i=0; i <= pp; i++) 
							{
								V(i,k)=-V(i,k);
							}
						}
					}
					// Order the singular values.
					while(k<pp) 
					{
						if(s[k]>=s[k+1]) 
						{
							break;
						}
						double t=s[k];
						s[k]=s[k+1];
						s[k+1]=t;
						if(wantv &&(k<n-1)) 
						{
							for(i=0; i<n; i++) 
							{
								t=V(i,k+1); V(i,k+1)=V(i,k); V(i,k)=t;
							}
						}
						if(wantu &&(k<m-1)) 
						{
							for(i=0; i<m; i++) 
							{
								t=U(i,k+1); U(i,k+1)=U(i,k); U(i,k)=t;
							}
						}
						k++;
					}
					iter=0;
					p--;
				}
				break;
			}
		}

		//set S
		for(i=0; i<int(s.size()); i++) 
		{
			S(i,i) = s[i];
		}
	}

	// find Pseudo inverse matrix by SVD
	void pinv(Matrix& inA, Matrix& A)
	{
		unsigned int i;
		Matrix S,U,V;
		
		SVD(U,S,V,A);
		
		Matrix S1(S);
		for(i=0; i<S1.ColSize(); i++) S1(i,i)=1.0/S1(i,i);
		
		Matrix trU(U); trU.Trans();

		Matrix tmp;
		V.Multiply(S1,tmp);
		tmp.Multiply(trU,inA);
	}

} //namespace alg

} //namespace az
