//Random.cpp

#include <ctime>
#include <cmath>

namespace az
{
namespace rnd
{
	#define IM1		2147483563L		//!< const value
	#define IM2		2147483399L		//!< const value
	#define AM		( 1.0 / IM1 )	//!< const value
	#define IMM1	( IM1 - 1 )		//!< const value
	#define IA1		40014L			//!< const value
	#define IA2		40692L			//!< const value
	#define IQ1		53668L			//!< const value
	#define IQ2		52774L			//!< const value
	#define IR1		12211L			//!< const value
	#define IR2		3791			//!< const value
	#define NTAB	32				//!< const value
	#define NDIV	(1+IMM1/NTAB)	//!< const value
	#define EPS		1.2e-7			//!< const value
	#define RNMX	(1.0 - EPS)		//!< const value

	static long idum2	= 123456789;	//!< temporal variable
	static long iy		= 0;			//!< temporal variable
	static long iv[NTAB];				//!< temporal variable
	static long idum	= 0;			//!< temporal variable

	//initialize the random seed
	void seed( long seeds = 0 )
	{
		int	j;
		long	k;

		idum = ( seeds == 0 ) ? ( ( long )time( NULL ) ) : seeds;

		if ( idum == 0 )	idum	= 1;
		if ( idum < 0 )		idum	= -idum;

		idum2	= idum;

		for(	j = NTAB + 7; j >= 0; j-- )
		{
			k	= idum / IQ1;
			idum= IA1 * ( idum - k * IQ1 ) - k * IR1;
			if( idum < 0 ) idum	+= IM1;
			if( j < NTAB ) iv[j] = idum;
		}

		iy=iv[0];
	}

	//create a real random in (0.0,1.0)
	double rand()
	{
		int		j;
		long	k;
		double	temp;

		k			= idum / IQ1;
		idum	= IA1 * ( idum - k * IQ1 ) - k * IR1;
		if( idum < 0 ) idum	+= IM1;
		k			= idum2 / IQ2;
		idum2	= IA2 * ( idum2 - k * IQ2 ) - k * IR2;
		if( idum2 < 0 )	idum2	+= IM2;

		j		= iy / NDIV;
		iy	= iv[j] - idum2;
		iv[j] = idum;
		if( iy < 1 ) iy += IMM1;
		if( ( temp = AM * iy ) > RNMX ) return RNMX;
		else return temp;
	}

	//create a real random in (low,up) for real number and [low, up) for integer number
	double rand(double low, double up)
	{
		return low + (up - low)*rand();
	}

	//create a real random in (low,up) for real number and [low, up) for integer number
	int rand(int low, int up)
	{
		return (int)( double( low ) +double( up - low )*rand() );
	}

	//create a real random in (low,up) for real number and [low, up) for integer number
	unsigned int rand(unsigned int low, unsigned int up)
	{
		return (unsigned int)( double( low ) +double( up - low )*rand() );
	}

	//create a real Gaussian random with distribution (0.0,1.0)
	double gaussian()
	{
		static bool cached=false;
		static double cachevalue;
		if(cached == true)
		{
			cached = false;
			return cachevalue;
		}

		double rsquare, factor, var1, var2;
		do
		{
			var1 = 2.0 * rand() - 1.0;
			var2 = 2.0 * rand() - 1.0;
			rsquare = var1*var1 + var2*var2;
		} while(rsquare >= 1.0 || rsquare == 0.0);

		double val = -2.0 * log(rsquare) / rsquare;
		if(val > 0.0)	
			factor = sqrt(val);
		else	
			factor = 0.0;

		cachevalue	= var1 * factor;
		cached		= true;

		return (var2 * factor);
	}

	// create a real number with Triangular distribution
	double triangular(double min, double max, double mode)
	{
		double uni = rand(0.0, 1.0);
		if(uni <= (mode-min)/(max-min)) 
			return min+sqrt(uni*(max-min)*(mode-min));
		else
			return max-sqrt((1-uni)*(max-min)*(max-mode));
	}
} //namespace rnd
} //namespace az
