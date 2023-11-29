//IndividualMO.cpp

#include <float.h>
#include "emo/IndividualMO.h"

#if defined(WIN32)
    #define wxFinite(n) _finite(n)
#elif defined(_LINUX)
    #define wxFinite(n) finite(n)
#else
    #define wxFinite(n) ((n)==(n))
#endif
namespace az
{
namespace mea
{

//Constructor, assign space to parameters
CIndividualMO::CIndividualMO(CParameter& par)
{
	mConstraint		= 0.0;
	pPar			= &par;
	mRank			= 0;
	mID				= 0;
	mOpt			= 0;
	mvX.resize(pPar->XSize());
	mvF.resize(pPar->FSize());
	mvEq.resize(pPar->ESize());
	mvIneq.resize(pPar->ISize());
}

//Constructor
CIndividualMO::CIndividualMO(CIndividualMO& ind)	
{
	*this = ind;
}

//dominance check
int CIndividualMO::Dominate(CIndividualMO& ind)
{
	unsigned int better = 0, wbetter = 0, worse = 0, wworse = 0;
	for(unsigned int i=0; i<P().FSize(); i++)
	{
		if(	F(i) <= ind.F(i) + P().TolF())
		{
			better++;
			if(F(i) > ind.F(i) - P().TolF()) wbetter ++;
		}
		
		if(ind.F(i) <= F(i) + P().TolF())
		{
			worse++;
			if(ind.F(i) > F(i) - P().TolF()) wworse ++;
		}
	}

	if(better == P().FSize() && wbetter < better)	return  1;
	else if(worse == P().FSize() && wworse < worse)	return -1;
	else return 0;
}

//e-dominance check
int CIndividualMO::CEDominate(CIndividualMO& ind, std::vector<double>& e)
{
	if(mConstraint < ind.mConstraint) return 1;
	if(mConstraint > ind.mConstraint) return -1;

	unsigned int better = 0, sbetter = 0, worse = 0, sworse = 0;
	for(unsigned int i=0; i<P().FSize(); i++)
	{
		if(	F(i) <= ind.F(i) )
		{
			better++;
			if(F(i) < ind.F(i) - e[i]) sbetter ++;
		}
		
		if(ind.F(i) <= F(i))
		{
			worse++;
			if(ind.F(i) > F(i) - e[i]) sworse ++;
		}
	}

	if(better == P().FSize() && sbetter > 0)	return  1;
	else if(worse == P().FSize() && sworse > 0)	return -1;
	else return 0;
}
//constratint-dominance check
int CIndividualMO::CDominate(CIndividualMO& ind)
{
	if(mConstraint < ind.mConstraint) return 1;
	if(mConstraint > ind.mConstraint) return -1;
	return Dominate(ind);
}

//Evaluate the individual
void CIndividualMO::Evaluate()
{
	unsigned int i;

	//evaluate
	// with coding
	if(P().XCoding())
	{
		std::vector<double> XX(mvX.size());
		for(i=0; i<mvX.size(); i++) XX[i] = P().XRealLow(i) + (mvX[i]-P().XLow(i))/(P().XUpp(i)-P().XLow(i))*(P().XRealUpp(i)-P().XRealLow(i));
		(*P().Evaluator())(mvF, mvEq, mvIneq, XX);
		XX.clear();
	}
	// without coding
	else
	{
		(*P().Evaluator())(mvF, mvEq, mvIneq, mvX);
	}

	for(i=0; i<P().FSize(); i++) if( wxFinite(mvF[i]) == 0 ) mvF[i] = 1.0E20;

	//Evaluate equalities
	mConstraint = 0.0;
	for(i=0; i<P().ESize(); i++) {if( wxFinite(mvEq[i]) == 0 ) mvEq[i] = 1.0E20; mConstraint += fabs(mvEq[i]);}

	//Evaluate inequalities
	for(i=0; i<P().ISize(); i++) {if( wxFinite(mvIneq[i]) == 0 ) mvIneq[i] = 1.0E20; mConstraint += mvIneq[i] > 0.0 ? mvIneq[i] : 0.0;}
}

//Check the bounds
void CIndividualMO::Check()
{
	for(unsigned int i=0; i<P().XSize(); i++)
	{
		if(X(i) < P().XLow(i)) X(i) = rnd::rand(P().XLow(i), (P().XLow(i)+P().XUpp(i))/2.0);
		else if(X(i) > P().XUpp(i)) X(i) = rnd::rand((P().XLow(i)+P().XUpp(i))/2.0, P().XUpp(i));
	}
		//if(X(i) < P().XLow(i))		X(i) = P().XLow(i);
		//else if(X(i) > P().XUpp(i)) X(i) = P().XUpp(i);
}

CIndividualMO& CIndividualMO::operator= (CIndividualMO& ind)
{
	mRank	= ind.mRank;
	mvX		= ind.mvX;
	mvF		= ind.mvF;
	mvEq	= ind.mvEq;
	mvIneq	= ind.mvIneq;
	pPar	= ind.pPar;
	mID		= ind.mID;
	mOpt	= ind.mOpt;
	mConstraint	= ind.mConstraint;

	return *this;
}

bool CIndividualMO::operator==(CIndividualMO& ind)
{
	for(unsigned int i=0; i<P().XSize(); i++)
		if(fabs(X(i) - ind.X(i)) > P().TolX()) return false;
	return true;
}

//check to see which is better
bool CIndividualMO::operator<(CIndividualMO& ind)
{
	return (mRank < ind.mRank) || (mRank == ind.mRank && mConstraint < ind.mConstraint);
}

//Write to stream
std::ostream& operator<<(std::ostream& os, CIndividualMO& ind)
{
	unsigned int i;
	for(i=0; i<ind.P().FSize(); i++)		os<<ind.F(i)<<"\t";
	if(ind.P().XCoding()) // coding
		for(i=0; i<ind.P().XSize(); i++)	/*os<<ind.X(i)<<"\t";*/os<<ind.P().XRealLow(i) + (ind.X(i)-ind.P().XLow(i))/(ind.P().XUpp(i)-ind.P().XLow(i))*(ind.P().XRealUpp(i)-ind.P().XRealLow(i))<<"\t";
	else	// without coding
		for(i=0; i<ind.P().XSize(); i++)	os<<ind.X(i)<<"\t";
	os<<std::endl;
	return os;
}

//Read from stream
std::istream& operator>>(std::istream& is, CIndividualMO& ind)
{
	unsigned int i;
	for(i=0; i<ind.P().FSize(); i++)	is>>ind.F(i);
	for(i=0; i<ind.P().XSize(); i++)	is>>ind.X(i);
	return is;
}

} //namespace mea
} //namespace az
