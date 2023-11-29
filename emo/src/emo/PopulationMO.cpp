//PopulationMO.cpp

#include <fstream>
#include "emo/PopulationMO.h"

namespace az
{
namespace mea
{

//constructor
CPopulationMO::CPopulationMO(CParameter& par)
	:mbSort(true),pPar(&par)
{
	mvPop.empty();
}

//constructor
CPopulationMO::CPopulationMO(CPopulationMO& pop)
{
	*this = pop;
}

//deconstrutor
CPopulationMO::~CPopulationMO()
{
	Clear();
}

//set parameters
void CPopulationMO::P(CParameter& par)
{
	pPar = &par;
}

//see whether the individual is in the population
bool CPopulationMO::IsContain(IND_TYPE& ind)
{
	for(unsigned int i=0; i<Size(); i++) if(ind == In(i)) return true;
	return false;
}

//clear the population
void CPopulationMO::Clear()
{
	CLEAR(mvPop);
	mbSort = true;
}

//evaluate the population
void CPopulationMO::Evaluate()
{
	for(unsigned int i =0; i<Size(); i++) In(i).Evaluate();
	mbSort = false;
}

//shuffle the population
void CPopulationMO::Shuffle()
{
	for(unsigned int i=0; i<Size()-1; i++) Swap(i, rnd::rand(i, Size()));
	mbSort = false;
}

//resize the population
void CPopulationMO::Resize(unsigned int s)
{
	Clear();
	mvPop.resize(s);
	for(unsigned int i =0; i<s; i++) mvPop[i] = new IND_TYPE(P());
	mbSort = false;
}

//swap two individuals
void CPopulationMO::Swap(unsigned int i, unsigned int j)
{
	std::swap(mvPop[i], mvPop[j]);
	mbSort = false;
}

//erase individuals
void CPopulationMO::Erase(unsigned int i)
{
	for(unsigned int k=i; k<Size(); k++) if(mvPop[k]>0) delete mvPop[k];
	mvPop.erase(mvPop.begin()+i, mvPop.end());
}

//copy an individual to population
CPopulationMO& CPopulationMO::Copy(IND_TYPE*& pind)
{
	if(pind>0)
	{
		mvPop.push_back(new IND_TYPE(*pind));
		mbSort = false;
	}
	return *this;
}

//copy an individual to population
CPopulationMO& CPopulationMO::Copy(IND_TYPE& ind)
{
	mvPop.push_back(new IND_TYPE(ind));
	mbSort = false;
	return *this;
}

//copy a population to population
CPopulationMO& CPopulationMO::Copy(CPopulationMO& pop)
{
	for(unsigned int i=0; i<pop.Size(); i++) Copy(pop[i]);
	return *this;
}

//combine an individual to population
CPopulationMO& CPopulationMO::Combine(IND_TYPE*& pind)
{
	if(pind>0)
	{
		mvPop.push_back(pind);
		pind = 0;
		mbSort = false;
	}
	return *this;
}

//combine an individual to population
CPopulationMO& CPopulationMO::Combine(IND_TYPE& ind)
{
	mvPop.push_back(new IND_TYPE(ind));
	mbSort = false;
	return *this;
}

//combine a population to population
CPopulationMO& CPopulationMO::Combine(CPopulationMO& pop)
{
	for(unsigned int i=0; i<pop.Size(); i++) Combine(pop.At(i));
	pop.Clear();
	return *this;
}

//find the union of two sets (populations)
CPopulationMO& CPopulationMO::Unite(CPopulationMO& popA, CPopulationMO& popB)
{
	*this = popA;
	for(unsigned int i=0; i<popB.Size(); i++) if(!popA.IsContain(popB[i])) Copy(popB[i]);
	return *this;
}

CPopulationMO& CPopulationMO::Minus(CPopulationMO& popA, CPopulationMO& popB)
{
	Clear();
	for(unsigned int i=0; i<popA.Size(); i++) if(!popB.IsContain(popA[i])) Copy(popA[i]);
	return *this;
}

//assign a population 
CPopulationMO& CPopulationMO::operator=(CPopulationMO& pop)
{
	Clear();
	mbSort	= pop.mbSort;
	pPar	= pop.pPar;
	Resize(pop.Size());
	for(unsigned int i=0; i<Size(); i++) In(i) = pop[i];
	return *this;
}

//get a sub-population with Rank = r
CPopulationMO& CPopulationMO::RankSub(CPopulationMO& pop, unsigned int r)
{
	unsigned int i,start,end;
	RankSize(r, start, end);
	pop.Clear();
	if(start>Size()) return pop;
	for(i=start; i<end; i++)
		pop.Copy(In(i));
	pop.IsSort(true);
	return pop;
}

//get the maximum rank value
unsigned int CPopulationMO::RankMax()
{
	RankSort();
	return In(Size()-1).Rank();
}

//get the sub-population size with Rank=r
unsigned int CPopulationMO::RankSize(unsigned int r, unsigned int& start, unsigned int& end)
{
	RankSort();
	start = 0;
	while(start<Size() && In(start).Rank()<r) start++;
	end = start;
	while(end<Size() && In(end).Rank() == r) end++;
	if(start>=Size()) return 0;
	else return end - start;
}

//assign rank value and sort population
void CPopulationMO::RankSort()
{
	RankSort(false);
}

void CPopulationMO::ERankSort()
{
	RankSort(true);
}

void CPopulationMO::RankSort(bool esort)
{
	//has been sorted before
	if(IsSort()) return;

	//only has one individual
	if(Size()<2) {IsSort(true);return;}

	int s, t, better; unsigned int i, j, minRank, curRank = 0, noAssign, size;

	//Step 0: Move all infeasible solutions to the tail of the sequence
	s=0; t=Size()-1;
	while(s<t)
	{
		while(s<int(Size()) && In(s).IsFeasible()) s++;
		while(t>=0 && !In(t).IsFeasible()) t--;
		if(s<t)	{Swap(s,t); s++; t--;}
	}
	size = s;

	std::vector<bool>						vAssign(size);	
	std::list<unsigned int>					v2Assign;
	std::vector<unsigned int>				vBeDom(size);
	std::vector< std::list<unsigned int> >	vDom(size); 
	std::list<unsigned int>::iterator it1,it2;	
	std::vector<double>						vE(P().FSize()),vFMax(P().FSize()),vFMin(P().FSize());

	if(esort)
	{
		for(i=0; i<P().FSize(); i++)
		{
			vFMax[i] = -1.0E100; 
			vFMin[i] =  1.0E100;
		}
		for(i=0; i<Size(); i++)
		{
			for(j=0; j<P().FSize(); j++)
			{
				if(In(i).F(j)<vFMin[j]) vFMin[j] = In(i).F(j);
				if(In(i).F(j)>vFMax[j]) vFMax[j] = In(i).F(j);
			}
		}
		for(i=0; i<P().FSize(); i++) vE[i] = 0.001*(vFMax[i]-vFMin[i]);
	}

	//Step 1: Initialize
	for(i=0; i<size; i++) vAssign[i] = false;

	//Step 2: MOGAFonseca flow
	for(i=0; i<size; i++)
	{
		for (j=i+1 ; j<size; j++)
		{
			better = esort ? In(i).CEDominate(In(j),vE) : In(i).CDominate(In(j)) ;
			//i is dominated by j
			if(better<0) { vBeDom[i]++; vDom[j].push_back(i); }
			//j is dominated by i
			else if(better>0) { vBeDom[j]++; vDom[i].push_back(j); }
		}	//end for
	}	// end for

	//Step 3: Assign rank to feasible solutions
	noAssign = size;
	while(noAssign > 0)
	{
		curRank++;

		minRank = size;
		//find the cluster to assign a rank value
		for (i=0; i<size; i++) if(!vAssign[i])
		{
			if(vBeDom[i]<minRank) 
			{
				minRank = vBeDom[i];
				v2Assign.clear(); v2Assign.push_back(i);
			}
			else if(vBeDom[i]==minRank) v2Assign.push_back(i);
		}

		//CHECK( v2Assign.size()>0, "CRankSort::Sort()" );

		//assign rank
		it1 = v2Assign.begin();
		while(it1!=v2Assign.end())
		{
			vAssign[*it1] = true; 
			In(*it1).Rank(curRank);
			it2 = vDom[*it1].begin(); 
			while(it2!=vDom[*it1].end()) vBeDom[*it2++]--;
			it1++; noAssign--;
		}//for
	}//end while

	//Step 4: Assign rank to infeasible solutions
	curRank++; for(i=size; i<Size(); i++) In(i).Rank(curRank);

	//Step 5: Sort the population by rank and constraint
	for(s=0; s<int(Size()-1); s++) for(t=s+1; t<int(Size()); t++) if(In(t)<In(s)) Swap(s,t);

	IsSort(true);
}

unsigned int CPopulationMO::DominateSort()
{
	int s, t, better; unsigned int i, j;

	std::vector<bool> state(Size()); for(i=0; i<Size(); i++) state[i] = true;
	
	// find the domination state
	for(i=0; i<Size(); i++)
	{
		for(j=i+1; j<Size(); j++)
		{
			better = In(i).CDominate(In(j));
			if(better<0)		state[i] = false;
			else if(better>0)	state[j] = false;
		}
	}
	
	// assign rank
	for(i=0; i<Size(); i++) In(i).Rank(state[i] ? 1 : Size());

	// sort the population by rank
	for(s=0; s<int(Size()-1); s++) for(t=s+1; t<int(Size()); t++) if(In(t)<In(s)) Swap(s,t);

	unsigned int num = 0; for(i=0; i<Size(); i++) if(state[i]) num++;

	return num;
}

//write the population to I/O stream
std::ostream& CPopulationMO::Write(std::ostream& os)
{
	//os<<std::scientific<<std::setprecision(10);
	//os<<Size()<<std::endl<<P().FSize()<<"\t"<<P().XSize()<<std::endl;
	//for(unsigned int i=0; i<Size(); i++) os<<In(i);
	//return os;

	unsigned int i,s = DominateSort();
	os<<std::scientific<<std::setprecision(10);
	os<<s<<"\t"<<Size()-s<<"\t"<<P().FSize()<<"\t"<<P().XSize()<<std::endl;
	for( i=0; i<Size(); i++ ) os<<In(i);
	return os;
}

//write the population to a file
void CPopulationMO::Write(std::string name)
{
	std::ofstream file(name.c_str());	
	file<<*this;	
	file.close();
}

//write the population to a file
void CPopulationMO::Write(const char *name)
{
	std::ofstream file(name);	
	file<<*this;	
	file.close();
}

//write the population to I/O stream
std::ostream& operator<<(std::ostream& os, CPopulationMO& pop)
{
	return pop.Write(os);
}

//read population from I/O stream
std::istream& CPopulationMO::Read(std::istream& is)
{
	unsigned int size, size1, size2, x, y;
	is>>size1>>size2>>y>>x; size = size1+size2;
	(*this).Resize(size);
	for(unsigned int i=0; i<Size(); i++) is >> In(i);
	return is;
}

//read population from a file
void CPopulationMO::Read(std::string name)
{	
	std::ifstream file(name.c_str());	
	if(file != 0)	file>>*this;	
	file.close();
}

//read population from a file
void CPopulationMO::Read(const char *name)
{
	std::ifstream file(name);	
	if(file != 0)	file>>*this;	
	file.close();
}

//read population from I/O stream
std::istream& operator>>(std::istream& is, CPopulationMO& pop)
{
	return pop.Read(is);
}

} //namespace mea
} //namespace az
