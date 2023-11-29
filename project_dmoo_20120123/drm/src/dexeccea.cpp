#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "emo/Config.h"

#include "problem.h"
#include "COEA.h"

#include "stdlib.h"

az::mea::Config	config,	// statistic configure file
				confignew;	// configure for one run;
unsigned int sectionsize;	// section size
std::vector< std::string >	tmpstr(1), sections, items, names;
std::string pathname;

void scan_run( unsigned int index );
void execone();

int main(int argc, char* argv[])
{
	if(argc<2)
	{
		std::cout<<"need configure file"<<std::endl;
		return -1;
	}

	unsigned int i;
	std::string	 infile;	//config file name
	infile		= std::string(argv[1]);
	config.Load(infile, '#');
	confignew	= config;

	config.Get(std::string("COMMONA"), std::string("PATHNAME"),	0, pathname);
	sectionsize = config.GetSize(std::string("COMMONS"),std::string("SECTION"));
	sections.resize(sectionsize); items.resize(sectionsize); names.resize(sectionsize);
	for(i=0; i<sectionsize; i++)
	{
		config.Get(std::string("COMMONS"),std::string("SECTION"),i,sections[i]);
		config.Get(std::string("COMMONS"),std::string("ITEM")   ,i,items[i]);
	}

	scan_run(0);

	return 1;
}

void scan_run(unsigned int index)
{
	unsigned int i,k,isize;
	isize = config.GetSize(sections[index],items[index]);

	for(i = 0; i < isize; i++)
	{
		config.Get(sections[index],items[index],i,tmpstr[0]);
		confignew.Set(sections[index],items[index],tmpstr);
		names[index] = tmpstr[0];

		if(index < sectionsize-1)
		{
			scan_run(index+1);
		}
		else
		{
			tmpstr[0] = pathname + std::string("/")+names[0];
			for(k=1; k<sectionsize; k++) tmpstr[0] += std::string("_")+names[k];
			confignew.Set(std::string("COMMONA"),std::string("PATHNAME"),tmpstr);

			execone();
		}
	}
}

void execone()
{
	std::string method, problem, fname, path;
	unsigned int strategy, runs, generation, popsize, dimension, gen, ir, i, j;
	unsigned int taot, nt; 
	double		 t0;

	confignew.Get(std::string("COMMONA"), std::string("METHOD"),	0, method);
	confignew.Get(std::string("COMMONA"), std::string("RUNS"),		0, runs);
	confignew.Get(std::string("COMMONA"), std::string("POPSIZE"),	0, popsize);
	confignew.Get(std::string("COMMONA"), std::string("PROBLEM"),	0, problem);
	confignew.Get(std::string("COMMONA"), std::string("PATHNAME"),	0, path);
	confignew.Get(std::string("COMMONA"), std::string("STRATEGY"),	0, strategy);

	confignew.Get(std::string("DGTM"),	  std::string("TAOT"),		0, taot);
	confignew.Get(std::string("DGTM"),	  std::string("TINI"),      0, t0);
	confignew.Get(std::string("DGTM"),	  std::string("NT"),		0, nt);

	confignew.Get(std::string("COMMONP"),std::string("DIMENSION"),  0, dimension);

	std::vector< std::vector<double> > PF0;

	unsigned int ic,ipf0,ipf1;

	std::cout<<path<<std::endl;

	std::ofstream f;

	generation = 16*taot*nt + taot + 5;

	IND::evalFp = OBJ_FUN;

	RND.seed((long) time(NULL));

	unsigned int xdim = (dimension < 3) ? dimension:3;

	for(ir=0; ir<runs; ir++)
	{
		ic = ipf1 = ipf0 = 0;

		PF0.resize(popsize*20*nt);

		COEA* pEA = new COEA(taot, nt, dimension, generation, popsize, t0, problem);

		pEA->Init();

		//unsigned int cou = 0;
		while(!pEA->IsTerminate())
		{
			gen = pEA->Step();
			
			if(pEA->IsToChange())
			{
				for(i=0; i<(unsigned int)pEA->m_pop2.m_size; i++) 
				{
					PF0[ipf0+i].resize(pEA->DF()+xdim);
					for(j=0; j<pEA->DF(); j++)	PF0[ipf0+i][j]			 = pEA->m_pop2.m_pind[i].m_obj[j];
					for(j=0; j<xdim; j++)		PF0[ipf0+i][j+pEA->DF()] = pEA->m_pop2.m_pind[i].m_var[j];
				}
				for(i=(unsigned int)pEA->m_pop2.m_size; i<popsize; i++) 
				{
					PF0[ipf0+i].resize(pEA->DF()+xdim);
					for(j=0; j<pEA->DF(); j++)	PF0[ipf0+i][j]			 = pEA->m_pop2.m_pind[pEA->m_pop2.m_size-1].m_obj[j];
					for(j=0; j<xdim; j++)		PF0[ipf0+i][j+pEA->DF()] = pEA->m_pop2.m_pind[pEA->m_pop2.m_size-1].m_var[j];
				}
				ipf0 += popsize;
				//std::cout<<"T"<<cou<<"\t";
			}
			//cou++;
		}
		std::cout<<ir<<" ";

		// save data
		std::stringstream ss;
		ss<<path<<"_"<<ir<<".pop";
		f.open(ss.str().c_str());
		f<<std::scientific<<std::setprecision(5);
		for(i=0; i<ipf0; i++)
		{
			for(j=0; j<pEA->DF()+xdim; j++) f<<PF0[i][j]<<"\t";
			f<<std::endl;
		}
		f.close();

		PF0.clear();
		
		delete pEA;
	}

	std::cout<<std::endl;
}
