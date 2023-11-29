#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "AlgD.h"
#include "emo/Parameter.h"
#include "emo/Config.h"
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

			//tmpstr[0] = pathname + std::string("/");
			//for(k=0; k<sectionsize-1; k++) tmpstr[0] += names[k] + std::string("_");
			//tmpstr[0] += names[sectionsize-1];
			//confignew.Set(std::string("COMMONA"),std::string("PATHNAME"),tmpstr);
			//tmpstr[0] += std::string(".set");
			//confignew.Write(tmpstr[0]);
		}
	}
}

void execone()
{
	std::string method, problem, fname, path;
	unsigned int strategy, runs, generation, popsize, dimension, gen, ir, i, j;
	unsigned int torder, taot, nt; 
	double		 t0, alpha;
	az::mea::CParameter mPar;

	confignew.Get(std::string("COMMONA"), std::string("METHOD"),	0, method);
	confignew.Get(std::string("COMMONA"), std::string("RUNS"),		0, runs);
	confignew.Get(std::string("COMMONA"), std::string("POPSIZE"),	0, popsize);
	confignew.Get(std::string("COMMONA"), std::string("PROBLEM"),	0, problem);
	confignew.Get(std::string("COMMONA"), std::string("PATHNAME"),	0, path);
	confignew.Get(std::string("COMMONA"), std::string("STRATEGY"),	0, strategy);

	confignew.Get(std::string("DGTM"),	  std::string("TAOT"),		0, taot);
	confignew.Get(std::string("DGTM"),	  std::string("TORDER"),    0, torder);
	confignew.Get(std::string("DGTM"),	  std::string("TINI"),      0, t0);
	confignew.Get(std::string("DGTM"),	  std::string("NT"),		0, nt);
	confignew.Get(std::string("DGTM"),	  std::string("ALPHA"),		0, alpha);

	confignew.Get(std::string("COMMONP"),std::string("TOLERANCEF"), 0, mPar.TolF());
	confignew.Get(std::string("COMMONP"),std::string("TOLERANCEX"), 0, mPar.TolX());
	confignew.Get(std::string("COMMONP"),std::string("TOLERANCEC"), 0, mPar.TolC());
	confignew.Get(std::string("COMMONP"),std::string("DIMENSION"),  0, dimension);
	mPar.XSize(dimension);
	mPar.Problem(problem);
	mPar.XCoding() = false;

	std::vector< std::vector<double> > PF0, PF1;

	unsigned int ic,ipf0,ipf1;

	bool justinit = true;

	std::cout<<path<<std::endl;

	std::ofstream f0,f1;

	generation = 16*taot*nt + taot + 5;

	az::mea::dea::DMOO* pEA = new az::mea::dea::DMOO(strategy, method, popsize, generation, taot, nt, torder, t0, alpha, mPar);

	unsigned int xdim = (dimension < 3) ? dimension:3;

	for(ir=0; ir<runs; ir++)
	{
		ic = ipf1 = ipf0 = 0;

		PF0.resize(popsize*20*nt);
		PF1.resize(popsize*20*nt);

		pEA->Reset();

		//unsigned int cou = 0;
		while(!pEA->IsTerminate())
		{
			gen = pEA->Step();

			if(justinit)
			{
				justinit = false;
				for(i=0; i<pEA->Population().Size(); i++) 
				{
					PF1[ipf1+i].resize(pEA->P().FSize()+xdim);
					for(j=0; j<pEA->P().FSize(); j++) PF1[ipf1+i][j]		= pEA->Population()[i].F(j);
					for(j=0; j<xdim; j++) PF1[ipf1+i][j+pEA->P().FSize()]	= pEA->Population()[i][j];
				}
				ipf1 += pEA->Population().Size();
			}
			
			if(pEA->IsToChange())
			{
				for(i=0; i<pEA->Population().Size(); i++) 
				{
					PF0[ipf0+i].resize(pEA->P().FSize()+xdim);
					for(j=0; j<pEA->P().FSize(); j++) PF0[ipf0+i][j]		= pEA->Population()[i].F(j);
					for(j=0; j<xdim; j++) PF0[ipf0+i][j+pEA->P().FSize()]	= pEA->Population()[i][j];
				}
				ipf0 += pEA->Population().Size();
				justinit = true;
				//std::cout<<"T"<<cou<<"\t";
			}
			//cou++;
		}
		std::cout<<ir<<" ";

		// save data
		std::stringstream ss0,ss1;
		ss1<<path<<"_"<<ir<<".pin";
		//std::cout<<ss1.str()<<std::endl;
		f1.open(ss1.str().c_str());
		f1<<std::scientific<<std::setprecision(5);
		for(i=0; i<ipf1; i++)
		{
			for(j=0; j<pEA->P().FSize()+xdim; j++) f1<<PF1[i][j]<<"\t";
			f1<<std::endl;
		}
		f1.close();

		ss0<<path<<"_"<<ir<<".pop";
		f0.open(ss0.str().c_str());
		//std::cout<<ss0.str()<<std::endl;
		f0<<std::scientific<<std::setprecision(5);
		for(i=0; i<ipf0; i++)
		{
			for(j=0; j<pEA->P().FSize()+xdim; j++) f0<<PF0[i][j]<<"\t";
			f0<<std::endl;
		}
		f0.close();

		PF0.clear();
		PF1.clear();
	}
	delete pEA;

	std::cout<<std::endl;
}

//==================================================================================================
//int main(int argc, char* argv[])
//{
//	if(argc<2)
//	{
//		std::cout<<"need configure file"<<std::endl;
//		return -1;
//	}
//
//	az::mea::Config	confignew;	// configure for one run;
//	std::string	 infile;	//config file name
//	std::string method, problem, fname, path, str1, str2, str3;
//	unsigned int strategy, runs, generation, popsize, objs, dimension, gen, ir, i, j, sectionsize;
//	unsigned int torder, taot, nt; 
//	double		 t0;
//	az::mea::CParameter mPar;
//
//	infile		= std::string(argv[1]);
//	confignew.Load(infile, '#');
//
//	confignew.Get(std::string("COMMONA"), std::string("METHOD"),	0, method);
//	confignew.Get(std::string("COMMONA"), std::string("RUNS"),		0, runs);
//	confignew.Get(std::string("COMMONA"), std::string("POPSIZE"),	0, popsize);
//	confignew.Get(std::string("COMMONA"), std::string("PROBLEM"),	0, problem);
//	confignew.Get(std::string("COMMONA"), std::string("PATHNAME"),	0, path);
//	confignew.Get(std::string("COMMONA"), std::string("STRATEGY"),	0, strategy);
//
//	confignew.Get(std::string("DGTM"),	  std::string("TAOT"),		0, taot);
//	confignew.Get(std::string("DGTM"),	  std::string("TORDER"),    0, torder);
//	confignew.Get(std::string("DGTM"),	  std::string("TINI"),      0, t0);
//	confignew.Get(std::string("DGTM"),	  std::string("NT"),		0, nt);
//
//	confignew.Get(std::string("COMMONP"),std::string("TOLERANCEF"), 0, mPar.TolF());
//	confignew.Get(std::string("COMMONP"),std::string("TOLERANCEX"), 0, mPar.TolX());
//	confignew.Get(std::string("COMMONP"),std::string("TOLERANCEC"), 0, mPar.TolC());
//	confignew.Get(std::string("COMMONP"),std::string("DIMENSION"),  0, dimension);
//	mPar.XSize(dimension);
//	mPar.Problem(problem);
//	mPar.XCoding() = false;
//
//	std::vector< std::vector<double> > PCS, CCS, PF0, PF1;
//
//	unsigned int ic,ipf0,ipf1;
//
//	std::cout<<path<<std::endl;
//
//	std::ofstream f;
//
//	generation = 16*taot*nt + 5;
//
//	az::mea::dea::DMOO* pEA = new az::mea::dea::DMOO(strategy, method, popsize, generation, taot, nt, torder, t0, mPar);
//
//	for(ir=0; ir<runs; ir++)
//	{
//		ic = ipf1 = ipf0 = 0;
//
//		PCS.resize(500); CCS.resize(500); PF0.resize(popsize*20*nt); PF1.resize(popsize*20*nt);
//
//		pEA->Reset();
//
//		//unsigned int cou = 0;
//		while(!pEA->IsTerminate())
//		{
//			gen = pEA->Step();
//			
//			if(pEA->IsToChange())
//			{
//				for(i=0; i<pEA->Population().Size(); i++) 
//				{
//					PF0[ipf0+i].resize(pEA->P().FSize());
//					for(j=0; j<pEA->P().FSize(); j++) PF0[ipf0+i][j] = pEA->Population()[i].F(j);
//				}
//				ipf0 += pEA->Population().Size();
//
//				//std::cout<<"T"<<cou<<"\t";
//			}
//			else if(pEA->IsChange())
//			{
//				if(pEA->Centers(CCS[ic],PCS[ic])) ic++;
//
//				for(i=0; i<pEA->Population().Size(); i++)
//				{
//					PF1[ipf1+i].resize(pEA->P().FSize());
//					for(j=0; j<pEA->P().FSize(); j++) PF1[ipf1+i][j] = pEA->Population()[i].F(j);
//				}
//				ipf1 += pEA->Population().Size();
//
//				//std::cout<<"C"<<cou<<"\t";
//			}
//
//			//cou++;
//		}
//		std::cout<<ir<<" ";
//
//		// save data
//		char chr[100];
//		sprintf(chr,"_%d_0.pop",ir);
//		fname = path + std::string(chr);
//		f.open(fname.c_str());
//		f<<std::scientific<<std::setprecision(5);
//		for(i=0; i<ipf0; i++)
//		{
//			for(j=0; j<pEA->P().FSize(); j++) f<<PF0[i][j]<<"\t";
//			f<<std::endl;
//		}
//		f.close();
//
//		sprintf(chr,"_%d_1.pop",ir);
//		fname = path + std::string(chr);
//		f.open(fname.c_str());
//		f<<std::scientific<<std::setprecision(5);
//		for(i=0; i<ipf1; i++)
//		{
//			for(j=0; j<pEA->P().FSize(); j++) f<<PF1[i][j]<<"\t";
//			f<<std::endl;
//		}
//		f.close();
//
//		sprintf(chr,"_%d.cen",ir);
//		fname = path + std::string(chr);
//		f.open(fname.c_str());
//		f<<std::scientific<<std::setprecision(5);
//		for(i=0; i<ic; i++)
//		{
//			for(j=0; j<CCS[0].size(); j++) f<<CCS[i][j]<<"\t";
//			f<<std::endl;
//		}
//		for(i=0; i<ic; i++)
//		{
//			for(j=0; j<PCS[0].size(); j++) f<<PCS[i][j]<<"\t";
//			f<<std::endl;
//		}
//		f.close();
//
//		PCS.clear(); CCS.clear(); PF0.clear(); PF1.clear(); 
//	}
//	delete pEA;
//
//	std::cout<<std::endl;
//}