/*! \file	LogFile.h
	
	\brief	writing log file

	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk

	\date	May.29 2006 create
*/

#ifndef AZ_LOGFILE_H
#define	AZ_LOGFILE_H

#include <ctime>
#include <fstream>

//!\brief	az namespace, the top namespace
namespace az
{

//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
namespace mea
{

//\brief Log File class
class LogFile
{
protected:
	std::ofstream	fhand;
	std::string		fname;
	char			fmode;
public:
	//!\brief	constructor
	//!\param	mode write option: a or A - append, c or C - create
	//!\return	void
	LogFile(const char mode = 'A')
	{
		if(mode == 'C' || mode == 'c')	fmode = 'C';
		else							fmode = 'A';

		fname = "meda.log";

		if(fmode == 'A')fhand.open(fname.c_str(),std::ios::app);
		else			fhand.open(fname.c_str());
	}

	//!\brief	destructor
	//!\return	void
	~LogFile()
	{
		if(fhand != 0) fhand.close();
	}

	//!\brief	clear logfile
	//!\return	void
	void clear()
	{
		if(fhand != 0)	fhand.close();
		if(fmode == 'A')fhand.open(fname.c_str(),std::ios::app);
		else			fhand.open(fname.c_str());
	}

	//!\brief	write time lable in log file
	//!\return	void
	void label()
	{
		time_t t;
		time(&t);
		struct tm *times= localtime(&t);
		fhand<<times->tm_hour<<":"<<times->tm_min<<":"<<times->tm_sec<<","<<times->tm_mday<<"/"<<times->tm_mon+1<<"/"<<times->tm_year+1900<<std::endl;
	}

	//!\brief	write data
	//!\return	output stream
	template< class T >
		std::ofstream & operator<<(const T& value)
	{
		fhand<<value;
		return fhand;
	}
};//class LogFile

} //namespace mea

} //namespace az

#endif	//AZ_LOGFILE_H
