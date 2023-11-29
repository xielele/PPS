/*! \file	Config.h
	
	\brief	read parameters from config file
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
		
	\date	Apr.12 2006 create
*/

#ifndef AZ_CONFIG_H
#define	AZ_CONFIG_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

//!\brief namespace of config
namespace CFG
{
	//!\brief definition of table
	typedef std::map< std::string, std::vector< std::string > > SECTION;
	typedef std::map< std::string, SECTION > TABLE;

	class Config
	{
	protected:
		TABLE	mTable;			//!< table to store parameters
		char	mComment;		//!< comment symbol
		std::string mSection;	//!< current section
	public:
		Config() {}
		
		void Load(const std::string& filename, const char& comment = '#')
		{
			mTable.clear();
			mComment = comment;
			mSection = "default";

			std::string sline;

			std::ifstream file(filename.c_str());
			
			while(std::getline(file, sline)) 
				PhraseLine(sline);
			
			file.close();
		}

		void Write(const std::string& filename)
		{
			std::ofstream file(filename.c_str());
			TABLE::iterator itt;
			SECTION::iterator its;
			for(itt=mTable.begin(); itt!=mTable.end(); itt++)
			{
				file<<"["<<itt->first<<"]"<<std::endl;
				for(its=itt->second.begin(); its!=itt->second.end(); its++)
				{
					file<<std::left<<std::setw(15)<<its->first<<"\t";
					for(unsigned int i=0; i<(its->second).size(); i++) file<<(its->second)[i]<<" ";
					file<<std::endl;
				}
				file<<std::endl;
			}
			file.close();
		}
		
		//bool&			Get(const std::string& section, const std::string& parameter, unsigned int index, bool& value);
		//int&			Get(const std::string& section, const std::string& parameter, unsigned int index, int& value);
		//unsigned int&	Get(const std::string& section, const std::string& parameter, unsigned int index, unsigned int& value);
		//double&			Get(const std::string& section, const std::string& parameter, unsigned int index, double& value);
		//std::string&	Get(const std::string& section, const std::string& parameter, unsigned int index, std::string& value);

		template<class T>
			T& Get(const std::string& section, const std::string& parameter, unsigned int index, T& value)
		{
			std::string value_string = mTable[section][parameter][index];
			std::istringstream istr(value_string);
			istr >> value;
			return value;
		}

		template<class T>
			T& Get(const std::string& section, const std::string& parameter, unsigned int index)
		{
			T value;
			std::string value_string = mTable[section][parameter][index];
			std::istringstream istr(value_string);
			istr >> value;
			return value;
		}

		unsigned int GetSize(const std::string& section, const std::string& parameter)
		{
			return (unsigned int)mTable[section][parameter].size();
		}

		void Set(const std::string& section, const std::string& parameter, const std::vector< std::string >& value)
		{
			mTable[section][parameter] = value;
		}
	protected:
		void PhraseLine(std::string& line)
		{
			unsigned int pos = 0;
			std::string	word,parameter;
			std::vector< std::string > value;
			
			SkipSpaces(line, pos);

			//blank line
			if(pos >= line.length()) return;
			
			PhraseWord(word, line, pos);

			//the line is comments 
			if(word[0] == mComment) return;
			
			//section
			if(word == std::string("["))
			{
				SkipSpaces(line, pos);
				PhraseWord(mSection, line, pos);
				SkipSpaces(line, pos);
				PhraseWord(word, line, pos);
				//if(word != std::string("]")) wrong!!!!
				return;
			}
		
			//parameters
			parameter = word;
			SkipSpaces(line, pos);
			value.clear();
			while(pos < line.length())
			{
				PhraseWord(word, line, pos);
				SkipSpaces(line, pos);
				value.push_back(word);
			}
			if(parameter.size()>0) mTable[mSection][parameter] = value;
		}

		void PhraseWord(std::string& word, std::string& line, unsigned int& pos)
		{
			unsigned int start = pos;
			while(pos < line.length() && !isspace(line[pos]) && line[pos] != '[' && line[pos] !=']' && line[pos] != mComment ) 
					pos++;
			if(pos == start && (line[pos] == '[' || line[pos] ==']' || line[pos] == mComment)) pos++;
			word = line.substr(start, pos-start);
		}

		void SkipSpaces(std::string& line, unsigned int& pos)
		{
			while(pos < line.length() && isspace(line[pos]) )	pos++;
		}
	}; //class Config

} //namespace CFG

#endif //AZ_CONFIG_H
