/*! \file	Phrase.h
	
	\brief	phrase parameters from a file
	
	\author Aimin ZHOU
	\author Department of Computer Science,
	\author University of Essex, 
	\author Colchester, CO4 3SQ, U.K
	\author azhou@essex.ac.uk
	
	\date	Jun.22 2005 create
	\date	Sep.28 2005	redesign and rewrite
*/

#ifndef AZ_PHRASE_H
#define AZ_PHRASE_H

#include <fstream>
#include <vector>
#include <map>
#include <string>

//!\brief namespace PHR: phrase 
namespace PHR
{
	//!\brief definition of table
	typedef std::map< std::string, std::vector< double > > TABLE;

	//!\brief	skip spaces
	//!\param	line a line
	//!\param	start start phrase point
	//!\return	next start point
	int SkipSpaces(std::string& line, int start);

	//!\brief	phrase a word
	//!\param	word a word
	//!\param	line a line
	//!\param	start phrase start point
	//!\return	next start point
	int PhraseWord(std::string& word, std::string& line, int start);

	//!\brief	phrase a float value
	//!\param	val float value
	//!\param	line a line
	//!\param	start start phrase point
	//!\return	next phrase point
	int PhraseFloat(double& val, std::string& line, int start);

	//!\brief	phrase a line
	//!\param	table pharse table
	//!\param	line a lint
	//!\return	table
	TABLE& PhraseLine(TABLE& table, std::string& line);

	//!\brief	phrase a file
	//!\param	table table contains words and values
	//!\param	filename file name
	//!\return	table
	TABLE& PhraseFile(TABLE& table, std::string filename);
} //namespace PHR

#endif//AZ_PHRASE_H
