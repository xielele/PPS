// Phrase.h

#include "Phrase.h"

//namespace PHR: phrase 
namespace PHR
{

	//skip spaces
	int SkipSpaces(std::string& line, int start)
	{
		int pos = start;
		while (pos < (int)line.length() && (line[pos]=='\t' || line[pos]==' ' || line[pos]=='\r' || line[pos]=='\n'))
			pos++;
		return pos;
	}

	//phrase a word
	int PhraseWord(std::string& word, std::string& line, int start)
	{
		int pos = start;
		while(pos < (int)line.length() && (isalnum(line[pos]) || line[pos] == '#')) 
				pos++;
		word = line.substr(start, pos);

		return pos;
	}	

	//phrase a float value
	int PhraseFloat(double& val, std::string& line, int start)
	{
		int pos = start;
		while(pos < (int)line.length() && (line[pos] == '+' || line[pos] == '-' || line[pos] == '.' || line[pos] == 'E'|| line[pos] == 'e'|| isdigit(line[pos]))) 
				pos++;
		
		val = atof(line.substr(start, pos).c_str());

		return pos;
	}

	//phrase a line
	TABLE& PhraseLine(TABLE& table, std::string& line)
	{
		int		pos = 0;
		double	val;
		std::string	key;
		std::vector< double > value;
		
		pos = SkipSpaces(line, pos);

		//blank line
		if(pos >= (int)line.length()) return table;

		pos = PhraseWord(key, line, pos);

		//comments
		if(key == std::string("#")) return table;

		pos = SkipSpaces(line, pos);
		
		//float parameters
		while(pos < (int)line.length())
		{
			pos = PhraseFloat(val, line, pos);
			pos = SkipSpaces(line, pos);
			value.push_back(val);
		}
	
		table[ key ] = value;

		return table;
	}

	//phrase a file
	TABLE& PhraseFile(TABLE& table, std::string filename)
	{
		std::string sline;

		table.clear();

		std::ifstream file(filename.c_str());
		
		while(std::getline(file, sline)) 
		{
			PhraseLine(table, sline);
		}
		file.close();
		return table;
	}
} //namespace PHR
