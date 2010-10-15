/** @file fileparser.cpp
*
*   Copyright (C) 2006 Frederic Guillaume    <guillaum@zoology.ubc.ca>
*   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
*
*   quantiNEMO:
*   quantiNEMO is an individual-based, genetically explicit stochastic
*   simulation program. It was developed to investigate the effects of
*   selection, mutation, recombination, and drift on quantitative traits
*   with varying architectures in structured populations connected by
*   migration and located in a heterogeneous habitat.
*
*   quantiNEMO is built on the evolutionary and population genetics
*   programming framework NEMO (Guillaume and Rougemont, 2006, Bioinformatics).
*
*
*   Licensing:
*   This file is part of quantiNEMO.
*
*   quantiNEMO is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   quantiNEMO is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with quantiNEMO.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <sstream>
#include <fstream>
#include <errno.h>
#include "fileparser.h"
#include "output.h"

#include <streambuf>
#include <locale>

using namespace std;

//------------------------------------------------------------------------------
/** adds the new parameter with its argument to the map.
  * if the parameter has previously be set a warning is plotted
  */
void
FileParser::add_parsedParam (string& param, const vector<string>& argvect)
{
  map<string,vector<string> >::iterator pos = _parsedParams.find(param);
  if(pos != _parsedParams.end()) warning("Parameter '%s' has previously been set!\n", param.c_str());
  _parsedParams[param] = argvect;

}

//------------------------------------------------------------------------------
// FileParser::read
// ----------------------------------------------------------------------------------------
map< string, vector<string> >&
FileParser::read(string name, string dir)
{
	_dir = dir;
  // open the settings file
	ifstream FILE((dir+name).c_str());
	if(!FILE) fatal("Settings file '%s' could not be opened!\n", (dir+name).c_str());

  // read the settings file
  if(!read(FILE)) fatal("Reading settings file '%s' failed!\n",(dir+name).c_str());

	FILE.close();
  return _parsedParams;
}

//----------------------------------------------------------------------------------------
// FileParser::read
// ----------------------------------------------------------------------------------------
/** read the settings file */
bool
FileParser::read(istream & IN)
{
  int linecnt = 0;
  string key;                 // string to store parameter name
  string args;                // output string to collect parameter arguments
  string strEOL = "1";   strEOL[0] = EOL;
  char c;                     // temporary character
  vector< string > argvect;   // vector to store sequential arguments
  bool otherParams;           // true if there are more arguments ot be read for the same parameter
  _parsedParams.clear();      // remove previous parameters if present


  //--------------------------------------------------------------------------------
  //read the file parameter by parameter
	while(IN.good() && !IN.eof()) {
    linecnt++;
    argvect.clear();
    otherParams = true;

    // remove the forgoing space of the line
    if(!STRING::removeCommentAndSpace(IN,linecnt,0,EOL)) continue;

    //read the parameter name:
    key="";
    while(IN.get(c) && IN.good() && !IN.eof() && !isspace(c)){
      key += c;
    }
    if(c==EOL) IN.putback(c);

		#ifdef _DEBUG
			message("\n  line %i: %s", linecnt, key.c_str());
		#endif

    // remove space between parameter name and argument
    if(!STRING::removeCommentAndSpace(IN,linecnt,0,EOL)) continue;

    // get the arguments
    try{
      // check if the argument is given in an external file (argument starts with '$')
      if(IN.peek() == '$'){
        IN.get(c);      // remove the $
        string file;
        while(IN.get(c) && IN.good() && !IN.eof() && !isspace(c)){
          if(c=='#'){      // check if a comment is included
            IN.putback(c);                                       // comment
            if(!STRING::removeCommentAndSpace(IN, linecnt, 0, EOL)) break;
          }
          file += c;
        }
				if(c==EOL) IN.putback(c);

        // read the extern file
				ifstream EXTERN((_dir+file).c_str());
        if(!EXTERN) fatal("External file '%s' could not be found!\n", file.c_str());
        while(otherParams){
          args = readArguments(EXTERN, linecnt, c, strEOL, otherParams);
          if(!args.empty()) argvect.push_back(args);
        }
        EXTERN.close();
      }
      else{
        while(otherParams){
          args = readArguments(IN, linecnt, c, strEOL, otherParams);
          if(!args.empty()) argvect.push_back(args);
        }
      }
    }
    catch(const char* error){
      fatal("Parameter '%s' in line %i could not be read: %s\n", key.c_str(), linecnt, error);
    }

    #ifdef _DEBUG
	   message(":");
     unsigned int nb=argvect.size();
	   for(unsigned int i=0 ;i<nb; ++i){
      if(i) message(" |");
       message(" %s", argvect[i].c_str());
     }
     message(" (%i args)", nb);
    #endif

   //remove newline if EOL = \r in DOS files:
    if(EOL == '\r' && IN.peek() == '\n') IN.get(c);

		add_parsedParam(key, argvect);

	}//__END__WHILE__

	return true;
}

//------------------------------------------------------------------------------
/** read the arguments character by character until the "end" character.
  */
string
FileParser::readArguments(istream & IN, int& linecnt, char& c, string end, bool& otherParams){
	string args, temp;
  unsigned int i;

  // character by character
	while(IN.get(c) && IN.good() && !IN.eof()){
    // test for all ending characters
    for(i=0; i<end.size(); ++i){
     if(c != end[i]) continue;
     otherParams = false;
     return args;
    }

    switch(c){
      default  : // read
                 args += c;
                 break;

      case '\\': // line continuation
                 while(IN.get(c) && IN.good() && !IN.eof() && c != EOL){}   // remove rest of line
                 linecnt++;
                 break;

      case '{' : IN.putback(c);
                 args += STRING::readBrackets2String(IN, linecnt, '}');
                 break;

      case '(' : IN.putback(c);
                 args += STRING::readBrackets2String(IN, linecnt, ')');
                 break;

      case '[' : IN.putback(c);
                 args += STRING::readBrackets2String(IN, linecnt, ']');
                 break;

      case '\"': IN.putback(c);
                 args += STRING::readUntilCharacter(IN, linecnt, '\"', true);
                 break;

      case EOL : ++linecnt;
                 break;

      case '#' : IN.putback(c);                                       // comment
      case ' ' :                                                      // with space
			case '\t': if(!STRING::removeCommentAndSpace(IN, linecnt, 0, EOL)){    // tab
										otherParams = false; // that was the last argument
								 }

								 // remove trailing space
								 int i = (int)args.length();
								 while(isspace(args[i-1])) {
										--i;
								 }
								// if(i!=(int)args.length()) args = args.substr(0, i);
								 return args;
		}

	}//end while read args

  if(IN.eof()) otherParams = false;
  return args;
}


