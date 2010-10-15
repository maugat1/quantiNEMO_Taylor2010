/** @file tstring.cpp
*
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

#include "tstring.h"
#include "output.h"
#include <cerrno> 
using namespace std;

//------------------------------------------------------------------------------
/** read white space and comment unitl the end of the file, or the specified character "end".
  * returns FALSE if the "end "character is first reached return a false
  * if not commented text is reached first return true
  * #: out commented line (until the end of the line)
  * #/ ... /# commented parts (may be among several lines)
  */
bool
STRING::removeCommentAndSpace(istream & IN, int& line, int comment, char end, bool removed){
  char c;

  switch(comment){
    case 0: // no comment
        while(IN.get(c) && c != end){
          if(isspace(c)){
            if(c == EOL) ++line;
            continue;
          }
          if(c == '#'){      // start of a comment
            if(IN.peek() == '/'){
              IN.get();
              return removeCommentAndSpace(IN, line, 2, end);
            }
            else                 return removeCommentAndSpace(IN, line, 1, end);
          }

          // no comment
          IN.putback(c);  // put back the character and return
          if(removed) IN.putback(' ');    // if it was a comment 2 control that there is at least a space
          return true;
        }
        return false;

    case 1: // comment #
        while(IN.get(c)){
          if(c == EOL){
            if(c == end) return false;
            ++line;
            return removeCommentAndSpace(IN, line, 0, end);
          }
        }
        return false;

    case 2: // comment #/ ... /#
        while(IN.get(c)){
          if(c == EOL) ++line;
          if(c == '/' && IN.peek() == '#'){
            IN.get();     // remove #
            return removeCommentAndSpace(IN, line, 0, end, true);
          }
        }
        return false;
  }
  return false;     // will never be used, but the compiler preferes it...
}

//------------------------------------------------------------------------------
/** file information starts with the key word "[FILE_INFO]" and the infomration
  * is enclosed by {...}
  */
map<string, int>
STRING::readFileInfo(istream & IN, int& line){
  char c;
  string key, keyWord="[FILE_INFO]";
  int value, i;
  map<string, int> info;

  // check if it is really a file info
  for(i=0; i<(int)keyWord.size(); ++i){
    IN.get(c);
    if(c != keyWord[i]){
      // it is not a file info: put all read characters back
      IN.putback(c);     // current character
      while(i >= 0){     // all previous read characters
        IN.putback(keyWord[i]);
        --i;
      }
      return info;
    }
  }

  // get the opening bracket
  if(!removeCommentAndSpace(IN, line)) fatal("File information is not complete!\n");
  IN.get(c);
  if(c != '{') fatal("File info could not be properly read: misplaced '{'!\n");
  if(!removeCommentAndSpace(IN, line)) fatal("File information is not complete!\n");

  // read the info line by line
  while(IN.good() && !IN.eof()) {
    if(!removeCommentAndSpace(IN, line)) fatal("File information is not complete1!\n");
    if(IN.peek() == '}') break;

    //read the parameter name:
    key="";
    while(IN.get(c) && IN.good() && !IN.eof() && !isspace(c) && c != EOL){
      if(c == '}') fatal("File info could not be read properly!\n");
      if(c == '#'){
        IN.putback(c);
        if(!removeCommentAndSpace(IN, line)) fatal("File information is not complete!\n");
      }
      key += c;
    }
    if(c==EOL) fatal("File info could not be properly read!\n");

    // read the argument (is an integer!)
    removeCommentAndSpace(IN, line);
    IN >> value;

    // save the file info to a map
    if(info.find(key) != info.end()) warning("Column '%s' has previously been set!\n", key.c_str());
    info[key] = value;
  }

  // remove the }
  IN.get(c);
  assert(c == '}');

  return info;
}

//------------------------------------------------------------------------------
/** read a matrix and return it as a string
  * a matrix is read unitl the number of opening gand closing brackets is the same
  * comments are removed
  * no control if numbers are read, or text is between rows */
string
STRING::readBrackets2String(istream & IN, int& linecnt, char dim){
  char c;
	int whiteSpace = -1;  // -1: ignore coming ws, 0: last char was important, 1: ws needed
  string args;

  // read the first char
  IN.get(c);
  args += c;
  assert(c=='{' || c=='(' || c=='[');

	// read brackets
	while(IN.get(c) && IN.good() && !IN.eof()){
    switch(c){
      case '{': IN.putback(c);  args += readBrackets2String(IN, linecnt, '}'); continue;
      case '(': IN.putback(c);  args += readBrackets2String(IN, linecnt, ')'); continue;
      case '[': IN.putback(c);  args += readBrackets2String(IN, linecnt, ']'); continue;
    }

    if(c == dim){
      args +=c;
      return args;
    }

		if(isspace(c)){
      if(c == '\n') ++linecnt;
			if(whiteSpace != -1) whiteSpace = 1;
			continue;
		}

		if(c == '#'){
      IN.putback(c);
			if(!removeCommentAndSpace(IN, linecnt)) throw "Could not read text between brackets!";
			if(whiteSpace != -1) whiteSpace = 1;
			continue;
    }

		// that character is needed, add it to the string
		if(whiteSpace == 1) args += " ";
		whiteSpace = 0;
		args += c;
	}

	throw "Could not read text between brackets: missing closing bracket!";
}

//------------------------------------------------------------------------------
/** read string until the specified character occurs. The starting and ending
  * character may be included or excluded.
  */
string
STRING::readUntilCharacter(istream & IN, int& linecnt, char item, bool include){
  char c;
  string args;
  int whiteSpace = 0;

  IN.get(c);
  if(include) args = c;

  // read brackets
  do{
    if(!IN.get(c)) throw "Could not read between two characters!";
    if(c == item) break;
    if(c == '#'){
      IN.putback(c);
      if(!removeCommentAndSpace(IN, linecnt)) throw "Could not read between two characters!";
      if(!whiteSpace) args += " ";
      continue;
    }

    if(isspace(c)){
      if(!whiteSpace){
        ++whiteSpace;
        args += " ";
      }
    }
	  else{
      whiteSpace = 0;
      args += c;
    }
  }while(1);

  //remove the trailling space
  int i;
  for(i=(int)args.size(); i>0; --i){
    if(!isspace(args[i-1])) break;
  }
  args = args.substr(0,i);

  if(include) args += item;
  return args;
}

//------------------------------------------------------------------------------
bool
STRING::file_exists(const string& name){
  ifstream ifs(name.c_str());
  if(!ifs) return false;
  ifs.close();
  return true;
}
// ----------------------------------------------------------------------------------------
// is_number
// ----------------------------------------------------------------------------------------
bool
STRING::is_number(const string& str)
{
  unsigned int i = 0;
  while(i < str.size()) {
	  if(!isdigit(str[i])){
	    if(str[i] != '.' && str[i] != 'e' && str[i] != '-') return false;
    }
	  i++;
  }

  return true;
}
// ----------------------------------------------------------------------------------------
// is_integer
// ----------------------------------------------------------------------------------------
bool
STRING::is_integer(const string& str)
{
  unsigned int i = 0;
  while(i < str.size()) {
	  if(!isdigit(str[i])){
	    if(str[i] != 'e' && str[i] != '-') return false;
    }
	  i++;
  }

  return true;
}

// ----------------------------------------------------------------------------------------
// seq
// ----------------------------------------------------------------------------------------
  /** function seq adapted to the function seq in R:
    * format: seq(from, to, by)
    * input: the entire format including seq and the parantheses  "seq(1,10,2)"
    * output: the sequence specified  "1 3 5 7 9"
    */
string
STRING::seq(const string& t){
  assert(t.substr(0,4)=="seq(");

  istringstream IN;                             // make a new stream
  IN.str(t.substr(4,t.length()-5));             // allocate the line to the new stream
  double from1, from2=my_NAN, to1, to2;
  int by;                                       // number of steps
  char c1, c2;                                  // comma 1, comma 2
  IN >> from1 >> ws;                            // read the first value

  if(IN.peek() == ','){                         // if simple sequence
    IN >> c1 >> to1 >> c2 >> by;                // read the rest
  }
  else{                                         // if temporal sequence
    IN >> from2 >> c1 >> to1 >> to2 >> c2 >> by;// read the rest
  }

  // test if the seq was correctly set
  if(c1!=',' || c2!=',') fatal("Sequence '%s': The parameters have to be separated by a comma ','!\n", t.c_str());

  // create the sequence
  string text = STRING::int2str(from1);
  if(from2 == my_NAN){                          // if simple sequence
    double step = (to1-from1)/(by-1);
    for(int i=1; i<by; ++i){
      from1 += step;
      text += " " + STRING::int2str(from1);
    }
  }
  else{                                         // if temporal sequence
    double step1 = (to1-from1)/(by-1);
    double step2 = (to2-from2)/(by-1);
    text += " " + STRING::int2str(from2);
    for(int i=1; i<by; ++i){
      from1 += step1;
      from2 += step2;
      text  += ", " + STRING::int2str(my_round(from1)) + " " + STRING::int2str(from2);
    }
  }
  return text;
}

// ----------------------------------------------------------------------------------------
// replace_seq
// ----------------------------------------------------------------------------------------
/** replaces in the text string all seq() by its sequence.
  * Start to search fromt he back
  */
void
STRING::replace_seq(string& text){
  int cur = text.length(), nb;
	string::size_type next;
  next = text.rfind("seq(", cur);
	while(next != string::npos){
    nb  = text.find(')', next)+1-next;
    text.replace(next, nb, seq(text.substr(next, nb)));
    cur = next -1;
    next = text.rfind("seq(", cur);
  }
}

// ----------------------------------------------------------------------------------------
// rep
// ----------------------------------------------------------------------------------------
  /** function rep adapted to the function rep in R:
    * format: seq(val, nb)
    * input: the entire format including seq and the parantheses  "rep(1,10)"
    * output: the sequence specified  "1 1 1 1 1 1 1 1 1 1"
    */
string
STRING::rep(const string& t){
  assert(t.substr(0,4)=="rep(");

  istringstream IN;                             // make a new stream
  IN.str(t.substr(4,t.length()-5));             // allocate the line to the new stream
  double val;
  int nb;                                       // number of steps
  char c;                                       // comma 1, comma 2
  IN >> val >> c >> nb;                         // read the first value

  // test if the rep was correctly set
  if(c!=',') fatal("Repetition '%s': The parameters have to be separated by a comma ','!\n", t.c_str());

  // create the sequence
  string valStr = STRING::int2str(val);
  string text = valStr;
  for(int i=1; i<nb; ++i){
    text += " " + valStr;
  }
  return text;
}
// ----------------------------------------------------------------------------------------
// replace_rep
// ----------------------------------------------------------------------------------------
/** replaces in the text string all rep() by its sequence.
  * Start to search fromt he back
  */
void
STRING::replace_rep(string& text){
  int cur = text.length(), nb;
	string::size_type next;
  next = text.rfind("rep(", cur);
	while(next != string::npos){
    nb  = text.find(')', next)+1-next;
    text.replace(next, nb, rep(text.substr(next, nb)));
    cur = next -1;
    next = text.rfind("rep(", cur);
  }
}


