/** @file tstring.h
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
//---------------------------------------------------------------------------

#ifndef tstringH
#define tstringH

#include "types.h"
#include <iomanip>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>

using namespace std;
//---------------------------------------------------------------------------
class STRING {
private:


public:
  static int getNbDigits(const int& val){
  	if(val>0) return (int)log10((double)val ) + 1;
  	if(val<0) return (int)log10((double)-val) + 1;
  	return 1;
	}
  static int getNbDigits(const unsigned int& val){
  	if(val>0) return (int)log10((double)val ) + 1;
  	return 1;
	}
  // -----------------------------------------------------------------------------
  // transform a number to a string
	// -----------------------------------------------------------------------------
	// number to string or string to number conversions

  template <typename T> static string int2str(const T& value){
    ostringstream o;
  	if (!(o << value)) fatal("Could not convert number '%i' to string!\n", value);
  	return o.str();
	}

  template <typename T> static string int2str(const T& value, const T& max){
    ostringstream o;
		if(max && max>value) o << setfill('0') << setw(STRING::getNbDigits(max)) << right;
  	if (!(o << value)) fatal("Could not convert number '%i' to string!\n", value);
  	return o.str();
	}

	template <typename T> static string int2str_(const T& value){  			// _02
  	return ("_" + int2str(value));
  }
	template <typename T> static string int2str_(const T& value, const T& max, string s=""){  			// _t02
  	return ("_" + s + int2str(value, max));
  }

  // -----------------------------------------------------------------------------
  // transform a string to a number
	// -----------------------------------------------------------------------------
	template <typename T> static T str2int (const string& str) {
  	istringstream ss(str);
   	T x;
   	char c;
   	if (!(ss >> x) || ss.get(c)) return my_NAN;
   	return x;
	}

  /** file reading routines */
	static bool removeCommentAndSpace(istream & IN, int& line, int comment=0, char end='\0', bool removed = false);
  static string readBrackets2String(istream & IN, int& linecnt, char dim=')');

	static string readUntilCharacter(istream & IN, int& linecnt, char item, bool include);
  static map<string, int> readFileInfo(istream & IN, int& line);

  static bool file_exists(const string& name);
	static bool is_number(const string& str);
  static bool is_integer(const string& str);

// -----------------------------------------------------------------------------
	template <typename T>
	vector<T> static strMatrix2vector(string m){
		vector<T> vec;
		T val;

		istringstream LINE;                   // make a new stream
		// assert(m[0] == '{' && m[m.length()-1] == '}');
		string s = m.substr(1, m.length()-2); // removing the brackets
		LINE.str(s);  // allocate the matrix to the new stream

		while(LINE.good() && !LINE.eof()){
			LINE >> val;
			vec.push_back(val);
		}
		return vec;
	}

//------------------------------------------------------------------------------
	/** transformes the string matrix (1 dimensional) to an array of doubles, the
			array has to be created before and be passed  */
	template <typename T>
	void static strMatrix2array(const string& m, T* array, const int& size){
		istringstream LINE;                   // make a new stream
		// assert(m[0] == '{' && m[m.length()-1] == '}');
		LINE.str(m.substr(1, m.length()-2));  // allocate the matrix to the new stream (removing the brackets)

		for(int i=0; i<size; ++i){
			LINE >> array[i];
		}
		// assert(LINE.eof());
	}

//------------------------------------------------------------------------------
  string static seq(const string& t);
  void   static replace_seq(string& text);
  string static rep(const string& t);
  void   static replace_rep(string& text);
};
#endif
