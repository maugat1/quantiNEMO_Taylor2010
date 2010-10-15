/** @file fileparser.h
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

#ifndef fileparserH
#define fileparserH

#include "output.h"


/**Text input parameter file parser.
 * This class provides the StreamParser with the whole content of the input file.
 */
class FileParser {

public:
  FileParser() {}
  virtual ~FileParser(){}

  map<string, vector<string> >&  read(string name, string dir="");
  bool read(istream& stream);
  map<string, vector< string> >& get_parsedParams() {return _parsedParams;}


private:
  map< string, vector< string > > _parsedParams;
  string readArguments(istream & IN, int& linecnt, char& c, string end, bool&);
	void add_parsedParam (string& param, const vector<string>& argvect);

	string _dir;

};
#endif
