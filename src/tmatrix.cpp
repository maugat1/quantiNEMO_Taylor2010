/** @file tmatrix.cpp
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

#include "tmatrix.h"
// ----------------------------------------------------------------------------------------
// read
// ----------------------------------------------------------------------------------------
void
TMatrix::read(string text)
{
  reset();
	unsigned int i, j, size;

  // read it a as a TMatrixVar object
	TMatrixVar<double> m;
	try{
		m.read(text);
	}
  catch(const char* error){
		fatal("ParameterMatrix '%s': %s\n", text.c_str(), error);
	}

	_rows = m.getNbRows();

  // find the number of columns (it could have empty columns)
  for(i = 0; i < _rows; i++) {
    _cols = m.getNbCols(i);
    if(_cols) break;
  }

  //check for matrix coherence:
  for(i = 0; i < _rows; i++) {
    size = m.getNbCols(i);
    if(!size){                        // if it is a empty column, fill it with my_NANs
      for(j = 0; j < _cols; j++) {
        m(i)->push_back(my_NAN);
      }
    }
	  else if(size != _cols){
	    throw "Could not read matrix: not same number of columns per row!";
    }
  }

  //copy to input TMatrix:
  if(_val) delete _val;
  _val = new double[_cols*_rows];
  for(i = 0; i < _rows; ++i){
	  for(j = 0; j < _cols; ++j){
      _val[i*_cols + j] = m.get(i,j);
    }
  }
}

// ----------------------------------------------------------------------------------------
// read_matrix
// ----------------------------------------------------------------------------------------
// read a matrix file and return the file infomration if available
map<string, int>
TMatrix::read_matrix (string filename)
{
  ifstream IN(filename.c_str());
  if(!IN) fatal("File '%s' could not be opened!", filename.c_str());

  int line=0;
  map<string, int> fileInfo;

  // read the file info if available
  STRING::removeCommentAndSpace(IN, line);
  fileInfo = STRING::readFileInfo(IN, line);
  if(fileInfo.empty()) fatal("The file info box of file '%s' could not be read!\n", filename.c_str());


  STRING::removeCommentAndSpace(IN, line);

  // read the matrix (check if the matrix has brackets)
  try{
    //if(IN.peek()=='{') read_matrix_bracket(IN);
    //else               read_matrix_none(IN);
    read_matrix_none(IN);
  }
  catch(const char* error){
    fatal("File '%s': \n", filename.c_str());
  }

  IN.close();

  return  fileInfo;
}

  // check if the matrix is
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
// read_matrix_bracket
// ----------------------------------------------------------------------------------------
/** reading a matrix from a stream. The matrix has to be specified by brackets */
void
TMatrix::read_matrix_bracket (istream & IN)
{
  vector< vector<double> > tmpMat;

  int line;
  double elmnt;
  char c;
  int openBrack = 1,  // the one which was just removed
      closeBrack = 0;

  //remove the first enclosing bracket
  STRING::removeCommentAndSpace(IN, line);
  IN.get(c);
  assert(c=='{');

  //then read the rows
  for(int i=0; openBrack != closeBrack && IN; ++i){
    tmpMat.push_back( vector<double>());  // add a new row

    // go to next row
    STRING::removeCommentAndSpace(IN, line);
    IN.get(c);             //first character:
    if(c == '{') ++openBrack;
	  else throw "Could not read matrix: brackets are misplaced!";

    //read a row enclosed by {...}:
    while(IN) {
      STRING::removeCommentAndSpace(IN, line);
      c=IN.peek();
      if(c == ',' || c == ';') continue;
      if(c == '}') {
        ++closeBrack;
        IN.get(c);
        break;
      }    //go to next row

      //read a row element:
      IN >> elmnt;
      tmpMat[i].push_back(elmnt);
	  }
  }

  //check for matrix coherence:
  unsigned int rows = tmpMat.size();       // get the number of rows
  unsigned int cols = tmpMat[0].size();   // get the size of the first row
  unsigned int i, j;
  for(i = 1; i < rows; i++) {
	  if(tmpMat[i].size() != cols){
	    throw "Could not read matrix: not same number of elements in all rows!";
    }
  }

  //copy to input TMatrix:
  reset(rows, cols);
  for(i = 0; i < rows; ++i){
	  for(j = 0; j < cols; ++j){
	    set(i,j,tmpMat[i][j]);
    }
  }
}

// ----------------------------------------------------------------------------------------
// read_matrix_none
// ----------------------------------------------------------------------------------------
/** reading a matrix from a stream. The matrix has to be specified WIHTOUT brackets,
  * each line corresponds to a row
  */
void
TMatrix::read_matrix_none (istream & IN)
{
  vector< vector<double>* > tmpMat;
  int line, lineTemp;
  double elmnt;
  string lineStr;
  char c;
  string text;


    //remove the first enclosing bracket
	STRING::removeCommentAndSpace(IN, line, 0, '\n', false);

	// remove the comments in the table
	while(IN.get(c)){
    if(c == '#'){
      IN.putback(c);
			if(!STRING::removeCommentAndSpace(IN, lineTemp, 0, '\n')){
				text += '\n';           // end of line is reached (but we need the return...)
				continue;
			}
		}
    else text += c;
  }

  // for each line
  istringstream ALL_LINES;             // make a new stream
  ALL_LINES.str(text);              // allocate the line to the new stream
  unsigned int i=0;
  while(ALL_LINES.good() && !ALL_LINES.eof()){
    getline(ALL_LINES, lineStr);           // read the line
    if(lineStr.empty()) continue;
    istringstream LINE;             // make a new stream
    LINE.str(lineStr);              // allocate the line to the new stream
    LINE >> ws;


    //read a row :
    tmpMat.push_back( new vector<double>);  // add a new row
    while(LINE.good() && !LINE.eof()) {
      //read a row element:
      if(LINE.peek() == '{'){
        if(!_strMatrix) _strMatrix = new vector< map<int, string> >;
        if(_strMatrix->size() <= i) _strMatrix->push_back(map<int, string>());  // add a new row
        (*_strMatrix)[i][tmpMat[i]->size()] = STRING::readBrackets2String(LINE, lineTemp, '}');
        tmpMat[i]->push_back(my_NAN);
      }
      else{
        LINE >> elmnt;
        tmpMat[i]->push_back(elmnt);
      }
      LINE >> ws;
	  }
    ++i;
  }

  //check for matrix coherence:
  unsigned int rows = tmpMat.size();       // get the number of rows
  unsigned int cols = tmpMat[0]->size();   // get the size of the first row
  unsigned int j;
  for(i = 1; i < rows; i++) {
	  if(tmpMat[i]->size() != cols){
      throw "Could not read matrix: not same number of elements in all rows of the matrix!";
    }
  }

  //copy to input TMatrix:
  reset(rows, cols);
  for(i = 0; i < rows; ++i){
	  for(j = 0; j < cols; ++j){
	    set(i,j,(*tmpMat[i])[j]);
    }
    delete tmpMat[i];
  }
}


