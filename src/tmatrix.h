/** @file tmatrix.h
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
//---------------------------------------------------------------------------

#ifndef tmatrixH
#define tmatrixH
//---------------------------------------------------------------------------

#include <map>
#include "types.h"
#include "output.h"
#include <vector>

using namespace std;



/** This matrix is variable in both dimensions and allows to read any matrix input
	* The variable matrix may be of any kind. Values are stored in a double vector.
	* If the number of column is not identical for all rows, the function get_nbCols returns NaN.
	**/
template <class T> class TMatrixVar {
private:

  unsigned int _rows, _cols;   // if the rows have not the same number of elements: _cols = my_NAN

	vector< vector<T> > _vals;

public:

	TMatrixVar() : _rows(0), _cols(0){ }
	TMatrixVar(string text){read(text);}

	TMatrixVar  (const TMatrixVar<T>& mat) {
		_rows   = mat._rows;
		_cols   = mat._cols;
		_vals   = mat._vals;
	}

	~TMatrixVar () {	}

	/**@brief Accessor to element at row i and column j**/
	T operator()(const int& i, const int& j) {return get(i,j);}
  T get (unsigned int i, unsigned int j) {
		assert(i<_rows && j<_vals[i].size());
    return _vals[i][j];
  }

  /**@brief Accessor to row i**/
	vector<T>* operator()(const int& i) {return get(i);}
  vector<T>* get (unsigned int i) {
	  assert(i < _rows);
    return &_vals[i];
  }
  /**@brief Accessor to the whole double vector**/
	vector<vector<T> >* get() {return &_vals;}

  /**@brief Accessor to the matrix dimensions
    *@param dims an array of at least 2 elements to store the row [0] and column [1] numbers. May be NULL.
    *@return the total size of the matrix
   **/
  unsigned int get_dims (unsigned int* dims) {
	  if(dims) {
      dims[0] = _rows;
      dims[1] = _cols;
    }
	  return _rows*_cols;
  }


	/** set the dimension of the matrix. If not all the rows have the same number
    * of elements, _cols will be set to my_NAN.
    */
  void set_dimensions(){
		
		// get the number of rows
		_rows = _vals.size();
		if(!_rows) {_cols = my_NAN;  return;}    // check if the matrix is empty

    // check if all rows have the same number of elements
    _cols = _vals[0].size();
    for(unsigned int i=1; i<_rows; ++i){
      if(_vals[i].size() != _cols){
				_cols = my_NAN;
        return;
      }
    }
  }

  /**Gives the number of rows.*/
  unsigned int getNbRows ( ) {return _rows;}

	/**Gives the number of columns. (my_NAN if the columns have different sizes)*/
	unsigned int getNbCols ( ) {return _cols;}

	/**Gives the number of columns of the row i */
	unsigned int getNbCols (int i) {
		assert(i<_rows);
		return _vals[i].size();
	}

  /**Reset members to zero state.*/
	void reset ( ) {
		_rows = 0;
		_cols = 0;
    _vals.clear();
	}

	/** extracts from a string the matrix (1D or 2D) */
	void read(string text){
		reset();
		istringstream IN;
		IN.str(text);
		string rowStr;

		T elmnt;
		char c;
		bool oneDims = false;
		int line=0;

		//remove the first enclosing bracket
		IN >> c;
		assert(c=='{');

		//then read the rows
		while(IN.peek() != '}' && !oneDims){
			IN >> ws;
			if(IN.peek()!='{'){
				if(!_vals.empty()) throw "Could not read matrix: text between rows!";

				// it is a one dimensional array
				IN.putback('{');
				oneDims = true;
			}

			_vals.push_back(vector<T>());                    // add a new row

			//read a row enclosed by {...}:
			rowStr = STRING::readUntilCharacter(IN, line, '}', false);

			// check if it has a row indicator and add the missing rows if necessary
			string::size_type pos = rowStr.find_first_of(':');
			if(pos != string::npos){
				unsigned int nb = STRING::str2int<unsigned int>(rowStr.substr(0, pos));
				rowStr = rowStr.substr(pos+1, rowStr.size());
				if(nb <_vals.size()) throw "Could not read matrix: row indicator must exceed the number of the previous row!";
				while(_vals.size() < nb){
					_vals.push_back(vector<T>());               // add a new row
				}
			}
	
			// read the row
			istringstream ROW;
			ROW.str(rowStr);
			ROW >> ws;
			while(!ROW.eof()) {
				ROW >> elmnt;
				_vals.back().push_back(elmnt);                  // add the element
				ROW >> ws;
			}
		}

		set_dimensions();
	}
};
	
	
////////////////////////////////////////////////////////////////////////////////
	/**A class to handle matrix in params, coerces matrix into a vector of same total size**/
class TMatrix {
	private:
	
  unsigned int _rows, _cols;
	
	  double* _val;
	
	  vector<map<int, string> >* _strMatrix;     // _strMatrix[line][col]
	
public:

  TMatrix  () : _rows(0), _cols(0), _val(0), _strMatrix(0) { }

  /**copy constructor.**/
  TMatrix  (const TMatrix& mat)
  {
    _rows = mat._rows;
    _cols = mat._cols;
    _strMatrix = NULL;
    _val = new double [_rows*_cols];
    memcpy(_val,mat._val,_rows*_cols*sizeof(double));
  }

  /**@brief Creates an array of doubles of size = rows*cols**/
  TMatrix ( unsigned int rows, unsigned int cols )
  {
    _val = new double [rows*cols];
    _strMatrix = NULL;
    _rows = rows;
    _cols = cols;
  }

  ~TMatrix () {
    if(_val) delete [] _val;
    if(_strMatrix) delete _strMatrix;
  }

  vector<map<int, string> >* get_strMatrix(){return _strMatrix;}

  /**Sets element at row i and column j to value val**/
  void set (unsigned int i, unsigned int j, double val) {
	  assert(i*j < _rows*_cols);
    _val[i*_cols + j] = val;
  }

  void assign (int val){
    memset(_val, val, _rows * _cols * sizeof(double));
  }

  /**Re-allocate the existing matrix with assigned rows and cols dimensions**/
  void reset (unsigned int rows, unsigned int cols) {
	  if(_val) delete [] _val;
    _val = new double [rows*cols];
	  _rows = rows;
    _cols = cols;
  }

  /**Reset the existing matrix to the new dimensions and copies the array**/
  void reset (unsigned int rows, unsigned int cols, double* array) {
	  reset(rows, cols);
    memcpy(_val, array, rows * cols * sizeof(double));
  }
  /**Reset members to zero state.*/
  void reset ( ){
    _rows = 0;
    _cols = 0;
    if(_val) delete [] _val;
    _val = NULL;
  }

  /**@brief Accessor to element at row i and column j**/
  double get (unsigned int i, unsigned int j) {
	  assert(i*_cols + j < _cols*_rows);
	  return _val[i*_cols + j];
  }
  /**Accessor to the whole array.*/
  double* get () const {return _val;}

  /**Accessor to the matrix dimensions.
    *@param dims an array of at least 2 elements to store the row [0] and column [1] numbers. May be NULL.
    *@return the total size of the matrix
   **/
  unsigned int get_dims (unsigned int* dims) {
	  if(dims) { dims[0] = _rows; dims[1] = _cols; }
	  return _rows*_cols;
  }

  /**Gives the number of rows.*/
  unsigned int getNbRows ( ) {return _rows;}

  /**Gives the number of columns.*/
  unsigned int getNbCols ( ) {return _cols;}

  /**Returns the number of elements in the matrix.*/
  unsigned int length    ( ) {return _rows*_cols;}

  /**Gives access to a column of the matrix.*/
  void getColumnView (unsigned int col, unsigned int n, double* array){
    assert(col < _cols);
    assert(n == _rows);
    for(unsigned int i = 0; i < _rows; ++i){
      array[i] = _val[i*_cols + col];
    }
  }

  /**Gives access to a row of the matrix.*/
  void getRowView (unsigned int row, unsigned int n, double* array){
    assert(row < _rows);
    assert(n == _cols);
    for(unsigned int i = 0, stride = row*_cols; i < _cols; ++i){
      array[i] = _val[stride + i];
    }
  }

  /**Adds a value to an element of the matrix.*/
  void plus (unsigned int i, unsigned int j, double value){
	  assert(i*j < _rows*_cols);
    _val[i*_cols + j] += value;
  }

  /**Substracts a value from an element of the matrix.*/
  void minus (unsigned int i, unsigned int j, double value){
	  assert(i*j < _rows*_cols);
    _val[i*_cols + j] -= value;
  }

  /**Multiply an element of the matrix by a value.*/
  void multi (unsigned int i, unsigned int j, double value){
	  assert(i*j < _rows*_cols);
    _val[i*_cols + j] *= value;
  }

  /**Divide an element of the matrix by a value.*/
  void divide (unsigned int i, unsigned int j, double value){
	  assert(i*j < _rows*_cols);
    if(value) _val[i*_cols + j] /= value;
    else      _val[i*_cols + j] = my_NAN;
  }

  void show_up(){
    message("TMatrix dimensions: rows = %i, columns = %i\n",_rows,_cols);
    for(unsigned int i = 0; i < _rows; i++) {
      for(unsigned int j = 0; j < _cols; j++){
        message("%.3f ",_val[i*_cols + j]);
      }
      message("\n");
    }
  }

  // functions to read matrixes from files or streams
  map<string, int> read_matrix(string filename);
  void             read_matrix_bracket (istream & IN);
  void             read_matrix_none (istream & IN);

  void read(string text);
};
#endif

