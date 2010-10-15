/** @file stat_rec_base.h
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

#ifndef stat_rec_baseH
#define stat_rec_baseH

#include <vector>
#include <list>
#include <functional>
#include "types.h"
#include "output.h"


/**Base class for the StatRecorder's, stores the stat values in a matrix.
 * The way the stat values are stored in the _val matrix is set by the _ordering flag (see StatRecorder::setVal()).
 * one stat per object
 * **/
class StatRecBase {

private:
  /**The title of the stat recorder, longer and more explicite than the name.*/
  string  _title;

  /**Name of the stat, should be short (20 char) and R compliant (no '-', '+', ' ')*/
  string  _name;

  /**Dimensions of the recording matrix*/
  unsigned int _rows,    // 1 dimension: replicate
               _cols;    // 2 dimension: generation

  /**A argument to be passed to one of the function variable stored in the StatRecorder structure.*/
  unsigned int _arg;

  /**The age class for which this stat applies.*/
  age_t        _age;

  /**A flag specifying the way the stat values are recorded in the matrix, is initialized to FLAT by default.*/
  st_order     _ordering;

  /**The recording matrix, stores the stat values in a [_rows X _cols] matrix, i.e. first indice for row, second for column.*/
  double**     _val;
  unsigned int _cur_row;     // 1 dimension: replicate
  unsigned int _cur_col;     // 2 dimension: generation

  /**The values of the index value of each row, might be the generation number (with GEN and FLAT ordering) or the replicate number.*/
  unsigned int* _index;    // don't delete it here

public:

  StatRecBase  ( ) : _title(), _name(), _rows(0), _cols(0), _arg(0), _age(ALL),
                     _ordering(FLAT), _val(0), _index(0){ }

  StatRecBase  (unsigned int n_rows, unsigned int n_cols) : _rows(n_rows), _cols(n_cols),
    _val(0), _index(0)
  { init(); }

  ~StatRecBase ( );

  /**Creates the _val matrix if its dimensions have been set and sets each elements to "NaN".*/
  void init    ( );

  /**Sets the recorder attributes.
	* @param T the stat title
	* @param N the stat name (headers in the output file)
	* @param Order stat table ordering flag
	* @param AGE the age class for which the stat will be recorded
	* @param ARG the argument to pass to the S function
  **/
  void set     (string T, string N, st_order Order, age_t AGE, unsigned int ARG);

  ///@name Accessors
  ///@{
  void setName                      (string N)            {_name = N;}

  /**@param i the number of rows
     @param j the number of columns */
  void setTableDims                 (int i, int j)             {_rows = i; _cols = j;}
  void setIndex                     (unsigned int* i)           {_index = i;}
  st_order     getOrdering           ( )                        {return _ordering;}
  string       getTitle              ( )                        {return _title;}
  string       getName               ( )                        {return _name;}
  age_t        getAge                ( )                        {return _age;}
  unsigned int getArg                ( )                        {return _arg;}
  unsigned int getRows               ( )                        {return _rows;}
  unsigned int getCols               ( )                        {return _cols;}

  /**Returns the jth element of the ith row of the _val matrix.*/
  double getVal                     (unsigned int i, unsigned int j) {
    assert(i<_rows && j<_cols);
    if(!_val[i]) return my_NAN;    // no data available
    return _val[i][j];
  }

  double   getMean(const unsigned int& index) {return (this->*getMean_func_ptr)(index);}
  double   (StatRecBase::*getMean_func_ptr)(const unsigned int& index);
  double   getVar(const unsigned int& index) {return (this->*getVar_func_ptr)(index);}
  double   (StatRecBase::*getVar_func_ptr)(const unsigned int& index);
  double   getMean_FLAT             (const unsigned int& gen_index);
  double   getVar_FLAT              (const unsigned int& gen_index);
  double   get_GEN                  (const unsigned int& gen_index);
  double   get_RPL                  (const unsigned int& rpl_index);

  unsigned int getIndex             (unsigned int i)           {return _index[i];}

  /**Returns the number of elements in _index (should be equal to _rows!).*/
  unsigned int getIndexSize         ()                         {return _rows;}

  ///@}
  bool   (StatRecBase::*setVal_func_ptr)(int crnt_gen, int rpl_cntr, double value);
  bool   setVal_FLAT                (int crnt_gen, int rpl_cntr, double value);
  bool   setVal_GEN                 (int crnt_gen, int rpl_cntr, double value);
  bool   setVal_RPL                 (int crnt_gen, int rpl_cntr, double value);
};

/**Stores the pointers to the StatHandler's stat functions.*/
template <class S>
class StatRecorder : public StatRecBase
{
private:
  /**Pointer to a 'stat getter' function of S using no argument.*/
  double (S::* _getStat)     		(void);
  double (S::* _getStatAGE)     (const age_t& AGE);

  /**Pointer to a 'stat getter' function of S using a bool argument.*/
  double (S::* _getStatBool) 		(bool);
	double (S::* _getStatBoolAGE) (bool, const age_t& AGE);

  /**Pointer to a 'stat getter' function of S using a unsigned int argument.**/
  double (S::* _getStatUI)   		(unsigned int);
  double (S::* _getStatUIAGE) 	(unsigned int, const age_t& AGE);

public:


  StatRecorder() : _getStat(0),    _getStatBool(0),    _getStatUI(0),
                   _getStatAGE(0), _getStatBoolAGE(0), _getStatUIAGE(0) {}

  /**
   * @param n_rows number of rows of the stat table (usually nbr of records per replicates)
   * @param n_cols number of columns of the stat table (usually nbr of replicates)
  **/
  StatRecorder(unsigned int n_rows,unsigned int n_cols, unsigned int* index);

//  ~StatRecorder() {}
  /**@brief sets the recorder attributes
   * @param Title 						the stat title
   * @param Name 							the stat name (headers in the output file)
   * @param Order 						stat table ordering flag
   * @param AGE 							age on which the stat should be processed
   * @param ARG 							the argument to pass to the S function
   * @param getSt 						function ptr to a S getter
   * @param getStBoolArg 			function ptr to a S getter with boolean argument
   * @param getStUintArg 			function ptr to a S getter with unsigned int argument
   * @param getStAGE 					function ptr to a S getter with age argument
   * @param getStBoolArgAGE 	function ptr to a S getter with boolean and age argument
   * @param getStUintArgAGE		function ptr to a S getter with unsigned int and age argument
   **/
  void set(string Title,string Name,st_order Order,age_t AGE,unsigned int ARG,
  				double(S::* getSt)(void),
          double(S::* getStBoolArg)(bool)=0,
          double(S::* getStUintArg)(unsigned int)=0,
  				double(S::* getStAGE)(const age_t&)=0,
          double(S::* getStBoolArgAGE)(bool, const age_t&)=0,
          double(S::* getStUintArgAGE)(unsigned int, const age_t&)=0);

//  void setGetStat                   (double(S::* gst)(void))   {_getStat = gst;}
  /**@brief launch the linked stat function and store the result in the stat table
   * @return true if at least one getter ptr is != 0 and can be executed
   * @return false otherwise
   * @param AGE age on which the stat should be processed
   * @param crnt_gen the current running generation
   * @param rpl_cntr the current running replicate
   * @param StatHandler an instance of S
  **/
	bool setVal (age_t AGE,int crnt_gen,int rpl_cntr, S* StatHandler);
	bool setValDirectly (age_t AGE,int crnt_gen,int rpl_cntr, S* StatHandler,  ostream& FH);

};
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
//
//                               ****** StatRecorder ******

// ----------------------------------------------------------------------------------------
// StatRecorder::StatRecorder()
// ----------------------------------------------------------------------------------------
template <class S> StatRecorder<S>::StatRecorder(unsigned int n_rows,unsigned int n_cols, unsigned int* index)
{
  setTableDims(n_rows,n_cols);
  setIndex(index);

  //init the params, set all table cells to "NaN"
  init();
}
// ----------------------------------------------------------------------------------------
// StatRecorder::set()
// ----------------------------------------------------------------------------------------
template<class S>
void StatRecorder<S>::set(string T, string N, st_order Order, age_t AGE,unsigned int ARG,
											double (S::* getSt) (void),
                      double (S::* getStBoolArg) (bool),
											double (S::* getStUintArg)(unsigned int),
											double (S::* getStAGE) (const age_t&),
                      double (S::* getStBoolArgAGE) (bool, const age_t&),
											double (S::* getStUintArgAGE)(unsigned int, const age_t&))
{
  StatRecBase::set(T, N, Order, AGE, ARG);

  _getStat        = getSt;
  _getStatBool    = getStBoolArg;
  _getStatUI      = getStUintArg;
  _getStatAGE     = getStAGE;
  _getStatBoolAGE = getStBoolArgAGE;
  _getStatUIAGE   = getStUintArgAGE;
}

// ----------------------------------------------------------------------------------------
// StatRecorder::setVal
// ----------------------------------------------------------------------------------------
/** store te sta in the db */
template<class S>
bool StatRecorder<S>::setVal(age_t AGE, int crnt_gen, int rpl_cntr, S* StatHandler)
{
  double statValue;
	age_t age = getAge();

  // test if the stat has to be computed for this age
  if(age & AGE) {
    #ifdef _DEBUG
      message(" %s\n",getName().c_str());
    #endif

    //get the value:
		if(_getStat)           		statValue = (StatHandler->*_getStat)       ();
		else if(_getStatBool)  		statValue = (StatHandler->*_getStatBool)   ((bool)getArg());
		else if(_getStatUI)    		statValue = (StatHandler->*_getStatUI)     (getArg());
    else if(_getStatAGE)      statValue = (StatHandler->*_getStatAGE)    (age);
		else if(_getStatBoolAGE)  statValue = (StatHandler->*_getStatBoolAGE)((bool)getArg(), age);
    else if(_getStatUIAGE)    statValue = (StatHandler->*_getStatUIAGE)  (getArg(), age);
    else{
      error("StatRecorder::setVal: no _getStat fct ptr !!\n");
      return false;
    }
	  return (this->*setVal_func_ptr)(crnt_gen, rpl_cntr, statValue);
  }

  return true; //get out not complaining
}

// ----------------------------------------------------------------------------------------
// StatRecorder::setVal
// ----------------------------------------------------------------------------------------
/** this function writes the stt directly to the file (file is just opended fo rhte output and then closwed again */
template<class S>
bool StatRecorder<S>::setValDirectly(age_t AGE, int crnt_gen, int rpl_cntr, S* StatHandler, ostream& FH)
{
	double statValue;
	age_t age = getAge();

	// test if the stat has to be computed for this age
	if(age & AGE) {
		#ifdef _DEBUG
			message(" %s\n",getName().c_str());
		#endif


		//get the value:
		if(_getStat)           		statValue = (StatHandler->*_getStat)       ();
		else if(_getStatBool)  		statValue = (StatHandler->*_getStatBool)   ((bool)getArg());
		else if(_getStatUI)    		statValue = (StatHandler->*_getStatUI)     (getArg());
		else if(_getStatAGE)      statValue = (StatHandler->*_getStatAGE)    (age);
		else if(_getStatBoolAGE)  statValue = (StatHandler->*_getStatBoolAGE)((bool)getArg(), age);
		else if(_getStatUIAGE)    statValue = (StatHandler->*_getStatUIAGE)  (getArg(), age);
		else{
			error("StatRecorder::setValDirectly: no _getStat fct ptr !!\n");
			return false;
		}

		// write the stat directly to the file
		if(statValue == my_NAN) FH << "\t" << "NaN";
		else                    FH << "\t" << statValue;
	}

	return true; //get out not complaining
}


#endif

