/** @file stat_rec_base.cpp
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


#include <cmath>
#include "stat_rec_base.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/
//
//                               ****** StatRecBase ******
// ----------------------------------------------------------------------------------------
// StatRecBase::~StatRecBase()
// ----------------------------------------------------------------------------------------
StatRecBase::~StatRecBase()
{
  if(_val){
    for(unsigned int i = 0; i < _rows; ++i){
      if(_val[i]) delete [] _val[i];
    }
    delete [] _val;
  }
}

// ----------------------------------------------------------------------------------------
// StatRecBase::init()
// ----------------------------------------------------------------------------------------
void StatRecBase::init()
{
  _title = "";
  _name = "";
  _age = ALL;
  _ordering = FLAT;
  _cur_row = 0;
  _cur_col = 0;
  _arg = 0;

  if(_rows) ARRAY::create_1D(_val, _rows, (double*)NULL);
}

// ----------------------------------------------------------------------------------------
// StatRecBase::set()
// ----------------------------------------------------------------------------------------
void StatRecBase::set(string T, string N, st_order Order, age_t AGE, unsigned int ARG)
{
  _title = T;
  _name = N;
  _ordering = Order;
  _age = AGE;
  _arg = ARG;

  switch(Order){
    case FLAT:  // values are stored for each generation and replicate
                setVal_func_ptr     = &StatRecBase::setVal_FLAT;
                getMean_func_ptr    = &StatRecBase::getMean_FLAT;    // across replciates
                getVar_func_ptr     = &StatRecBase::getVar_FLAT;     // across replciates
                break;
    case GEN:   // values are stored for each generation: replicate values are added
                setVal_func_ptr     = &StatRecBase::setVal_GEN;
                getMean_func_ptr    = &StatRecBase::get_GEN;
                getVar_func_ptr     = &StatRecBase::get_GEN;
                break;
    case RPL:   // values are stored for each replicate: generation values are added
                setVal_func_ptr     = &StatRecBase::setVal_RPL;
                getMean_func_ptr    = &StatRecBase::get_RPL;
                getVar_func_ptr     = &StatRecBase::get_RPL;
                break;
  }
}

// ----------------------------------------------------------------------------------------
// StatRecorder::setVal_FLAT
// ----------------------------------------------------------------------------------------
/** stats are stored for each generation and replicate separately
  * crnt_gen: is the real generation (starting at 1)
  * rpl_cntr: is the real replciate (starting at 1)
  * _val[gen][rpl]: the replicate arrays are created when needed
  */
bool StatRecBase::setVal_FLAT(int crnt_gen, int rpl_cntr, double value)
{
  _cur_col = rpl_cntr-1;

  // scroll to the correct position
  while(_index[_cur_row] != (unsigned int)crnt_gen){
    ++_cur_row;                     // it could be that a value is missing
    if(_cur_row >= _rows){          // do we have to change the column (new replicate)?
      _cur_row = 0;                 // reset the current column
      if(!_val[_cur_row]) ARRAY::create_1D(_val[_cur_row], _cols, (double)my_NAN);
    }
  }

  if(!_val[_cur_row]) ARRAY::create_1D(_val[_cur_row], _cols, (double)my_NAN);

  // gen 0 is the index of the inital population value (generation 0)
	//remember that the replicate and generation counters start with 1, not 0!!!
	_val[_cur_row][_cur_col] = value;

  ++_cur_col;

  return true;
}

// ----------------------------------------------------------------------------------------
// StatRecorder::setVal_GEN
// ----------------------------------------------------------------------------------------
/** Stats across replicates are summed up, i.e. there is a single stat for each generation.
  * crnt_gen: is the real generation (starting at 1)
  * rpl_cntr: is the real replciate (starting at 1)
  * that is the only place where rows and cols are inverted:
  * rows: replicate
  * cols: not used
  */
bool StatRecBase::setVal_GEN(int crnt_gen, int rpl_cntr, double value)
{
  while(_index[_cur_row] != (unsigned int)crnt_gen){
    ++_cur_row;                     // it could be that a value is missing
    if(_cur_row >= _rows){          // do we have to change the column (new replicate)?
      _cur_row = 0;                 // reset the current column
    }
  }

   // does a new row has to be added
  if(!_val[_cur_row]){
    _val[_cur_row] = new double[1];  // first time to be passed (either first replicate or previous replicates have gone extincted)
    _val[_cur_row][0] = value;
  }
  else{
    //add the replicate value to the generation cell of the _val matrix:
    _val[_cur_row][0] += value;
  }

  ++_cur_row;
  if(_cur_row >= _rows){          // do we have to change the column (new replicate)?
    _cur_row = 0;                 // reset the current column
  }

  return true;
}

// ----------------------------------------------------------------------------------------
// StatRecorder::setVal_RPL
// ----------------------------------------------------------------------------------------
/** Stats across generation are summed up, i.e. there is a sinlge stat for each replicate.
  * crnt_gen: is the real generation (starting at 1)
  * rpl_cntr: is the real replciate (starting at 1)
  */
bool StatRecBase::setVal_RPL(int crnt_gen, int rpl_cntr, double value)
{
  while(_cur_row != (unsigned int)rpl_cntr-1){
    ++_cur_row;                     // it could be that a value is missing
  }

   // does a new row has to be added
  if(!_val[_cur_row]){
    _val[_cur_row] = new double[1];  // first time to be passed (either first replicate or previous replicates have gone extincted)
    _val[_cur_row][0] = value;
  }
  else{
    //add the replicate value to the generation cell of the _val matrix:
    _val[_cur_row][0] += value;
  }

  ++_cur_row;
  if(_cur_row >= _rows){          // do we have to change the column (new replicate)?
    _cur_row = 0;                 // reset the current column
  }

  return true;
}

// ----------------------------------------------------------------------------------------
// StatRecorder::get
// ----------------------------------------------------------------------------------------
/** values are present for each generation and replicate
  * stats across replicates, i.e. per generation
  * gen_ind: generation index and NOT the generation itself!
  */
double StatRecBase::getMean_FLAT(const unsigned int& gen_index){
  if(!_val[gen_index]) return my_NAN;     // there are no values stored
  return ARRAY::mean(_val[gen_index], _cols);
}

double StatRecBase::getVar_FLAT(const unsigned int& gen_index){
  if(!_val[gen_index]) return my_NAN;     // there are no values stored
  return ARRAY::var(_val[gen_index], _cols);
}

// ----------------------------------------------------------------------------------------
// StatRecorder::get_GEN
// ----------------------------------------------------------------------------------------
/** stats of different replicats are summed up in the first position
  * stats summed up across replicates, i.e. per generation
	* gen_index: generation index and NOT the generation itself!
	* example: alive.rpl
  */
double StatRecBase::get_GEN(const unsigned int& gen_index){
  assert(gen_index<_rows);
  if(!_val[gen_index]) return 0;
  return *_val[gen_index];    // first position
}

// ----------------------------------------------------------------------------------------
// StatRecorder::get_RPL
// ----------------------------------------------------------------------------------------
/** stats of different generations are summed up in the first position
  * stats summed up across generations, i.e. per replicate
  * rpl_index: is not the replicate, but cur_repl - 1
  */
double StatRecBase::get_RPL(const unsigned int& rpl_index){
  assert(rpl_index<_cols);
  if(!_val[rpl_index]) return my_NAN;
  return *_val[rpl_index];    // first position
}


