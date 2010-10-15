/** @file stathandler.h
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

#ifndef stathandlerbaseH
#define stathandlerbaseH

#include <list>
#include "handler.h"
#include "statservices.h"
#include "tmatrix.h"

class Metapop;
class Patch;
//#include "metapop.h"


class StatRecBase;
/**Base class of the StatHandler class, implements the Handler interface.
 * This class stores the Handler state and the list of links to the StatRecBase.
 * This list is duplicated in the StatHandler as a list of the StatRecorder templates.
 * It allows the StatService to access the stat recorders without knowledge of the
 * actual StatHandler and StatRecorder template type.
 */
class StatHandlerBase : public Handler {

private:
  /**Link to the StatRecorder list elements in the StatHandler derived class.*/
  list<StatRecBase*> _stats;

  /**Structure of the stats table in the stat recorders.*/
  unsigned int       _nrows, _ncols;

  unsigned int*       _index;        // don't delete it

  /**Link to the StatService.*/
  StatServices*      _service;

  /**Occurence of the stats recording as set by the user (see parameter "stat_log_time").*/
  unsigned int       _GenerationOccurrence;

  typedef list< StatRecBase* >::iterator STAT_IT;

protected:
  /**Link to the current population, set through the link to the StatService.*/
  Metapop* _pop;

  /**Replicate state of this Handler.*/
	static unsigned int       _current_replicate;

	/**Generation state of this Handler.*/
	static unsigned int       _current_generation;

	/** the following vectors allow to make a subselection of the patches
			it is possible to sample a subnumber of patches (all the other patches are never considered for stats) */
	static Patch**            _sample_pops;              // array of the patches to sample  (do not delete the array here!!!)
	static unsigned int       _sample_pops_size;         // number of patches to sample ( = _stat_pops.size())

	// the following vectors are in relation to the parameter _sample_pops!!! (OFFSx=0, PDISPx=1, ADLTx=2)
	static vector<Patch*>     _current_empty_pops[NB_AGE_CLASSES];  // current un-populated pops
	static vector<Patch*>     _current_pops[NB_AGE_CLASSES];        // current populated pops
	static vector<Patch*>     _current_popsS[2][NB_AGE_CLASSES];    // current populated pops seperated by sex and age: [sex][age]

	static unsigned int       _current_nbInds[NB_AGE_CLASSES];      // total number of individuals in the sampled patches
	static unsigned int       _current_nbIndsS[2][NB_AGE_CLASSES];  // total number of sampled patches separated by sex and age: [sex][age]
	static unsigned int       _current_nbPops[NB_AGE_CLASSES];      // current number of populated pops ( = _current_pops.size())

	inline unsigned int get_current_nbPops(age_idx i)          {return _current_nbPops[i];}
	inline double       get_current_nbPopsStat(age_idx i)      {return _current_nbPops[i];}  // return must be double for a stat...
	inline unsigned int get_current_nbPops(age_t i)            {return _current_nbPops[i-1];}

	inline unsigned int get_current_nbInds(age_t i)            {return _current_nbInds[i-1];}
	inline unsigned int get_current_nbIndsS(age_t i, sex_t SEX){return _current_nbIndsS[SEX][i-1];}

public:
	StatHandlerBase( ) { }

	virtual           ~StatHandlerBase     ( ) {  }

	///@name Accessors
  ///@{
  Metapop*          get_pop_ptr      ( )                 {return _service->get_pop_ptr();}

  void              set_service      (StatServices* srv) {_service = srv;}

  StatServices*     get_service      ( )                 {return _service;}

  unsigned int      getNbRecorders   ( )                 {return _stats.size();}

  unsigned int      get_nrows        ( )                 {return _nrows;}

  unsigned int      get_ncols        ( )                 {return _ncols;}

  unsigned int*     get_index        ( )                 {return _index;}

  list<StatRecBase*>& getStats       ( )                 {return _stats;}

  virtual void      add              (StatRecBase* rec)  {_stats.push_back(rec);}
  ///@}
  /**Empties the _stats list and calls clear() (defined in the derived class).*/
  virtual void      reset            ( );

  ///@name Stat recorder interface
  ///@{
  unsigned int      getStatRecIndex  (unsigned int i);

  void              setMean       ( );

	void              print_headers    (ostream& FH, unsigned int order);
	void              print_legend     (ostream& FH, unsigned int order);
	void              print_value      (ostream& FH, unsigned int i, unsigned int j);
	void              print_mean       (ostream& FH, unsigned int i);
	void              print_variance   (ostream& FH, unsigned int i);

  ///@}
  ///@name Handler implementation
  ///@{
  virtual bool      init             ( );
	virtual void      update           ( );
	virtual void      update_patch_states ( );
  ///@}
  ///@name StatHandler interface declaration
  ///@{
  virtual void      execute          ( ) = 0;

  virtual bool      setStatRecorders (const string& token) = 0;

  virtual string getName() = 0;

  virtual void      clear            ( ) = 0;
  ///@}
};

#endif //STATHANDLERBASE_H

