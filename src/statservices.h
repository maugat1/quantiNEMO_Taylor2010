/** @file statservices.h
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

#ifndef statservicesH
#define statservicesH

#include <list>
#include "service.h"
#include "stat_rec_base.h"

class Metapop;

class StatHandlerBase;
/**The Service class used to manage the StatHandler objects.*/
class StatServices : public Service {

private:

  Metapop* _popPtr;
  list< StatHandlerBase* > _children;
  string _statArg;
	unsigned int 	_occurrence;
	unsigned int 	_tot_occurrence;
	unsigned int* _index;       // pointer to an array in the class LCE_StatServicesNotifier (don't delete it!)

	int           _save_choice; // 0: all; 1: each rep; 2: mean and var, 3: only mean; 4: only var; 5: nothing

	string 				_file_stats;  // file name for the stats (each generation and replicate)
	string				_file_mean;   // file name for the mean across replicates
	string			  _file_var;    // file name for the variance across replciates

public:

	void set_file_name_stats(string s) {_file_stats = s;}
	void set_file_name_mean(string s)  {_file_mean = s;}
	void set_file_name_var(string s)   {_file_var = s;}
	void set_save_choice(const int& i) {_save_choice = i;}

	string get_file_name_stats( ) const {return _file_stats;}
	string get_file_name_mean( )  const {return _file_mean;}
	string get_file_name_var( )   const {return _file_var;}
	int    get_save_choice( )     const {return _save_choice;}

	typedef list< StatHandlerBase* >::const_iterator stat_it;

  StatServices ( ) : _popPtr(0), _occurrence(0), _tot_occurrence(0){ }

  virtual ~StatServices ( ) { }

  virtual bool init ( );

  Metapop* get_pop_ptr          ( )      {return _popPtr;}

  void set_pop_ptr              (Metapop* pop) {_popPtr=pop;}

  void set                      (string& str, unsigned int occ, unsigned int tot_occ, unsigned int* index) {
    _statArg = str;
    _occurrence = occ;
    _tot_occurrence = tot_occ;
    _index = index;
  }

  unsigned int getOccurrence    ( ) {return _occurrence;}
  unsigned int getTotOccurrence ( ) {return _tot_occurrence;}
  unsigned int* getIndex        ( ) {return _index;}
  void         setOccurrence    (const unsigned int& occ) {_occurrence = occ;}
  void         setTotOccurrence (const unsigned int& occ) {_tot_occurrence = occ;}

	void printStatHeaders         (ostream& FH, unsigned int order);
  void printStatLegend          (ostream& FH, unsigned int order);

	void printStatValue           (ostream& FH, unsigned int i, unsigned int j);
	void printStatMean            (ostream& FH, unsigned int i);
	void printStatVariance        (ostream& FH, unsigned int i);

  stat_it getFirst () {return _children.begin();}

  stat_it getLast () {return _children.end();}

  unsigned int getStatRecIndex  (unsigned int i);
  /*
   *
   */
  virtual void notify ();

  /**
	*  tell the SimComponent to load its stat handlers
	*  @param sc the SimComponent
	*/
  virtual void load ( SimComponent* sc );

  /**
	*  attach the StatHandler to the current list (_children) of the StatServices
	*  @param SH the StatHandler
	*/
  virtual void attach ( StatHandlerBase* SH);

  /**
	*  clear the list of StatHandler
	*/
  virtual void reset ( );
};
#endif //STATSERVICES_H

