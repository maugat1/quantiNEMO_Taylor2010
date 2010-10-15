/** @file simulation.h
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


#ifndef simulationH
#define simulationH

#include <list>
#include <map>
#include <vector>
#include <time.h>
#include "basicsimulation.h"
#include "fileservices.h"
#include "statservices.h"
#include "metapop.h"


/**Performs the setup of the Metapop and SimComponents and runs the simulations.
 * A SimRunner brings together the basic simulation components and a metapopulation on which they act.
 * It perfoms the setups necessary to have a Metapop ready for the simulation and runs the different simulations
 * stored in its ParamManager base class. Also hosts the file and stat services.
 */
class SimRunner: public SimBuilder {
private:
  Metapop*                    _thePop;
  string                      _logfile;
  string                      _postexec_script;
  bool                        _do_postexec;
	char                        _startTime[20],
                              _endTime[20];
  string                      _simElapsedTime;

public:

  FileServices					   _FileServices;
  StatServices					   _StatServices;

	static RAND r;           // random number generator

  /**Cstor.
    @callgraph **/
  SimRunner   (){
    Metapop() ;   // no idea why this is needed: If removed Borland compiler returns a linker error, but for gcc this would work as well
  }

  /**Dstror.
    @callgraph **/
                                   ~SimRunner            ( );

  list< map<string, string> >& readInputFile( int ARGC, char **ARGV );
  void set_seed(map< string,string >& simparams);

  void loadDefaultTemplates(map<string, string>& inputParams);
  void loadDefaultTemplates(map<string, vector<string> >& inputParams);

  /**Checks simulation parameters and init the FileServices with the base filename. @callgraph**/
  bool                             init                   ( );

  /**Performs the setup of the Metapop and the services.
    *Builds the list of the simulation parameters and load the components that have their ParamSet in the set state.
    *Builds the population with the selected TraitPrototype and LifeCycleEvent, register the various services and init
    *the StatServices.
    *@param simparams the hashtable containing the parameters and their arguments parsed from the init file
		*@callgraph
		*/
	bool                             setup                  (map< string,string >& simparams);

	void                             print_info               ( );

  /**Resets all the parameters to the unset state, resets the services.
    @callgraph**/
  void                             reset                  ( );

  /**Register the FileHandler and StatHandler attached to the SimComponent. @param cmpt a SimComponent**/
  void                             register_services      (SimComponent* cmpt);

  /**Resets the FileServices and StatServices.
    @callgraph**/
  void                             reset_services         ( );

  /**Returns the FileServices. **/
  FileServices*					   get_FileServices       ( )                      {return &_FileServices;}

  /**Returns the StatServices. **/
  StatServices*					   get_StatServices       ( )                      {return &_StatServices;}

		/**First loop of the simulation, performs the simulations stored in the ParamManager base class.
    @callgraph **/
  bool                             run                    ( );

  /**not yet fully functional.
    @callgraph **/
  bool                             run_event              (string& name);

  /**Accessor to the pop ptr.
   * @return the pop ptr
   **/
  Metapop*                         get_pop                ( )                      {return _thePop;}

  /**Calls the Metapop init procedure with current traits and LCEs.
    @callgraph **/
  bool                             build_pop              ( );

  /**Calls the the funstions to initialzing the genetic map.
    @callgraph **/
  void                             initGeneticMap         ( );

  /**Detach a pop from the simulation. Removes it from the components list.
   * @param pop ptr to the pop object   Samuel
   **/
  void                             detach_pop             () {
    if(_thePop) delete _thePop;
    this->_components.pop_front();
  }

  /**Attach a pop to the simulation. Adds it to the components list.
   * @param pop ptr to the pop object
   **/
  void                             attach_pop             (Metapop* pop=NULL) {
   if(!pop) pop= new Metapop();
/*     _thePop = pop;
    this->_components.push_back(_thePop);
*/  }

  void printLogHeader ();
  void printLog ();
  void runPostExec();
  int get_nbTraits(string name, map<string, string>& inputParams);
  int get_nbTraits(string name, map<string, vector<string> >& inputParams);
};

#endif

