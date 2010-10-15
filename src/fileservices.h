/** @file fileservices.h
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

#ifndef fileservicesH
#define fileservicesH

#include <list>
#include "service.h"
#include "param.h"

class Metapop;

class FileHandler;

/**A class to manage the files associated with each components of the simulation.

   Implements the Observer design pattern (is the concrete subject), stores the base filename of the simulation and
   updates the replicate filenames. It also performs files checking and saves the simulation parameters on init.
*/

class FileServices : public Service {

private:
  /**a pointer to the current Metapop*/
  Metapop*    _popPtr;

  /**the list of the FileHandler's registered by the SimComponent*/
  list< FileHandler* > _children;

  /**the file name associated with the current simulation replicate*/
  string _rep_filename;

  /**the base file name of the simulation, read from the init file (param "filename") */
  string _basename;

  /** folder name contining all the simulation outputs */
  static string _simfolder;

  /**the list of the simulation parameters*/
  list< ParamSet* > _params;

  /** should existing files be overwriten without asking? */
  bool      _overwriteFiles;
  int       _logfile_type;       // 0: as input, 1: minimal, 2: maximal
  static int _genetic_map_output; // 0: yes, if used, 2: no

  /** time parameters */
  char      _startTime[20],
            _endTime[20];
  clock_t   _startClock,
            _endClock;
  string    _elapsedTime;


public:

    typedef list< FileHandler* >::const_iterator file_it;

  FileServices ( ) : _popPtr(0), _rep_filename(""), _basename("") { }

  virtual ~FileServices ( ) { }

  virtual bool init ( ) {return false;}

  /**Checks if files with _basename already exist and save the simulation parameters in log files.
   * @param params a ref to the list of the current parameters of the simulation
   * @return true if the files check is ok
   * @return false if the user wants to skip this simulation
   */
  bool init (list< ParamSet* >&  params);
  /**
   * @return the pointer to current Metapop
   */
  virtual Metapop*   get_pop_ptr ( )      {return _popPtr;}

  /**Sets the Metapop reference.*/
  virtual void set_pop_ptr (Metapop* pop) {_popPtr=pop;}

  /**Sets the base file name of the simulation*/
  void set_basename (string& name) {_basename = name;}

  /**Sets the folder name of the simulation*/
  static void set_simfolder (string& name);

  /**Saves the current simulation parameters in log files.
   * @param params a ref to the list of the current parameters of the simulation
   * @param FH the handle to the log file
   */
  void save_simparams(list< ParamSet* >&  params, ostream& FH);
  void end_logfile();

  /**
   * @return the list of the current parameters of the simulation
   */
  list< ParamSet* >& get_params() {return _params;};

  file_it getFirst() {return _children.begin();}
  file_it getLast() {return _children.end();}
 // int     getNbFileHandler() {return _children.size();}

  /**
   * @return the base file name
   */
  string&  getBaseFileName ();

  /**
   * @return the simulation folder
   */
  static string&  getSimfolder ();

  /**
   * @return the current replicate file name
   */
  string&  getReplicateFileName ();

  bool getOverwriteFiles()             {return _overwriteFiles;}
  void setOverwriteFiles(bool b)       {_overwriteFiles=b;}
  void set_logfile_type(int val)       {_logfile_type=val;}
  static void set_genetic_map_output(int val) {_genetic_map_output=val;}
  static int  get_genetic_map_output()        {return _genetic_map_output;}

  /**
   * @return the current file name with generation and replicate counters
   */
  string  getGenerationReplicateFileName ();

  /**Tells the SimComponent to load its file handlers.
   *  @param sc the SimComponent
   */
  virtual void load ( SimComponent* sc );

  /**Attaches the FileHandler to the current list (_children) of the FileServices.
   *  @param FH the FileHandler
   */
  virtual void attach ( FileHandler* FH );

  /**Clears the list of FileHandlers. */
  virtual void reset ( ) {
 //   for(list< FileHandler* >::iterator HIT  = _children.begin(); HIT != _children.end(); ++HIT) {
 //     (*HIT)->reset();
 //   }
    Service::reset();
    _children.clear();
  }
};
#endif //FILESERVICES_H

