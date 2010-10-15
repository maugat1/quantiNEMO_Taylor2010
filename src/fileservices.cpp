/** @file fileservices.cpp
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


#include <iomanip>
#include <sstream>
#include <time.h>
#include <cmath>
#include "simcomponent.h"
#include "fileservices.h"
#include "filehandler.h"
#include "version.h"
#include "metapop.h"
#include "output.h"


using namespace std;
string FileServices::_simfolder;
int FileServices::_genetic_map_output;

// ----------------------------------------------------------------------------------------
// attach
// ----------------------------------------------------------------------------------------
void FileServices::attach ( FileHandler* FH )
{
  Service::attach(FH);
  _children.push_back(FH);
  FH->set_service(this);
  _popPtr->set_service(this);
}
// ----------------------------------------------------------------------------------------
// load
// ----------------------------------------------------------------------------------------
void FileServices::load ( SimComponent* sc ) {
  sc->loadFileServices(this);
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool FileServices::init (list< ParamSet* >&  params)
{
  char yn;
  bool ok;
  list< FileHandler* >::iterator HIT;

  _params = params;

  // while the filename is not properly set
  do{
    ok = true;
    for(HIT = _children.begin(); HIT != _children.end(); ++HIT) {
      ok &= (*HIT)->init();
    }

    if(!ok && !_overwriteFiles) {
			message(" Do you want to overwrite all the files that use it ? (y/n/s(kip)): \n");
			cin>>yn;

			switch(yn) {
				case 'y':
          ok = true;
          break;
        case 's':
          message(" Skipping this simulation\n");
          return false;
        default: {
          message(" Please give a new output filename: ");
          cin>>_basename;
        }
      }
    }
  } while(!ok && !_overwriteFiles);

  // write log file
  string logfile = _simfolder + _basename + ".log";
  FileHandler::check_path(_simfolder);
  ofstream FH(logfile.c_str(),ios::out);
  if(!FH){
    error("FileServices::init: could not open log file '%s'!\n",logfile.c_str());
  }
  else {
    save_simparams(params, FH);
    FH.close();
	}

  return true;
}


// ----------------------------------------------------------------------------------------
// save_simparams
// ----------------------------------------------------------------------------------------
/**Sets the folder name of the simulation*/
void FileServices::set_simfolder (string& name) {
  _simfolder = FileHandler::get_base_path() + name;
  if(!_simfolder.empty() && _simfolder[_simfolder.length()-1] != SEP) _simfolder += SEP;
}

// ----------------------------------------------------------------------------------------
// save_simparams
// ----------------------------------------------------------------------------------------
void FileServices::save_simparams(list< ParamSet* >&  params, ostream& FH)
{
  time_t t = time(NULL);
  strftime(_startTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&t));
  _startClock = clock();

	FH << "# quantiNemo    v[" << VERSION_DATE << "; " << VERSION_TIME <<"]\n#\n"
     << "# simulation started on " << _startTime << "\n\n";

  switch(_logfile_type){
    case 0:  // as input file
        for(list< ParamSet* >::iterator Pit = params.begin(); Pit != params.end(); ++Pit) {
          FH << "\n";
          (*Pit)->print(FH);
        }
        break;
    case 1: // minimal
        FH << "# minimal log\n\n";
        for(list< ParamSet* >::iterator Pit = params.begin(); Pit != params.end(); ++Pit) {
          FH << "\n";
          (*Pit)->print_minimal(FH);
        }
        break;
    case 2: // maximal
        FH << "# maximal log\n\n";
        for(list< ParamSet* >::iterator Pit = params.begin(); Pit != params.end(); ++Pit) {
          FH << "\n";
          (*Pit)->print_maximal(FH);
        }
        break;
  }
  FH<<endl;
}

// ----------------------------------------------------------------------------------------
// save_simparams
// ----------------------------------------------------------------------------------------
void FileServices::end_logfile()
{
  string logfile = _simfolder + _basename + ".log";
	ofstream FH(logfile.c_str(),ios::app);
  if(!FH){
    error("FileServices::init: could not output sim parameters to\"%s\"\n",logfile.c_str());
  }
  else {
    time_t t = time(NULL);
    strftime(_endTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&t));
    _endClock = clock();
    _elapsedTime = getElapsedTime(_endClock - _startClock);

    FH << "\n# simulation ended on " << _endTime
       << "\n#\n# duration of the simulation was " << _elapsedTime;

    FH.close();
  }
}

// ----------------------------------------------------------------------------------------
// getReplicateFileName
// ----------------------------------------------------------------------------------------
string& FileServices::getReplicateFileName ()
{
  _rep_filename = _basename + _popPtr->getReplicateCounter_();

  return _rep_filename;
}
// ----------------------------------------------------------------------------------------
// getGenerationReplicateFileName
// ----------------------------------------------------------------------------------------
string FileServices::getGenerationReplicateFileName ()
{
  string name = _basename + _popPtr->getGenerationCounter_g() + _popPtr->getReplicateCounter_r();
  return name;
}
// ----------------------------------------------------------------------------------------
// getBaseFileName
// ----------------------------------------------------------------------------------------
string& FileServices::getBaseFileName (){
  return _basename;
}
// ----------------------------------------------------------------------------------------
// getSimfolder
// ----------------------------------------------------------------------------------------
string& FileServices::getSimfolder (){
  return _simfolder;
}
// ----------------------------------------------------------------------------------------



