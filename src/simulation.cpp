/** @file simulation.cpp
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
#include <cerrno>
#include "metapop.h"
#include "fileservices.h"
#include "statservices.h"
#include "lce_misc.h"
#include "lce_breed.h"
#include "lce_disperse.h"
#include "ttneutral.h"
#include "ttquanti.h"
#include "random.h"
#include "output.h"
#include "version.h"
#include "fileparser.h"
#include "simulation.h"

#ifndef LONG_MAX
	#define LONG_MAX 2147483647L
#endif
#ifndef ULONG_MAX
	#define ULONG_MAX (LONG_MAX * 2UL + 1)
#endif


RAND SimRunner::r;

//----------------------------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------------------------
/*
SimRunner::SimRunner   (): _thePop(0) {
  attach_pop();
}
*/
//----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool SimRunner::init ( )
{
	_FileServices.set_basename(this->_paramSet.getArg("filename"));

	if(!this->_paramSet.isSet("folder")){
    // generat the default folder name with current date
    char datetime[20];
    time_t t = time(NULL);
    strftime(datetime, 20, "%Y-%m-%d_%H-%M-%S", localtime(&t));
		string name = "simulation_";
    name += datetime;
    _FileServices.set_simfolder(name);
  }
  else _FileServices.set_simfolder(this->_paramSet.getArg("folder"));     // folder name is set


  _logfile      = FileHandler::get_base_path() + _paramSet.getArg("logfile");
  _FileServices.set_logfile_type((int)_paramSet.getValue("logfile_type"));
  _FileServices.setOverwriteFiles((bool)this->_paramSet.getValue("overwrite"));
  _FileServices.set_genetic_map_output((int)_paramSet.getValue("genetic_map_output"));

  printLogHeader();

  _postexec_script = _paramSet.getArg("postexec_script");
  _do_postexec = (_postexec_script.empty() ? false : true);

  return true;
}

//----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
SimRunner::~SimRunner ()
{
 //  message("SimRunner::~SimRunner\n");
  reset();
  detach_pop();    // somwhere the metapop has to be deleted!
}

//----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void SimRunner::reset ( )
{
  //reset all the parameter to the "unset" state
  list<ParamSet*>::iterator current_paramset = this->_allParams.begin();
  while(current_paramset != this->_allParams.end()) {
    (*current_paramset)->reset();
    current_paramset++;
  }

  reset_services();
}
//----------------------------------------------------------------------------------------
// reset_services
// ----------------------------------------------------------------------------------------
void SimRunner::reset_services ( )
{
  _FileServices.reset();
  _StatServices.reset();
}
//----------------------------------------------------------------------------------------
// register_services
// ----------------------------------------------------------------------------------------
void SimRunner::register_services (SimComponent* cmpt)
{
  _FileServices.load(cmpt);
  _StatServices.load(cmpt);
}
//----------------------------------------------------------------------------------------
// build_pop
// ----------------------------------------------------------------------------------------
bool SimRunner::build_pop ()
{
  if(!_thePop->init(this->build_currentTraits(), this->build_currentLifeCycle())) return false;

  initGeneticMap();
	return true;
}

//----------------------------------------------------------------------------------------
// init geneti map
// ----------------------------------------------------------------------------------------
/** initialises the genetic map */
void SimRunner::initGeneticMap()
{
  // create the static super chromosome
  if(!TTraitProto::initStaticGeneticMap()) return;

  // and the positions for each trait
  map<string, TTraitProto* >::iterator trait = this->_currentTraits.begin();
  for(;trait != this->_currentTraits.end(); ++trait) {
    trait->second->initTraitGeneticMap();
  }
}

//----------------------------------------------------------------------------------------
// setup
// ----------------------------------------------------------------------------------------
bool SimRunner::setup (map< string,string >& simparams)
{
  //first reset all paramSets and services (clear the stat handlers)
  reset();

  _FileServices.set_pop_ptr(_thePop);
  _StatServices.set_pop_ptr(_thePop);

  //build the list of parameters from the record:
  #ifdef _DEBUG
   message("SimRunner::run:building current params\n");
  #endif

	if(!this->build_currentParams(simparams)){
		error("SimRunner::run:couldn't build current params\n");
		return false;
	}

	//init the sim and pop params
	_thePop->set_replicates((unsigned int)_paramSet.getValue("replicates"));
	_thePop->set_generations((unsigned int)_paramSet.getValue("generations"));
	if( !(init()) || !(build_pop()) ) return false;

	//load the stats and files handlers of the Metapop:
	register_services(_thePop);

	// register the traits
	map<string, TTraitProto* >::iterator trait = this->_currentTraits.begin();
	for(; trait != this->_currentTraits.end(); ++trait) {
		register_services(trait->second);
	}

	// check if the sequence of life cycle events makes sense
	checkLCEconsistency();

	// register the life cycle events
	map< int, LCE* >::iterator iterLCE = this->_currentLifeCycle.begin();
	for(; iterLCE != this->_currentLifeCycle.end(); ++iterLCE) {
		register_services(iterLCE->second);
	}

	//StatServices: build the lists of stat recorders:
	if( !(_StatServices.init()) ) return false;

	print_info();

	return true;
}

//----------------------------------------------------------------------------------------
// print_info
// ----------------------------------------------------------------------------------------
void SimRunner::print_info(){
	message("\nSETTINGS");

	// simulation settings
	message("\n  Simulation:");
	message("\n    %i generations", _thePop->getGenerations());
	message("\n    %i replicates\n", _thePop->getReplicates());

	// traits
	message("\n  Loaded traits: ");
//	map<string, TTraitProto* >::iterator trait = this->_currentTraits.begin();
//	for(int i=0; trait != this->_currentTraits.end(); ++trait, ++i) {
//		message("\n    %s", trait->second->get_info().c_str());    // problem: alphabetically ordered
//	}

  list< TTraitProto* >::iterator qtrait = this->_TTrait_Templates.begin();
  for(;qtrait != this->_TTrait_Templates.end(); ++qtrait) {
    if((*qtrait)->get_paramset()->isSet() ){       // _currentTraits cannot be used since there the traits are sorted alphabetically
      message("\n    %s", (*qtrait)->get_info().c_str());
    }
	}

	message("\n");

	// life cycle events
	message("\n  Life cycle sequence:");
	map< int, LCE* >::iterator  iterLCE = this->_currentLifeCycle.begin();
	LCE_Disperse* pLCE;
	for(int i=1; iterLCE != this->_currentLifeCycle.end(); ++iterLCE, ++i) {
		message("\n    %i. %s", i, iterLCE->second->get_event_name().c_str());

		// find the migrate LCE
		pLCE = dynamic_cast<LCE_Disperse*>(iterLCE->second);
		if(pLCE){
			_thePop->set_pDisperse_LCE(pLCE);
		}
	}
	assert(_thePop->get_pDisperse_LCE());       // the disperse LCE must be present!
	message("\n");

	// metapopulation
	message("\n  Metapopulation:");
	message("\n    %i populations", _thePop->getPatchNbr());
	message("\n    Migration model: %s", _thePop->get_pDisperse_LCE()->get_disp_model_str().c_str());
	message("\n    Mating system: %s", _thePop->get_pBreed_LCE()->getMatingSystem_str().c_str());

  // genetic map
	message("\n\n  Genetic map:");
  TTraitProto::printGeneticMapInfo();


	message("\n");


}

//----------------------------------------------------------------------------------------
// read input file ("quantiNemo.ini")
// ----------------------------------------------------------------------------------------
/** the param ARGV[0] was changed such as to be the folder name of the executable (working directory) */
list< map<string, string> >&
SimRunner::readInputFile( int ARGC, char **ARGV ){
	FileParser Reader;
  map< string, vector<string> >* pParams;

  // check if the name of the settings file is passed
  if (ARGC == 1){  // if no argument is passed use the default file name
    // is there a default settings file?
    ARGV[1] = "quantiNemo.ini";
		ifstream FILE(((string)ARGV[0]+ARGV[1]).c_str());
    if(!FILE){    // no default settings file
      message("\nPlease enter the settings file name: ");
      cin >> ARGV[1];
    }
    else FILE.close();
    ARGC = 2; // now we have a settings file name
  }

	// read the settings files
  for (int i = 1; i < ARGC; ++i){
	  message("\nReading settings file '%s' ...", ARGV[i]);
		pParams = &Reader.read(ARGV[i], ARGV[0]);         // read the settings file
		loadDefaultTemplates(*pParams);                   // create the objects in order to check the names
		check_parameter_names(*pParams);                  // consistency check of the input parameters
		build_records(*pParams);                          // build the different sim records
	}
  if(ARGC == 2) message("\nReading settings file done!\n", ARGC, get_simRecords().size());
  else          message("\nReading settings file done (%i file(s) and %i batch simulations)!\n", ARGC-1, get_simRecords().size());

  return get_simRecords();
}

//----------------------------------------------------------------------------------------
// set_seed
// ----------------------------------------------------------------------------------------
/** set the seed (_paramSet is not yet set!)
	* seed may be a single unsigned int or an array of unsigned ints
	*/
void SimRunner::set_seed (map< string,string >& simparams)
{
	unsigned long* seed=NULL;
	unsigned int nbSeed=0;
  bool seedFromTime = true;

	// check if the seed is set
	map<string, string>::iterator iter = simparams.find("seed");
	if(iter != simparams.end()){                   // if seed is given by input
    string text = iter->second;
    if(text != "-1"){                            // if default seed (-1) was given stop
      // seed is given by input
      seedFromTime = false;

		  if(text[0] == '{'){       		// seed is given as an array
			  vector<unsigned long> m = STRING::strMatrix2vector<unsigned long>(text);
			  nbSeed = m.size();
			  seed = new unsigned long[nbSeed];
			  for(unsigned int i=0; i<nbSeed; ++i){
				  seed[i] = m[i];
			  }
		  }
		  else{                     		// seed is given as a single value
			  nbSeed = 1;
			  seed = new unsigned long[nbSeed];
			  seed[0] = STRING::str2int<unsigned long>(iter->second);
		  }
		  #ifdef _DEBUG
		   message("Random generator seed taken from the setting file: %s\n", text.c_str());
		  #endif
    }
	}

	if(seedFromTime) {             							// if initialized by time
		nbSeed = 2;
		seed = new unsigned long[nbSeed];
		seed[0] = time(NULL) % ULONG_MAX; // seconds since January 1, 1970.
		seed[1] = clock() % ULONG_MAX;    // number of clock ticks elapsed since the program started.
		simparams["seed"] = "{" + STRING::int2str(seed[0]) + " " + STRING::int2str(seed[1]) + "}";
		#ifdef _DEBUG
		 message("Random generator is initialized with a time seed: %i %i\n", seed[0], seed[1]);
		#endif
	}

	// initialize random generator
	r.set_seed(seed, nbSeed);

	if(seed) delete[] seed;
}

//----------------------------------------------------------------------------------------
// run
// ----------------------------------------------------------------------------------------
bool SimRunner::run ( )
{

	time_t t;
	string text;
	unsigned int sim, simnbre = _simRecords.size();


  //first loop: perform all simulations contained in _simRecords:
  //---------------------------------------------------------------------------------------
  list< map<string, string> >::iterator currentSim = _simRecords.begin();
  for(sim=1; currentSim != _simRecords.end(); ++sim, ++currentSim) {
		if(simnbre>1) message("\n\n----- SIMULATION %i/%i -----\n",sim,simnbre,_simElapsedTime.c_str());
		t = time(NULL);
		strftime(_startTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&t));

    loadDefaultTemplates(*currentSim);
		set_seed(*currentSim);
    //clear and build the lists of traits, LCEs, stats and files handlers, etc.
		if(!setup(*currentSim))  return false;

    //init the file services:
    //->files check, false means the user wants to skip this simulation.
		//->save simparameters in log files
		if( !(_FileServices.init(this->_currentParams)) ) {
      continue;
    }
		clock_t start = clock();

    //run the simulation
    //-------------------------------------------------------------------------------------
    _thePop->Replicate_LOOP();
    //-------------------------------------------------------------------------------------

		clock_t stop = clock();
    t = time(NULL);
    strftime(_endTime, 20, "%d-%m-%Y %H:%M:%S", localtime(&t));
    _simElapsedTime = getElapsedTime(stop - start);
		printLog();
		_FileServices.end_logfile();     // wirte the duration of the simulation to the log file

		if(simnbre==1) message("\n\n----- SIMULATION done (CPU time: %ss) -----\n",_simElapsedTime.c_str());
		else message("\n\n----- SIMULATION %i/%i done (CPU time: %ss) -----\n",sim,simnbre,_simElapsedTime.c_str());
	}
	//---------------------------------------------------------------------------------------
  if(_do_postexec) runPostExec();

  return true;
}

//----------------------------------------------------------------------------------------
// run_event
// ----------------------------------------------------------------------------------------
bool SimRunner::run_event (string& name)
{
  //before calling that function, first set the params and call setup()!!

  LCE* event = this->get_current_event(name);

  if(!event) {
    error("SimRunner::run_event:event \"%s\" not found in set events!\n",name.c_str());
    return false;
  }

  if( !(event->get_paramset()->isSet()) ) return false;//highly unlikely!!

  event->execute();

  return true;
}

// ----------------------------------------------------------------------------------------
// SimRunner::printLogHeader()
// ----------------------------------------------------------------------------------------
void SimRunner::printLogHeader()
{
  ofstream FH;
  ifstream IF;

  IF.open(_logfile.c_str(),ios::in);
  if(IF){
    IF.close();
    return;
  }

  FH.open(_logfile.c_str(),ios::out);
  if(!FH) {
    error("could not create simulation logfile \"%s\"\n",_logfile.c_str());
    return;
  }

  FH<<"--- N E M O ---\n"
    <<"    LOGFILE\n\n\n";
  FH<<"| basename                                | simfolder                               |      start time     |      stop time      | e-time CPU |"
	<<" repl done | rpl e-time | mean gen | version                 | hostname             | output files \n";

  FH.close();
}
// ----------------------------------------------------------------------------------------
// SimRunner::printLog()
// ----------------------------------------------------------------------------------------
void SimRunner::printLog()
{
  ofstream FH(_logfile.c_str(),ios::app);

  if(!FH.is_open()){
    error("could not open simulation logfile \"%s\"\n",_logfile.c_str());
    return;
  }

  string simfolder = _FileServices.getBaseFileName();

  FH<<"| ";
  FH.width(40);
  FH.setf(ios::left,ios::adjustfield);
  if(simfolder.empty()) FH << "-";
  else                  FH << simfolder;

  FH<<"| ";
  FH.width(40);
  FH.setf(ios::left,ios::adjustfield);
  FH << _FileServices.getSimfolder( );

  FH <<"| "<< _startTime <<" | "<< _endTime <<" | ";

  FH.width(10);
  FH.setf(ios::right,ios::adjustfield);
  FH<<_simElapsedTime<<" | ";

  FH.width(9);
  FH << _thePop->getReplicates( ) <<" | ";

  FH.width(10);
  FH << getElapsedTime( _thePop->getMeanReplElapsedTime() ) << " | ";

  FH.width(8);
  FH << _thePop->getMeanGenLength( ) <<" | ";

//  FH<<" "<<MAIN_VERSION<<"."<<MINOR_VERSION<<"."<<REVISION<<RELEASE
//    <<" "<<VERSION_DATE;
  FH<<"["<<VERSION_DATE<<"; "<<VERSION_TIME<<"]";

  FH<<" | ";
  FH.width(20);
  FH.setf(ios::left,ios::adjustfield);
  char* host;
  if( (host = getenv("HOST")) != NULL )
    FH << host << " |";
  else if ( (host = getenv("HOSTNAME")) != NULL )
    FH << host << " |";
  else
    FH << "-" << " |";

  FileServices::file_it file = _FileServices.getFirst(), last = _FileServices.getLast() ;

  for(;file != last; ++file)
    FH << " \"" << (*file)->get_extension() << "\":" << (*file)->get_short_path();

  FH << "\n";

  FH.close();
}

// ----------------------------------------------------------------------------------------
/** launch the script at the end of all simulations */
void SimRunner::runPostExec ( )
{
  ifstream script(_postexec_script.c_str(),ios::in);
  string cmd;

  if(!script.is_open()) {
	  error("could not open post simulation shell script!\n");
	  return;
  }

  message("Executing shell script \"%s\" ",_postexec_script.c_str());
  fflush(stdout);

  cmd = "sh " + _postexec_script;

  if(system(cmd.c_str()) < 0){
	  error("execution of `sh %s' failed: %s\n",_postexec_script.c_str(),strerror(errno));
	  return;
  }

  message("...done\n");
}

// ----------------------------------------------------------------------------------------
/** get the maximal number of traits */
int SimRunner::get_nbTraits(string name, map<string, vector<string> >& inputParams){

  // if the number of loci is not specified no trait will be inizialized
  map<string, vector<string> >::iterator pos = inputParams.find(name+"_loci"); // general param
  if(pos == inputParams.end()) pos = inputParams.find(name+"_loci_1");         // or at least for the frist trait
  if(pos == inputParams.end()) return 0;     // the number of loci is not specified
	if(!STRING::str2int<unsigned int>(pos->second[0])) return 0;

  // get the number of traits
  pos = inputParams.find(name+"_nb_trait");
  if(pos == inputParams.end()) return 1;     // the number of traits is not specified
  int nb, max = 0;
  for(unsigned int i=0; i<pos->second.size(); ++i){
    nb = STRING::str2int<int>(pos->second[i]);
    if(nb < 0) fatal("The parameter '%s_nb_trait' must be positive or zero!\n", name.c_str());
    if(max < nb) max = nb;
  }
  return max;                               // return the max number of traits
}

// ----------------------------------------------------------------------------------------
/** get the number of traits of the current simulation */
int SimRunner::get_nbTraits(string name, map<string, string>& inputParams){

  // if the number of loci is not specified no trait will be inizialized
  map<string, string>::iterator pos = inputParams.find(name+"_loci");          // general param
  if(pos == inputParams.end()) pos = inputParams.find(name+"_loci_1");         // or at least for the frist trait
  if(pos == inputParams.end()) return 0;     // the number of loci is not specified

  // get the number of traits
  pos = inputParams.find(name+"_nb_trait");
  if(pos == inputParams.end()) return 1;     // the number of traits is not specified
  int nb = STRING::str2int<int>(pos->second);
  if(nb < 0) fatal("The parameter '%s_nb_trait' must be positive or zero!\n", name.c_str());
  return nb;                                 // return the specified number of traits
}

// ----------------------------------------------------------------------------------------
// loadDefaultTemplates
// ----------------------------------------------------------------------------------------
/** load all traits and LCEs for the maximal simulation */
void SimRunner::loadDefaultTemplates(map<string, vector<string> >& inputParams)
{
  #ifdef _DEBUG
   message("\nSimManager::loadDefaultTemplates:\n");
  #endif

  // reset the previous templates
  reset_sim_components();
  get_allParams().clear();
  reset_services();

  //attach_pop(new Metapop());
  _thePop = new Metapop();
  _components.push_back(_thePop);


  //add TTraitProto
  int nb = get_nbTraits("quanti", inputParams);
  switch(nb){
    case 0: break;
    case 1: add_trait_template(new TTQuantiProto());   break;
    default:
      for(int i=0; i<nb; ++i){
        add_trait_template(new TTQuantiProto(i+1));
      }
      break;
  }

  nb = get_nbTraits("ntrl", inputParams);
  switch(nb){
    case 0: break;
    case 1: add_trait_template(new TTNeutralProto());  break;
    default:
      for(int i=0; i<nb; ++i){
        add_trait_template(new TTNeutralProto(i+1));
      }
      break;
  }

  //add LifeCycleEvents
  add_LCE_template(new LCE_Breed(1));
  add_LCE_template(new LCE_StatServiceNotifier(2));
  add_LCE_template(new LCE_FileServicesNotifier(3));
  add_LCE_template(new LCE_Aging(4));
	add_LCE_template(new LCE_Regulation(OFFSx, 5));
	add_LCE_template(new LCE_Disperse(6));
	add_LCE_template(new LCE_Regulation(ADLTx, 7));
	add_LCE_template(new LCE_Patch_Extinction(8));

	build_allParams();
}

// ----------------------------------------------------------------------------------------
// loadDefaultTemplates
// ----------------------------------------------------------------------------------------
/** load all traits and LCEs of a certain simualtion*/
void SimRunner::loadDefaultTemplates(map<string, string>& inputParams)
{
  #ifdef _DEBUG
   message("\nSimManager::loadDefaultTemplates:\n");
  #endif

  // reset the previous templates
  reset_sim_components();

  //attach_pop(new Metapop());
  _thePop = new Metapop();
  _components.push_back(_thePop);
  reset_services();

  //add TTraitProto
  int nb = get_nbTraits("quanti", inputParams);
  switch(nb){
    case 0: break;
    case 1: add_trait_template(new TTQuantiProto());   break;
    default:
      for(int i=0; i<nb; ++i){
        add_trait_template(new TTQuantiProto(i+1));
      }
      break;
  }

  nb = get_nbTraits("ntrl", inputParams);
  switch(nb){
    case 0: break;
    case 1: add_trait_template(new TTNeutralProto());  break;
    default:
      for(int i=0; i<nb; ++i){
        add_trait_template(new TTNeutralProto(i+1));
      }
      break;
  }

  //add LifeCycleEvents
  add_LCE_template(new LCE_Breed(1));
  add_LCE_template(new LCE_StatServiceNotifier(2));
  add_LCE_template(new LCE_FileServicesNotifier(3));
  add_LCE_template(new LCE_Aging(4));
	add_LCE_template(new LCE_Regulation(OFFSx, 5));
	add_LCE_template(new LCE_Disperse(6));
	add_LCE_template(new LCE_Regulation(ADLTx, 7));
	add_LCE_template(new LCE_Patch_Extinction(8));

  build_allParams();

}
