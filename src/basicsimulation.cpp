/** @file basicsimulation.cpp
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


#include "basicsimulation.h"
#include "output.h"

//----------------------------------------------------------------------------------------
// build_component_list
// ----------------------------------------------------------------------------------------
void ComponentManager::build_component_list ()
{
  _components.clear();

  list< TTraitProto* >::iterator TT_it = _TTrait_Templates.begin();

  while(TT_it != _TTrait_Templates.end()) {
    _components.push_back( (*TT_it) );
    TT_it++;
  }

  list< LCE* >::iterator LCE_it = _LCE_Templates.begin();

  while(LCE_it != _LCE_Templates.end()) {
    _components.push_back( (*LCE_it) );
    LCE_it++;
  }
}
//----------------------------------------------------------------------------------------
// get_trait_template
// ----------------------------------------------------------------------------------------
TTraitProto* ComponentManager::get_trait_template (string& name)
{
  list< TTraitProto* >::iterator TT = _TTrait_Templates.begin();

  while(TT != _TTrait_Templates.end()) {
	  if( (*TT)->get_paramset()->getName() == name) return (*TT);
	  TT++;
  }

  return NULL;
}
//----------------------------------------------------------------------------------------
// get_LCE_template
// ----------------------------------------------------------------------------------------
LCE* ComponentManager::get_LCE_template (string& name)
{
  list< LCE* >::iterator LCE = _LCE_Templates.begin();

  while(LCE != _LCE_Templates.end()) {
    if( (*LCE)->get_event_name() == name) return (*LCE);
    LCE++;
  }

  return NULL;
}
/*/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******ParamManager*******

//----------------------------------------------------------------------------------------
// cstor
// ---------------------------------------------------------------------------------------
ParamManager::ParamManager()
{
  #ifdef _DEBUG
   message(" ParamManager::ParamManager()\n");
  #endif
  _paramSet.set_name("simulation");
  _paramSet.set_isRequired(true);
  _paramSet.add_param("filename",STR,false,false,0,0,"simulation");
  _paramSet.add_param("folder",STR,false,false,0,0,"simulation");
  _paramSet.add_param("logfile",STR,false,false,0,0,"quantinemo.log");
  _paramSet.add_param("overwrite",INT2,false,true,0,1,"0");
  _paramSet.add_param("postexec_script",STR,false,false,0,0,"");
  _paramSet.add_param("seed",INT_MAT,false,false,0,0,"-1");
  _paramSet.add_param("logfile_type",INT2,false,true,0,2,"0");
  _paramSet.add_param("genetic_map_output",INT2,false,true,0,1,"0");

  // simulation
  _paramSet.add_param("replicates",INT2,false,false,0,0,"1");
  _paramSet.add_param("generations",INT2,true,false,0,0,"0");
}

// ---------------------------------------------------------------------------------------
/** destructor */
ParamManager::~ParamManager ( ) {
  #ifdef _DEBUG
   message(" ParamManager::~ParamManager (%i ParamSets)\n", _allParams.size());
  #endif
}

//----------------------------------------------------------------------------------------
// build_allParams
// ----------------------------------------------------------------------------------------
void ParamManager::build_allParams ()
{
  _allParams.clear();

  _allParams.push_back(&_paramSet);

  list< SimComponent* >::iterator cmpt = this->_components.begin();
  while(cmpt != this->_components.end()) {
    _allParams.push_back( (*cmpt)->get_paramset() );
    cmpt++;
  }
}

//----------------------------------------------------------------------------------------
// reset_temporalParams
// ----------------------------------------------------------------------------------------
void ParamManager::reset_params ()
{
  list<ParamSet*>::iterator pos = _allParams.begin();
  for(; pos != _allParams.end(); ++pos){
    (*pos)->reset();
  }
}

//----------------------------------------------------------------------------------------
// param_consistency_check
// ----------------------------------------------------------------------------------------
bool ParamManager::param_consistency_check ()
{
  list<ParamSet*>::iterator current_paramset = _allParams.begin();

  bool check = true;

  while(current_paramset != _allParams.end()){
    if(!(*current_paramset)->check_consistency()){
      error("ParamManager::param_consistency_check::consistency not satisfied for \"%s\"\n",
            (*current_paramset)->getName().c_str());
      check = false;
    }
    current_paramset++;
  }
  return check;
}
//----------------------------------------------------------------------------------------
// set_parameters
// ----------------------------------------------------------------------------------------
/** check if the input parameters are really parameters
  * same function as below, but for another type of container
  */
void
ParamManager::check_parameter_names (map<string, vector<string> >& inputParams)
{
  list<ParamSet*>::iterator cur_paramset;
  map<string, vector<string> >::iterator input;
  map<string, Param*>::iterator cur_param;
  bool set;
  vector< string > notSetParams;
  vector< string >::iterator iterVec;
  string generalName, token;
  int cur_number, number;
  string::size_type pos;

  reset_params();     // all parameters are resetted ("set" is set to false)

  // 1: set all "precice" parameters
  for(input = inputParams.begin(); input != inputParams.end(); ++input) {
    set = false;

    // find the corresponding parameter set
		for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
			//cout << "\n" << (string&)input->first;
			set |= (*cur_paramset)->set_param((string&)input->first, input->second[0]);
      //is true if at least one paramset could set the param
    }

    // add it to the vector if the parameter was never used
    if(!set) notSetParams.push_back(input->first);
  }

  // 2: for each multiple trait find the corresponding parameter:
  //   multiple traits possibilities (example parameter: quanti_nb_all_5):
  //   - for each trait a specific paramter is passed in the settings file (eg "quanti_nb_all_5")
  //   - if such a parameter is not set it is verified if a parameter with a lower number is set (eg "quanti_nb_all_3")
  //   - if such a parameter is not set it is verified if a general paramter is set (eg "quanti_nb_all"),
  //   - general parameters without numbers and paramters with a 0 are identical (eg quanti_nb_all = quanti_nb_all_0)

  for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
    // check if this parameter set is a multiple trait
    token = (*cur_paramset)->getName();
    pos = token.rfind("_");
    if(pos == string::npos) continue;                              // "_" did not exist
		number = STRING::str2int<int>(token.substr(pos+1).c_str());    // check if the last substring is a number
    if(number == my_NAN) continue;                                 // if the substring is not a number

    // it is a multiple trait:
    // scroll through all the trait parameters and search the correct input parameter
    cur_param = (*cur_paramset)->getAllParams().begin();
    for(; cur_param != (*cur_paramset)->getAllParams().end(); ++cur_param){
      if(cur_param->second->isSet()) continue;  // if the parameter is already set continue

      // check if the parameter is grouped (i.e. the parameter has a lower number)
      // or finally no number (general parameter)
      generalName = cur_param->first.substr(0, cur_param->first.rfind("_"));
      cur_number = number;
      while(cur_number >= -1){
        if(cur_number == -1)  token = generalName;                        // general parameter
        else                  token = generalName + STRING::int2str_(cur_number); // number between 0 and n
        input = inputParams.find(token);
        if(input != inputParams.end()) break;
        --cur_number;
      }

      // if no input was found continue
      if(input == inputParams.end()) continue;

      // set the general parameter
      cur_param->second->set(input->second[0], *cur_paramset);

      // remove this parameter from the not set vector
      for(iterVec = notSetParams.begin(); iterVec != notSetParams.end(); ++iterVec){
        if(*iterVec== generalName){
          notSetParams.erase(iterVec);
          break;
        }
      }
    }
  }

  // output all not used input parameters and remove them from the map
  for(iterVec = notSetParams.begin(); iterVec != notSetParams.end(); ++iterVec){
      warning("Settings file: parameter '%s' could not be set -> parameter is skipped!",(*iterVec).c_str());
      assert(inputParams.find(*iterVec) != inputParams.end());
      inputParams.erase(inputParams.find(*iterVec));
  }

  param_consistency_check();
}

//----------------------------------------------------------------------------------------
// set_parameters
// ----------------------------------------------------------------------------------------
/**  attribute the input to the parameters */
bool ParamManager::set_parameters (map< string,string >& simparams)
{
  _inputParams = simparams;
  list<ParamSet*>::iterator cur_paramset;
  map< string,string >::iterator input;
  map<string, Param*>::iterator cur_param;
  bool set;
  vector< string > notSetParams;
  vector< string >::iterator iterVec;
  string generalName, token, paramName;
  int cur_number, number;
  string::size_type pos;

  // 1: set all "precice" parameters
  for(input = _inputParams.begin(); input != _inputParams.end(); ++input) {
    set = false;

    // find the corresponding parameter set
    for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
      set |= (*cur_paramset)->set_param((string&)input->first,input->second);
      //is true if at least one paramset could set the param
    }

    // add it to the vector if the parameter was never used
    if(!set) notSetParams.push_back(input->first);
  }

  // 2: for each multiple trait find the corresponding parameter:
  //   multiple traits possibilities (example parameter: quanti_nb_all_5):
  //   - for each trait a specific paramter is passed in the settings file (eg "quanti_nb_all_5")
  //   - if such a parameter is not set it is verified if a parameter with a lower number is set (eg "quanti_nb_all_3")
  //   - if such a parameter is not set it is verified if a general paramter is set (eg "quanti_nb_all"),
  //   - general parameters without numbers and paramters with a 0 are identical (eg quanti_nb_all = quanti_nb_all_0)

  for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
    // check if this parameter set is a multiple trait
    paramName = (*cur_paramset)->getName();
    pos = paramName.rfind("_");
    if(pos == string::npos) continue;                                  // "_" did not exist
		number = STRING::str2int<int>(paramName.substr(pos+1).c_str());    // check if the last substring is a number
    if(number == my_NAN) continue;                                     // if the substring is not a number

    // it is a multiple trait:
    // scroll through all the trait parameters and search the correct input parameter
    cur_param = (*cur_paramset)->getAllParams().begin();
    for(; cur_param != (*cur_paramset)->getAllParams().end(); ++cur_param){
      if(cur_param->second->isSet()) continue;  // if the parameter is already set continue

      // check if the parameter is grouped (i.e. the parameter has a lower number)
      // or finally no number (general parameter)
      generalName = cur_param->first.substr(0, cur_param->first.rfind("_"));
      cur_number = number;
      while(cur_number >= -1){
        if(cur_number == -1)  token = generalName;                        // general parameter
        else                  token = generalName + STRING::int2str_(cur_number); // number between 0 and n
        input = _inputParams.find(token);
        if(input != _inputParams.end()) break;
        --cur_number;
      }

      // if no input was found continue
      if(input == _inputParams.end()) continue;

      // set the general parameter
      cur_param->second->set(input->second, *cur_paramset);

      // remove this parameter from the not set vector
      for(iterVec = notSetParams.begin(); iterVec != notSetParams.end(); ++iterVec){
        if(*iterVec== generalName){
          notSetParams.erase(iterVec);
          break;
        }
      }
    }
  }

  // output all not used input parameters
  for(iterVec = notSetParams.begin(); iterVec != notSetParams.end(); ++iterVec){
      warning("ParamManager::could not set param '%s'\n",(*iterVec).c_str());
  }

  // check if parameters change over time
  _temporalParameters = false;
  for(cur_paramset = _allParams.begin(); cur_paramset != _allParams.end(); ++cur_paramset){
    _temporalParameters |= (*cur_paramset)->isTemporalParamSet();
    //is true if at least one paramset could set the param
  }

  return param_consistency_check();
}

//----------------------------------------------------------------------------------------
// get_paramset
// ----------------------------------------------------------------------------------------
ParamSet* ParamManager::get_paramset (string& name)
{
  list<ParamSet*>::iterator pset = _allParams.begin();

  while(pset != _allParams.end()) {
    if( (*pset)->getName().compare(name) == 0) return (*pset);
    pset++;
  }

  return NULL;
}

//----------------------------------------------------------------------------------------
// build_records
// ----------------------------------------------------------------------------------------
/** generates the individual simulation parameter sets if sequential parameters are used */
void ParamManager::build_records (map< string, vector<string> >& initParams)
{
  map<string,string> params;
	unsigned int RecNb = 1, ArgNb, SeqParam, BlockSize;
  vector<unsigned int> sequence; //stores the number of args of each sequence parameters
  vector<string> currSeqArg;
  string NAME = "simulation", arg;

	//then process the InitParams map to find sequence parameters:
  map< string,vector<string> >::iterator Pit;
  for(Pit = initParams.begin(); Pit != initParams.end(); ++Pit) {
	  if(Pit->second.size() > 1) {
      sequence.push_back(Pit->second.size());           //fetch the number of the current sequence parameter:
      RecNb *= Pit->second.size();                      //increase the total number of simulations records:
    }
  }

  if(RecNb > 100) warning("Your settings file has %u sequential simulations specified! Is that your aim, or is there a problem in your settings file?\n", RecNb);

  // for each sequence
  for(unsigned int i=0; i<RecNb; ++i) {
    //now build the simulation records with the right params!
    //the map 'param' will get all the params used for one simulation
    //it is then added to the list of simulations' parameters map

		SeqParam = 0;//used as index of the 'sequence' vector

    // for each input parameter
    for(Pit = initParams.begin(); Pit != initParams.end(); ++Pit) {

			// if it is the filename
			if (Pit->first == "filename") {
        if(RecNb>1) NAME = Pit->second[0];       // sequence simulations
				else  params["filename"] = Pit->second[0];
        continue;
			}

      // for all other parameters
      //get the number of arguments for the current parameter:
      switch(Pit->second.size()){
        case 0:
          //current param has no argument (bool type param) we give it value 1 (true)
          params[Pit->first] = "1";
          break;

				case 1:
					//the current param has one argument
          params[Pit->first] = Pit->second[0];
          break;

        default:
          //the current param is a sequence param (ArgNb > 1)
          //increase the index of the sequence parameter
          SeqParam++;
          //then compute the right argument to give to the current simulation record:
          BlockSize = RecNb;
          ArgNb = Pit->second.size();

          //message("\n%s (%i): %s",Pit->first.c_str(), Pit->second.size(), Pit->second[ (i/BlockSize) % ArgNb ].c_str());

					for(unsigned int j=0;j < SeqParam;++j){
            BlockSize /= sequence[j];
          }
          arg = Pit->second[ (i/BlockSize) % ArgNb ];
          if(arg[0]=='{' || arg[0]=='('){         // if matrix or temporal -> argument is too big, replaced by the number
            currSeqArg.push_back(STRING::int2str(1+(i/BlockSize) % ArgNb, ArgNb));
          }
          else currSeqArg.push_back(arg);        // add argument
          params[Pit->first] = Pit->second[ (i/BlockSize) % ArgNb ];
          break;
      }
    }

		if(RecNb>1) params["filename"] = setFilename(NAME,i+1,RecNb,currSeqArg);   // sequence simulation

    //add all the params previously computed to the main simulation recorder list (of maps)
    _simRecords.push_back(params);

    currSeqArg.clear();
  } // end of each sequence of simulation
}

// ----------------------------------------------------------------------------------------
// setFilename
// ----------------------------------------------------------------------------------------
string ParamManager::setFilename(string& fstring,const unsigned int& sim,const unsigned& nbSims, vector<string>& args)
{
  string out, tail;
  int cur, digits;
	string::size_type next;
  unsigned int index, i, size = args.size();

  // control array that all parameters were used, i.e. that the file name is unique!
  int* usedParams = new int[size];
  for(i=0; i<size; ++i){
    usedParams[i] = 0;    // set it to false
  }

  // parse the name
  cur  = 0;
  next = fstring.find_first_of('%');
	while(next != string::npos){
    out += fstring.substr(cur, next-cur);
    tail = fstring.substr(next+1, string::npos);
    index = (unsigned int) strtol(tail.c_str(),NULL,10);          // get the param number
    if(index > args.size()) fatal("too much arguments in filename!\n");
    digits = STRING::getNbDigits(index);

    out += args[index-1];
    usedParams[index-1] = 1;    // this param was used
    cur = next + digits + 1;
    next = fstring.find('%',cur);  // get the next position
  }
  out += fstring.substr(cur, next-cur);

  // control if all params were used
  for(i=0; i<size; ++i){
    if(!usedParams[i]) break;
  }
  delete[] usedParams;

  if(i != size){    // add the simulation number if the name is not unique
    out += "-" + STRING::int2str(sim, nbSims);
  }

  return out;
}

// ----------------------------------------------------------------------------------------
// lowercase
// ----------------------------------------------------------------------------------------
string ParamManager::lowercase(string& input)
{
  for(unsigned int i=0;i<input.size();++i){
    input[i] = tolower(input[i]);
  }
  return input;
}

/*/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******SimBuilder*******

//----------------------------------------------------------------------------------------
// copy cstor
// ---------------------------------------------------------------------------------------
SimBuilder::SimBuilder (const SimBuilder& SB)
{
  fatal("SimBuilder::SimBuilder(const SimBuilder&) not used and thus not implemented yet!\n");
}
//----------------------------------------------------------------------------------------
// get_current_trait
// ----------------------------------------------------------------------------------------
TTraitProto* SimBuilder::get_current_trait (string type)
{
  map<string, TTraitProto* >::iterator trait;

  trait = _currentTraits.find(type);

  if(trait != _currentTraits.end())
    return trait->second;
  else
    return NULL;
}
//----------------------------------------------------------------------------------------
// get_current_event
// ----------------------------------------------------------------------------------------
LCE* SimBuilder::get_current_event (string& name)
{
  map< int, LCE* >::iterator LCE = _currentLifeCycle.begin();

  while(LCE != _currentLifeCycle.end()) {
    if( LCE->second->get_event_name() == name)
      return LCE->second;
    LCE++;
  }

  return NULL;
}
//----------------------------------------------------------------------------------------
// build_currentParams
// ----------------------------------------------------------------------------------------
bool SimBuilder::build_currentParams (map< string,string >& simparams)
{
  if(! this->set_parameters(simparams) ) return false;

  //empty the current parameters container:
  _currentParams.clear();

  list<ParamSet*>::iterator current_paramset = this->_allParams.begin();

  //add all set parameters to the list of the current params:
  while(current_paramset != this->_allParams.end()){
    if((*current_paramset)->isSet()) {
      #ifdef _DEBUG
       message("SimBuilder::build_currentParams:adding paramset %s\n",
              (*current_paramset)->getName().c_str());
      #endif
      _currentParams.push_back(*current_paramset);
    }
    current_paramset++;
  }

  return true;
}
//----------------------------------------------------------------------------------------
// build_currentTraits
// ----------------------------------------------------------------------------------------
map< string,TTraitProto* >& SimBuilder::build_currentTraits()
{
  //build the list of Traits for this simulation:
  list< TTraitProto* >::iterator trait = this->_TTrait_Templates.begin();
  _currentTraits.clear();

   while(trait != this->_TTrait_Templates.end()) {
    if( (*trait)->get_paramset()->isSet() ){
      _currentTraits[(*trait)->get_type()] = (*trait);
    }
    trait++;
	}

  return _currentTraits;
}
//----------------------------------------------------------------------------------------
// build_currentLifeCycle
// ----------------------------------------------------------------------------------------
map< int, LCE* >& SimBuilder::build_currentLifeCycle()
{
  //build the list of the life cycle events for this simulation:
  list< LCE* >::iterator LCEs = this->_LCE_Templates.begin();

  _currentLifeCycle.clear();

	for(;LCEs != this->_LCE_Templates.end();++LCEs){
    _currentLifeCycle[(int)(*LCEs)->get_rank()] = (*LCEs);
	}
	return _currentLifeCycle;
}


//----------------------------------------------------------------------------------------
// checkLCEconsistency
// ----------------------------------------------------------------------------------------
/** check if the sequence of LCE makes sense */
bool
SimBuilder::checkLCEconsistency(){

  // check if a LCE is really used
  map< int, LCE* >::iterator iterLCE = _currentLifeCycle.begin();
  for(; iterLCE != _currentLifeCycle.end();){
    string name = iterLCE->second->get_event_name();
    if(!iterLCE->second->get_paramset()->isSet()){
	    _currentLifeCycle.erase(iterLCE++);
	  }
	  else  ++iterLCE;
  }

  //find the age class required to start the life cycle with:
  age_t availableAge = NONE;
  iterLCE = _currentLifeCycle.begin();
  for(; iterLCE != _currentLifeCycle.end() && availableAge == NONE; ++iterLCE){
    availableAge = iterLCE->second->requiredAgeClass();
  }

  // check if the LCE makes sense
  LCE *curLCE=NULL, *lastLCE=NULL;
  bool stop = false;
  iterLCE = _currentLifeCycle.begin(); // start with the first LCE
  while(!stop){
    // check if the last and the first LCE correspond
    if(iterLCE == _currentLifeCycle.end()){
      stop = true;
      iterLCE = _currentLifeCycle.begin();
    }

    // get the LCE
    curLCE = iterLCE->second;

    // if a age class is required
    if(curLCE->requiredAgeClass()){
      // check if the required age class is available or if no age classes are needed
      if(!(availableAge & curLCE->requiredAgeClass())){
        string name = curLCE->get_event_name();
        fatal("LCE order does not make sense (%s -> %s)\n",
            lastLCE->get_event_name().c_str(), curLCE->get_event_name().c_str());
      }
      lastLCE = curLCE;
    }

    // update the current age class
    availableAge &= ~curLCE->removeAgeClass();
    availableAge |= curLCE->addAgeClass();
    ++iterLCE;
  }
  return true;
}
