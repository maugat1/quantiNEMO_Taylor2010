/** @file LCEmisc.cpp
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
#include "lce_misc.h"
#include "metapop.h"
#include "stat_rec_base.h"
#include "stathandler.cpp"
#include "simulation.h"
#include "output.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** LCE_FileServicesNotifier ********

// ----------------------------------------------------------------------------------------
// LCE_FileServicesNotifier::execute
// ----------------------------------------------------------------------------------------
void LCE_FileServicesNotifier::execute ()
{
  #ifdef _DEBUG
   message("  LCE_FileServicesNotifier ... ");
  #endif

  _service->notify();

  #ifdef _DEBUG
   message("done! (gen: %i rpl: %i)\n",
		  this->get_pop_ptr()->getCurrentGeneration(),
		  this->get_pop_ptr()->getCurrentReplicate());
  #endif
}


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** LCE_StatServiceNotifier ********

//-----------------------------------------------------------------------------
// LCE_StatServiceNotifier
//-----------------------------------------------------------------------------
LCE_StatServiceNotifier::LCE_StatServiceNotifier (int rank)
: LCE("save_stats","",rank), _service(0), _tot_occurrence(0), _index(0){

  this->add_parameter("stat",STR_MAT,true,false,0,0,"");
  this->add_parameter("stat_save",INT2,false,true,0,6,"0");
  this->add_parameter("stat_log_time",INT2,false,false,0,0,"1",true);
  this->add_parameter("stat_dir",STR,false,false,0,0,"");
}
//-----------------------------------------------------------------------------
// loadStatServices
//-----------------------------------------------------------------------------
void LCE_StatServiceNotifier::loadStatServices ( StatServices* loader )
{
  _service = loader;
	_fileHandler.set_statService(loader);
	_fileHandler.set_save_choice(_save_choice);
  _service->set(_arg,_occurrence, _tot_occurrence, _index);
  loader->attach(&_statHandler);
}
//-----------------------------------------------------------------------------
// LCE_StatServiceNotifier::execute
//-----------------------------------------------------------------------------
void LCE_StatServiceNotifier::execute ()
{
  #ifdef _DEBUG
   message("  LCE_StatServiceNotifier ... ");
  #endif

  if(_service) _service->notify();

  #ifdef _DEBUG
   if(_service) message("done! (occurrence %i)\n",_service->getOccurrence());
  #endif
}

//-----------------------------------------------------------------------------
// init
//-----------------------------------------------------------------------------
bool LCE_StatServiceNotifier::init (Metapop* popPtr)
{
  LCE::init(popPtr);

  _arg = _paramSet->getArg("stat");

  _save_choice = (unsigned int)_paramSet->getValue("stat_save");
  if(_save_choice == 5){      // no stats are written, LCE is not used
    _paramSet->set_isSet(false);
    return false;
  }

  unsigned int occ = _occurrence = (unsigned int)_paramSet->getValue("stat_log_time");
  if(occ < 1) fatal("Parameter 'stat_log_time' must be 1 or larger!\n");

  // set the index array
  vector<unsigned int> vIndex;   // get the indexes
  int gen = _popPtr->getGenerations();

  // if the log time changes over time
  if(_paramSet->get_param("stat_log_time")->isTemporalParam()){
    map<int, string>* pArgs = _paramSet->get_param("stat_log_time")->get_temporal_args();
    map<int, string>::iterator pos = pArgs->begin();  // set the iterator

    // set the starting conditions
    ++pos;        // go to the next one

    // loop through the generations and sum up the occurences
    for(int i=1; i<=gen; ++i){
      if(i == pos->first){                      // new temporal parameter
        occ = STRING::str2int<unsigned int>(pos->second);     // change occurence
        ++pos;                                  // go the the next position
      }
      if(!((i) % occ)) vIndex.push_back(i); //++ _tot_occurrence;
    }
  }
  else {     // no temporal change
    for(int i=1; i<=gen; ++i){
      if(!((i) % occ)) vIndex.push_back(i); //++ _tot_occurrence;
    }
  }
  // copy vector to array
  _tot_occurrence = vIndex.size();
  if(!_tot_occurrence){
    warning("The log time of the summary statistics (%i) is larger than the number of generations (%i): no summary statistics will be computed!\n",
      _occurrence, gen);
  }

  if(_index) delete[] _index;
  _index = new unsigned int[_tot_occurrence];
  for(unsigned int i=0; i<_tot_occurrence; ++i){
    _index[i] = vIndex[i];
  }

  _fileHandler.set(false, false,
                   popPtr->getReplicates(),
                   popPtr->getGenerations(),
                   this->get_rank(),
                   _paramSet->getArg("stat_dir"),
                   "", 0, 0, NULL, 1);
  return true;
}

//-----------------------------------------------------------------------------
// executeBeforeEachGeneration
//-----------------------------------------------------------------------------
void
LCE_StatServiceNotifier::executeBeforeEachGeneration(const int& gen)
{
  // temporal parameters
  map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);

  // if it is a temporal paramter
  if(pParam){
    // check if a change has to be made
    map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
    if(pMap){
      // iterate through the map and performe the updates
      map<string, Param*>::iterator pos = pMap->begin();
      for(; pos != pMap->end(); ++pos){
        if(pos->first == "stat_log_time" && _service){
          _occurrence = (unsigned int)pos->second->get_value();
          _service->setOccurrence(_occurrence);
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// get_isAlive
//-----------------------------------------------------------------------------
double LCE_StatSH::get_isAlive () {
  return (double)_pop->isAlive();
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                           ******** LCE_Aging ********

// ----------------------------------------------------------------------------------------
// LCE_Aging::execute
// ----------------------------------------------------------------------------------------
void LCE_Aging::execute ()
{
  #ifdef _DEBUG
   message("  LCE_Aging ... ");
  #endif

  for(unsigned int i = 0; i < _popPtr->getPatchNbr(); i++){
		_popPtr->getPatch(i)->flush(ADLTx);

	  //set the Patch extinction tag:
	  _popPtr->getPatch(i)->set_isExtinct(_popPtr->getPatch(i)->size(OFFSx)==0);
  }

  #ifdef _DEBUG
	 message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
    _popPtr->getPatchNbr(),
		_popPtr->size(FEM, OFFSPRG),  _popPtr->size(MAL, OFFSPRG),
		_popPtr->size(FEM, ADULTS),   _popPtr->size(MAL, ADULTS));
  #endif
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** LCE_Patch_Extinction ********

//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::execute
//-----------------------------------------------------------------------------
void LCE_Patch_Extinction::execute ()
{
  #ifdef _DEBUG
  message("  LCE_Patch_Extinction ... ");
  unsigned int cnt = 0;
  #endif

  int nbExt, nbPatch, i, j, rand;
  Patch *current_patch;

  // get the number of extinctions
  nbPatch = _popPtr->getPatchNbr();
	nbExt = (int)SimRunner::r.Binomial(_Xtion_rate, nbPatch);
  if(nbExt){
    int* extinctions = new int[nbExt];   // to store the patches of extinction in order to prevent to extinct twice a patch
    for(i = 0; i< nbExt;){
			rand = SimRunner::r.Uniform(nbPatch);    	// draw randomly a patch for extinction

      // check if this patch was already drawn
      for(j = 0; j<i; ++j){
        if(extinctions[j] == rand) break;
      }

      // if this patch is not yet drawn
      if(j>=i){
				extinctions[i] = rand;                 	// add the patch to the array
				_popPtr->getPatch(rand)->flush();  			// perform the extinction
        ++i;
      }
    }
    delete[] extinctions;
  }

  // loop through all the patches to check if they are extinct or not
  for(i = 0; i < nbPatch; ++i) {
    current_patch = _popPtr->getPatch(i);
    current_patch->set_isExtinct(current_patch->isEmpty());

    #ifdef _DEBUG
    cnt += (current_patch->get_isExtinct());
    #endif
  }
  #ifdef _DEBUG
  message("done! (%i extinct patches)\n",cnt);
  #endif
}

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// LCE_Patch_Extinction::init
//-----------------------------------------------------------------------------
bool LCE_Patch_Extinction::init (Metapop* popPtr)
{
  LCE::init(popPtr);
  _Xtion_rate = this->get_parameter_value("extinction_rate");
  if(!_Xtion_rate) _paramSet->set_isSet(false);
  return (_Xtion_rate!=0);  // false if no extinction occurs
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Regulation ********

// ----------------------------------------------------------------------------------------
// LCE_Regulation::init
// ----------------------------------------------------------------------------------------
bool LCE_Regulation::init(Metapop* popPtr){
  LCE::init(popPtr);

    // does selection acts at this stage?
  if( ( _age == OFFSx && _popPtr->get_selection_position() == 2) || 
  (_age == ADLTx && _popPtr->get_selection_position() == 3) )
    {

		switch(_popPtr->get_selection_level()){
			case 0: regulation = &LCE_Regulation::regulation_fitness_patch;    break;   // patch level
			case 1: regulation = &LCE_Regulation::regulation_fitness_metapop;           // metapop level
							_popPtr->set_total_carrying_capacity();
							break;
			case 2: regulation = &LCE_Regulation::regulation_fitness_hard;     break;   // hard selection
		}
  }
	else if((int)this->get_parameter_value("regulation_model_"+_ageStr)){           // neutral regulation
	  regulation = &LCE_Regulation::regulation_neutral;
	}
	else{ 					// no regulation
		_paramSet->set_isSet(false);
		return false;
	}
	return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::execute
// ----------------------------------------------------------------------------------------
void LCE_Regulation::execute (){
	#ifdef _DEBUG
	string reg;
	if(     regulation == &LCE_Regulation::regulation_fitness_patch)   reg = "patch selection";
	else if(regulation == &LCE_Regulation::regulation_fitness_metapop) reg = "metapop selection";
	else if(regulation == &LCE_Regulation::regulation_fitness_hard)    reg = "hard selection";
	else                                                               reg = "neutral";
	 message("  LCE_Regulation (%s, %s)... ", _ageStr.c_str(), reg.c_str());
  #endif

  (this->*regulation)();

  #ifdef _DEBUG
   message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
    _popPtr->getPatchNbr(),
    _popPtr->size(FEM, OFFSPRG),  _popPtr->size(MAL, OFFSPRG),
		_popPtr->size(FEM, ADULTS),   _popPtr->size(MAL, ADULTS));
	#endif
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::execute
// ----------------------------------------------------------------------------------------
void LCE_Regulation::regulation_neutral ()
{
	unsigned int N, K;
	Patch* current_patch;

	#ifdef _DEBUG
	 _ex_cnt = _col_cnt = _ph_cnt = 0;
	#endif

	for(unsigned int i = 0, size = _popPtr->getPatchNbr(); i < size; ++i) {
		current_patch = _popPtr->getPatch(i);

		// females
		K = current_patch->get_KFem();
		N = current_patch->size(FEM, _age);
		if(N > K){
			if(N < K/2) drawSuccessfullIndividuals(current_patch, K, FEM);
			else        drawUnSuccessfullIndividuals(current_patch, K, FEM);    // if(N>K/2)
		}

		// males
		K = current_patch->get_KMal();
		N = current_patch->size(MAL, _age);
		if(N > K){
			if(N < K/2) drawSuccessfullIndividuals(current_patch, K, MAL);
			else        drawUnSuccessfullIndividuals(current_patch, K, MAL);    // if(N>K/2)
		}
		
		#ifdef _DEBUG
		 _ex_cnt  += current_patch->get_isExtinct();
		 _col_cnt += (current_patch->get_isExtinct() ? current_patch->nbKolonisers : 0);
		#endif
	}//end for patches
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::drawSuccessfulIndividuals
// ----------------------------------------------------------------------------------------
/** regulates randomly the size of the pop, by randomly drawing survivors
	* (this needs a new array of the succeeders which will then be exchanged)
	*/
void
LCE_Regulation::drawSuccessfullIndividuals(Patch* curPatch,
																				const unsigned int& K, const sex_t& SEX)
{
	unsigned int ind, cur_size=curPatch->size(SEX, _age);

	deque<Individual*> array;
	array.assign(K, NULL);

	for(unsigned int i=0; i<K; ++i, --cur_size){
		ind = SimRunner::r.Uniform(cur_size);

		assert(!(curPatch->get_isExtinct()
					 && curPatch->get(SEX, _age, ind)->getCurrentPatch()->get_ID() == curPatch->get_ID()));

		//this one won a breeding place, congratulations!
		array[i] = curPatch->get(SEX, _age, ind);
		curPatch->remove(SEX, _age, ind);
	}

	// reassign the container
	deque<Individual*> &container = curPatch->get_containers()[SEX][_age];
	delete &container; 		 // delete the remianing individuals
	container = array;     // assign the new array to the container
	curPatch->get_sizes()[SEX][_age]      = K;
	curPatch->get_capacities()[SEX][_age] = K;

	assert(curPatch->size(SEX, _age) == K);
}

// ----------------------------------------------------------------------------------------
// LCE_Regulation::drawUnSuccessfulIndividuals
// ----------------------------------------------------------------------------------------
/** regulates randomly the pop size, by removing randomly supernumerous individuals */
void
LCE_Regulation::drawUnSuccessfullIndividuals(Patch* curPatch,
																				const unsigned int& K, const sex_t& SEX)
{
	unsigned int ind, nbInd;

	// remove randomly supernumerous individuals
	for(nbInd = curPatch->size(SEX, _age); nbInd > K; --nbInd) {
		ind = SimRunner::r.Uniform(curPatch->size(SEX, _age) );
		curPatch->recycle(SEX, _age, ind);
	}
	assert(curPatch->size(SEX, _age) == K);
}
// ----------------------------------------------------------------------------------------


// ----------------------------------------------------------------------------------------
//                               ******** LCE_StatFH ********

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
/** 0: all; 1: each rep; 2: mean and var, 3: only mean; 4: only var; 5: nothing */
void LCE_StatFH::FHwrite()
{
	switch(_statService->get_save_choice()){
    case 0:                              // all
						printStat_each();
            printStat_mean();
            printStat_variance();
            printStat_legend();
            break;
    case 1:                              // each replicate
						printStat_each();
            printStat_legend();
            break;
    case 2:                              // mean and var
            printStat_mean();
            printStat_variance();
            printStat_legend();
            break;
    case 3:                              // only mean
            printStat_mean();
            printStat_legend();
            break;
    case 4:                              // only var
            printStat_variance();
            printStat_legend();
            break;
		case 5: break;                        // no ouput
		case 6: break;                        // output without db

	}
}

// ----------------------------------------------------------------------------------------
// printStat_legend
// ----------------------------------------------------------------------------------------
/** prints the legend of the statistics (only once at generation 0) */
void LCE_StatFH::printStat_legend()
{
  ofstream FH;

  string filename = get_path() + get_service()->getBaseFileName() + "_legend.txt";

  #ifdef _DEBUG
   message("  LCE_StatFH::PrintStat_legend (%s)\n",filename.c_str());
  #endif

  FH.open(filename.c_str(),ios::out);

  if(!FH)  fatal("Could not open stat legend file '%s'!\n", filename.c_str());

  // print the heading
  FH << "Legend of the used statistics\n"
     << "*****************************\n\n";
  FH.width(20);
  FH.setf(ios::left,ios::adjustfield);
  FH << "Statistic" << "  Legend\n";
  FH << "--------------------------------------------------------------\n";

  //print the stat names:
  _statService->printStatLegend(FH,(GEN | FLAT));

  FH.close();
}

// ----------------------------------------------------------------------------------------
// printStat_each
// ----------------------------------------------------------------------------------------
/** prints the stats for each replicate separately */
void
LCE_StatFH::printStat_each_header(ostream& FH)
{
	//print the first row with stat headers:
	FH.width(12);
	FH.setf(ios::left,ios::adjustfield);
  FH<<"replicate"<<"\t";
  FH.width(12);
  FH<<"generation";

	//print the stat names:
	_statService->printStatHeaders(FH, FLAT);

	FH<<"\n";   //next rows:
}

// ----------------------------------------------------------------------------------------
// printStat_each
// ----------------------------------------------------------------------------------------
/** prints the stats for each replicate separately */
void
LCE_StatFH::printStat_each()
{
	#ifdef _DEBUG
	 message("  LCE_StatFH::FHwrite (%s)\n",_statService->get_file_name_stats().c_str());
	#endif

	ofstream FH(_statService->get_file_name_stats().c_str(),ios::out);
	if(!FH) fatal("Could not open stat output file \"%s\"!\n",_statService->get_file_name_stats().c_str());

  //print the first row with stat headers:
  FH.width(12);
  FH.setf(ios::left,ios::adjustfield);
  FH<<"replicate"<<"\t";
  FH.width(12);
  FH<<"generation";

  //print the stat names:
	_statService->printStatHeaders(FH, FLAT);

  //next rows:
  FH<<"\n";

	unsigned int gen_index, i, j;
  unsigned int ncols = _statService->get_pop_ptr()->getReplicates();
  unsigned int nrows = _statService->getTotOccurrence();

  for(i = 0; i < ncols; i++) {
    for(j = 0; j < nrows; j++) {
      //val is ordered rows x cols, row = generation, col = replicate
      //we write: for 1st rpl all the stat for all the gen recorded, then we move to the next rpl.
      //print generation number:
      gen_index = _statService->getStatRecIndex(j);
      //if the generation at index j is null, all replicates crashed before that index

      if(gen_index != my_NAN) {
        FH.width(12);
        FH.setf(ios::left,ios::adjustfield);
        FH<<setprecision(4);
        FH<<(i+1);
        FH<<"\t";
        FH.width(12);
        FH<<gen_index;

				_statService->printStatValue(FH,j,i);

        FH<<"\n";
      }
      else break; //break the stat loop, switch to the next column/replicate
    }
  }
  FH.close();
}

// ----------------------------------------------------------------------------------------
// PrintStat_mean
// ----------------------------------------------------------------------------------------
/** merges the stats across replicates */
void LCE_StatFH::printStat_mean()
{
	#ifdef _DEBUG
	 message("  LCE_StatFH::PrintStat_mean (%s)\n",_statService->get_file_name_mean().c_str());
	#endif

	ofstream FH(_statService->get_file_name_mean().c_str(),ios::out);
	if(!FH)  fatal("Could not open stat mean file '%s'!\n", _statService->get_file_name_mean().c_str());

  //print the first row with stats headers:
  FH.width(12);
  FH.setf(ios::left,ios::adjustfield);
  FH<<"generation";

  //print the stat names:
  _statService->printStatHeaders(FH,(GEN | FLAT));

  FH<<"\n";

  //next rows:
  unsigned int gen_index;
  unsigned int nrows    = _statService->getTotOccurrence();

  for(unsigned int i = 0; i < nrows; i++) {
    //print generation number:
    gen_index = _statService->getStatRecIndex(i);

    //if the index is my_NAN, all replicates crashed before that index
    if(gen_index != my_NAN) {
      FH.width(12);
      FH.setf(ios::left,ios::adjustfield);
      FH<<setprecision(4);
      FH<<gen_index;

			_statService->printStatMean(FH,i);

      FH<<"\n";
    }
    else break;
  }
  FH.close();
}

// ----------------------------------------------------------------------------------------
// PrintStat_var
// ----------------------------------------------------------------------------------------
/** merges the stats across replicates */
void LCE_StatFH::printStat_variance()
{
	#ifdef _DEBUG
	 message("  LCE_StatFH::PrintStat_var (%s)\n",_statService->get_file_name_var().c_str());
  #endif

	ofstream FH(_statService->get_file_name_var().c_str(),ios::out);
	if(!FH)  fatal("Could not open stat var file '%s'!\n", _statService->get_file_name_var().c_str());

  //print the first row with stats headers:
  FH.width(12);
  FH.setf(ios::left,ios::adjustfield);
  FH<<"generation";

  //print the stat names:
  _statService->printStatHeaders(FH,(GEN | FLAT));

  FH<<"\n";

  //next rows:
  unsigned int gen_index;
  unsigned int nrows    = _statService->getTotOccurrence();

  //val is ordered rows x cols, row = generation, col = replicate
  //we take the variance over rows
  //_statService->setStatVariance();

  for(unsigned int i = 0; i < nrows; i++) {
    //print generation number:
    gen_index = _statService->getStatRecIndex(i);

    //if the index is my_NAN, all replicates crashed before that index
    if(gen_index != my_NAN) {
      FH.width(12);
      FH.setf(ios::left,ios::adjustfield);
      FH<<setprecision(4);
      FH<<gen_index;

      _statService->printStatVariance(FH,i);

      FH<<"\n";
    }
    else break;
  }
  FH.close();
}



