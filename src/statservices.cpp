/** @file statservices.cpp
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


#include <sstream>
#include <list>
#include "statservices.h"
#include "stathandler.h"
#include "metapop.h"
#include "output.h"

//-----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool StatServices::init(){

	if(!_tot_occurrence) return true;

  list< StatHandlerBase* >::iterator HIT;

  #ifdef _DEBUG
   message("Loaded stats:\n");
   message("  init StatHandler (");
   for(HIT = _children.begin(); HIT != _children.end(); ++HIT) {
      message(" %s", (*HIT)->getName().c_str());
   }
   message(")\n");
  #endif

  //first init the stat handlers -> will define the dims of the stat tables
  for(HIT = _children.begin(); HIT != _children.end(); ++HIT) {
	  (*HIT)->init();
  }

	// for each stat of the input file
	TMatrixVar<string> m(_statArg);
	unsigned int i, j, s1, s2;
	string statName;
	bool is_set;
	for(i=0, s1=m.getNbRows(); i<s1; ++i){
		for(j=0, s2=m.get(i)->size(); j<s2; ++j){
			is_set = false;
			statName = m.get(i,j);      // get the next stat key word

			#ifdef _DEBUG
			 // find the corresponding stat handler
			 message("  %s (", statName.c_str());
			#endif

			for(HIT = _children.begin(); HIT != _children.end(); ++HIT) {
				is_set |= (*HIT)->setStatRecorders(statName);
			}

			#ifdef _DEBUG
			 message(")\n");
			#endif

			if(!is_set) {
				warning("The string \"%s\" is not a valid statistic option and is therfore not considered!\n", statName.c_str());
			}
		}
	}


	if(_save_choice == 1){
		#ifdef _DEBUG
		 message("  LCE_StatFH::FHwrite (%s)\n",_file_stats.c_str());
		#endif

		ofstream FH(_file_stats.c_str(),ios::out);
		if(!FH) fatal("Could not open stat output file \"%s\"!\n",_file_stats.c_str());

		//print the first row with stat headers:
		FH.width(12);
		FH.setf(ios::left,ios::adjustfield);
		FH<<"replicate"<<"\t";
		FH.width(12);
		FH<<"generation";

		//print the stat names:
		printStatHeaders(FH, FLAT);

		//next rows:
		FH<<"\n";
	}

	return true;
}

//-----------------------------------------------------------------------------------------
// load
// ----------------------------------------------------------------------------------------
void StatServices::load ( SimComponent* sc ){
  sc->loadStatServices(this);
}

//-----------------------------------------------------------------------------------------
// attach
// ----------------------------------------------------------------------------------------
void StatServices::attach ( StatHandlerBase* SH){
  Service::attach(SH);
  _children.push_back(SH);
  SH->set_service(this);
}

//-----------------------------------------------------------------------------------------
// notify
// ----------------------------------------------------------------------------------------
void StatServices::notify ()
{
	if(!(_popPtr->getCurrentGeneration() % _occurrence)){
		Service::get_observers()->front()->update_patch_states();
		Service::notify();
	}
}

//-----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void StatServices::reset (){
  Service::reset();

  // StatHandlerBase themself are deleted by the corresponding traits
  _children.clear();
}

//-----------------------------------------------------------------------------------------
// print_headers
// ----------------------------------------------------------------------------------------
void StatServices::printStatHeaders(ostream& FH,unsigned int order){
  list< StatHandlerBase* >::iterator HIT = _children.begin();
  for(; HIT != _children.end(); ++HIT) {
		(*HIT)->print_headers(FH,order);
	}
}

//-----------------------------------------------------------------------------------------
// print_legend
// ----------------------------------------------------------------------------------------
void StatServices::printStatLegend(ostream& FH,unsigned int order){
  list< StatHandlerBase* >::iterator HIT = _children.begin();
	for(; HIT != _children.end(); ++HIT) {
	  (*HIT)->print_legend(FH,order);
  }
}

//-----------------------------------------------------------------------------------------
// printStatValue
// ----------------------------------------------------------------------------------------
void StatServices::printStatValue(ostream& FH, unsigned int i, unsigned int j){
  list< StatHandlerBase* >::iterator HIT = _children.begin();
  for(; HIT != _children.end(); ++HIT) {
	  (*HIT)->print_value(FH,i,j);
  }
}

//-----------------------------------------------------------------------------------------
// printStatMean
// ----------------------------------------------------------------------------------------
void StatServices::printStatMean(ostream& FH, unsigned int i){
  list< StatHandlerBase* >::iterator HIT = _children.begin();
  for(; HIT != _children.end(); ++HIT) {
	  (*HIT)->print_mean(FH,i);
  }
}

//-----------------------------------------------------------------------------------------
// printStatVariance
// ----------------------------------------------------------------------------------------
void StatServices::printStatVariance(ostream& FH, unsigned int i){
  list< StatHandlerBase* >::iterator HIT = _children.begin();
  for(; HIT != _children.end(); ++HIT) {
	  (*HIT)->print_variance(FH,i);
  }
}

//-----------------------------------------------------------------------------------------
// getStatRecIndex
// ----------------------------------------------------------------------------------------
unsigned int StatServices::getStatRecIndex(unsigned int i){
  list< StatHandlerBase* >::iterator el;
  //find the first statHandler that has registered statRecorders and has computed some stats:
  for(el = _children.begin(); el != _children.end(); el++){
    if((*el)->getStatRecIndex(i) != my_NAN) return (*el)->getStatRecIndex(i);
  }

  return my_NAN;
}


