/** @file stathandlerbase.cpp
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

#include "stathandler.h"
#include "stat_rec_base.h"
#include "metapop.h"
using namespace std;


// ----------------------------------------------------------------------------------------
// static
// ----------------------------------------------------------------------------------------

unsigned int   StatHandlerBase::_current_replicate;
unsigned int   StatHandlerBase::_current_generation;
Patch**        StatHandlerBase::_sample_pops;              // array of the patches to sample  (do not delete the array here!!!)
unsigned int   StatHandlerBase::_sample_pops_size;         // number of patches to sample ( = _stat_pops.size())
vector<Patch*> StatHandlerBase::_current_empty_pops[NB_AGE_CLASSES];  // current un-populated pops
vector<Patch*> StatHandlerBase::_current_pops[NB_AGE_CLASSES];        // current populated pops
vector<Patch*> StatHandlerBase::_current_popsS[2][NB_AGE_CLASSES];    // current populated pops seperated by sex and age: [sex][age]
unsigned int   StatHandlerBase::_current_nbInds[NB_AGE_CLASSES];      // total number of individuals in the sampled patches
unsigned int   StatHandlerBase::_current_nbIndsS[2][NB_AGE_CLASSES];  // total number of sampled patches separated by sex and age: [sex][age]
unsigned int   StatHandlerBase::_current_nbPops[NB_AGE_CLASSES];      // current number of populated pops ( = _current_pops.size())
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool StatHandlerBase::init()
{
	_pop = get_pop_ptr();
  _nrows = _service->getTotOccurrence();
  _ncols = _pop->getReplicates();
  _index = _service->getIndex();    // overtake the array

  _sample_pops      = _pop->get_sample_pops();
  _sample_pops_size = _pop->get_sample_pops_size();

  return true;
}

// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
/** updates the patche vectors of the colonized patches
	* must be called once before stats are computed
	*/
void StatHandlerBase::update_patch_states()
{
	_current_replicate = _pop->getCurrentReplicate();
  _current_generation = _pop->getCurrentGeneration();

  // reset vectors
  int a, s;
  for(a=0; a<2; ++a){    // for each age (caution i has to be 0 (OFFSx) and 2 (ADLTx))
		_current_pops[a].clear();             // pops with age class
    _current_empty_pops[a].clear();       // pops without any offsprings

    for(s = 0; s<2; ++s){    // for each age (0: MAL, 1: FEM)
      _current_popsS[s][a].clear();       // pops with sex and age class
      _current_nbIndsS[s][a]=0;           // total number of females
    }
  }

  unsigned int sizeM, sizeF;
  for(Patch **curPatch=_sample_pops, **end=curPatch+_sample_pops_size; curPatch!=end; ++curPatch){
    // offspring related vectors
    _current_nbIndsS[FEM][OFFSx] += sizeF = (*curPatch)->size(FEM, OFFSx);
    _current_nbIndsS[MAL][OFFSx] += sizeM = (*curPatch)->size(MAL, OFFSx);

    if(sizeF)          _current_popsS[FEM][OFFSx].push_back(*curPatch);      // female offspring
    if(sizeM)          _current_popsS[MAL][OFFSx].push_back(*curPatch);      // male offsrping
    if(sizeF || sizeM) _current_pops[OFFSx].push_back(*curPatch);            // any offspring
    else               _current_empty_pops[OFFSx].push_back(*curPatch);      // without any offspring

    // adult related vectors
    _current_nbIndsS[FEM][ADLTx] += sizeF = (*curPatch)->size(FEM, ADLTx);
    _current_nbIndsS[MAL][ADLTx] += sizeM = (*curPatch)->size(MAL, ADLTx);
    if(sizeF)          _current_popsS[FEM][ADLTx].push_back(*curPatch);      // female adults
    if(sizeM)          _current_popsS[MAL][ADLTx].push_back(*curPatch);      // male adults
    if(sizeF || sizeM) _current_pops[ADLTx].push_back(*curPatch);            // any adults
    else               _current_empty_pops[ADLTx].push_back(*curPatch);      // without any adults
  }

  // get the total number of inhabited pops
  _current_nbPops[OFFSx] = _current_pops[OFFSx].size();
  _current_nbPops[ADLTx] = _current_pops[ADLTx].size();

  _current_nbInds[OFFSx] = _current_nbIndsS[FEM][OFFSx] + _current_nbIndsS[MAL][OFFSx];
  _current_nbInds[ADLTx] = _current_nbIndsS[FEM][ADLTx] + _current_nbIndsS[MAL][ADLTx];
}

// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
void StatHandlerBase::update()
{
	execute();
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void StatHandlerBase::reset ( )
{
  for(STAT_IT IT = _stats.begin(); IT != _stats.end(); ++IT) {
    delete (*IT);
  }

  _stats.clear();
  _nrows = _ncols = 0;
  clear();
}

// ----------------------------------------------------------------------------------------
// getStatRecIndex
// ----------------------------------------------------------------------------------------
unsigned int StatHandlerBase::getStatRecIndex (unsigned int i)
{
  STAT_IT el = _stats.begin();

  if(el != _stats.end()) {
    if(i > (*el)->getIndexSize()){
      fatal("StatHandlerBase::getStatRecIndex: index vector overflow !!\n");
    }

    //find the first statRecorder that has computed some stats
    for(; el != _stats.end(); el++){
      if((*el)->getIndex(i) != my_NAN) return (*el)->getIndex(i);
    }
  }

  return my_NAN;
}
// ----------------------------------------------------------------------------------------
// print_headers
// ----------------------------------------------------------------------------------------
void StatHandlerBase::print_headers(ostream& FH, unsigned int order)
{
  STAT_IT IT = _stats.begin();

  while(IT != _stats.end()) {
    if( (*IT)->getOrdering() & order) {
      FH<<"\t";
      FH.width(12);
      FH<<(*IT)->getName();
    }
    IT++;

  }
}

// ----------------------------------------------------------------------------------------
// print_legend
// ----------------------------------------------------------------------------------------
void StatHandlerBase::print_legend(ostream& FH, unsigned int order)
{
  STAT_IT IT = _stats.begin();
  for(int i=1; IT != _stats.end(); ++i, ++IT) {
    if( (*IT)->getOrdering() & order) {
      FH.width(20);
      FH.setf(ios::left,ios::adjustfield);
      FH<<(*IT)->getName();
      FH<<": " << (*IT)->getTitle() << "\n";
    }
  }
}

// ----------------------------------------------------------------------------------------
// print_value
// ----------------------------------------------------------------------------------------
void StatHandlerBase::print_value(ostream& FH,unsigned int i, unsigned int j)
{
  STAT_IT IT = _stats.begin();
  double value;

  while(IT != _stats.end()) {
    if((*IT)->getOrdering() == FLAT) {
      FH<<"\t";
      FH.width(12);
      value = (*IT)->getVal(i,j);
      if(value == my_NAN) FH << "NaN";
      else             FH << value;
    }
    IT++;
  }
}

// ----------------------------------------------------------------------------------------
// print_mean
// ----------------------------------------------------------------------------------------
void StatHandlerBase::print_mean(ostream& FH, unsigned int i)
{
  double value;

  for(STAT_IT IT = _stats.begin(); IT != _stats.end(); ++IT) {
    FH<<"\t";
    FH.width(12);

    value = (*IT)->getMean(i);
    if(value == my_NAN) FH << "NaN";
    else             FH << value;
  }
}

// ----------------------------------------------------------------------------------------
// print_variance
// ----------------------------------------------------------------------------------------
void StatHandlerBase::print_variance(ostream& FH, unsigned int i)
{
  double value;

  for(STAT_IT IT = _stats.begin(); IT != _stats.end(); ++IT) {
    FH<<"\t";
    FH.width(12);

    value = (*IT)->getVar(i);
    if(value == my_NAN) FH << "NaN";
    else             FH << value;
  }
}

