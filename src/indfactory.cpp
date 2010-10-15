/** @file indfactory.cpp
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

#include "indfactory.h"
#include "random.h"
#include "output.h"
#include "patch.h"
#include "tselectiontype.h"
#include "tselection.h"

// ----------------------------------------------------------------------------------------
// makePrototype
// ----------------------------------------------------------------------------------------
void IndFactory::makePrototype(map< string,TTraitProto* >& TTlist, Metapop* pMetapop){
  #ifdef _DEBUG
   message("IndFactory::makePrototype\n");
  #endif
  _popPtr = pMetapop;

  //and clear the ttraits:
  _protoIndividual.clearTraits();

  //store the traits list for future use:
  _protoTraits = TTlist;

  //then add the traits:
  _TraitsIndex.clear();
  int i = 0;
  map< string, TTraitProto* >::iterator trait;
  for(trait = TTlist.begin(); trait != TTlist.end(); ++trait) {
    #ifdef _DEBUG
     message("IndFactory::makePrototype::addTrait: %s\n",trait->first.c_str());
    #endif

    trait->second->init(pMetapop);
    trait->second->set_absolute_index(i);
	  _protoIndividual.addTrait(trait->second->hatch(),i);
    _TraitsIndex[trait->first] = i++;
    _protoTraitsVector.push_back(trait->second);
  }

}

// ----------------------------------------------------------------------------------------
// resetPrototype
// ----------------------------------------------------------------------------------------
void IndFactory::resetPrototype(){
  #ifdef _DEBUG
   message("IndFactory::resetPrototype\n");
  #endif

  //first, reset the ID counters of the indivuals in the patch container:
  for(int i=0, size=_popPtr->getPatchNbr(); i<size; ++i){
    _popPtr->getPatch(i)->reset_ID_individual();
  }

  vector<TTraitProto* >::iterator trait;
  for(trait = _protoTraitsVector.begin(); trait != _protoTraitsVector.end(); ++trait) {
    (*trait)->reset();
  }
}

// ----------------------------------------------------------------------------------------
// resetPrototype
// ----------------------------------------------------------------------------------------
void IndFactory::resetPrototypeTotal(){
  #ifdef _DEBUG
   message("IndFactory::resetPrototypeTotal\n");
  #endif

  vector<TTraitProto* >::iterator trait;
  for(trait = _protoTraitsVector.begin(); trait != _protoTraitsVector.end(); ++trait) {
    (*trait)->resetTotal();
  }
}

// ----------------------------------------------------------------------------------------
/** returns a vector containing all indexes of the linked traits  */
vector<int> IndFactory::getTraitIndex (string type){
  vector<int> index;
  map< string, int >::iterator cur_trait = _TraitsIndex.begin();
  for(; cur_trait != _TraitsIndex.end(); ++cur_trait){
    if(cur_trait->first.find(type) != string::npos){
      index.push_back(cur_trait->second);
    }
  }
  return index;
}

// ----------------------------------------------------------------------------------------
// makeNewIndividual
// ----------------------------------------------------------------------------------------
Individual*
IndFactory::makeNewIndividual(Individual* mother, Individual* father, sex_t sex, Patch* homepatch)
{
  Individual* newind;

  if(RecyclingPOOL.empty()) {       //create new Individual
	  newind = _protoIndividual.clone();
	  newind->init();                 //allocate memory for the traits' sequences:
  }
  else {                            //recycle an Individual from the POOL
	  newind = RecyclingPOOL[0];
		RecyclingPOOL.pop_front();      //homepatch->get_ID()
		newind->reset();
		newind->reset_counters();       //to reset the matings and fecundity counters
  }

  newind->setSex(sex);
  newind->setNatalPatch(homepatch);
  newind->setID(homepatch->get_next_IDofIndividual());
  newind->setCurrentPatch(homepatch);
  newind->setMother(mother);
	newind->setFather(father);
  return newind;
}

// ----------------------------------------------------------------------------------------
// makeOffsprg
// ----------------------------------------------------------------------------------------
Individual* IndFactory::makeOffsprg(Individual* mother, Individual* father, sex_t sex, Patch* homepatch)
{
  Individual* NewOffsprg = makeNewIndividual(mother,father,sex,homepatch);

  NewOffsprg->create(mother,father);  //inherit and mutate the traits:
  NewOffsprg->setIsSelfed(mother == father);
  bool cat = (mother->getNatalPatch() == father->getNatalPatch());
  mother->DidHaveABaby(cat);
  father->DidHaveABaby(cat);
  homepatch->add(sex, OFFSx, NewOffsprg);   // add individual to the patch
  return NewOffsprg;
}

