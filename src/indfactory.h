/** @file indfactory.h
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

#ifndef indfactoryH
#define indfactoryH

#include <map>
#include <deque>
#include "individual.h"
#include "types.h"
#include "ttrait.h"

class TSelection;

/**Factory of Individual, stores the individual prototype and the trait prototypes, manages the individual garbage collector.
 * Provides methods to generate new individuals within the Metapop. Each new individual is created by cloning a prototype itself
 * created at the simulation setup. New individuals are decorated with the appropriate traits as set by the trait prototypes and
 * receives a unique ID (unique within a simulation).
 **/
class IndFactory {
protected:
  /**Map of the trait prototypes.*/
  map< string,TTraitProto* > _protoTraits;

  Metapop* _popPtr;


  /**Table containing the index of each trait.*/
  map< string, int > _TraitsIndex;

  /** vector containing a pointer to the TraitPrototype (allows to access the prototypes directly) */
  vector<TTraitProto*> _protoTraitsVector;


  /**The individuals prototype used to create any new individual in a simulation.*/
  Individual _protoIndividual;

  /**Garbage collector for unused Individual's.*/
  deque<Individual*> RecyclingPOOL;

public:

  IndFactory ( ):_popPtr(0) {
  }
  virtual ~IndFactory ( ) { }

  /**Put an individual in the recycling pool.*/
  void recycle(Individual* ind) {
    assert(ind);
    RecyclingPOOL.push_back(ind);
  }

  /**Creates the individuals prototype from the selected trait prototypes.
    Resets the individual's ID counter to 0 and sets the traits index table.
    @callgraph
    @param TTlist the list of the current trait prototype selected from the current simulation parameters.
    **/
  void                    makePrototype               (map< string,TTraitProto* >& TTlist, Metapop* pMetapop);

  /**Resets the individuals prototype between replicates from the selected trait prototypes.
    Resets the individual's ID counter to 0 and sets the traits index table.
    @callgraph
    @param TTlist the list of the current trait prototype selected from the current simulation parameters.
    **/
  void                    resetPrototype               ( );
  void                    resetPrototypeTotal          ( );


  /**Creates a blank individual which has to be decorated.
   * ID is set and new traits are allocated but no genetic data is created. Sex has to be set too.
   * @callgraph
   **/
  Individual*             getNewIndividual() {return makeNewIndividual(NULL,NULL,MAL,NULL);}

  /**Creates an individual with pointers to parents, sex and home ID set but no genetic data.
   * No inheritance or mutations on the trait sequences are done.
   * @callgraph
   * @param mother ptr to the mother
   * @param father ptr to the father
   * @param sex gender of the individual
   * @param homepatch ID of the Patch where this individual is born, usually the current position in the Patch array
   **/
  Individual*    makeNewIndividual           (Individual* mother, Individual* father, sex_t sex, Patch* homepatch);

  /**Completely creates an individual with inheritance and mutations on all traits.
   * @callgraph
   * @param mother ptr to the mother
   * @param father ptr to the father
   * @param sex gender of the individual
   * @param homepatch ID of the Patch where this individual is born, usually the current position in the Patch array
   **/
  Individual*     makeOffsprg (Individual* mother, Individual* father, sex_t sex, Patch* homepatch);

  /**Individual prototype accessor.*/
  const Individual*       getIndividualProtoype       ( )   {return &_protoIndividual;}

  /**Accessor to the list of TTraitProto's.*/
  map< string,TTraitProto* >& getTraitPrototypes  ( )  {return _protoTraits;}

  TTraitProto& getTraitPrototype(int i) {
    assert(i>=0 && i<(int)_protoTraitsVector.size());
    return *_protoTraitsVector[i];
  }

  unsigned int getTraitPrototypeSize() {return _protoTraitsVector.size();}

  /**Gives the index of trait with \a type.
    @param type the type of the trait*/
  vector<int> getTraitIndex (string type);
};

#endif

