/** @file individual.h
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


#ifndef individualH
#define individualH

#include <map>
#include <deque>
#include "types.h"
#include "ttrait.h"

/**This class contains traits along with other individual information (sex, pedigree, etc.\ ).
 * The Individual class can be view as a trait container. It encapsulates the basic interface to manipulate
 * the traits an individual is carrying like initialization, inheritance, mutation and phenotype getter. It also
 * stores the basic individual info like its sex, pedigree (home Patch, mother and father id, etc.) and its
 * mating and fecundity values. It is not a StorableComponent but declares a similar interface to store and retrieve
 * all the previous info to/from a binary file.
 * All individuals in a simulation are instantiated by a call to the IndFactory class (through its derived class Metapop)
 * which contains the individual and traits prototypes.
 */
  class Patch;
  class LCE_Breed_fitness;
 class Individual {
private:
  /**ID tag, unique for one simulation.*/
  string _id;           // "123_1": individual 123 of patch 1

  /**Sex tag.*/
  sex_t _sex;

  /**Parents pointers.*/
  Individual *_mother, *_father;

  /** ID tags of the parents */
  string _motherID, _fatherID;    // "123_1": individual 123 of patch 1

  /**Patch tag.*/
  Patch *_natalPatch, *_currentPatch;

  /**Selfing flag.*/
	bool _isSelfed;

  /**Assigned fecundity*/
	double _fecundity;

  /** current fitness of the individual */
  double _fitness;

	/**Number of times the female mated with [0]: a male from a different Patch or [1]: from the same Patch.*/
	unsigned short _matings[2];

	/**Number of surviving offspring from the different mating categories (see matings).*/
	unsigned short _realizedFecundity[2];

  /**Number of traits in the table*/
  unsigned int _trait_nb;

public:
  typedef int IDX;

  /**The traits table.*/
  deque<TTrait*> Traits;

  Individual ();
  ~Individual () {clearTraits();}

  /**Inits parameters and traits.
    @callgraph
  **/
  Individual*     init                    ();

  /**Resets parameters and traits values.
    @callgraph
  **/
  void            reset                   ();
  void            reset_traits            ();

  ///@name Setters
  ///@{
  void            setID                  (string value)         {_id = value;}
  void            setFather              (Individual* f)        {_father = f; if(f){_fatherID = f->getID();}}
  void            setMother              (Individual* m)        {_mother = m; if(m){_motherID = m->getID();}}
  void            setFatherID            (string f)             {_fatherID = f;}
  void            setMotherID            (string m)             {_motherID = m;}
  void            setNatalPatch          (Patch* p)             {_natalPatch = p;}
  void            setCurrentPatch        (Patch* p)             {_currentPatch = p;}
  void            setSex                 (sex_t sex)            {_sex = sex;}
  void            setIsSelfed            (bool s)               {_isSelfed = s;}
  void            setFitness             (double value)         {_fitness = value;}
  ///@}

  ///@name Getters
  ///@{
  string          getID                  ()                      {return _id;}
  unsigned int    getID_individual      ()                      {return STRING::str2int<unsigned int>(_id.substr(0,_id.find('_')));}
  Individual*     getFather              ()                      {return _father;}
  Individual*     getMother              ()                      {return _mother;}
  string          getFatherID            ()                      {return _fatherID;}
  string          getMotherID            ()                      {return _motherID;}
  Patch*          getNatalPatch          ()                      {return _natalPatch;}
  Patch*          getCurrentPatch        ()                      {return _currentPatch;}
  sex_t           getSex                 ()                      {return _sex;}
  bool            isFemale               ()                      {return (_sex == FEM);}
  /**@param cat the mating category, 0 = between Patch, 1 = within Patch*/
  unsigned short  getMatings             (unsigned int cat)      {return _matings[cat];}
  /**@param cat the mating category, 0 = between Patch, 1 = within Patch*/
  unsigned short  getRealizedFecundity   (unsigned int cat)      {return _realizedFecundity[cat];}
  unsigned int    getTotRealizedFecundity()                      {return _realizedFecundity[0] + _realizedFecundity[1];}
  unsigned int    getTotMatings          ( )                     {return _matings[0] + _matings[1];}
  bool            getIsSelfed            ( )                     {return _isSelfed;}
  double          getFitness             ( )                     {return _fitness;}
  ///@}

  /**Resets the mating and fecundity counters.*/
  void            reset_counters         ( )                     {
   _matings[0] = _matings[1] = 0;
   _realizedFecundity[0] = _realizedFecundity[1] = 0;
  }


  void switch_sex(age_idx AGE, const int& i);

  /**Adds a mating event to the right mating category.
   * @param category the mating category, 0 = between Patch, 1 = within Patch
   */
  void   addMating  (unsigned int category) {_matings[category]++;}

  /**Increments the fecundity counter following the mating category
   * @param category the mating category, 0 = between Patch, 1 = within Patch
   */
  void   DidHaveABaby  (unsigned int category) {
   _matings[category]++;
   _realizedFecundity[category]++;
  }

  /**Returns the proportion of succesfull local matings.
    * @return _realizedFecundity[1]/_matings[1]
   */
  double getFecWithHomePatchMate  () {return (_matings[1] != 0 ? (double) _realizedFecundity[1]/_matings[1] : 0.0);}
  /**Returns the proportion of successfull remote matings.
    *@return _realizedFecundity[0]/_matings[0]
   */
  double getFecWithOtherPatchMate () {return (_matings[0] != 0 ? (double) _realizedFecundity[0]/_matings[0] : 0.0);}
  ///@}

  ///@name Trait interface
  ///@{
  unsigned int getTraitNumber() { return _trait_nb;}

  /**Traits table accessor.
    * @return the traits table
   */
  deque<TTrait *>& getTraits() {return Traits;}

  /**Sets the phenotype/value of a trait to a particular value.
    * @param T the trait index
   * @param value the value passed to the trait
   */
  void*  setTrait (IDX i, void* value){
   getTrait(i)->set_trait(value);
   return value;
  }

  /**Calls the value setting procedure of a particular trait.
    * @param T the trait index
   */
  void   setTraitValue (IDX T){
    getTrait(T)->set_value();
  }
  void   setTraitValue (IDX T, double values){
    getTrait(T)->set_value(values);
  }

  /**Gets the value (phenotype) of a particular trait.
    * @param T the trait type
   */
  double  getTraitValue (IDX T){
    return getTrait(T)->get_value();
  }

  /**Gets the value (genotype) of a particular trait.
    * @param T the trait type
   */
  double  getTraitGenotype (IDX T){
		return getTrait(T)->get_genotype();
	}

  /**Gets the value (phenotype) of a particular trait.
    * @param T the trait type
   */
  double  getTraitPhenotype (IDX T){
		return getTrait(T)->get_phenotype();
  }

	/**Gets the fitness factor of a particular trait (only used for quantitative traits).
    * @param T the trait type
   */
	double  getFitnessFactor (IDX T){
		return getTrait(T)->get_fitnessFactor();
	}

	/**Cpmutes the fitness factor of a particular trait (only used for quantitative traits).
		* @param T the trait type
	 */
	void  setFitnessFactor (IDX T){
		getTrait(T)->set_fitness_factor();
  }

	/**Trait accessor.
    * @param T the trait index
   */
  TTrait*  getTrait (IDX T)
  {
    assert(T >= 0 && T < (int)_trait_nb);
    return Traits[T];
  }

  /**Adds a trait to the table.
   * @param theTrait pointer to the trait to add
     @param pos the position where the trait should be added
   */
  void  addTrait (TTrait* theTrait, IDX pos){
    if((int)_trait_nb != pos) fatal("Individual::adding a trait to the wrong position (%i) (size %i)!\n",pos,_trait_nb);
    Traits.push_back(theTrait);
    _trait_nb++;
  }

  /**Removes a trait from the table.
    * @param T the trait index
   */
  void  removeTrait (IDX T){
   delete Traits[T]; Traits.erase(Traits.begin() + T); _trait_nb--;
  }

  /**Clears the traits container.**/
  void  clearTraits ( ){
    TTraitProto::deleteGeneticMap();

    for(unsigned int i = 0; i < _trait_nb; ++i){
      delete Traits[i];
    }
    Traits.clear();
    _trait_nb = 0;
  }

  /**Calls the inheritance procedure of a particular trait.
    * @param T the trait index
   * @param mother the mother
   * @param father the father
   */
  void  inheritTrait (IDX T, Individual* mother, Individual* father){
    getTrait(T)->inherit(mother->getTrait(T), father->getTrait(T));
  }

  /**Calls the recombination procedure of the traits.
    *This function has to be called only once when a new offspring is created.
    */
  void  recombine (){
    TTraitProto::recombine();
  }

  /**Calls the mutation procedure of a particular trait.
    * @param T the trait index
   */
  void  mutateTrait (IDX T){
    getTrait(T)->mutate();
  }

  /**The creation of a trait includes its inheritance, mutation and value setting.
    @param i the index of the trait to create
    @param mother the mother
    @param father the father
    **/
  void create (IDX i, Individual* mother, Individual* father){
    assert(mother && father);
    TTrait* T = Traits[i];
    T->inherit(mother->getTrait(i),father->getTrait(i));
    T->mutate();
    T->set_value();
  }

  /**Creates an individual's genotypes and phenotypes for first generation.**/
  void  create (){
    for(unsigned int i = 0; i < _trait_nb; i++) {
      Traits[i]->ini_sequence(_natalPatch);
      Traits[i]->set_value();
    }
  }
  /**Creates an individual, inherit, mutate and set all its traits.
    * @param mother the mother
    * @param father the father
    **/
  void  create (Individual* mother, Individual* father)
  {
    assert(mother && father);

    recombine();  // make a recombination if the genetic map exists

    TTrait *TT;
    for(unsigned int i = 0; i < _trait_nb; i++) {
      TT = Traits[i];
      TT->inherit(mother->getTrait(i), father->getTrait(i));
      TT->mutate();
      TT->set_value();
    }
  }

  /**Calls the inheritance procedure of all the traits present in the individual.
    * @param mother the mother
    * @param father the father
    **/
  void  inherit (Individual* mother, Individual* father);

  /**Calls the mutation procedure of all the traits present in the individual. */
  void mutate (){
    for(unsigned int i = 0; i < _trait_nb; i++){
      Traits[i]->mutate();
    }
  }

  ///@}
  /**Write some info to stdout.*/
  void  show_up();
  /**Cloning procedure, clones all the traits present in the individual.*/

  Individual*  clone ();
  ///@name Operators
  ///@{
  /**Assignment, make a deep copy of the parameter values and traits.*/
  Individual& operator=(const Individual& i);
  /**Only checks for traits equivalence. */
  bool operator==(const Individual& i);
  bool operator!=(const Individual& i);
  ///@}
};

#endif //INDIVIDUAL_H

