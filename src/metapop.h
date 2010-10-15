/** @file metapop.h
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

#ifndef metapopH
#define metapopH


#include "patch.h"


class LCE;
class LCE_Breed;
class LCE_Disperse;

//CLASS METAPOP

/**Top class of the metapopulation structure, contains patches, implements the replicate and generation loops.
 * This class implements the two main loops of a simulation, the replicate and the generation loops. The replicate
 * loop iterates the generation loop which itself iterates the life cycle loop composed of the LCEs selected by the user.

 * The basic design for the metapopulation structure is a top-down chain of responsibility where the Metapop class
 * takes care of the patches it contains which are themselves concerned by the management of their individual containers.
 * The Individual class is only concerned by the management of its traits. Thereby, a metapopulation can be viewed
 * as an interleaving of containers where one container class takes care of its directly contained class only,
 * without knowledge of the upward container state.

 * The Metapop class thus implements methods used to manage and get information from the patches, to manage the parameters
 * necessary to build a population and to get the population state information.

 * <b>Population states:</b> given by the number and the position of the individuals (both spatially in demes and temporally
 * in age class containers). The Metapop::_currentAge flag is set according to the age state of the metapopulation.
 * The life cycle events modify that state by moving individuals between individuals containers within the metatpopulation.
 * The Metapop::_currentAge flag can contain the following age class bits as defined in types.h.
 * The OFFSPRNG (=1) age class bit is set whenever the offspring containers are not empty. The ADULTS (=2) age class bit is set
 * whenever the adult containers are not empty. The ALL (=3) age class is the addition of the previous tags and NONE (=0) is the negation of them.
 * The individual containers are stored in the patches and handled through the Patch class interface. Each age class is represented
 * by two containers, one for the males (index 0) and the other for the females (index 1). These containers are store in a table
 * and are accessed through their age class index as defined by the age_idx enum (see types.h). These indexes are as follows:
 * The OFFSPRNG age class has index OFFSx = 0, and the ADULTS age class has index ADLTx = 2.
 *
*/
class Metapop: public SimComponent, public IndFactory
{
private:

	/**The life cycle events*/
	list < LCE* > _theCycle;

	/**The stat handler for the population stats*/
	MetapopSH _statHandler;

	/**The Patch container*/
	deque<Patch*> _vPatch;

	/** pointer to the breed LCE */
	LCE_Breed* _pBreed_LCE;

	/** pointer to the disperse LCE */
	LCE_Disperse* _pDisperse_LCE;

  Patch** _sample_pops;           // array with pointer to patches used for stats and other outputs (do not delete the patches here)
  unsigned int _sample_pops_size; // size of the array
  void set_sample_pops();         // creates the array

	//parameters:
	/**Number of patches in the population.*/
  unsigned int _patchNbr;

  /**Patch carrying capacity.*/
  unsigned int _patchK;

  /**Sex specific carrying capacities.*/
  unsigned int _patchKfem, _patchKmal;

  /**Patch population size.*/
  unsigned int _patchN;

  /**Sex specific population sizes.*/
  unsigned int _patchNfem, _patchNmal;

  /**Number of generations to iterate.*/
  unsigned int _generations;

  /**Number of replicates to iterate.*/
  unsigned int _replicates;
  //counters:

  /**The current generation in the generation loop, starts at 1.*/
  unsigned int _currentGeneration;
  /**The current replicate in the replicate loop, starts at 1.*/
  unsigned int _currentReplicate;
  //  unsigned int Current_LC_Rank;

  /**The current age class, might be changed by the LCEs.*/
  age_t        _currentAge;

  /**Clock counter, for logging.**/
  clock_t _meanReplElapsedTime;

  /**Generation counter, for logging.**/
  unsigned int _meanGenLength;

  /** inital sex ratio */
  double  _sexInitRatio;   // 0: only females => hermaphrodite; 0.5: f = m

	FileServices* _service;

	TSelection* 	_pSelection;   // pointer to the selection stuff if needed
	int           _selection_position; // when does selectiona acts?:
																		 //      0: reproductive success (fitness of adults, nbOffs in relation to adults)
																		 //      1: reproductive success (fitness of offsprings, nbOffs in relation to adults))
																		 //      2: pre-dispersal regulation (offspring)
																		 //      3: post-dispersal regulation (adults)
																		 // 		 4: no selection at all (or when no quantitative traits are specified)
	int           _selection_level;    // 0: soft selection (patch level)
																		 // 1: soft selection (metapopulation level)
																		 // 2: hard selection (fitness is directly translated to survival/reproductive success) 

	unsigned int  _total_carrying_capacity;

	double** _density_threshold;					// allows to change the disp rate depending on the density (0: patch, 1: density, 2: change)
	Param**  _density_threshold_param;    // pointer to the param
	unsigned int _density_threshold_nbChanges;

public:

  Metapop();
  virtual ~Metapop();

  /**Inits the population parameters from the ParamSet and builds the pop (adds patches), the prototypes and the life cycle.
    Called at the start of each simulation, resets the individual garbage collector.
    @param traits the map of the selected trait prototypes from the SimBuilder
    @param LCEs the map of the selected LCEs from the SimBuilder
    @callgraph
    */
  bool init( map< string, TTraitProto* >& traits, map< int, LCE* >& LCEs );

  /** function which changes parameters over time */
  void temporal_change(const int& gen);

  /**Called to empty the patches, individuals are move to the garbage collector.*/
  void reset();

  /**Called at the end of each simulation, empties the pop and the garbage collector, Individuals are destroyed.*/
  void clear();

  /**Creates the list of LCEs from the selected ones in the SimBuilder instance.
    \b Note: the LCEs inserted in the list are the LCE templates stored in the SimBuilder. They aren't copied or cloned.
    **/
  void setLifeCycle(map< int, LCE* >& lifeCycle);

  ///@name Main loops
  ///@{
  /**Replicate loop, iterates the life cycle \a _replicates times.
    @callgraph
    */
  void Replicate_LOOP();

  /** function which is executed before/after each replicate */
  void executeBeforeEachReplicate(const int& rep);
  void executeAfterEachReplicate(const int& rep);

  /** function which is executed before/after each generation */
  void executeBeforeEachGeneration(const int& gen);
  void executeAfterEachGeneration(const int& gen);

  /** function which prints the entire genetic map to a file named "genetic_map.txt" */
  void printEntireGeneticMap();

  /**Life cycle loop, executes the list of LCEs \a _generations times.
    @param startTime the starting time of the current replicate.
    @callgraph
    */
  void Cycle(clock_t& startTime);

  ///@}

  void createPopulations();

  ///@name Population builders
  ///@{

  /** function to set the inital population sizes */
  void setInitPopulationSizes();
  void setInitPopulationSizes_1sex();        // selfing/cloning
  void setInitPopulationSizes_2sexes();      // for all other mating systems
  void setInitPopulationSizesOfPatches();
  void setInitPopulationSizesOfPatches(TMatrix* popNfem);
  void setInitPopulationSizesOfPatches(sex_t SEX, TMatrix* popNmal);
  void setInitPopulationSizesOfPatches(TMatrix* popNfem, TMatrix* popNmal);

  /** function to set the carrying capacities */
  void setCarryingCapacities();
  void setCarryingCapacities_1sex();        // selfing/cloning
  void setCarryingCapacities_2sexes();      // for all other mating systems
  void setCarryingCapacitiesOfPatches();
  void setCarryingCapacitiesOfPatches(TMatrix* popNfem);
  void setCarryingCapacitiesOfPatches(sex_t SEX, TMatrix* popNmal);
  void setCarryingCapacitiesOfPatches(TMatrix* popNfem, TMatrix* popNmal);

  /**Sets the first generation of each replicates.
    @callgraph
    */
	void setPopulation ();
	void setPopulation_FSTAT ();  // if genotypes are given by FSTAT  file

	void set_pDisperse_LCE    (LCE_Disperse* l) {_pDisperse_LCE = l;}
	void set_pBreed_LCE       (LCE_Breed* l)    {_pBreed_LCE = l;}

	///@}

  /**Sets the Patch trait optima.
    @param optima a matrix containing the patches optima and intensity (see manual). */
  void set_patch_parameter(const unsigned int& nbTrait, const string& name, const string& name_full,
                                  void (Patch::*pt2Func)(double*, sex_t));
  void set_patch_value(const unsigned int& nbTrait, const double& value, sex_t SEX,
                                  void (Patch::*pt2Func)(double*, sex_t));
  void set_patch_matrix(const unsigned int& nbTrait, TMatrix* m, sex_t SEX, const string& name, 
                                  void (Patch::*pt2Func)(double*, sex_t));

  ///@name Implementations
  ///@{
  //SimComponent implementation:
  virtual void loadFileServices ( FileServices* loader ) {}

  virtual void loadStatServices ( StatServices* loader ) {loader->attach(&_statHandler);}

  /**Patch accessor, return the ith+1 patch in the metapop.*/
  Patch*       getPatch              (unsigned int i) {return _vPatch[i];}

  /**A secure version of the getPatch() method.*/
  Patch*       getPatchPtr           (unsigned int patch);
	deque<Individual*>* getAllIndividuals();
  unsigned int getGenerations        ( ) {return _generations;}
  unsigned int getReplicates         ( ) {return _replicates;}
  unsigned int getPatchNbr           ( ) {return _patchNbr;}
  unsigned int getPatchKFem          ( ) {return _patchKfem;}
  unsigned int getPatchKMal          ( ) {return _patchKmal;}
  unsigned int getPatchCapacity      ( ) {return _patchK;}
  Patch**      get_sample_pops       ( ) {return _sample_pops;}
  unsigned int get_sample_pops_size  ( ) {return _sample_pops_size;}
  /**Returns the mean processor time per replicate in seconds.*/
  clock_t  getMeanReplElapsedTime    ( ) {return _meanReplElapsedTime;}
  /**Returns the mean number of generations performed per replicate.*/
	unsigned int getMeanGenLength      ( ) {return _meanGenLength;}
	LCE_Disperse* get_pDisperse_LCE    ( ) {return _pDisperse_LCE;}
	LCE_Breed*    get_pBreed_LCE       ( ) {return _pBreed_LCE;}
  double      get_sexInitRatio       ( ) {return _sexInitRatio;}

	///@}

  ///@name Population state interface
  ///@{
  unsigned int getCurrentReplicate   ( ) {return _currentReplicate;}
  unsigned int getCurrentGeneration  ( ) {return _currentGeneration;}
//  unsigned int getCurrentRank        ( ) {return Current_LC_Rank;}
  age_t        getCurrentAge         ( ) {return _currentAge;}

  /**Sets the age flag.
    @param age the current age. */
  void setCurrentAge(age_t age) {_currentAge = age;}

  /**Checks if the population still contains at least one individual in any sex or age class.*/
  bool isAlive( ) {return size() != 0;}

  /**Get the total number of individuals present in the population, all sex and age classes together.*/
  unsigned int size( ) {return size(ALL);}

  /**Interface to get the size of a praticular age and sex class(es).
    @param AGE age class flags
    @param SEX sex class
    */
  unsigned int size ( sex_t SEX, age_t AGE );

  /**Interface to get the size of a praticular age and sex class within a patch.
    @param AGE age class flags
    @param SEX sex class
    @param deme the focal patch
    */
  unsigned int size (sex_t SEX, age_t AGE, unsigned int deme);

  /**Simplified interface to get the size of both sexes of the appropriate age class(es) in the whole population.
    @param AGE age class flags
  */
  unsigned int size ( age_t AGE ){
    return size( FEM, AGE ) + size( MAL, AGE );
  }

  /**Simplified interface to get the size of both sexes of the appropriate age class(es) in one patch.
    @param AGE age class flags
    @param deme the focal deme
    */
  unsigned int size ( age_t AGE, unsigned int deme ){
    return size( FEM, AGE, deme ) + size( MAL, AGE, deme );
  }

  /**Returns a pointer to the appropriate individual.
    @param SEX sex class container index
    @param AGE age class container index
    @param at the index of the individual in its container
    @param deme the patch where to grab the individual*/
  Individual* get (sex_t SEX, age_idx AGE, unsigned int at, unsigned int deme);

  /**Moves an individual from a deme to an other one, both demes sizes are modified.
    @param SEX sex class container index
    @param from_age age class container index in the deme of origin
    @param from_deme index of the deme of origin
    @param to_age age class container index in the destination deme
    @param to_deme index of the destination deme
    @param at index of the focal individual in the 'from' deme
    */
  void move (sex_t SEX, age_idx from_age, unsigned int from_deme,
             age_idx to_age, unsigned int to_deme, unsigned int at);

    /**
   * @return the current replicate counter string
   */
  string  getReplicateCounter ();
  string  getReplicateCounter_ ();
  string  getReplicateCounter_r ();

  /**
   * @return the current generation counter string
   */
  string  getGenerationCounter ();
  string  getGenerationCounter_ ();
  string  getGenerationCounter_g ();

  /** returns tre if only one sex is used for the simulations */
  void    set_SexInitRatio(map< int,LCE* >& LCEs);
  void    resize_nbPopulations(unsigned int nbPatch);

  void            set_service(FileServices* ptr){_service=ptr;}
  FileServices*   get_service(){return _service;}

  void    set_replicates(unsigned int i)  {_replicates = i;}
	void    set_generations(unsigned int i) {_generations = i;}

	void    regulate_selection_fitness_patch(age_idx AGE, unsigned int* Kmal=NULL, unsigned int* Kfem=NULL);
	void    regulate_selection_fitness_metapop(age_idx AGE);
	void    regulate_selection_fitness_hard(age_idx AGE);

	void 		set_total_carrying_capacity();
	unsigned int get_total_carrying_capacity(){return _total_carrying_capacity;}
	TSelection* get_pSelection()              {return _pSelection;}
	int     get_selection_position()          {return _selection_position;}
	int     get_selection_level()             {return _selection_level;}
	void    set_selection_level(const int& i) {_selection_level = i;}

	/** functions for changing dispersal rate following pop density of a certain patch */
	void change_disp_rate_after_density(const int& gen);
	void set_change_disp_rate_after_density(); 
	void reset_dynamic_params(); 

	///@}
};


inline unsigned int Metapop::size ( sex_t SEX, age_t AGE ){
  unsigned int s = 0;
  for(unsigned int i = 0; i < _patchNbr; i++){
    s += _vPatch[i]->size(SEX, AGE);
  }
  return s;
}

inline unsigned int Metapop::size (sex_t SEX, age_t AGE, unsigned int deme){
  return getPatchPtr(deme)->size(SEX, AGE);
}

inline Individual* Metapop::get (sex_t SEX, age_idx AGE, unsigned int at, unsigned int deme){
  return getPatchPtr(deme)->get(SEX, AGE, at);
}

inline void Metapop::move (sex_t SEX, age_idx from_age, unsigned int from_deme,
                           age_idx to_age, unsigned int to_deme, unsigned int at){
  _vPatch[to_deme]->add(SEX, to_age, get(SEX, from_age, at, from_deme));
  _vPatch[from_deme]->remove(SEX, from_age, at);

}
#endif
