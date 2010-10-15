/** @file patch.h
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
//---------------------------------------------------------------------------

#ifndef patchH
#define patchH

#include <list>
#include <deque>
#include <map>
#include <time.h>
#include "types.h"
#include "indfactory.h"
#include "metapop_sh.h"
//---------------------------------------------------------------------------
//CLASS PATCH

/**Second class in the metapopulation design structure, between the Metapop and Individual classes.
 * The Patch class is an abstraction of a sub-population or patch concept (also called a deme) in a metapopulation context.
 * It contains the individual containers for the different age classes. Three main age classes are currently implemented,
 * the offspring, post-dispersal and adult classes (see the age_t enum) which are all subdivided into male and female individuals.
 * These containers are accessed using the interface defined here or through the Metapop class interface.
 * The different LCE will use these interfaces to handle the individuals. They are also responsible to change the age flag of the population (see Metapop).
 *
 * The individuals are accessed using their age index value and not the age flag value (e.g., using the Patch::get() method).
 * These indexes are defined in the age_idx enum (see type.h) and differ from the age class flag values. For instance the adults'
 * containers have the index 2 (ADLTx = 2) whereas their age flag is 4 (ADULTS = 4). This might be confusing but it saves a lot
 * of checks at runtime! It also allows to give a flag containing several class bits set instead of a unique index value when needed
 * (see the Patch::size() and Metapop::size() suite of functions).
*/

class TSelectionTrait;
class Patch
{
  unsigned int _ID;               // id of the patch (starting at 0)
  string       _IDstr;            // same as _ID (but starting at 1), but stored as a string   int2str(_ID+1)
	static unsigned int _nbPatch;   // total number of patches
	static Metapop* _popPtr;


  /** counter for the individual id (starts at 0 for each simulation) */
  unsigned long _ID_individual;

  /**Carrying capacity for males and females*/
  unsigned int _K;          // _KFem + _KMal

  /**Sex specific carrying capacity*/
  unsigned int _KFem, _KMal;

  /** Sex specific inital population sizes */
  unsigned int _N_ini, _NFem_ini, _NMal_ini;

  // selection pressure  selection[sex][trait][param]
  double*** _localSelection;         // if stabilizing selection:     _localSelection[sex][trait][pos]
                                     //   1. pos: optima       (default: 0)
                                     //   2. pos: intensity    (default: 1)
                                     // if directional selection:
                                     //   1. pos: min          (min fitness,                 default: 0)
                                     //   2. pos: max          (max fitness,                 default: 1)
																		 //   3. pos: max growth   (phenotype of maximal growth, default: 1)
																		 //   4. pos: growth rate  (steepness of the slope,      default: 1)
                                     //   5. pos: symmetry     (symmetry of the slope,       default: 1)
  TSelectionTrait** _pSelectionType; // pointer to the selection trait neutral, stabilizing, or directional

  double*   _meanVe[2];              // mean environmental Factor  _meanVe[sex][trait]
  double*   _sdVe[2];                // environmental sd computed not input   _sdVe[sex][trait]
  double*   _h2[2];                  // Ve or h2: stored input     _sdVe[sex][trait]
  unsigned int _nbLinkedTraits;



  /**Extinction flag*/
  bool _isExtinct;

  /**Containers size counters, sex X age.*/
	unsigned int** _sizes;               //  _sizes[sex][age]

	/**Total size of the containers, amount of allocated memory.*/
	unsigned int** _capacities;          //  _capacities[sex][age]

	/**Individuals containers, sex X age.*/
	deque <Individual*>** _containers;   //  _containers[sex][age]

 public:
	//counters:
  unsigned int nbEmigrant, nbImmigrant, nbPhilopat, nbKolonisers;

  //methods:
	Patch() : _ID(0), _ID_individual(0), _K(0), _KFem(0), _KMal(0), _N_ini(0), _NFem_ini(0), _NMal_ini(0),
    _localSelection(0), _pSelectionType(0), _nbLinkedTraits(0), _isExtinct(0)
	{
		_sizes      = ARRAY::new_2D<unsigned int>(2, NB_AGE_CLASSES, (unsigned int)0);
		_capacities = ARRAY::new_2D<unsigned int>(2, NB_AGE_CLASSES, (unsigned int)0);
		_containers = ARRAY::new_2D<deque <Individual*> >(2, NB_AGE_CLASSES);

		for(unsigned int i = 0; i < 2; i++) {  // for each sex
      _meanVe[i]= NULL;
      _sdVe[i]  = NULL;
      _h2[i]    = NULL;
    }
  }

  ~Patch();

  Patch*        init                      (unsigned int id);
  inline void   set_ID_individual         (unsigned int i) {_ID_individual = i;}   // used when genotypes are passed
  void          set_PopSizes_ini          (unsigned int nbfem, unsigned int nbmal);
  void          set_PopSizes_ini_carrying_capacity();

  ///@name Setters
  ///@{
  void          set_isExtinct             (bool status)    {_isExtinct = status;}
  void          set_localMeanVe		        (double *array, sex_t SEX);
  void          set_localh2Ve		          (double *array, sex_t SEX);
  void          set_localOptima		        (double *array, sex_t SEX);
  void          set_localIntensity		    (double *array, sex_t SEX);
  void          set_localGrowthRate		    (double *array, sex_t SEX);
  void          set_localMaxGrowth		    (double *array, sex_t SEX);
  void          set_localSymmetry 		    (double *array, sex_t SEX);
  void          set_localMin 		          (double *array, sex_t SEX);
  void          set_localMax      		    (double *array, sex_t SEX);
  bool          set_Ve                    (const int& model, const int& trait, sex_t sex);
	void          set_h2                    (const double& v, const int& i){   // input of the h2
		if(_h2[0]) _h2[0][i] = v;   // check if there are males present
		_h2[1][i] = v;              // females are always present
  }
  static void   set_nbPatch               (const unsigned int& n){_nbPatch=n;}
	static void   set_pMetapop              (Metapop* p)           {_popPtr=p;}

	void          set_localParameter(double* array, sex_t sex,
                                  void (Patch::*pt2Func)(double*, sex_t));
  void          reset_ID_individual       ()               {_ID_individual=0;}


  ///@}
  ///@name Getters
  ///@{
  inline unsigned int  get_ID             ()               {return _ID;}
  string        get_next_IDofIndividual   ()               {return STRING::int2str(_ID_individual++)+"_"+_IDstr;}  // "234_1": individual 234 of patch 1
  unsigned int  get_K                     ()               {return _K;}
  unsigned int  get_KFem                  ()               {return _KFem;}
  unsigned int  get_KMal                  ()               {return _KMal;}
  unsigned int  get_N_ini                 ()               {return _N_ini;}
  unsigned int  get_NFem_ini              ()               {return _NFem_ini;}
  unsigned int  get_NMal_ini              ()               {return _NMal_ini;}
  double        get_density               (age_t AGE)      {return (double)size(AGE)/_K;}
  double        get_densityF              (age_t AGE)      {return (double)size(FEM, AGE)/_KFem;}
  double        get_densityM              (age_t AGE)      {return (double)size(MAL, AGE)/_KMal;}
  double        get_density               (age_idx AGE)    {return (double)size(AGE)/_K;}
  double        get_densityF              (age_idx AGE)    {return (double)size(FEM, AGE)/_KFem;}
	double        get_densityM              (age_idx AGE)    {return (double)size(MAL, AGE)/_KMal;}

	unsigned int**        &get_sizes        ()               {return _sizes;}
	unsigned int**        &get_capacities   ()               {return _capacities;}
	deque <Individual*>** &get_containers   ()               {return _containers;}



  bool          get_isExtinct             ()               {return _isExtinct;}
  double**      get_localSelection        (sex_t SEX)      {return _localSelection[SEX];}
  double*       get_localSelection        (sex_t s, int t) {return _localSelection[s][t];}
  double        get_localSelection        (sex_t s, int t, int p) {return _localSelection[s][t][p];}
  double*       get_meanVe                (sex_t s)        {return _meanVe[s];}
  double        get_meanVe                (sex_t s,int pos){return _meanVe[s][pos];}
  double*       get_sdVe                  (sex_t s)        {return _sdVe[s];}
  double        get_sdVe                  (sex_t s,int pos){return _sdVe[s][pos];}
  double*       get_h2                    (sex_t s)        {return _h2[s];}
  double        get_h2                    (sex_t s,int pos){return _h2[s][pos];}
  int           get_nbLinkedTraits        ()               {return _nbLinkedTraits;}
  bool          isEmpty                   ()               {return (size(ALL) == 0);}
  unsigned int  getAdultsNumber           ()               {return size(ADLTx);}
  unsigned int  getOffspringNumber        ()               {return size(OFFSx);}
  unsigned int  getIndividualNumber       ()               {return size(ALL);}
  double        getPatchOffsprgSatLevel   ()               {return (double)size(OFFSx)/_K;}
  double        getPatchAdultsSatLevel    ()               {return (double)size(ADLTx)/_K;}

  unsigned int  size       (sex_t SEX, age_t AGE);
  unsigned int  size       (age_t AGE);
  unsigned int  size       (sex_t SEX, age_idx AGE);
  unsigned int  size       (age_idx AGE);
  Individual*   get        (sex_t SEX, age_idx AGE, unsigned int at);
  void          set        (sex_t SEX, age_idx AGE, unsigned int at, Individual* ind);
  void          add        (sex_t SEX, age_idx AGE, Individual* ind);
  void          assign     (sex_t SEX, age_idx AGE, unsigned int n);
  void          remove     (sex_t SEX, age_idx AGE, unsigned int at);
	void          remove     (sex_t SEX, age_idx AGE, Individual* ind);
	void          recycle    (sex_t SEX, age_idx AGE, unsigned int at);
	void          recycle    (sex_t SEX, age_idx AGE, Individual* ind);
	void          move       (sex_t SEX, age_idx from, age_idx to, unsigned int at);
  void          move       (sex_t SEX, age_idx from, age_idx to);
  void          move       (age_idx from, age_idx to);
  void          swap       (sex_t SEX, age_idx from, age_idx to);
  void          swap       (age_idx from, age_idx to);
  void          clear      (sex_t SEX, age_idx AGE);
	void          flush      (sex_t SEX, age_idx AGE);
	void          flush      (age_idx AGE);
	void          flush      (age_t AGE);
	void          flush      ();
  void reset_counters();

  void set_LinkedTraits(const unsigned int& size, TSelectionTrait** type, bool bothSexes);
  void reset_LinkedTraits();

  /** return an array with all phenotypes of the specific trait, sex and age */
	double* getPhenotypes(sex_t SEX, age_idx AGE, const int& trait);

	/** pos dize regulation based on their fitness */
	void regulate_selection_fitness(const unsigned int& K, TSelection* pSel, sex_t SEX, age_idx AGE);
	void regulatation_selection_draw_survivors(TSelection* pSel, const unsigned int& K, sex_t SEX, age_idx AGE);
	void regulatation_selection_draw_loosers(TSelection* pSel, const unsigned int& K, sex_t SEX, age_idx AGE);

  /**Fills the patch containers corresponding to the age flags passed, for both sexes.*/
  void setNewGeneration(age_t AGE);

  /**Fills the patch container corresponding to the age class index passed, for both sexes.*/
	void setNewGeneration(age_idx AGE);

  void set_carrying_capacities(const unsigned int& fem, const unsigned int& mal);

};
#endif
