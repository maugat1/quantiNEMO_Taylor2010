/** @file tselection.h
*
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

#ifndef tselectionH
#define tselectionH

#include "simcomponent.h"
#include "tpatchfitness.h"

class TTQuantiProto;
class LCE_Breed;
class LCE;
class Metapop;
class Patch;
class Individual;
class FileServices;
class StatServices;
class TSelectionTrait;


class TSelection: public SimComponent {
private:
  // object pointers
  Metapop*            _popPtr;
  TTQuantiProto**     _pQuantiProto;


  TPatchFitness*      _fit[2];          // pointer to obejcts (male/female) containg the fitness, and individuals

  double*             _aPheno;          // temp array to stoe the phenotype of all traits

  // linked traits
	vector<int>         _vTraits;         // vector of indexes to the qTraits
  unsigned int        _vTraitsSize;     // number of qTraits

  // number of populations (used for the destructor as metapop is already destroit)
  int                 _nbPops;

  TSelectionTrait**   _selTrait;        // array for each selective trait

  double              _Ve_prop;         // Proportion of the used Ve (0: only natal patch; 1: only mating patch(default))
         int     _Ve_model;            // environmental model: 0: Ve directly, 1: h2; 2: h2 at every gen; 3: H2; 4: H" at every gen

  double _optima_sd;
  double _intensity_sd;
  double _growth_rate_sd;
  double _max_growth_sd;
  double _symmetry_sd;
  double _min_sd;
  double _max_sd;

  bool  _Ve_mean_set;     // is there a constant environmental factor?
  bool  _Ve_var_set;      // is there an evnironmental effect set for each patch separately?

 //	double*             _Ve_h2;           // array of the heritability (for each trait)

  double _get_fitness_multiplicative();

  /** phenotype specific functions ********************************************/
	void    set_phenotype(Individual* ind);

public:
  // constructor
	TSelection(Metapop* popPtr, vector<int> t){init(popPtr,t);}

	void init(Metapop* popPtr, vector<int> linked_traits);

  // destructor
  ~TSelection ( );

  //////////////////////////////////////////////////////////////////////////////
  // functions to get an individual depending on its fitness
  // random not taking into account the fitness
  inline Individual* get_RAND_noFit(sex_t sex){
    return _fit[sex]->get_RAND_noFit();
  }
  inline Individual* get_RAND_noFit_index(sex_t sex, unsigned int& index){
    return _fit[sex]->get_RAND_noFit_index(index);
  }
  inline Individual* get_RAND_noFit(sex_t sex, unsigned int nb){
    return _fit[sex]->get_RAND_noFit_subset(nb);
  }

  // most fittest
  inline Individual* get_mostFit(sex_t sex){
    return _fit[sex]->get_mostFit();
  }
  inline Individual* get_RAND_mostFit(sex_t sex){
    return _fit[sex]->get_RAND_mostFit();
  }
  inline Individual* get_RAND_mostFit_index(sex_t sex, unsigned int& index){     // changes the index
    return _fit[sex]->get_RAND_mostFit_index(index);
  }
  inline Individual* get_RAND_mostFit_of_mostFit(sex_t sex, const unsigned int& nb){
    return _fit[sex]->get_RAND_mostFit_of_mostFit(nb);
  }
  inline Individual* get_RAND_mostFit_of_RAND_mostFit(sex_t sex, const unsigned int& nb){
    return _fit[sex]->get_RAND_mostFit_of_RAND_mostFit(nb);
  }

  // less fittest
  inline Individual* get_lessFit(sex_t sex){
    return _fit[sex]->get_lessFit();
  }
  inline Individual* get_RAND_lessFit(sex_t sex){
    return _fit[sex]->get_RAND_lessFit();
  }
  inline Individual* get_RAND_lessFit_index(sex_t sex, unsigned int& index){
    return _fit[sex]->get_RAND_lessFit_index(index);
  }
  inline Individual* get_RAND_lessFit_of_lessFit(sex_t sex, const unsigned int& nb){
    return _fit[sex]->get_RAND_lessFit_of_lessFit(nb);
  }
  inline Individual* get_RAND_lessFit_of_RAND_lessFit(sex_t sex, const unsigned int& nb){
    return _fit[sex]->get_RAND_lessFit_of_RAND_lessFit(nb);
  }

  //////////////////////////////////////////////////////////////////////////////

  // returns the object containg the fitness and the individuals
  // (the array has to be deleted afterwards if it is detached)
  TPatchFitness* get_FitnessObject(sex_t sex, bool detachObject = false){
    assert(_fit[sex]);
    TPatchFitness* temp = _fit[sex];
    if(detachObject) _fit[sex] = NULL;
    return temp;
  }

  // adds any fitness arrays, which can then be used. It will be deleted later!
  void set_FitnessObject(sex_t sex, TPatchFitness* obj){
    if(_fit[sex]) delete _fit[sex];
    _fit[sex] = obj;
  }

  double getMeanFitness(sex_t sex){
    assert(_fit[sex]);
    return _fit[sex]->getMeanFitness();
  }

  double getSumFitness(sex_t sex){
    if(!_fit[sex]) return 0;
    return _fit[sex]->getSumFitness();
  }

	Individual**  get_aInd(sex_t sex)  {return _fit[sex]->_aInd;}
	double*       get_aFit(sex_t sex)  {return _fit[sex]->_aFit;}
	unsigned int  get_nbInd(sex_t sex) {return _fit[sex]->_nbInd;}
	void          remove(sex_t sex, unsigned int i){_fit[sex]->remove(i);}    

  // function setting the fitness for individuals of an entire patch
  void set_fitness(Patch*,        age_idx, int* sort=NULL, int subset=0);
	void set_fitness(Patch*, sex_t, age_idx, int sort=0,     int subset=0); // only for a single sex

	void sort_fitness(sex_t SEX, int how, int subset=0);


  void executeBeforeEachReplicate(const int& rep);
	void executeBeforeEachGeneration(const int& rep);

	void      reset_Ve(bool all=true);
	void      set_ve_mean();
	void      set_ve_var();

  Metapop*  get_popPtr()                        {return _popPtr;}
  int       get_vTraits(const int& i) const     {return _vTraits[i];}
  TSelectionTrait* get_selTrait(const int& i)   {return _selTrait[i];}
  double    get_optima_sd()                     {return _optima_sd;}
	double    get_intensity_sd()                  {return _intensity_sd;}
	double    get_min_sd()                        {return _min_sd;}
  double    get_max_sd()                        {return _max_sd;}
  double    get_growth_rate_sd()                {return _growth_rate_sd;}
  double    get_max_growth_sd()                 {return _max_growth_sd;}
  double    get_symmetry_sd()                   {return _symmetry_sd;}
  bool      get_Ve_mean_set()                   {return _Ve_mean_set;}
  bool      get_Ve_var_set()                    {return _Ve_var_set;}

//  virtual TSelection* clone ( ) {return new TSelection();}

  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader ) {}
};

#endif
