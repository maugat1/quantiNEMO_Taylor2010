/** @file tselection.cpp
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
#include "tselection.h"
#include "tselectiontype.h"
#include "lce_breed.h"
#include "patch.h"
#include "ttquanti.h"

// -----------------------------------------------------------------------------
// destructor
// -----------------------------------------------------------------------------
TSelection::~TSelection ( )
{
  if(_selTrait){
    for(unsigned int i=0; i<_vTraitsSize; ++i){
      if(_selTrait[i]) delete _selTrait[i];
    }
    delete[] _selTrait;
  }

  if(_aPheno)   delete[] _aPheno;

  for(int i=0; i<2; ++i){
    if(_fit[i])   delete _fit[i];
  }
}

//-----------------------------------------------------------------------------
// set_phenotype
//-----------------------------------------------------------------------------
void
TSelection::set_phenotype(Individual* ind)
{
	for(unsigned int t = 0; t<_vTraitsSize; ++t){
		(_selTrait[t]->*(_selTrait[t]->func_ptr_set_phenotype))(ind);

	}
}

//-----------------------------------------------------------------------------
// LCE_Breed_fitness::
//-----------------------------------------------------------------------------
/** computes the fitness of the individuals of the patch and age. returns an array of the fitness of all individuals of the specified patch, sex and age
  * First the phenotypes will be computed and set.
  * patch:   individuals of this patch will be considered
  * age:     only this age class will be considered
  * sort:    if the selection is based on something like fittest or less fittest,
  *          it is neccessary to sort the fitness array. Therfore a two-dimensional
  *          array has to be passed which informes about the sort
  *          (0: no (default); 1: upwards; 2: downwards)?
  */
void
TSelection::set_fitness(Patch* patch, age_idx age, int* sort, int subset)
{
	if(sort){
		set_fitness(patch, FEM, age, sort[FEM], subset);
    set_fitness(patch, MAL, age, sort[MAL], subset);
  }
  else{
		set_fitness(patch, FEM, age, 10, subset);
		set_fitness(patch, MAL, age, 10, subset);
  }
}

//-----------------------------------------------------------------------------
// set_fitness::
//-----------------------------------------------------------------------------
/** returns an array of the fitness of all individuals of the specified patch, sex and age
  * First the phenotypes will be computed and set.
  * patch:   individuals of this patch will be considered
  * sex:     only this sex will be considered
  * age:     only this age class will be considered
  * isCumul: should the array contain cumulative fitnesses (0: no, 1: yes);
  * sort:    sould the array be sorted on the bases of the fitness (0: no (default); 1: upwards; 2: downwards)?
  */
void
TSelection::set_fitness(Patch* patch, sex_t sex, age_idx age, int sort, int subset)
{
  unsigned int i, newSize;
  Individual* ind;

  // if the patch is empty for this sex, return
	newSize = patch->size(sex, age);

	// resize TPatchFitness to correct number of individuals
	if(!_fit[sex]) _fit[sex] = new TPatchFitness(newSize);
	else _fit[sex]->resize(newSize);

	if(!newSize) return;

  // get the selection pressure of the patch
  for(i=0; i<_vTraitsSize; ++i){
    _selTrait[i]->get_selection_pressure(patch, sex);
  }


  // get the current arrays
  double*      aFit = _fit[sex]->_aFit;
  Individual** aInd = _fit[sex]->_aInd;

  // for each individual
  for(i=0; i<newSize; ++i){
    aInd[i] = ind = patch->get(sex, age, i);            // get the individual
		set_phenotype(ind);              									  // set the phenotypes of this individual at each trait
    aFit[i] = _get_fitness_multiplicative();            // compute overall fitness
    ind->setFitness(aFit[i]);                           // set the fitness to the individual
  }

  _fit[sex]->sort(sort, subset);                        // sort
}

//-----------------------------------------------------------------------------
// sort_fitness
//-----------------------------------------------------------------------------
/** sort the array if nbecessary and make it cumulative */
void
TSelection::sort_fitness(sex_t SEX, int how, int subset)
{
	assert(_fit[SEX]->_sort == 10);
	_fit[SEX]->sort(how, subset);
}

//-----------------------------------------------------------------------------
// LCE_Breed_fitness::
//-----------------------------------------------------------------------------
/** fitness are multiplicative */
double
TSelection::_get_fitness_multiplicative(){
  double w = 1;
  for(unsigned int t=0; t<_vTraitsSize; ++t){
		w *= _selTrait[t]->get_fitness();
	}
  return w;
}


//-----------------------------------------------------------------------------
// LCE_Breed_fitness::init
//-----------------------------------------------------------------------------
void
TSelection::init(Metapop* popPtr, vector<int> linked_traits)
{
	// for each sex
	for(int i=0; i<2; ++i){
		_fit[i]        = NULL;
	}
	_selTrait = NULL;

	// get the linked traits
  _popPtr      = popPtr;
	_vTraits     = linked_traits;
	_vTraitsSize = _vTraits.size();
  _aPheno      = new double[_vTraitsSize];
  _nbPops      = _popPtr->getPatchNbr();

	_optima_sd        = sqrt(_popPtr->get_paramset()->getValue("patch_stab_sel_optima_var"));
	_intensity_sd     = sqrt(_popPtr->get_paramset()->getValue("patch_stab_sel_intensity_var"));

	_min_sd           = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_min_var"));
	_max_sd           = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_max_var"));
	_growth_rate_sd   = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_growth_rate_var"));
	_max_growth_sd    = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_max_growth_var"));
	_symmetry_sd      = sqrt(_popPtr->get_paramset()->getValue("patch_dir_sel_symmetry_var"));

  // create the TSelectionTrait objects
  _selTrait = new TSelectionTrait*[_vTraitsSize];
  unsigned int selKind = 0;                 // bit: 1: neutral, 2: stabilizing, 4: directional
  for(unsigned int t = 0; t < _vTraitsSize; ++t) {
    // find the kind of selection of the trait
    assert((dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t]))));
    switch((dynamic_cast <TTQuantiProto*> (&_popPtr->getTraitPrototype(_vTraits[t])))->get_selection_model()){
      case 0: // stabilizing: 2. digit
              _selTrait[t] = new TSelectionStabilizing(this, t);
              selKind |= 2;
              break;
      case 1: // directional: 3. digit
              _selTrait[t] = new TSelectionDirectional(this, t);
              selKind |= 4;
              break;
      case 2: // neutral:     1. digit
              _selTrait[t] = new TSelectionNeutral(this, t);
              selKind |= 1;
              break;
    }
  }


	// set the nb linked traits in the pops
	for(int p=0; p< _nbPops; ++p){
    _popPtr->getPatch(p)->set_LinkedTraits(_vTraitsSize, _selTrait, (bool)_popPtr->get_sexInitRatio());
  }

	// environmental variance
	set_ve_mean();
	set_ve_var();

	// what must be read?
  if(selKind & 2){      // stabilizing selection
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_stab_sel_optima",    "optima",    &Patch::set_localOptima);    // default 0
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_stab_sel_intensity", "intensity", &Patch::set_localIntensity); // default 1
	}
	if(selKind & 4){      // directional selection
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_growth_rate", "growth_rate", &Patch::set_localGrowthRate);  // default 1
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_max_growth",  "max_growth",  &Patch::set_localMaxGrowth);   // default 1
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_symmetry",    "symmetry",    &Patch::set_localSymmetry);    // default 0.5
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_min",         "min",         &Patch::set_localMin);         // default 0
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_dir_sel_max",         "max",         &Patch::set_localMax);         // default 1
  }
}

// ----------------------------------------------------------------------------------------
// reset_Ve
// ----------------------------------------------------------------------------------------
	// set the environmental mean if used and the corresponding function pointer
void
TSelection::set_ve_mean(){
	_Ve_mean_set = (_popPtr->get_parameter("patch_ve_mean")->isSet() ||_popPtr->get_parameter("patch_ve_mean_fem")->isSet());
	if(_Ve_mean_set){
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_ve_mean", "mean_environmental_effect", &Patch::set_localMeanVe);
	}
	for(unsigned int t = 0; t < _vTraitsSize; ++t){    // set the function pointers
		_selTrait[t]->set_ve_mean_func_ptr();
	}
}

// ----------------------------------------------------------------------------------------
// reset_Ve
// ----------------------------------------------------------------------------------------
	// set the environmental variance
	// here only the input is passed to the patches,
	// Ve itself can only be computed when pops are populated if h2 or H2 is defined
	// therfore this has to be done in the function execute_before_each_replicate()
	// at this time the corresponding function pointer will also be set
void
TSelection::set_ve_var(){
	_Ve_var_set = (_popPtr->get_parameter("patch_ve_var")->isSet() || _popPtr->get_parameter("patch_ve_var_fem")->isSet());
	if(_Ve_var_set){            // Ve is set by the parameter patch_ve_var (has precedence)
		_popPtr->set_patch_parameter(_vTraitsSize, "patch_ve_var", "heritability", &Patch::set_localh2Ve);

		// heritability cannot change over time for environmental model 1 and 3
		if(_popPtr->get_parameter("patch_ve_var")->isTemporalParam()){
			for(unsigned int t = 0; t < _vTraitsSize; ++t) {
				if(_selTrait[t]->get_Ve_model() == 1 || _selTrait[t]->get_Ve_model() == 3){
					fatal("The heritability cannot change over time if the parameter 'quanti_environmental_model is set to 1 or 3!\n");
				}
			}
		}
	}
	else{                    	// Ve is set by the parameter quanti_heritability
		for(unsigned int t = 0; t < _vTraitsSize; ++t){
			_selTrait[t]->set_quantiHeritability();
		}
	}
}

// ----------------------------------------------------------------------------------------
// reset_Ve
// ----------------------------------------------------------------------------------------
/** compute the environmental variance for each patch and trait
	* all is set by default to true, i.e. all Ve are set
  * if all=False, the environment has to be set at each generation (model 3)
  */
void
TSelection::reset_Ve(bool all){
	// set the Ve in each patch and trait
  for(unsigned int t = 0; t < _vTraitsSize; ++t){
    if(all || _selTrait[t]->get_Ve_model() == 2 || _selTrait[t]->get_Ve_model() == 4){
      _selTrait[t]->set_Ve();
    }
	}
}

// ----------------------------------------------------------------------------------------
// executeBeforeEachReplicate::
// ----------------------------------------------------------------------------------------
/** compute the environmental variance if used */
void
TSelection::executeBeforeEachReplicate(const int& rep){
	reset_Ve(true); // it has to be done here, since the pops have to be populated

}

// ----------------------------------------------------------------------------------------
// executeBeforeEachGeneration::
// ----------------------------------------------------------------------------------------
/** compute the environmental variance if used */
void
TSelection::executeBeforeEachGeneration(const int& gen){
	reset_Ve(false); // check if the Ve has to be recomputed
}



// ----------------------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------------------








