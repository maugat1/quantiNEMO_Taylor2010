/** @file LCEbreed.cpp
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

#include <deque>
#include "lce_breed.h"
#include "metapop.h"
#include "random.h"
#include "ttrait.h"
#include "tselectiontrait.h"

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** LCE_Breed ********

// ----------------------------------------------------------------------------------------
// LCE_Breed_Base   NATING FUNCTIONS
// ----------------------------------------------------------------------------------------
// get a specific individual (no stochasity) /////////////////////////////////
Individual*
LCE_Breed::Index_MatingFunc   (Patch* thePatch, unsigned int& index, sex_t sex){
  return thePatch->get(sex, ADLTx, index);         // index has a meaning
}

// random mating /////////////////////////////////////////////////////////////
Individual*
LCE_Breed::Random_MatingFunc   (Patch* thePatch, unsigned int& index, sex_t sex){
	return thePatch->get(sex, ADLTx, SimRunner::r.Uniform(_nbIndividuals[sex]));
}
Individual*
LCE_Breed::Random_Index_MatingFunc   (Patch* thePatch, unsigned int& index, sex_t sex){
	index = SimRunner::r.Uniform(_nbIndividuals[sex]);
  return thePatch->get(sex, ADLTx, index);
}
Individual*
LCE_Breed::Random_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
	return _pSelection->get_RAND_mostFit(sex);
}
Individual*
LCE_Breed::Random_Index_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
	return _pSelection->get_RAND_mostFit_index(sex, index);   // index has a meaning
}

// full polygyny (one male) //////////////////////////////////////////////////
Individual*
LCE_Breed::fullPolygyny_oneMale_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
  return thePatch->get(sex, ADLTx, 0);    // return any individual (first one)  // to be checked !!!
}
Individual*
LCE_Breed::fullPolygyny_oneMale_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
	return _pSelection->get_mostFit(sex);         // return the fittest individual
}
Individual*
LCE_Breed::fullPolygyny_oneMale_S_MatingFunc2 (Patch* thePatch, unsigned int& index, sex_t sex){
	return _pSelection->get_RAND_mostFit(sex);    // return the randomly fittest individual
}

Individual*
LCE_Breed::fullPolygyny_manyMales_MatingFunc  (Patch* thePatch, unsigned int& index, sex_t sex){
  if(_nbIndividuals[sex]<_mating_males) return Random_MatingFunc(thePatch,index,sex);                    // random mating
	else return thePatch->get(sex, ADLTx, SimRunner::r.Uniform(_mating_males));  // return any of the first x individuals
}
Individual*
LCE_Breed::fullPolygyny_manyMales_S_MatingFunc  (Patch* thePatch, unsigned int& index, sex_t sex){
  if(_nbIndividuals[sex]<_mating_males) return Random_S_MatingFunc(thePatch,index,sex);   // random mating
	else return _pSelection->get_RAND_mostFit_of_mostFit(sex, _mating_males);
}
Individual*
LCE_Breed::fullPolygyny_manyMales_S_MatingFunc2  (Patch* thePatch, unsigned int& index, sex_t sex){
	if(_nbIndividuals[sex]<_mating_males) return Random_S_MatingFunc(thePatch,index,sex);   // random mating
	else return _pSelection->get_RAND_mostFit_of_RAND_mostFit(sex, _mating_males);
}

// partial polygyny (one male) //////////////////////////////////////////////
Individual*
LCE_Breed::partialPolygyny_oneMale_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
	if(SimRunner::r.Uniform()>_mating_proportion)return Random_MatingFunc(thePatch, index, sex);
	else return thePatch->get(sex, ADLTx, 0);   // to be checked!!!
}
Individual*
LCE_Breed::partialPolygyny_oneMale_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
	if(SimRunner::r.Uniform()>_mating_proportion)return Random_S_MatingFunc(thePatch, index, sex);
	else return _pSelection->get_mostFit(sex);             // return the fittest individual
}
Individual*
LCE_Breed::partialPolygyny_oneMale_S_MatingFunc2 (Patch* thePatch, unsigned int& index, sex_t sex){
	if(SimRunner::r.Uniform()>_mating_proportion)return Random_S_MatingFunc(thePatch, index, sex);
	else return _pSelection->get_RAND_mostFit(sex);        // return the ranom fittest individual
}

// partial polygyny (many males) /////////////////////////////////////////////
Individual*
LCE_Breed::partialPolygyny_manyMales_MatingFunc  (Patch* thePatch, unsigned int& index, sex_t sex){
	if(SimRunner::r.Uniform()>_mating_proportion)return Random_MatingFunc(thePatch, index, sex);
	else return fullPolygyny_manyMales_MatingFunc(thePatch, index, sex);
}
Individual*
LCE_Breed::partialPolygyny_manyMales_S_MatingFunc  (Patch* thePatch, unsigned int& index, sex_t sex){
	if(SimRunner::r.Uniform()>_mating_proportion)return Random_S_MatingFunc(thePatch, index, sex);
	else return fullPolygyny_manyMales_S_MatingFunc(thePatch, index, sex);
}
Individual*
LCE_Breed::partialPolygyny_manyMales_S_MatingFunc2  (Patch* thePatch, unsigned int& index, sex_t sex){
	if(SimRunner::r.Uniform()>_mating_proportion)return Random_S_MatingFunc(thePatch, index, sex);
	else return fullPolygyny_manyMales_S_MatingFunc2(thePatch, index, sex);
}

// monogamy //////////////////////////////////////////////////////////////////
Individual*
LCE_Breed::Monogyny_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
  assert(_aMatingPairs[sex] && index<_aMatingPairs_size);
	if(SimRunner::r.Uniform()>_mating_proportion || _nbIndividuals[sex]<index+1)
    return Random_MatingFunc(thePatch, index, sex);
  else
    return _aMatingPairs[sex][index];              // index has a meaning
}
Individual*
LCE_Breed::Monogyny_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
  assert(_aMatingPairs[sex] && index<_aMatingPairs_size);
	if(SimRunner::r.Uniform()>_mating_proportion || _nbIndividuals[sex]<index+1)
    return Random_S_MatingFunc(thePatch, index, sex);
  else
    return _aMatingPairs[sex][index];             // index has a meaning
}

// one sex ///////////////////////////////////////////////////////////////////
Individual*
LCE_Breed::oneSex_notSameIndex_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
  if(_nbIndividuals[sex]<2) return thePatch->get(sex, ADLTx, 0);    // if there is only a single individual in the patch
  unsigned int newIndex;
  do {
		newIndex = SimRunner::r.Uniform(_nbIndividuals[sex]);
  } while(index == newIndex);    // while it is the same index
  return thePatch->get(sex, ADLTx, newIndex);
}
Individual*
LCE_Breed::partialSelfing_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
	if(SimRunner::r.Uniform()>_mating_proportion) return oneSex_notSameIndex_MatingFunc(thePatch, index, sex);  // random mating
	else return thePatch->get(sex, ADLTx, index);         // return the same female (index has a meaning)
}
Individual*
LCE_Breed::oneSex_notSameIndex_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
  if(_nbIndividuals[sex]<2) return thePatch->get(sex, ADLTx, 0);    // if there is only a single individual in the patch
  unsigned int newIndex, i=0;
  Individual* ind;
  do {
		ind = _pSelection->get_RAND_mostFit_index(sex, newIndex);
    if(++i > 1e5){ // for security reasons pick any individual
			ind = _pSelection->get_RAND_noFit_index(sex, newIndex);
    }
  } while(index == newIndex);  // while it is the same index
  return ind;
}
Individual*
LCE_Breed::partialSelfing_S_MatingFunc (Patch* thePatch, unsigned int& index, sex_t sex){
	if(SimRunner::r.Uniform()>_mating_proportion) return oneSex_notSameIndex_S_MatingFunc(thePatch, index, sex); // random mating
	else return thePatch->get(sex, ADLTx, index);          // return the same female (index has a meaning)
}



// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
LCE_Breed::LCE_Breed(int rank): LCE("breed","",rank),
	_aMatingPairs_size(0),
  _mating_system(0), _mating_proportion(1), _mean_fecundity(0),
  _growth_rate(0),
	getMother_func_ptr(0), getFather_func_ptr(0), _maleSex(MAL), _pSelection(0){

  this->add_parameter("mating_system",INT2,false,true,0,4,"0");
  this->add_parameter("mating_proportion",DBL,false,true,0,1,"1");
  this->add_parameter("mating_males",INT2,false,false,0,0,"1");
  this->add_parameter("mean_fecundity",DBL,false,false,0,0,"0");
  this->add_parameter("growth_rate",DBL,false,false,0,0,"0");
  this->add_parameter("sex_ratio",DBL,false,false,0,0,"1");
  this->add_parameter("mating_nb_offspring_model",INT2,false,true,0,5,"0");

  // selection
  this->add_parameter("breed_model",INT2, false, true,0,2,"0");

  this->add_parameter("sex_ratio_threshold",DBL,false,false,0,0,"0");

  _aMatingPairs[MAL] = _aMatingPairs[FEM] = NULL;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::execute
// ----------------------------------------------------------------------------------------
void LCE_Breed::execute () {
  #ifdef _DEBUG
   message("  LCE_Breed ");
  #endif

  (this->*breed)();

   if(_threshold != my_NAN) reset_sex_after_phentoype(OFFSx);

  #ifdef _DEBUG
   message("done! (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
   _popPtr->getPatchNbr(),
   _popPtr->size(FEM, OFFSPRG),  _popPtr->size(MAL, OFFSPRG),
	 _popPtr->size(FEM, ADULTS),   _popPtr->size(MAL, ADULTS));
  #endif
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::init
// ----------------------------------------------------------------------------------------
bool LCE_Breed::init(Metapop* popPtr)
{
  LCE::init(popPtr);

  _mating_system      = (int)this->get_parameter_value("mating_system");
  _mating_proportion  = this->get_parameter_value("mating_proportion");
	_mean_fecundity     = this->get_parameter_value("mean_fecundity");
	_growth_rate        = this->get_parameter_value("growth_rate");
	_mating_males       = (int)_paramSet->getValue("mating_males");

	// breed_model (official) and selection_model (inofficial) are redundant
	if(_paramSet->get_param("breed_model")->isSet()){
		_popPtr->set_selection_level((int)_paramSet->getValue("breed_model"));
	}

	// how should the number of offspring be computed?
	_nbOffspring_model = (int)_paramSet->getValue("mating_nb_offspring_model");
	switch(_nbOffspring_model){
		case 0: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_CarryCapacity;  break;
		case 1: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_KeepNb;         break;
		case 2: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_Fecundity;      break;
		case 3: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_RandFecundity;  break;
		case 4: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_Logistic;       break;
		case 5: setNbOffspring_func_ptr = &LCE_Breed::setNbOffspring_RandLogistic;   break;
	}

	// does selection acts at the reproduction?
	switch(_popPtr->get_selection_position()){
		case 0: // selection acts at the reproductive success
			_pSelection = _popPtr->get_pSelection();
			switch(_popPtr->get_selection_level()){
				case 0: breed = &LCE_Breed::breed_selection_patch;    break;   // patch level
				case 1: breed = &LCE_Breed::breed_selection_metapop;           // metapop level
								_popPtr->set_total_carrying_capacity();
								break;
				case 2: breed = &LCE_Breed::breed_selection_hard;     break;   // hard selection
			}
			break;

		
		case 1: // selection acts at the survival of the newborn
			_pSelection = _popPtr->get_pSelection();
			switch(_popPtr->get_selection_level()){
				case 0: breed = &LCE_Breed::breed_selection_offspring_patch;    break;   // patch level
				case 1: breed = &LCE_Breed::breed_selection_offspring_metapop;           // metapop level
								_popPtr->set_total_carrying_capacity();
								break;
				case 2: breed = &LCE_Breed::breed_selection_offspring_hard;     break;   // hard selection
			}
			break;

		// no selection at this stage
		default:  breed = &LCE_Breed::breed_selection_neutral;              break;   // neutral
							
	}

  // one sex or two?
	if(_mating_system <2) isMatingPossible_func_ptr = &LCE_Breed::isMatingPossbile_1_sex;
	else                  isMatingPossible_func_ptr = &LCE_Breed::isMatingPossible_2_sex;

  //set the mating function ptr:
	_sort[MAL] = _sort[FEM] = -1;
	if(_popPtr->get_selection_position()==0){ // selection acts at the reproductive success of the parents
    switch(_mating_system) {
				// with selection
			case 0:    // random mating hermaphrodite
        _maleSex = FEM;    // only one sex
        getMother_func_ptr = &LCE_Breed::Random_S_MatingFunc;
        getFather_func_ptr = &LCE_Breed::Random_S_MatingFunc;
        break;
      case 1:    // selfing hermaphrodite
        _maleSex = FEM;    // only one sex
        getMother_func_ptr = &LCE_Breed::Random_Index_S_MatingFunc;
        if(_mating_proportion == 1){
          getFather_func_ptr = &LCE_Breed::Index_MatingFunc;        // get exactly the same individual
        }
        else {
          getFather_func_ptr = &LCE_Breed::partialSelfing_S_MatingFunc;
        }
        break;
      case 2:    // promiscuity/random mating.
        getMother_func_ptr = &LCE_Breed::Random_S_MatingFunc;
        getFather_func_ptr = &LCE_Breed::Random_S_MatingFunc;
        break;
      case 3:    // polygyny (random subset after the fitness of the males)
        _sort[MAL] = -3;
        getMother_func_ptr = &LCE_Breed::Random_S_MatingFunc;
        if(_mating_proportion == 1){
          if(_mating_males == 1)  getFather_func_ptr = &LCE_Breed::fullPolygyny_oneMale_S_MatingFunc2;
          else                    getFather_func_ptr = &LCE_Breed::fullPolygyny_manyMales_S_MatingFunc2;
        }
        else{
          if(_mating_males == 1)  getFather_func_ptr = &LCE_Breed::partialPolygyny_oneMale_S_MatingFunc2;
          else                    getFather_func_ptr = &LCE_Breed::partialPolygyny_manyMales_S_MatingFunc2;
        }
        break;
      case 4:    // monogamy
        getMother_func_ptr = &LCE_Breed::Random_Index_S_MatingFunc;
        getFather_func_ptr = &LCE_Breed::Monogyny_S_MatingFunc;
        break;
    }
  }
	else{        // neutral mating (selection may act at a different stage)
    switch(_mating_system) {
			case 0:    // random mating hermaphrodite
					_maleSex = FEM;    // only one sex
          getMother_func_ptr = &LCE_Breed::Random_MatingFunc;
          getFather_func_ptr = &LCE_Breed::Random_MatingFunc;
				break;
      case 1:    // selfing hermaphrodite
        _maleSex = FEM;    // only one sex
        getMother_func_ptr = &LCE_Breed::Random_Index_MatingFunc;
        if(_mating_proportion == 1){
          getFather_func_ptr = &LCE_Breed::Index_MatingFunc;        // get exactly the same individual
        }
        else {
          getFather_func_ptr = &LCE_Breed::partialSelfing_MatingFunc;
        }
        break;
      case 2:     // promiscuity/random mating
        getMother_func_ptr = &LCE_Breed::Random_MatingFunc;
        getFather_func_ptr = &LCE_Breed::Random_MatingFunc;
        break;
      case 3:     // polygyny
      {
        getMother_func_ptr = &LCE_Breed::Random_MatingFunc;
        if(_mating_proportion == 1){
          if(_mating_males == 1)  getFather_func_ptr = &LCE_Breed::fullPolygyny_oneMale_MatingFunc;
          else                    getFather_func_ptr = &LCE_Breed::fullPolygyny_manyMales_MatingFunc;
        }
        else{
          if(_mating_males == 1)  getFather_func_ptr = &LCE_Breed::partialPolygyny_oneMale_MatingFunc;
          else                    getFather_func_ptr = &LCE_Breed::partialPolygyny_manyMales_MatingFunc;
        }
        break;
      }
      case 4:     // monogamy
        getMother_func_ptr = &LCE_Breed::Random_Index_MatingFunc;
        getFather_func_ptr = &LCE_Breed::Monogyny_MatingFunc;
        break;
    }
  }

  // sex ratio
  if(_mating_system < 2){ // hermaphrodite
    setSexRatio_func_ptr  = &LCE_Breed::setSexRatio_Selfing;
    getRandomSex_func_ptr = &LCE_Breed::getRandomSex_Selfing;
  }
  else {
    // is a sex ratio specified?
		_sex_ratio = this->get_parameter_value("sex_ratio");        // input: males/females
		if(_sex_ratio<0){     // is the sex ratio within the limits?
      setSexRatio_func_ptr  = &LCE_Breed::setSexRatio_KeepSexRatio;
      getRandomSex_func_ptr = &LCE_Breed::getRandomSex_KeepSexRatio;
		}
		else{
		  _sex_ratio /= _sex_ratio+1;                            // now:   males/(males+females)
      setSexRatio_func_ptr  = &LCE_Breed::setSexRatio_NoSelfing;
      getRandomSex_func_ptr = &LCE_Breed::getRandomSex_NoSelfing;
    }
  }

	// if the sex is determined by the phenotype of the first quantitative trait
	if(this->get_parameter("sex_ratio_threshold")->isSet()){
    _threshold = this->get_parameter_value("sex_ratio_threshold");

		// we need two sexes
		if(isMatingPossible_func_ptr != &LCE_Breed::isMatingPossible_2_sex) fatal("Two sexes are requiered when using a sex threshold!\n");

    // the phenotype of the first trait has to be set just after the birth to set the sex
		if(!_pSelection) fatal("Sex-chromosome: At least a quantitative trait is needed to specify the sex!\n");
		if(_pSelection->get_selTrait(0)->get_Ve_prop()) fatal("Sex-chromosome: only the natal environment can be taken into account!\n");
  }
  else _threshold = my_NAN;       // not used, simulations are normal


  return true;
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed
// ----------------------------------------------------------------------------------------
/** This function is used for selection at the patch level.
  * First there is random mating depending on the fecundity of the females.
	* In a second step the popoulation size is downregulated depending on the fitness of the offspring,
	* but the number of indivduals depends ont he number of adults.
	*/
void LCE_Breed::breed_selection_offspring_patch ()
{
  #ifdef _DEBUG
	 message("(patch (selection on offspring)) ... ");
	#endif

	// get the number of survivors
	unsigned int nbBaby, nbSons, nbDaughters, nbPairs, nbPatch = _popPtr->getPatchNbr();
	Patch* cur_patch;
	for(unsigned int home = 0; home < nbPatch; ++home){
		cur_patch = _popPtr->getPatch(home);

		// check if the mating requirements are met
		if(!(this->*isMatingPossible_func_ptr)(cur_patch)) continue;

		// create the mating pairs
		if(_mating_system == 4) create_mating_pairs(cur_patch, nbPairs);

		// create the offspring with neutral mating
		nbBaby = my_round(_mean_fecundity*_nbIndividuals[FEM]);
		(this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
		createOffspring(cur_patch, nbDaughters, nbSons);

		// check if the pop size has to be regulated
		nbBaby = (this->*setNbOffspring_func_ptr)(_nbIndividuals[MAL], _nbIndividuals[FEM], cur_patch->get_K());
		if(nbBaby<cur_patch->size(OFFSx)){
			(this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
			_pSelection->set_fitness(cur_patch, OFFSx); // compute the fitnesses, but don't sort or make yet the array cumulative
			cur_patch->regulate_selection_fitness(nbSons,      _pSelection, MAL, OFFSx);  // regulate the males
			cur_patch->regulate_selection_fitness(nbDaughters, _pSelection, FEM, OFFSx);  // regulate the females
		}
	}
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed_selection_metapop
// ----------------------------------------------------------------------------------------
/** This function is used for selection at the metapopulation level. First the
	* fitness of all individuals is computed. The number of total offspring (of the metapop)
	* is distributed among the patches according to the mean fitness.
	* It does not make sence to use this function for a neutral case.
	*/
void LCE_Breed::breed_selection_offspring_metapop ()
{
  #ifdef _DEBUG
   message("(metapop (selection on offspring)) ... ");
	#endif

  Patch* cur_patch;
	unsigned int home, nbBaby, i, nbPatches = _popPtr->getPatchNbr();
	unsigned int nbDaughters, nbSons, nbOff, nbPairs;

  // array to buffer the fitness arrays for each patch (female and male)
	double totFitness = 0;                             // total fitness
	double* sumFitness = new double[nbPatches];        // sum of the fitness of each patch
	unsigned nbInd[2];                                 // total number of adult females and males
  TPatchFitness** fitnessArrays[2];
  for(i=0; i<2; ++i){
    fitnessArrays[i] = new TPatchFitness*[nbPatches];
    nbInd[i] = 0;
  }

	// create the offspring with neutral mating and store their fitnesses
	for(unsigned int home = 0; home < nbPatches; home++) {
		cur_patch = _popPtr->getPatch(home);

		// check if the mating requirements are met
		if(!(this->*isMatingPossible_func_ptr)(cur_patch)){
			fitnessArrays[MAL][home] = NULL;
			fitnessArrays[FEM][home] = NULL;
			continue;
		}

		// create the mating pairs
		if(_mating_system == 4) create_mating_pairs(cur_patch, nbPairs);

		// create the offspring with neutral mating
		nbBaby = my_round(_mean_fecundity*_nbIndividuals[FEM]);
		(this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
		createOffspring(cur_patch, nbDaughters, nbSons);
		nbInd[MAL] += _nbIndividuals[MAL];      // get total number of male adults
		nbInd[FEM] += _nbIndividuals[FEM];      // get total number of female adults

		// compute and store the OFFSPRING fitness of the current patch (0: male, 1: female)
		_pSelection->set_fitness(cur_patch, OFFSx);  // do not sort or make the array cumulative
		sumFitness[home]   = _pSelection->getSumFitness(MAL) + _pSelection->getSumFitness(FEM);
		totFitness         += sumFitness[home];
		if(nbSons)      fitnessArrays[MAL][home] = _pSelection->get_FitnessObject(MAL, true);
		if(nbDaughters) fitnessArrays[FEM][home] = _pSelection->get_FitnessObject(FEM, true);
	}

	// compute the total number of offspring of the metapopulation (taking the adults as reference...)
	nbOff = (this->*setNbOffspring_func_ptr)(nbInd[MAL], nbInd[FEM], _popPtr->get_total_carrying_capacity());

	// population size regulation depending on the fitness of the offspring?
	for(home = 0; home < nbPatches; ++home) {
		if(!fitnessArrays[FEM][home]) continue;  // female missing -> no offspring
		cur_patch = _popPtr->getPatch(home);

		// get the number of sons/daughters to produce of this patch (in relation to the adults)
		nbBaby = my_round(sumFitness[home]/totFitness*nbOff);
		(this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, cur_patch->size(MAL, ADLTx), cur_patch->size(FEM, ADLTx));

		// perform population size regulation for the males
		if(nbSons<cur_patch->size(MAL, OFFSx)){              	// are there too many sons?
			_pSelection->set_FitnessObject(MAL, fitnessArrays[MAL][home]);  // set the fitness arrays of the current patch
			cur_patch->regulate_selection_fitness(nbSons, _pSelection, MAL, OFFSx);
		}

		// perform population size regulation for the females
		if(nbDaughters<cur_patch->size(FEM, OFFSx)){        	// are there too many daughters?
			_pSelection->set_FitnessObject(FEM, fitnessArrays[FEM][home]);  // set the fitness arrays of the current patch
			cur_patch->regulate_selection_fitness(nbDaughters, _pSelection, FEM, OFFSx);
		}
	}//end_for_nbPatch

	// delete the fitness arrays
  if(sumFitness)      delete[] sumFitness;
  for(i=0; i<2; ++i){
    if(fitnessArrays[i]) delete[] fitnessArrays[i];
  }
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed
// ----------------------------------------------------------------------------------------
/** This function is used for hard selection. The fitness is assumed to absolut
  * This function can also be used for neutral mating.
  */
void LCE_Breed::breed_selection_offspring_hard ()
{
  #ifdef _DEBUG
   message("(hard (selection on offspring)) ... ");
  #endif

  Patch* cur_patch;
	unsigned int nbBaby, nbSons, nbDaughters, nbPairs;
  double meanFitness;

	// for each patch
  for(unsigned int home = 0; home < _popPtr->getPatchNbr(); home++) {
    cur_patch = _popPtr->getPatch(home);

    // check if the mating requirements are met
    if(!(this->*isMatingPossible_func_ptr)(cur_patch)) continue;

    // create the mating pairs
    if(_mating_system == 4) create_mating_pairs(cur_patch, nbPairs);

		// create the offspring with neutral mating depending on the fecundity
		nbBaby = my_round(_mean_fecundity*_nbIndividuals[FEM]);
    (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);
    createOffspring(cur_patch, nbDaughters, nbSons);

		// compute the mean fitness of the current patch (offspring, 0: male, 1: female)
		_pSelection->set_fitness(cur_patch, OFFSx); 	// do not sort or make the array cumulative
		meanFitness = (_pSelection->getSumFitness(MAL) + _pSelection->getSumFitness(FEM))/(nbSons + nbDaughters);

		// get the number of sons/daughters to produce of this patch (in relation to adults)
		nbBaby = my_round(meanFitness * (this->*setNbOffspring_func_ptr)(_nbIndividuals[MAL], _nbIndividuals[FEM], cur_patch->get_K()));
		(this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);

		// regulate pop size
		cur_patch->regulate_selection_fitness(nbSons,      _pSelection, MAL, OFFSx);
		cur_patch->regulate_selection_fitness(nbDaughters, _pSelection, FEM, OFFSx);
	}//end_for_nbPatch
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed
// ----------------------------------------------------------------------------------------
/** This function is used for selection at the patch level. The other patches are
	* not taken into account at any time.
	* This function can also be used for neutral mating.
  */
void LCE_Breed::breed_selection_neutral ()
{
  #ifdef _DEBUG
   message("(neutral) ... ");
  #endif

  Patch* cur_patch;
  unsigned int nbBaby, nbSons, nbDaughters, nbPairs;

  // for each patch
  for(unsigned int home = 0; home < _popPtr->getPatchNbr(); home++) {
    cur_patch = _popPtr->getPatch(home);

    // check if the mating requirements are met
    if(!(this->*isMatingPossible_func_ptr)(cur_patch)) continue;

    // create the mating pairs
    if(_mating_system == 4) create_mating_pairs(cur_patch, nbPairs);

    // get the number of sons/daughters to produce of this patch
    nbBaby = (this->*setNbOffspring_func_ptr)(_nbIndividuals[MAL], _nbIndividuals[FEM], cur_patch->get_K());
    (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);

    createOffspring(cur_patch, nbDaughters, nbSons);
  }//end_for_nbPatch
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed
// ----------------------------------------------------------------------------------------
/** This function is used for selection at the patch level. The other patches are
  * not taken into account at any time.
  * This function can also be used for neutral mating.
  */
void LCE_Breed::breed_selection_patch ()
{
  #ifdef _DEBUG
   message("(patch) ... ");
	#endif

	assert(_pSelection);

  Patch* cur_patch;
  unsigned int nbBaby, nbSons, nbDaughters, nbPairs;

	// for each patch
  for(unsigned int home = 0; home < _popPtr->getPatchNbr(); home++) {
    cur_patch = _popPtr->getPatch(home);

    // check if the mating requirements are met
    if(!(this->*isMatingPossible_func_ptr)(cur_patch)) continue;

    // create the mating pairs
    if(_mating_system == 4) create_mating_pairs(cur_patch, nbPairs);

		// compute the fitness of the current patch (0: male, 1: female)
		_pSelection->set_fitness(cur_patch, ADLTx, _sort, _mating_males);

    // get the number of sons/daughters to produce of this patch
		nbBaby = (this->*setNbOffspring_func_ptr)(_nbIndividuals[MAL], _nbIndividuals[FEM], cur_patch->get_K());
    (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);

    createOffspring(cur_patch, nbDaughters, nbSons);
  }//end_for_nbPatch
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed_selection_metapop
// ----------------------------------------------------------------------------------------
/** This function is used for selection at the metapopulation level. First the
  * fitness of all individuals is computed. The number of total offspring (of the metapop)
  * is distributed among the patches according to the mean fitness.
  * It does not make sence to use this function for a neutral case.
  */
void LCE_Breed::breed_selection_metapop ()
{
  #ifdef _DEBUG
   message("(metapop) ... ");
  #endif

	assert(_pSelection);

	Patch* cur_patch;
  unsigned int home, nbBaby, i, nbPatches = _popPtr->getPatchNbr();
  unsigned int nbDaughters, nbSons, nbOff, nbPairs;

  // array to buffer the fitness arrays for each patch (female and male)
  double totFitness = 0;                              // total fitness
  double* sumFitness = new double[nbPatches];        // mean fitness of each patch
  unsigned nbInd[2];
  TPatchFitness** fitnessArrays[2];
  for(i=0; i<2; ++i){
    fitnessArrays[i] = new TPatchFitness*[nbPatches];
    nbInd[i] = 0;
  }

  // for each patch compute the fitness of the individuals
  for(home = 0; home < nbPatches; ++home) {
    cur_patch = _popPtr->getPatch(home);

    // check if the mating requirements are met
    if(!(this->*isMatingPossible_func_ptr)(cur_patch)){
      fitnessArrays[MAL][home] = NULL;
      fitnessArrays[FEM][home] = NULL;
      continue;
    }

    // compute and store the fitness of the current patch (0: male, 1: female)
    _pSelection->set_fitness(cur_patch, ADLTx, _sort, _mating_males);
		nbInd[MAL] += _nbIndividuals[MAL];      // get total number of males
		nbInd[FEM] += _nbIndividuals[FEM];      // get total number of females

		// sum of the fitnesses and store the fitness arrays
		sumFitness[home]   = _pSelection->getSumFitness(MAL) + _pSelection->getSumFitness(FEM);
		totFitness          += sumFitness[home];
		if(_nbIndividuals[MAL]) fitnessArrays[MAL][home] = _pSelection->get_FitnessObject(MAL, true);
		if(_nbIndividuals[FEM]) fitnessArrays[FEM][home] = _pSelection->get_FitnessObject(FEM, true);
	}

  // compute the total number of offspring
  nbOff = (this->*setNbOffspring_func_ptr)(nbInd[MAL], nbInd[FEM], _popPtr->get_total_carrying_capacity());

  // create the offspring for each patch separately
  for(home = 0; home < nbPatches; ++home) {
    if(!fitnessArrays[MAL][home]) continue;  // female or male missing -> no offspring
    cur_patch = _popPtr->getPatch(home);

    // (re)set the patch parameters
    _nbIndividuals[MAL] = cur_patch->size(MAL, ADLTx);                // reset number of males
    _nbIndividuals[FEM] = cur_patch->size(FEM, ADLTx);                // reset number of females
		if(_nbIndividuals[MAL]) _pSelection->set_FitnessObject(MAL, fitnessArrays[MAL][home]);    // set the fitness arrays of the current patch
		if(_nbIndividuals[FEM]) _pSelection->set_FitnessObject(FEM, fitnessArrays[FEM][home]);    // set the fitness arrays of the current patch

    // create the mating pairs
    if(_mating_system == 4) create_mating_pairs(cur_patch, nbPairs);

    // get the number of sons/daughters to produce of this patch
    nbBaby = my_round(sumFitness[home]/totFitness*nbOff);
    (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);

    createOffspring(cur_patch, nbDaughters, nbSons);
  }//end_for_nbPatch

  // delete the fitness arrays
  if(sumFitness)      delete[] sumFitness;
  for(i=0; i<2; ++i){
    if(fitnessArrays[i]) delete[] fitnessArrays[i];
  }
}

// ----------------------------------------------------------------------------------------
// LCE_Breed::breed
// ----------------------------------------------------------------------------------------
/** This function is used for hard selection. The fitness is assumed to absolut
  * This function can also be used for neutral mating.
  */
void LCE_Breed::breed_selection_hard ()
{
  #ifdef _DEBUG
   message("(hard) ... ");
  #endif

  Patch* cur_patch;
  unsigned int nbBaby, nbSons, nbDaughters, nbPairs;
  double meanFitness;

  // for each patch
  for(unsigned int home = 0; home < _popPtr->getPatchNbr(); home++) {
    cur_patch = _popPtr->getPatch(home);

    // check if the mating requirements are met
    if(!(this->*isMatingPossible_func_ptr)(cur_patch)) continue;

    // create the mating pairs
    if(_mating_system == 4) create_mating_pairs(cur_patch, nbPairs);

    // compute the fitness of the current patch (0: male, 1: female)
    if(_pSelection) _pSelection->set_fitness(cur_patch, ADLTx, _sort, _mating_males);
    meanFitness = (_pSelection->getSumFitness(MAL) + _pSelection->getSumFitness(FEM))/
                              (_nbIndividuals[MAL] + _nbIndividuals[FEM]);

    // get the number of sons/daughters to produce of this patch
    nbBaby = my_round(meanFitness * (this->*setNbOffspring_func_ptr)
                      (_nbIndividuals[MAL], _nbIndividuals[FEM], cur_patch->get_K()));
    (this->*setSexRatio_func_ptr)(nbBaby, nbSons, nbDaughters, _nbIndividuals[MAL], _nbIndividuals[FEM]);

    createOffspring(cur_patch, nbDaughters, nbSons);
  }//end_for_nbPatch
}

// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
/** create the daugthers and sons */
void
LCE_Breed::createOffspring(Patch* cur_patch, unsigned int nbDaughters, unsigned int nbSons){
  Individual *MotherPtr, *FatherPtr;
	unsigned int index;

  // daughters: mate randomly a female and a male following their fitness
  while(nbDaughters>0) {
    MotherPtr = getMotherPtr(cur_patch, index);        // get a random male and female (according to their fitness)
    FatherPtr = getFatherPtr(cur_patch, index);
    _popPtr->makeOffsprg(MotherPtr, FatherPtr, FEM, cur_patch); // create the baby
    --nbDaughters;
  }//end_for_nbDaughters

  // sons: mate randomly a female and a male following their fitness
  while(nbSons>0) {
    MotherPtr = getMotherPtr(cur_patch, index);        // get a random male and female (according to their fitness)
    FatherPtr = getFatherPtr(cur_patch, index);
    _popPtr->makeOffsprg(MotherPtr, FatherPtr, MAL, cur_patch);   // create the baby
    --nbSons;
  }//end_for_nbSons
}

// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
/** Create mating pairs for the monogamy mating system.
  * The number of paris corresponds to the number of females.
  * if(nbFemales > nbMales) males may mate with several females
  * if(nbFemales < nbMales) not all males may mate
  */
void
LCE_Breed::create_mating_pairs(Patch* cur_patch, unsigned int& nbPairs){
  int i, pos, nbMales;

  nbPairs = cur_patch->size(FEM, ADLTx);    //number of females/ number of mating paris
  if(nbPairs > _aMatingPairs_size){
    _aMatingPairs_size = nbPairs;
    if(_aMatingPairs[MAL]) delete[] _aMatingPairs[MAL];
    if(_aMatingPairs[FEM]) delete[] _aMatingPairs[FEM];
    _aMatingPairs[MAL] = new Individual*[_aMatingPairs_size];
    _aMatingPairs[FEM] = new Individual*[_aMatingPairs_size];
  }

  nbMales = cur_patch->size(MAL, ADLTx);    // number of males

  // while we do not have made all the pairs
  while(nbPairs>0){

    // put the males into a vector
    vector<Individual*> vecMales;
    for(i = 0; i<nbMales; ++i){
      vecMales.push_back(cur_patch->get(MAL, ADLTx, i));
    }

    // create the pairs (choose randomly a male for each female)
    for(i=nbMales-1; i>=0 && nbPairs>0; --i){
      --nbPairs;
      _aMatingPairs[FEM][nbPairs] = cur_patch->get(FEM, ADLTx, nbPairs); // get the female
      pos = SimRunner::r.Uniform(i);                                     // select randomly a male
      _aMatingPairs[MAL][nbPairs] = vecMales[pos];                       // put this male to the pair
      vecMales.erase(vecMales.begin()+pos);                       // remove the male from the vector
    }
  }
}


// ----------------------------------------------------------------------------------------
// LCE_Breed
// ----------------------------------------------------------------------------------------
/** This function to set the sexes according to the phenotype of the first trait.
	* female <= threshold < male
  */
void
LCE_Breed::reset_sex_after_phentoype(age_idx AGE)
{
  Patch* cur_patch;
  Individual* ind;
  int i, sizeF, sizeM;
	unsigned int p;

  for(p = 0; p < _popPtr->getPatchNbr(); p++) {            // for each patch
    cur_patch = _popPtr->getPatch(p);
    sizeF = (int)cur_patch->size(FEM, AGE);
    sizeM = (int)cur_patch->size(MAL, AGE);

    // check all females
    for (i = sizeF-1; i >= 0; --i) {
      ind = cur_patch->get(FEM, AGE, i);
			_pSelection->get_selTrait(0)->set_phenotype(ind);  // set the phenotype
      if(_pSelection->get_selTrait(0)->get_phenotype() > _threshold){    // above threshold it should be a male
        ind->switch_sex(AGE,i);
      }
    }

    // check all males
    for (i = sizeM-1; i >= 0; --i) {      // note that here sizeM is not anymore cur_patch->size(MAL, ADLTx)
      ind = cur_patch->get(MAL, AGE, i);
      _pSelection->get_selTrait(0)->set_phenotype(ind);  // set the phenotype
      if(_pSelection->get_selTrait(0)->get_phenotype() <= _threshold){    // above threshold it should be a male
        ind->switch_sex(AGE,i);
      }
    }
  }
}


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

