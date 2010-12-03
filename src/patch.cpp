/** @file patch.cpp
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

#include "metapop.h"
#include "random.h"
#include "output.h"
#include "tselectiontrait.h"
#include "ttquanti.h"
#include "tselection.h"


unsigned int Patch::_nbPatch = 0;
Metapop* Patch::_popPtr=NULL;

using namespace std;

// ----------------------------------------------------------------------------------------
// reset_LinkedTraitsParams
// ----------------------------------------------------------------------------------------
/** creates the space for parameters linked to the quantitative traits */
void
Patch::reset_LinkedTraits()
{
  if(_pSelectionType) {delete[] _pSelectionType; _pSelectionType = NULL;}
  ARRAY::delete_3D(_localSelection, 2, _nbLinkedTraits);
  _nbLinkedTraits = 0;

  for(int i=0; i<2; ++i){      // for each sex
    if(_meanVe[i]){delete[] _meanVe[i]; _meanVe[i]= NULL;}
    if(_sdVe[i])  {delete[] _sdVe[i];   _sdVe[i]= NULL;}
    if(_h2[i])    {delete[] _h2[i];     _h2[i]= NULL;}
  }
}

// ----------------------------------------------------------------------------------------
// set_LinkedTraits
// ----------------------------------------------------------------------------------------
/** creates the space for parameters linked to the quantitative traits */
void
Patch::set_LinkedTraits(const unsigned int& size, TSelectionTrait** type, bool bothSexes)
{

  reset_LinkedTraits();

  _nbLinkedTraits = size;
  ARRAY::create_2D(_localSelection, 2, _nbLinkedTraits);
  _pSelectionType = new TSelectionTrait*[_nbLinkedTraits];
  for(unsigned int i=0; i<_nbLinkedTraits; ++i){
    _pSelectionType[i] = type[i];
    if(bothSexes) _localSelection[MAL][i] = new double[type[i]->get_nb_selection_params()];
    else          _localSelection[MAL][i] = NULL;
    _localSelection[FEM][i] = new double[type[i]->get_nb_selection_params()];
  }

  if(bothSexes){
    _meanVe[MAL] = new double [_nbLinkedTraits];
    _sdVe[MAL]   = new double [_nbLinkedTraits];
    _h2[MAL]     = new double [_nbLinkedTraits];
  }
  _meanVe[FEM] = new double [_nbLinkedTraits];
  _sdVe[FEM]   = new double [_nbLinkedTraits];
  _h2[FEM]     = new double [_nbLinkedTraits];
}

// ----------------------------------------------------------------------------------------
// set_Ve
// ----------------------------------------------------------------------------------------
/** transfers the input _h2 to _sdVe depending on the selected model
  * computes the Ve for the specified trait and both sexes
  * returns true if Ve is larger than zero for males and/or females
  */
bool
Patch::set_Ve(const int& model, const int& trait, sex_t sex){
  // for each quantitative trait
	double h2 = _h2[sex][trait];

  // is Ve directly set?
  if(model==0){
    if(h2<0) fatal("The environmental variance has to be positive (%f)!\n", h2);
    _sdVe[sex][trait] = sqrt(h2);    // input is the variance, here we need the sd!!!
  }

  // if heritability is 1 => Ve=0
  else if(h2==1){ // fully heritable
    _sdVe[sex][trait] = 0;
  }

  // Ve is set by the heritability
  else{
    if(h2<=0 || h2>1) fatal("The heritability is out of range %f!\n", h2); // check the range of h2

    // check if the pop is populated in case where Ve is only set at the beginning of the sim
    if((model==1 || model==3) && !size(ADULTS)){
      fatal("The environmental variance cannot be computed since all patches have to be populated at the start of a simulation (true for quanti_environmental_model set to 1 or 3)!\n");
		}

		// narrow or broad-sense heritability?
		double var, mean;
		TTQuantiSH* pStats = _pSelectionType[trait]->get_pQuantiProto()->_stats;
    if(model==1 || model==2){ // narrow-sense heritability
      unsigned int nb_locus  = _pSelectionType[trait]->get_pQuantiProto()->get_nb_locus();
      unsigned int nb_allele = _pSelectionType[trait]->get_pQuantiProto()->get_nb_allele();
			double** freqs = ARRAY::new_2D<double>(nb_locus, nb_allele);
      unsigned int** counts = ARRAY::new_2D<unsigned int>(nb_locus, nb_allele);
			pStats->set_alleleCountTable_local_ofPatch(this, ADULTS, freqs, counts);
			pStats->get_Va_ofPatch(this, ADULTS, mean, var, freqs, counts);  // compute Va
      if(var == my_NAN){      // if it could not be computed -> use the method for random mating to compute Va
        warning("Ve of patch %i could not be computed since Va could not be estimated -> using Va computed for random mating to set Ve (see manual parameter 'quanti_va_model')!\n",_ID+1);
        pStats->get_Va_ofPatch_random_mating(this, ADULTS, mean, var, freqs, counts);
      }
      ARRAY::delete_2D(freqs,  nb_locus);
      ARRAY::delete_2D(counts, nb_locus);
    }
    else{   // broad-sense heritability
			pStats->setMeanAndVar_Vg_ofPatch(this, ADULTS, mean, var); // compute Vg
    }

    _sdVe[sex][trait] = var * (1-h2)/h2;    // set the Ve
  }

  return (_sdVe[sex][trait] != 0);  // return true if sdVe is not zero
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void
Patch::set_PopSizes_ini(unsigned int nbfem, unsigned int nbmal){
  _NFem_ini = nbfem;
  _NMal_ini = nbmal;
  _N_ini    = nbfem + nbmal;
  _isExtinct = (_N_ini==0);
}

// ----------------------------------------------------------------------------------------
// set_PopSizes_ini_carrying_capacity
// ----------------------------------------------------------------------------------------
void
Patch::set_PopSizes_ini_carrying_capacity(){
  _NFem_ini = _KFem;
  _NMal_ini = _KMal;
  _N_ini    = _K;
  _isExtinct = (_N_ini==0);
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
Patch*
Patch::init(unsigned int id){
  _ID = id;
  _IDstr = STRING::int2str(id+1);
  reset_counters();

  for(unsigned int i=0; i < NB_AGE_CLASSES; i++) {
    _containers[MAL][i].assign( _KMal, 0 );
    _containers[FEM][i].assign( _KFem, 0 );
    _sizes[MAL][i] = 0;
    _sizes[FEM][i] = 0;
    _capacities[MAL][i] = _KMal;
    _capacities[FEM][i] = _KFem;
  }

  return this;
}

// ----------------------------------------------------------------------------------------
// set_carrying_capacities
// ----------------------------------------------------------------------------------------
void
Patch::set_carrying_capacities(const unsigned int& nbfem, const unsigned int& nbmal){
  _KFem = nbfem;
  _KMal = nbmal;
  _K = _KFem + _KMal;
}

// ----------------------------------------------------------------------------------------
// environmental variance
// ----------------------------------------------------------------------------------------
void Patch::set_localMeanVe (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    _meanVe[SEX][i] = array[i];
  }
}

void Patch::set_localh2Ve (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    _h2[SEX][i] = array[i];
  }
}

// ----------------------------------------------------------------------------------------
// stabilizing selection pressure
// ----------------------------------------------------------------------------------------
void Patch::set_localOptima (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    if(_pSelectionType[i]->get_nb_selection_params() == 2){
      _localSelection[SEX][i][0] = array[i];     // 1. pos
    }
  }
}

void Patch::set_localIntensity (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    if(_pSelectionType[i]->get_nb_selection_params() == 2){
      _localSelection[SEX][i][1] = array[i];     // 2. pos
    }
  }
}

// ----------------------------------------------------------------------------------------
// directional selection pressure
// ----------------------------------------------------------------------------------------
void Patch::set_localMin (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    if(_pSelectionType[i]->get_nb_selection_params() == 5){
      _localSelection[SEX][i][0] = array[i];     // 1. pos
    }
  }
}

void Patch::set_localMax (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    if(_pSelectionType[i]->get_nb_selection_params() == 5){
      _localSelection[SEX][i][1] = array[i];     // 2. pos
    }
  }
}

void Patch::set_localMaxGrowth (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    if(_pSelectionType[i]->get_nb_selection_params() == 5){
			_localSelection[SEX][i][2] = array[i];     // 4. pos
		}
  }
}

void Patch::set_localGrowthRate (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    if(_pSelectionType[i]->get_nb_selection_params() == 5){
      _localSelection[SEX][i][3] = array[i];     // 3. pos
    }
  }
}

void Patch::set_localSymmetry (double *array, sex_t SEX)
{
  for(unsigned int i = 0; i<_nbLinkedTraits; i++) {
    if(_pSelectionType[i]->get_nb_selection_params() == 5){
      _localSelection[SEX][i][4] = array[i];     // 5. pos
    }
  }
}

// ----------------------------------------------------------------------------------------
// set_localParameter
// ----------------------------------------------------------------------------------------
void Patch::set_localParameter(double* array, sex_t sex,
                                  void (Patch::*pt2Func)(double*, sex_t))
{
  (this->*pt2Func)(array, sex);
}

// ----------------------------------------------------------------------------------------
// reset_counters
// ----------------------------------------------------------------------------------------
void Patch::reset_counters()
{
  nbEmigrant = 0;
  nbImmigrant = 0;
  nbPhilopat = 0;
  nbKolonisers = 0;
  colnAge = 0;
}
// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
void Patch::setNewGeneration(age_t AGE)
{
  unsigned int mask = 1;

  for(unsigned int i = 0; i < NB_AGE_CLASSES; i++, mask<<=1) {
		if(mask & AGE) setNewGeneration(static_cast<age_idx>(i));
  }
}
// ----------------------------------------------------------------------------------------
// setNewGeneration
// ----------------------------------------------------------------------------------------
void Patch::setNewGeneration(age_idx AGE)
{
	Individual *new_ind;

	//if too much females in the Patch, flush them into the RecyclingPOOL
	if(size(FEM, AGE)) flush(FEM, AGE);
	for(unsigned int i = 0; i < _NFem_ini; i++) {
		new_ind = _popPtr->makeNewIndividual(0,0,FEM,this);
		new_ind->create();
		add(FEM, AGE, new_ind);
	}

	//males: same as for the females....
	if(size(MAL, AGE)) flush(MAL, AGE);
	for(unsigned int i = 0; i < _NMal_ini; i++) {
		new_ind = _popPtr->makeNewIndividual(0,0,MAL,this);
		new_ind->create();
		add(MAL, AGE, new_ind);
	}
}

// ----------------------------------------------------------------------------------------
// ~Patch
// ----------------------------------------------------------------------------------------
Patch::~Patch()
{
//#ifdef _DEBUG
//  message("Patch::~Patch\n");
//#endif

	for (unsigned int i = 0; i < 2; ++i){
		for(unsigned int j = 0; j < NB_AGE_CLASSES; ++j){
			for(unsigned int k = 0; k < _sizes[i][j] ; ++k){
				delete _containers[i][j][k];
			}
		}
	}

	ARRAY::delete_2D(_sizes, 2);
	ARRAY::delete_2D(_capacities, 2);
	ARRAY::delete_2D(_containers, 2);

	reset_LinkedTraits();
}

//---------------------------------------------------------------------------
#include "patch.h"
#include "metapop.h"
//---------------------------------------------------------------------------
///@}
/**Returns the size of the container for the appropriate sex and age classes present in the age flag.
	@param SEX the sex class
	@param AGE the flag value of the age class
	*/
unsigned int
Patch::size(sex_t SEX, age_t AGE){
	unsigned int mask = 1, s = 0;
	for(unsigned int i = 0; i < NB_AGE_CLASSES; i++, mask <<= 1) {
		if(mask & AGE) s += _sizes[SEX][i];
	}
	return s;
}

//---------------------------------------------------------------------------
/**Returns the size of the container of the appropriate age class(es) for both sexes.
	@param AGE the flag value of the age class
	*/
unsigned int
Patch::size(age_t AGE){
	return size(MAL,AGE) + size(FEM,AGE);
}

//---------------------------------------------------------------------------
/**Returns the size of the container for the appropriate sex and age class.
	@param SEX the sex class
	@param AGE the index of the age class
*/
unsigned int
Patch::size(sex_t SEX, age_idx AGE){
	return _sizes[SEX][AGE];
}

//---------------------------------------------------------------------------
/**Returns the size of the container for the appropriate age class for both sexes.
	@param AGE the index of the age class
	*/
unsigned int
Patch::size(age_idx AGE){
	return _sizes[0][AGE] + _sizes[1][AGE];
}

//---------------------------------------------------------------------------
/**Returns a pointer to the individual sitting at the index passed.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param at the index of the individual in the container
	*/
Individual*
Patch::get(sex_t SEX, age_idx AGE, unsigned int at){
	return _containers[SEX][AGE][at];
}

//---------------------------------------------------------------------------
/**Modifies the appropriate container with value of the pointer given.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param at the index of the individual in the container
	@param ind the pointer to the individual
	*/
void
Patch::set(sex_t SEX, age_idx AGE, unsigned int at, Individual* ind){
	_containers[SEX][AGE][at] = ind;
}

//---------------------------------------------------------------------------
/**Adds an individual to the appropriate container, increments its size, eventually resizing it.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param ind the pointer to the individual
	*/
void
Patch::add(sex_t SEX, age_idx AGE, Individual* ind){

	while(_sizes[SEX][AGE] + 1 > _capacities[SEX][AGE]) {
		if(_K){
			_containers[SEX][AGE].resize(_capacities[SEX][AGE] + _K);
			_capacities[SEX][AGE] += _K;
		}
		else{
			_containers[SEX][AGE].resize(_capacities[SEX][AGE] + 1);
			_capacities[SEX][AGE] += 1;
		}
	}

	_containers[SEX][AGE][_sizes[SEX][AGE]] = ind;
	++_sizes[SEX][AGE];
}

//---------------------------------------------------------------------------
/**Assigns a new container of given size for the sex and age class passed, sets all values to NULL.*/
void
Patch::assign(sex_t SEX, age_idx AGE, unsigned int n){
	_containers[SEX][AGE].assign(n,0);
	_sizes[SEX][AGE] = 0;
	_capacities[SEX][AGE] = n;
}

//---------------------------------------------------------------------------
/**Removes the individual sitting at the given index in the appropriate container.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param at the index of the individual in the container
	*/
void
Patch::remove(sex_t SEX, age_idx AGE, unsigned int at){
	assert(at<_sizes[SEX][AGE]);
	_containers[SEX][AGE][at] = _containers[SEX][AGE][_sizes[SEX][AGE] - 1];  // swap the last individual to the position at 'at'
	_sizes[SEX][AGE]--;
}

//---------------------------------------------------------------------------
/**Removes the individual in the appropriate container (function is not very efficient).
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param ind pointer of an individual to remove
	*/
void
Patch::remove(sex_t SEX, age_idx AGE, Individual* ind){
	for(unsigned int i=0; i<_sizes[SEX][AGE]; ++i){
		if(ind == _containers[SEX][AGE][i]){
			remove(SEX, AGE, i);
			return;
		}
	}
	assert(1!=1); 	// function should never pass here
}

//---------------------------------------------------------------------------
/**Removes and deletes the individual sitting at the given index in the appropriate container.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param at the index of the individual in the container
	*/
void
Patch::recycle(sex_t SEX, age_idx AGE, unsigned int at){
	assert(at<_sizes[SEX][AGE]);
	_popPtr->recycle(_containers[SEX][AGE][at]);
	_containers[SEX][AGE][at] = _containers[SEX][AGE][_sizes[SEX][AGE] - 1];  // swap the last individual to the position at 'at'
	_sizes[SEX][AGE]--;
}

//---------------------------------------------------------------------------
/**Removes and deletes the individual in the appropriate container (function is not very efficient).
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param ind pointer of an individual to remove
	*/
void
Patch::recycle(sex_t SEX, age_idx AGE, Individual* ind){
	for(unsigned int i=0; i<_sizes[SEX][AGE]; ++i){
		if(ind == _containers[SEX][AGE][i]){
			_popPtr->recycle(ind);
			remove(SEX, AGE, i);
			return;
		}
	}
	assert(1!=1); 	// function should never pass here
}

//---------------------------------------------------------------------------
// MOVE: does not overwrite the individuals in the to container
//---------------------------------------------------------------------------
/**Moves all individual from an age class to an other one.
	\b Note: both containers are transformed by this operation. The 'from'
					 container size is set to zero while the 'to' container size
					 is increased by the number of the 'from' container.
	@param from the original age class of the individual
	@param to the destination age class of the individual
	*/
void
Patch::move(age_idx from, age_idx to){
	move(FEM, from, to);
	move(MAL, from, to);
}

//---------------------------------------------------------------------------
/**Moves all individual from an age class to another one.
	\b Note: both containers are transformed by this operation. The 'from'
					 container size is set to zero while the 'to' container size
					 is increased by the number of the 'from' container.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
	*/
void
Patch::move(sex_t SEX, age_idx from, age_idx to){
	for(int i=0, size=_sizes[SEX][from]; i<size; ++i){
		move(SEX, from, to, 0);
	}
}

//---------------------------------------------------------------------------
/**Moves an individual from an age class to an other one.
	\b Note: both containers are transformed by this operation. The 'from'
					 container size is reduced by one while the 'to' container size
					 is increased by one.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
	@param at the index of the individual in the container
	*/
void
Patch::move(sex_t SEX, age_idx from, age_idx to, unsigned int at){
	add(SEX, to, _containers[SEX][from][at]);
	remove(SEX, from, at);
}

//---------------------------------------------------------------------------
// SWAP: overwrites the individuals in the to container, thus the to container has to be empty
//---------------------------------------------------------------------------
/**Copies all elements in the 'from' age-class container to the 'to' age-class container.
	The previous elements of the 'to' container are overwritten, it is thus worth considering using flush()
	before swaping to avoid memory leaks! The size of the 'from' container is set to 0.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
*/
void
Patch::swap(age_idx from, age_idx to){
	swap(FEM, from, to);
	swap(MAL, from, to);
}

//---------------------------------------------------------------------------
/**Copies all elements in the 'from' age-class container to the 'to' age-class container of the same sex.
	The previous elements of the 'to' container are overwritten, it is thus worth considering using flush()
	before swaping to avoid memory leaks! The size of the 'from' container is set to 0.
	@param SEX the sex class of the individual
	@param from the original age class of the individual
	@param to the destination age class of the individual
*/
void
Patch::swap(sex_t SEX, age_idx from, age_idx to){
	assert(!_sizes[SEX][to]);

	//store values of the "to" container temporarily
	unsigned int       temp_capacity  = _capacities[SEX][to];
	deque<Individual*> temp_container = _containers[SEX][to];

	// reassign the "from" container
	_containers[SEX][to]    = _containers[SEX][from];
	_sizes[SEX][to]         = _sizes[SEX][from];
	_capacities[SEX][to]    = _capacities[SEX][from];

	// emtpy the "to" container
	_containers[SEX][from]  = temp_container;
	_capacities[SEX][from]  = temp_capacity;
	_sizes[SEX][from]       = 0;
}

//---------------------------------------------------------------------------
/**Sets the size of the appropriate container to zero.
	\b Note: no memory operation is performed, the capacity of the container is thus not affected.
					The individual pointers are not flushed to the recycling pool. It is thus a good idea
					to consider using Patch::flush to be sure no pointers remained in the container.
	@see flush()
	@param SEX the sex class
	@param AGE the index of the age class
	*/
void
Patch::clear(sex_t SEX, age_idx AGE){
	_sizes[SEX][AGE] = 0;
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
	\b Note: not memory operation is performed, the total amount of memory allocated is
	left untouched.
	@param SEX the sex class of the individual
	@param AGE the index of the age class
	@param pop the pointer to the metapop for access to the recycling pool
*/
void
Patch::flush(sex_t SEX, age_idx AGE){
	for (unsigned int i = 0; i < _sizes[SEX][AGE]; ++i) {
		_popPtr->recycle(_containers[SEX][AGE][i]);
	}
	_sizes[SEX][AGE] = 0;
}

//---------------------------------------------------------------------------
void
Patch::flush(age_idx AGE){
	flush(FEM, AGE);
	flush(MAL, AGE);
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of the appropriate sex and age class and flush them into the recycling pool.
	@param AGE an unsigned int containing the flags of the age classes to flush
	@param pop the pointer to the metapop for access to the recycling pool
	@see flush()
	*/
void
Patch::flush(age_t AGE){
	unsigned int mask = 1;

	for(unsigned int i = 0; i < NB_AGE_CLASSES; i++, mask <<= 1) {
		if(mask & AGE) {
			flush(MAL, static_cast<age_idx>(i));
			flush(FEM, static_cast<age_idx>(i));
		}
  }
}

//---------------------------------------------------------------------------
/**Removes all individual pointers of all sex and age classes and flush them into the recycling pool.*/
void
Patch::flush(){
  for(unsigned int i = 0; i < NB_AGE_CLASSES; i++) {
		flush(MAL, static_cast<age_idx>(i));
    flush(FEM, static_cast<age_idx>(i));
  }
}

//---------------------------------------------------------------------------
/** return an array with all phenotypes of the specific trait.
  * The array has to be deleted after use.
  */
double*
Patch::getPhenotypes(sex_t SEX, age_idx AGE, const int& trait){
  unsigned int i, size=_sizes[SEX][AGE];
  double* pheno = new double[size];
  for (i = 0; i < size; ++i) {
    pheno[i] = _containers[SEX][AGE][i]->getTraitPhenotype(trait);
  }
  return pheno;
}

//---------------------------------------------------------------------------
/** downregulates the pop size to the given K if necessary. The survivors are drawn
	* randomly in relation to their fitness
	*/
void
Patch::regulate_selection_fitness(const unsigned int& K, TSelection* pSel, sex_t SEX, age_idx AGE)
{
	// if K is zero: remove all individuals
	if(!K){
		flush(SEX, AGE);
		return;
	}

	unsigned int nbInd = pSel->get_nbInd(SEX);
	double* aFit       = pSel->get_aFit(SEX);
	Individual** aInd  = pSel->get_aInd(SEX);

	// remove frist all individuals with a fitness of 0.0
	for(int i=(int)nbInd-1; i>=0; --i){
		if(!aFit[i]) {
			recycle(SEX, AGE, aInd[i]);
			pSel->remove(SEX, i);
		}
	}

	nbInd = pSel->get_nbInd(SEX);
	if(nbInd>K){
		if((nbInd-K)>K/2) regulatation_selection_draw_survivors(pSel, K, SEX, AGE);
		else              regulatation_selection_draw_loosers(pSel, nbInd-K, SEX, AGE);
	}
}

//---------------------------------------------------------------------------
void
Patch::regulatation_selection_draw_survivors(TSelection* pSelection, const unsigned int& K, sex_t SEX, age_idx AGE)
{
	assert(K && K<=size(SEX, AGE));
	// make the array cumulative to draw the MOST fittest
	pSelection->sort_fitness(SEX, -3, K);
	Individual** aInd = pSelection->get_aInd(SEX);

	// remove the loosers (second part of the array)
	for(unsigned int i = K, N = size(SEX, AGE); i<N; ++i){
		recycle(SEX, AGE, aInd[i]);     // remove the "unsorted" individuals
	}
}

//---------------------------------------------------------------------------
void
Patch::regulatation_selection_draw_loosers(TSelection* pSelection, const unsigned int& K, sex_t SEX, age_idx AGE)
{

	assert(K && K<=size(SEX, AGE));
	// make the array cumulative to draw the MOST fittest
	pSelection->sort_fitness(SEX, 3, K);
	Individual** aInd = pSelection->get_aInd(SEX);

	// remove the loosers (first sorted part of the array)
	for(unsigned int i = 0; i<K; ++i){
		recycle(SEX, AGE, aInd[i]);     // remove the "sorted" individuals
	}
}

//---------------------------------------------------------------------------




