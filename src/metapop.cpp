/** @file metapop.cpp
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

#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <dirent.h>
#include "random.h"
#include "output.h"
#include "metapop.h"
#include "individual.h"
#include "lifecycleevent.h"
#include "simulation.h"
#include "lce_breed.h"
#include "lce_disperse.h"

#include "ttquanti.h"
#include "ttneutral.h"

using namespace std;

// ----------------------------------------------------------------------------------------
// constructor
// ----------------------------------------------------------------------------------------
Metapop::Metapop() : _statHandler(), _sample_pops(0), _patchNbr(0), _patchK(0),
    _patchKfem(0), _patchKmal(0), _generations(0), _replicates(0), _currentGeneration(1), _currentReplicate(1),
    _currentAge(NONE), _meanReplElapsedTime(0), _meanGenLength(0), _sexInitRatio(0.5), _service(0), _pSelection(0),
    _total_carrying_capacity(my_NAN), _density_threshold(0), _density_threshold_param(0){

#ifdef _DEBUG
        message(" Metapop::Metapop()\n");
#endif
        set_paramset("population", true);

        // carrying capacities
        _paramSet->add_param("patch_number",INT2,false,false,0,0,"1");
        _paramSet->add_param("patch_capacity",INT_MAT,false,false,0,0,"",true);
        _paramSet->add_param("patch_capacity_fem",INT_MAT,false,false,0,0,"",true);
        _paramSet->add_param("patch_capacity_mal",INT_MAT,false,false,0,0,"",true);

        // population inital sizes
        _paramSet->add_param("patch_ini_size",INT_MAT,false,false,0,0,"");
        _paramSet->add_param("patch_ini_size_fem",INT_MAT,false,false,0,0,"");
        _paramSet->add_param("patch_ini_size_mal",INT_MAT,false,false,0,0,"");

        // environmental influence
        _paramSet->add_param("patch_ve_mean",DBL_MAT,false,false,0,0,"0", true);
        _paramSet->add_param("patch_ve_mean_fem",DBL_MAT,false,false,0,0,"0", true);
        _paramSet->add_param("patch_ve_mean_mal",DBL_MAT,false,false,0,0,"0", true);
        _paramSet->add_param("patch_ve_var",DBL_MAT,false,false,0,0,"0", true);
        _paramSet->add_param("patch_ve_var_fem",DBL_MAT,false,false,0,0,"0", true);
        _paramSet->add_param("patch_ve_var_mal",DBL_MAT,false,false,0,0,"0", true);

        // stabilizing selection
        _paramSet->add_param("patch_stab_sel_optima",DBL_MAT,false,false,0,0,"0",true);
        _paramSet->add_param("patch_stab_sel_optima_fem",DBL_MAT,false,false,0,0,"0",true);
        _paramSet->add_param("patch_stab_sel_optima_mal",DBL_MAT,false,false,0,0,"0",true);
        _paramSet->add_param("patch_stab_sel_optima_var",DBL_MAT,false,false,0,0,"0");
        _paramSet->add_param("patch_stab_sel_intensity",DBL_MAT,false,false,0,0,"1",true);
        _paramSet->add_param("patch_stab_sel_intensity_fem",DBL_MAT,false,false,0,0,"1",true);
        _paramSet->add_param("patch_stab_sel_intensity_mal",DBL_MAT,false,false,0,0,"1",true);
        _paramSet->add_param("patch_stab_sel_intensity_var",DBL_MAT,false,false,0,0,"0");

        // directional selection
        _paramSet->add_param("patch_dir_sel_min",DBL_MAT,false,true,0,1,"0",true);
        _paramSet->add_param("patch_dir_sel_min_fem",DBL_MAT,false,true,0,1,"0",true);
        _paramSet->add_param("patch_dir_sel_min_mal",DBL_MAT,false,true,0,1,"0",true);
        _paramSet->add_param("patch_dir_sel_min_var",DBL_MAT,false,false,0,0,"0");
        _paramSet->add_param("patch_dir_sel_max",DBL_MAT,false,true,0,1,"1",true);
        _paramSet->add_param("patch_dir_sel_max_fem",DBL_MAT,false,true,0,1,"1",true);
        _paramSet->add_param("patch_dir_sel_max_mal",DBL_MAT,false,true,0,1,"1",true);
        _paramSet->add_param("patch_dir_sel_max_var",DBL_MAT,false,false,0,0,"0");
        _paramSet->add_param("patch_dir_sel_max_growth",DBL_MAT,false,false,0,0,"0",true);
        _paramSet->add_param("patch_dir_sel_max_growth_fem",DBL_MAT,false,false,0,0,"0",true);
        _paramSet->add_param("patch_dir_sel_max_growth_mal",DBL_MAT,false,false,0,0,"0",true);
        _paramSet->add_param("patch_dir_sel_max_growth_var",DBL_MAT,false,false,0,0,"0");
        _paramSet->add_param("patch_dir_sel_growth_rate",DBL_MAT,false,false,0,0,"1",true);
        _paramSet->add_param("patch_dir_sel_growth_rate_fem",DBL_MAT,false,false,0,0,"1",true);
        _paramSet->add_param("patch_dir_sel_growth_rate_mal",DBL_MAT,false,false,0,0,"1",true);
        _paramSet->add_param("patch_dir_sel_growth_rate_var",DBL_MAT,false,false,0,0,"0");
        _paramSet->add_param("patch_dir_sel_symmetry",DBL_MAT,false,true,0,1,"1",true);
        _paramSet->add_param("patch_dir_sel_symmetry_fem",DBL_MAT,false,true,0,1,"1",true);
        _paramSet->add_param("patch_dir_sel_symmetry_mal",DBL_MAT,false,true,0,1,"1",true);
        _paramSet->add_param("patch_dir_sel_symmetry_var",DBL_MAT,false,false,0,0,"0");

        _paramSet->add_param("sampled_patches",INT_MAT,false,false,0,0,"0");

        _paramSet->add_param("selection_position",INT2,false,true,0,5,"0");
        _paramSet->add_param("selection_level",INT2,false,true,0,2,"0");

        _paramSet->add_param("temporal_change_following_density",STR_MAT,false,false,0,0,"0");
    }

// ----------------------------------------------------------------------------------------
// destructor
// ----------------------------------------------------------------------------------------
Metapop::~Metapop()
{
    clear();
#ifdef _DEBUG
    message(" Metapop::~Metapop()\n");
#endif
}

// ----------------------------------------------------------------------------------------
// setCarryingCapacities
// ----------------------------------------------------------------------------------------
/** sets the carrying capacities
*/
void Metapop::setCarryingCapacities()
{
    if(!_sexInitRatio) setCarryingCapacities_1sex();
    else               setCarryingCapacities_2sexes();
}

// ----------------------------------------------------------------------------------------
// setCarryingCapacities_1sex
// ----------------------------------------------------------------------------------------
/** sets the carrying capacities of the populations for one sex
 * (mating systems selfing and cloning)
 */
void Metapop::setCarryingCapacities_1sex()
{
    _patchNmal = 0;   // we don't use males...

    // sex specific parameter
    if(_paramSet->isSet("patch_capacity_fem")) {
        // matrix
        if(_paramSet->isMatrix("patch_capacity_fem") ) {
            setCarryingCapacitiesOfPatches(FEM,_paramSet->getMatrix("patch_capacity_fem"));
        }
        // not matrix
        else {
            _patchN = _patchNfem = (unsigned int)_paramSet->getValue("patch_capacity_fem");
            setCarryingCapacitiesOfPatches();
        }
    }

    // general parameter
    else if(_paramSet->isSet("patch_capacity")) {
        // matrix
        if( _paramSet->isMatrix("patch_capacity") ) {
            setCarryingCapacitiesOfPatches(FEM, _paramSet->getMatrix("patch_capacity"));
        }
        // not matrix
        else {
            _patchN = _patchNfem = (unsigned int)_paramSet->getValue("patch_capacity");
            setCarryingCapacitiesOfPatches();
        }
    }
}

// ----------------------------------------------------------------------------------------
// setCarryingCapacities_2sexes
// ----------------------------------------------------------------------------------------
/** sets the carrying capacities of the populations for two sex
 * (all mating systems expect selfing and cloning)
 */
void Metapop::setCarryingCapacities_2sexes()
{
    // if sex specific parameters are passed
    if(_paramSet->isSet("patch_capacity_fem")){
        if(_paramSet->isSet("patch_capacity_mal")) {
            // both matrix
            if(_paramSet->isMatrix("patch_capacity_fem") && _paramSet->isMatrix("patch_capacity_mal")) {
                setCarryingCapacitiesOfPatches(_paramSet->getMatrix("patch_capacity_fem"),_paramSet->getMatrix("patch_capacity_mal"));
            }
            // only male matrix
            else if(!_paramSet->isMatrix("patch_capacity_fem") && _paramSet->isMatrix("patch_capacity_mal")) {
                _patchNfem = (unsigned int)_paramSet->getValue("patch_capacity_fem");
                setCarryingCapacitiesOfPatches(MAL,_paramSet->getMatrix("patch_capacity_mal"));
            }
            // only female matrix
            else if(_paramSet->isMatrix("patch_capacity_fem") && !_paramSet->isMatrix("patch_capacity_mal")) {
                _patchNmal = (unsigned int)_paramSet->getValue("patch_capacity_mal");
                setCarryingCapacitiesOfPatches(FEM,_paramSet->getMatrix("patch_capacity_fem"));
            }
            // no matrixes
            else{
                _patchNfem = (unsigned int)_paramSet->getValue("patch_capacity_fem");
                _patchNmal = (unsigned int)_paramSet->getValue("patch_capacity_mal");
                _patchN = _patchNfem + _patchNmal;
                setCarryingCapacitiesOfPatches();
            }
        }
        else fatal("Only one sex specific carring capacity is specified: both are required!\n");
    }

    // general parameter
    else if(_paramSet->isSet("patch_capacity")) {
        // matrix
        if( _paramSet->isMatrix("patch_capacity") ) {
            setCarryingCapacitiesOfPatches(_paramSet->getMatrix("patch_capacity"));
        }
        else{
            // split the total carrying capacity according to the sex ratio
            _patchN = (unsigned int)_paramSet->getValue("patch_capacity");
            _patchNfem = my_round(_sexInitRatio * _patchN);  // set the number of females
            _patchNmal = _patchN - _patchNfem;  // set the number of males
            setCarryingCapacitiesOfPatches();
        }
    }
}

// ----------------------------------------------------------------------------------------
// setPopulationSizes
// ----------------------------------------------------------------------------------------
/** set the inital population sizes
*/
void Metapop::setInitPopulationSizes()
{
    if(!_sexInitRatio) setInitPopulationSizes_1sex();
    else               setInitPopulationSizes_2sexes();
}

// ----------------------------------------------------------------------------------------
// setPopulationSizes_1sex
// ----------------------------------------------------------------------------------------
/** sets the population sizes for one sex
 * (mating systems selfing and cloning)
 */
void Metapop::setInitPopulationSizes_1sex()
{
    _patchNmal = 0;   // we don't use males...

    // sex specific parameter
    if(_paramSet->isSet("patch_ini_size_fem")) {
        // matrix
        if(_paramSet->isMatrix("patch_ini_size_fem") ) {
            setInitPopulationSizesOfPatches(FEM,_paramSet->getMatrix("patch_ini_size_fem"));
        }
        // not matrix
        else {
            _patchN = _patchNfem = (unsigned int)_paramSet->getValue("patch_ini_size_fem");
            setInitPopulationSizesOfPatches();
        }
    }

    // general parameter
    else if(_paramSet->isSet("patch_ini_size")) {
        // matrix
        if( _paramSet->isMatrix("patch_ini_size") ) {
            setInitPopulationSizesOfPatches(FEM, _paramSet->getMatrix("patch_ini_size"));
        }
        // not matrix
        else {
            _patchN = _patchNfem = (unsigned int)_paramSet->getValue("patch_ini_size");
            setInitPopulationSizesOfPatches();
        }
    }

    // instantaneous colonization
    else{
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_PopSizes_ini_carrying_capacity();
        }
    }
}

// ----------------------------------------------------------------------------------------
// setInitPopulationSizes_2sexes
// ----------------------------------------------------------------------------------------
/** sets the population sizes for two sex
 * (all mating systems expect selfing and cloning)
 */
void Metapop::setInitPopulationSizes_2sexes()
{
    // both sex specific parameters are passed
    if(_paramSet->isSet("patch_ini_size_fem")){
        if(_paramSet->isSet("patch_ini_size_mal")) {
            // both matrix
            if(_paramSet->isMatrix("patch_ini_size_fem") && _paramSet->isMatrix("patch_ini_size_mal")) {
                setInitPopulationSizesOfPatches(_paramSet->getMatrix("patch_ini_size_fem"),_paramSet->getMatrix("patch_ini_size_mal"));
            }
            // only male matrix
            else if(!_paramSet->isMatrix("patch_ini_size_fem") && _paramSet->isMatrix("patch_ini_size_mal")) {
                _patchNfem = (unsigned int)_paramSet->getValue("patch_ini_size_fem");
                setInitPopulationSizesOfPatches(MAL,_paramSet->getMatrix("patch_ini_size_mal"));
            }
            // only female matrix
            else if(_paramSet->isMatrix("patch_density_fem") && !_paramSet->isMatrix("patch_ini_size_mal")) {
                _patchNmal = (unsigned int)_paramSet->getValue("patch_ini_size_mal");
                setInitPopulationSizesOfPatches(FEM,_paramSet->getMatrix("patch_ini_size_fem"));
            }
            // no matrixes
            else{
                _patchNfem = (unsigned int)_paramSet->getValue("patch_ini_size_fem");
                _patchNmal = (unsigned int)_paramSet->getValue("patch_ini_size_mal");
                _patchN = _patchNfem + _patchNmal;
                setInitPopulationSizesOfPatches();
            }
        }
        fatal("Only one sex specific inital population size is specified: both are required!\n");
    }

    // general parameter
    else if(_paramSet->isSet("patch_ini_size")) {
        // matrix
        if( _paramSet->isMatrix("patch_ini_size") ) {
            setInitPopulationSizesOfPatches(_paramSet->getMatrix("patch_ini_size"));
        }
        else{
            _patchN    = (unsigned int)_paramSet->getValue("patch_ini_size");

            if(_paramSet->isSet("patch_ini_size_fem")) {
                _patchNfem = (unsigned int)_paramSet->getValue("patch_ini_size_fem");
                _patchNmal = _patchN - _patchNfem;
            }
            else if(_paramSet->isSet("patch_ini_size_mal")) {
                _patchNmal = (unsigned int)_paramSet->getValue("patch_ini_size_mal");
                _patchNfem = _patchN - _patchNmal;
            }
            else{
                // split the total population size according to the sex ratio
                _patchNfem = (unsigned int)my_round(_sexInitRatio * _patchN);  // set the number of females
                _patchNmal = _patchN - _patchNfem;  // set the number of males
            }
            setInitPopulationSizesOfPatches();
        }
    }

    // instantaneous colonization
    else{
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_PopSizes_ini_carrying_capacity();
        }
    }
}

// ----------------------------------------------------------------------------------------
// temporal_change
// ----------------------------------------------------------------------------------------
void
Metapop::temporal_change(const int& gen){
    map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);

    // if it is a temporal paramter
    if(pParam){
        // check if a change has to be made
        map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
        if(pMap){
            // iterate through the map and performe the updates
            bool optima = false, intensity = false, carCap = false,
                 growth_rate = false, max_growth = false, symmetry = false,
                 min = false, max = false, meanVe = false, h2 = false;
            map<string, Param*>::iterator pos = pMap->begin();
            for(; pos != pMap->end(); ++pos){
                if(pos->first.find("patch_capacity") != string::npos)            carCap =      true;

                else if(pos->first.find("patch_ve_mean") != string::npos)             meanVe =      true;
                else if(pos->first.find("patch_ve_var") != string::npos)              h2 =          true;

                else if(pos->first.find("patch_stab_sel_optima") != string::npos)     optima =      true;
                else if(pos->first.find("patch_stab_sel_intensity") != string::npos)  intensity =   true;

                else if(pos->first.find("patch_dir_sel_growth_rate") != string::npos) growth_rate = true;
                else if(pos->first.find("patch_dir_sel_max_growth") != string::npos)  max_growth =  true;
                else if(pos->first.find("patch_dir_sel_symmetry") != string::npos)    symmetry =    true;
                else if(pos->first.find("patch_dir_sel_min") != string::npos)         min =         true;
                else if(pos->first.find("patch_dir_sel_max") != string::npos)         max =         true;
            }

            // make the changes
            if(carCap){
                setCarryingCapacities();

                // check if the total carrying capacity has to be recomputed (soft/hard selection)
                if(_total_carrying_capacity != my_NAN) set_total_carrying_capacity();
            }

            unsigned int linkedTraits = _vPatch[0]->get_nbLinkedTraits();
            if(meanVe)      {
                _pSelection->set_ve_mean();
            }
            if(h2){
                set_patch_parameter(linkedTraits, "patch_ve_var", "heritability", &Patch::set_localh2Ve);
                _pSelection->reset_Ve();
            }

            if(optima)      set_patch_parameter(linkedTraits, "patch_stab_sel_optima", "optima", &Patch::set_localOptima);
            if(intensity)   set_patch_parameter(linkedTraits, "patch_stab_sel_intensity", "intensity", &Patch::set_localIntensity);

            if(min)         set_patch_parameter(linkedTraits, "patch_dir_sel_min", "min", &Patch::set_localMin);
            if(max)         set_patch_parameter(linkedTraits, "patch_dir_sel_max", "max", &Patch::set_localMax);
            if(growth_rate) set_patch_parameter(linkedTraits, "patch_dir_sel_growth_rate", "growth rate", &Patch::set_localGrowthRate);
            if(max_growth)  set_patch_parameter(linkedTraits, "patch_dir_sel_max_growth", "max growth", &Patch::set_localMaxGrowth);
            if(symmetry)    set_patch_parameter(linkedTraits, "patch_dir_sel_symmetry", "symmetry", &Patch::set_localSymmetry);
        }
    }
}

// ----------------------------------------------------------------------------------------
// ini_globs
// ----------------------------------------------------------------------------------------
bool Metapop::init(map< string,TTraitProto* >& traits, map< int,LCE* >& LCEs)
{
    if(!_paramSet->isSet()) fatal("parameters in 'population' are not properly set!\n");

    set_SexInitRatio(LCEs);

    createPopulations();
    setCarryingCapacities();

    setInitPopulationSizes();
    set_sample_pops();
    makePrototype(traits, this);

    // where does selection acts?
    _selection_position = (int)_paramSet->getValue("selection_position");
    _selection_level    = (int)_paramSet->getValue("selection_level");
    if(_selection_position != 4){    // 4: no selection at all
        // check if quantitatvie gtraits are specified
        vector<int> vQuanti = getTraitIndex("quanti");  // get the indexes to the quanti traits
        if(vQuanti.empty()) _selection_position = 4;    // if no quantitative traits are simulated => no selection at all
        else {
            if(_pSelection) delete _pSelection;
            _pSelection = new TSelection(this, vQuanti);
        }
    }

    setLifeCycle(LCEs);

    set_change_disp_rate_after_density();

    //empty and clean the RecyclingPOOL, safer...
    for(unsigned int i = 0; i < RecyclingPOOL.size(); ++i) {
        if(!RecyclingPOOL[i]) {
            error("Metapop::init: found null ptr in RecyclingPOOL!!\n");
            RecyclingPOOL.erase(RecyclingPOOL.begin() + i);
            continue;
        }
        delete RecyclingPOOL[i];
    }

    RecyclingPOOL.clear();

    return true;
}

// ----------------------------------------------------------------------------------------
// set_sample_pops
// ----------------------------------------------------------------------------------------
/** creates the vector of the pops used for the outputs
 * if matrix:       the numbers correspond to the ID of the patches to sample
 * if single value: the number correponds to the number of patches to randomly sample
 * if 0 (default):  all patches are sampled
 * Caution: input: id starts at 1 - in quantiNEMO id starts at 0
 */
void Metapop::set_sample_pops()
{
    if(_sample_pops) delete[] _sample_pops;

    // if the patches are explicitly defined
    if(_paramSet->get_param("sampled_patches")->is_matrix()){
        TMatrix* m = _paramSet->get_param("sampled_patches")->get_matrix();
        _sample_pops_size = m->get_dims(NULL);
        double* vec = m->get();

        _sample_pops = new Patch*[_sample_pops_size];
        int val;
        for(unsigned int i=0; i<_sample_pops_size; ++i){
            val = (int)vec[i]-1;
            if(val < 0)              fatal("Parameter 'sampled_patches': patch index (%i) has to be positive!\n", val+1);
            if(val > (int)_patchNbr) fatal("Parameter 'sampled_patches': patch index (%i) is exceeding the number of patches (%i)!\n", val+1, _patchNbr);
            _sample_pops[i] = getPatch(val);
        }
        return;
    }

    // it is a single value
    if(_paramSet->getValue("sampled_patches") < 0) fatal("Parameter 'sampled_patches' (%i) has to be positive!\n", _sample_pops_size);
    _sample_pops_size = (int)_paramSet->getValue("sampled_patches");
    if(!_sample_pops_size || _sample_pops_size == _patchNbr){ // default: if all patches are sampled
        _sample_pops_size = _patchNbr;
        _sample_pops = new Patch*[_sample_pops_size];
        for(unsigned int i=0; i<_sample_pops_size; ++i){
            _sample_pops[i] = getPatch(i);
        }
    }
    else{                                                     // n patches are randomly selected
        if(_sample_pops_size > _patchNbr) fatal("Parameter 'sampled_patches' (%i) has to be less or equal the total number of patches (%i)!\n", _sample_pops_size, _patchNbr);
        _sample_pops = new Patch*[_sample_pops_size];
        unsigned int* vSample = sample(_sample_pops_size, _patchNbr);   // draw n UNIQUE numbers
        for(unsigned int i=0; i<_sample_pops_size; ++i){
            _sample_pops[i] = getPatch(vSample[i]);
        }
        delete[] vSample;
    }
}

// ----------------------------------------------------------------------------------------
// buildPopulation
// ----------------------------------------------------------------------------------------
/** sets the inital sex ratio. For hemaphrodites it is 0, if both sexes are equaly set it is 1 */
void Metapop::set_SexInitRatio(map< int,LCE* >& LCEs)
{
    // find the breeding LCE
    map<int, LCE* >::iterator LCE = LCEs.begin();
    for(; LCE != LCEs.end(); ++LCE) {
        _pBreed_LCE = dynamic_cast<LCE_Breed*>(LCE->second);
        if(_pBreed_LCE) break;
    }
    assert(LCE != LCEs.end());       // the breed LCE must be present!

    if(_pBreed_LCE->get_parameter_value("mating_system") < 2){
        _sexInitRatio = 0;                      // hemaphrodites
    }
    else {
        _sexInitRatio = _pBreed_LCE->get_parameter_value("sex_ratio");
        if(!_sexInitRatio) fatal("Two sexes with a sex ratio of 0 makes no sense!\n");
        else if(_sexInitRatio < 0) _sexInitRatio = 0.5;  // if sex ratio should be kept start with equal number of males and females
        else _sexInitRatio /= _sexInitRatio+1;           // transform the sex ratio
    }
}

// ----------------------------------------------------------------------------------------
// buildPopulation
// ----------------------------------------------------------------------------------------
/** reset the right number of patches for the new simulation */
void Metapop::resize_nbPopulations(unsigned int nbPatch)
{
    if(_vPatch.size() == nbPatch) return;

    if(_vPatch.size() > nbPatch) {
        while(_vPatch.size() > nbPatch) {
            delete _vPatch[0];
            _vPatch.pop_front();
        }
    }
    else{
        while(_vPatch.size() < nbPatch) {
            _vPatch.push_back(new Patch());
        }
    }

    // set the number of pops (static variable)
    _vPatch[0]->set_nbPatch(nbPatch);
    _vPatch[0]->set_pMetapop(this);
}

// ----------------------------------------------------------------------------------------
// buildPopulation
// ----------------------------------------------------------------------------------------
void Metapop::createPopulations()
{
    // get the number of populations
    _patchNbr = (unsigned int)_paramSet->getValue("patch_number");
    resize_nbPopulations(_patchNbr);

    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->init(i);
    }
}

// ----------------------------------------------------------------------------------------
// setInitPopulationSizes
// ----------------------------------------------------------------------------------------
void Metapop::setInitPopulationSizesOfPatches()
{
    //set the population sizes:
    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->set_PopSizes_ini(_patchNfem,_patchNmal);
    }
}

// ----------------------------------------------------------------------------------------
// setInitPopulationSizes
// ----------------------------------------------------------------------------------------
void Metapop::setInitPopulationSizesOfPatches(TMatrix* popK)
{
    // get the vector and its size
    unsigned int nbTot, nbMal, nbK = popK->get_dims(NULL);
    double* vec = popK->get();

    // check the number of carrying capacities with the number of patches
    if(nbK>_patchNbr) warning("There are more inital population sizes than patches! Only a part of the inital population sizes is considered!\n");
    else if(_patchNbr % nbK) warning("The number of inital population sizes is not a entire subset of the number of patches!\n");

    for(unsigned int i = 0; i < _patchNbr; ++i){
        nbTot = (unsigned int) vec[i % nbK];  	// total number of individuals
        nbMal = my_round(_sexInitRatio*nbTot);  // number of males due to the sex ratio
        _vPatch[i]->set_PopSizes_ini(nbTot-nbMal, nbMal);
    }
}

// ----------------------------------------------------------------------------------------
// setInitPopulationSizes
// ----------------------------------------------------------------------------------------
void Metapop::setInitPopulationSizesOfPatches(sex_t SEX, TMatrix* popK)
{
    // get the vector and its size
    unsigned int nbK = popK->get_dims(NULL);
    double* vec = popK->get();

    // check the number of carrying capacities with the number of patches
    if(nbK>_patchNbr) warning("There are more inital population sizes than patches! Only a part of the inital population sizes is considered!\n");
    else if(_patchNbr % nbK) warning("The number of inital population sizes is not a entire subset of the number of patches!\n");

    if(SEX == FEM){
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_PopSizes_ini((unsigned int)vec[i % nbK],_patchNmal);
        }
    }
    else{
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_PopSizes_ini(_patchNfem,(unsigned int)vec[i % nbK]);
        }
    }
}

// ----------------------------------------------------------------------------------------
// setInitPopulationSizes
// ----------------------------------------------------------------------------------------
void Metapop::setInitPopulationSizesOfPatches(TMatrix* popKfem, TMatrix* popKmal)
{
    // get the vector and its size
    unsigned int nbKmal = popKmal->get_dims(NULL);
    double* vec_mal = popKmal->get();
    unsigned int nbKfem = popKfem->get_dims(NULL);
    double* vec_fem = popKfem->get();

    // check the number of carrying capacities with the number of patches
    if(nbKmal>_patchNbr) warning("There are more inital population sizes (males) than patches! Only a part of the inital population sizes is considered!\n");
    else if(_patchNbr % nbKmal) warning("The number of inital population sizes (males) is not a entire subset of the number of patches!\n");
    if(nbKfem>_patchNbr) warning("There are more inital population sizes (females) than patches! Only a part of the inital population sizes is considered!\n");
    else if(_patchNbr % nbKfem) warning("The number of inital population sizes (females) is not a entire subset of the number of patches!\n");

    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->set_PopSizes_ini((unsigned int)vec_fem[i % nbKfem],(unsigned int)vec_mal[i % nbKmal]);
    }
}

// ----------------------------------------------------------------------------------------
// setCarringCapacities
// ----------------------------------------------------------------------------------------
void Metapop::setCarryingCapacitiesOfPatches()
{
    //set the population sizes:
    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->set_carrying_capacities(_patchNfem,_patchNmal);
    }
}

// ----------------------------------------------------------------------------------------
// setCarringCapacities
// ----------------------------------------------------------------------------------------
void Metapop::setCarryingCapacitiesOfPatches(TMatrix* popK)
{
    // get the vector and its size
    unsigned int nbK = popK->get_dims(NULL);
    double* vec = popK->get();

    // check the number of carrying capacities with the number of patches
    if(nbK>_patchNbr) warning("There are more carrying capacities than patches! Only a part of the carrying capacities is considered!\n");
    else if(_patchNbr % nbK) warning("The number of carrying capacities is not a entire subset of the number of patches!\n");

    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->set_carrying_capacities((unsigned int)vec[i % nbK]/2,(unsigned int)vec[i % nbK]/2);
    }
}

// ----------------------------------------------------------------------------------------
// setCarringCapacities
// ----------------------------------------------------------------------------------------
void Metapop::setCarryingCapacitiesOfPatches(sex_t SEX, TMatrix* popK)
{
    // get the vector and its size
    unsigned int nbK = popK->get_dims(NULL);
    double* vec = popK->get();

    // check the number of carrying capacities with the number of patches
    if(nbK>_patchNbr) warning("There are more carrying capacities than patches! Only a part of the carrying capacities is considered!\n");
    else if(_patchNbr % nbK) warning("The number of carrying capacities is not a entire subset of the number of patches!\n");

    if(SEX == FEM){
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_carrying_capacities((unsigned int)vec[i % nbK],_patchNmal);
        }
    }
    else{
        for(unsigned int i = 0; i < _patchNbr; ++i){
            _vPatch[i]->set_carrying_capacities(_patchNfem,(unsigned int)vec[i % nbK]);
        }
    }
}

// ----------------------------------------------------------------------------------------
// setCarringCapacities
// ----------------------------------------------------------------------------------------
void Metapop::setCarryingCapacitiesOfPatches(TMatrix* popKfem, TMatrix* popKmal)
{
    // get the vector and its size
    unsigned int nbKmal = popKmal->get_dims(NULL);
    double* vec_mal = popKmal->get();
    unsigned int nbKfem = popKfem->get_dims(NULL);
    double* vec_fem = popKfem->get();

    // check the number of carrying capacities with the number of patches
    if(nbKmal>_patchNbr) warning("There are more carrying capacities (males) than patches! Only a part of the carrying capacities is considered!\n");
    else if(_patchNbr % nbKmal) warning("The number of carrying capacities (males) is not a entire subset of the number of patches!\n");
    if(nbKfem>_patchNbr) warning("There are more carrying capacities (females) than patches! Only a part of the carrying capacities is considered!\n");
    else if(_patchNbr % nbKfem) warning("The number of carrying capacities (females) is not a entire subset of the number of patches!\n");

    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->set_carrying_capacities((unsigned int)vec_fem[i % nbKfem],(unsigned int)vec_mal[i % nbKmal]);
    }
}

// ----------------------------------------------------------------------------------------
// set_patch_selection_slope
// ----------------------------------------------------------------------------------------
/** sets patch sepecific parameters for each quantitative trait. They may be passed as:
 *  - sex specific (in this case both sexes (if there are two sexes) have to be specified)
 *  - value: same vlaue for all patches
 *  - matrix: matrixes may be expanded to meet the number of patches and traits
 */
void Metapop::set_patch_parameter(const unsigned int& nbTrait,
        const string& name,             // "patch_dir_sel_slope"
        const string& name_full,        // "slope"
        void (Patch::*pt2Func)(double*, sex_t))
{
    // sex specific settings have precedence
    if(_paramSet->isSet(name+"_fem")){
        if(_paramSet->isSet(name+"_mal") || !_sexInitRatio){
            // if matrix
            if(_paramSet->isMatrix(name+"_fem")){
                set_patch_matrix(nbTrait, _paramSet->getMatrix(name+"_fem"), FEM, name_full, pt2Func);
            }
            else{
                set_patch_value(nbTrait, _paramSet->getValue(name+"_fem"), FEM, pt2Func);
            }
            if(_sexInitRatio){
                if(_paramSet->isMatrix(name+"_mal")){
                    set_patch_matrix(nbTrait, _paramSet->getMatrix(name+"_mal"), MAL, name_full, pt2Func);
                }
                else{
                    set_patch_value(nbTrait, _paramSet->getValue(name+"_mal"), MAL, pt2Func);
                }
            }
        }
        else fatal("Only one sex specific %s is specified: both are required!\n", name_full.c_str());
    }
    else{                         // general settings
        if(_paramSet->isMatrix(name)){
            set_patch_matrix(nbTrait, _paramSet->getMatrix(name), FEM, name_full, pt2Func);
            if(_sexInitRatio) set_patch_matrix(nbTrait, _paramSet->getMatrix(name), MAL, name_full, pt2Func);
        }
        else{
            set_patch_value(nbTrait, _paramSet->getValue(name), FEM, pt2Func);
            if(_sexInitRatio) set_patch_value(nbTrait, _paramSet->getValue(name), MAL, pt2Func);
        }
    }
}

// ----------------------------------------------------------------------------------------
// set_patch_value
// ----------------------------------------------------------------------------------------
/** this function sets the same value for each patch and trait*/
void Metapop::set_patch_value(const unsigned int& nbTrait,
        const double& value,
        sex_t SEX,
        void (Patch::*pt2Func)(double*, sex_t))
{
    unsigned int p, t;
    double *array = new double[nbTrait];

    // get the slope for each patch/trait
    for(p = 0; p < _patchNbr; ++p) {
        for(t = 0; t < nbTrait; ++t) {
            array[t] = value;
        }
        _vPatch[p]->set_localParameter(array, SEX, pt2Func);
    }
    delete[] array;
}

// ----------------------------------------------------------------------------------------
// set_patch_matrix
// ----------------------------------------------------------------------------------------
/** this function sets the value to each patch based on a matrix */
void Metapop::set_patch_matrix(const unsigned int& nbTrait,
        TMatrix* m,
        sex_t SEX,
        const string& name_full,        // "slope"
        void (Patch::*pt2Func)(double*, sex_t))
{
    unsigned int p, t;

    unsigned int dims[2];
    m->get_dims(&dims[0]);
    unsigned int count_trait = dims[0];   //nbr trait slope per patch
    unsigned int count_patch = dims[1];   //nbr of patches slope

    // check the number of slope per number of patches
    if(count_patch>_patchNbr) warning("There are more %s specified than patches! Only a part of the %s is considered!\n", name_full.c_str(), name_full.c_str());
    else if(_patchNbr % count_patch) warning("The number of %s is not a entire subset of the number of patches!\n", name_full.c_str());

    // check the number of slope per number of traits
    if(count_trait>nbTrait) warning("There are more %s specified than traits! Only a part of the %s is considered!\n", name_full.c_str(), name_full.c_str());
    else if(nbTrait % count_trait) warning("The number of %s is not a entire subset of the number of traits!\n", name_full.c_str());

    // create the arrays
    double* array = new double[nbTrait];

    // get the optima for each patch/trait
    for(p = 0; p < _patchNbr; ++p) {
        for(t = 0; t < nbTrait; ++t) {
            array[t] = m->get(t % count_trait, p % count_patch);
        }
        _vPatch[p]->set_localParameter(array, SEX, pt2Func);
    }

    delete[] array;
}

// ----------------------------------------------------------------------------------------
// setLifeCycle
// ----------------------------------------------------------------------------------------
void Metapop::setLifeCycle(map< int, LCE*>& lifeCycle)
{
    string name;
    _theCycle.clear();

    for(map< int,LCE* >::iterator LCE = lifeCycle.begin() ; LCE != lifeCycle.end(); ++LCE) {
        if(LCE->second->init(this)){
            _theCycle.push_back(LCE->second);
        }
    }
}

// ----------------------------------------------------------------------------------------
// Replicate_LOOP
// ----------------------------------------------------------------------------------------
void Metapop::Replicate_LOOP()
{
    clock_t start;
    clock_t stop;

    _meanReplElapsedTime = 0;
    _meanGenLength = 0;

    message("\n\nSIMULATION");

    // for each replicate
    for(_currentReplicate = 1; _currentReplicate <= _replicates; ++_currentReplicate) {
#ifdef _DEBUG
        message("\n\r**** replicate %i/%i [%s] ****           \n",_currentReplicate,_replicates,getElapsedTime(0).c_str());
#else
        message("\n\r    replicate %i/%i [%s] %i/%i",_currentReplicate,_replicates,getElapsedTime(0).c_str(),0,_generations);
#endif
        fflush(stdout);

        if(_currentReplicate!=1){
            reset_dynamic_params();
            resetPrototype();
        }

        //build metapopulation for the current replicate (new adult generation)
        if(TTQuantiProto::_ini_fstat_file!="" || TTNeutralProto::_ini_fstat_file!=""){
            setPopulation_FSTAT();
        }
        else setPopulation();

        executeBeforeEachReplicate(_currentReplicate);
        start = clock();
        //--------------------------- GENERATION LOOP ---------------------------

        Cycle(start);

        //-----------------------------------------------------------------------
        stop = clock();

        executeAfterEachReplicate(_currentReplicate);

        _meanReplElapsedTime += (stop - start);

        _meanGenLength += _currentGeneration - 1;

        if(!isAlive() && _currentGeneration != _generations) {
            _currentGeneration = _generations;
            if(get_service()) get_service()->notify();
        }
    } // end of for each replicate


    _meanReplElapsedTime /= _replicates;
    _meanGenLength /= _replicates;

    clear();

}
/***********************************************************************
 * Metapop::Cycle
 *
 * performs the generation loop for current replicate
 *
 */
void Metapop::Cycle(clock_t& startTime)
{
    list< LCE* >::iterator LCE;

#ifdef _DEBUG
    message("\n");
    message("Starting point (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            getPatchNbr(),
            size(FEM, OFFSPRG),  size(MAL, OFFSPRG),
            size(FEM, ADULTS),   size(MAL, ADULTS));
#endif

    int nbOutput = _generations/10;
    if(!nbOutput) nbOutput = 1;   // if less than 10 generations are simulated


    // ------------------------------ CYCLE --------------------------------
    for(_currentGeneration = 1; _currentGeneration <= _generations; _currentGeneration++) {
        executeBeforeEachGeneration(_currentGeneration);

        if(_currentGeneration<=10 || !(_currentGeneration % nbOutput)){
            message("\r    replicate %i/%i [%s] %i/%i",_currentReplicate,_replicates
                    ,getElapsedTime(clock() - startTime).c_str(),_currentGeneration,_generations);
            fflush(stdout);
        }

#ifdef _DEBUG
        message("\n%i. generation (%i. replicate) ***********\n", _currentGeneration, _currentReplicate);
#endif

        LCE = _theCycle.begin();
        while(LCE != _theCycle.end() && isAlive()) {
            (*LCE)->execute();
            _currentAge ^= (*LCE)->removeAgeClass();
            _currentAge |= (*LCE)->addAgeClass();
            LCE++;
        }
        if(!isAlive()) {
            message("\r    replicate %i/%i [%s] %i/%i -> Pop extinction !\n",_currentReplicate,_replicates
                    ,getElapsedTime(clock() - startTime).c_str(),_currentGeneration,_generations);
            _currentGeneration++;
            break;
        }

        executeAfterEachGeneration(_currentGeneration);

    }

    // --------------------------- END OF CYCLE --------------------------

}

//------------------------------------------------------------------------------
    void
Metapop::executeBeforeEachGeneration(const int& gen)
{
    // density dependent temporal change of the dispersal rate
    change_disp_rate_after_density(gen);

    // metapop
    temporal_change(gen);

    // proto traits
    map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
    for(; iter != getTraitPrototypes().end(); ++iter){
        iter->second->executeBeforeEachGeneration(gen);
    }

    // LCEs
    list<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->executeBeforeEachGeneration(gen);
    }

    if(_pSelection) _pSelection->executeBeforeEachGeneration(gen);
}

//------------------------------------------------------------------------------
    void
Metapop::executeAfterEachGeneration(const int& gen)
{
    // proto traits
    map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
    for(; iter != getTraitPrototypes().end(); ++iter){
        iter->second->executeAfterEachGeneration(gen);
    }

    // LCEs
    list<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->executeAfterEachGeneration(gen);
    }
}

//------------------------------------------------------------------------------
    void
Metapop::executeBeforeEachReplicate(const int& rep)
{
    // traits
    map< string,TTraitProto* >::iterator iterTrait = getTraitPrototypes().begin();
    for(; iterTrait != getTraitPrototypes().end(); ++iterTrait){
        iterTrait->second->executeBeforeEachReplicate(rep);
    }

    // selection
    if(_pSelection) _pSelection->executeBeforeEachReplicate(rep);
    if(_total_carrying_capacity != my_NAN){
        set_total_carrying_capacity();   // soft selection at the metapopualtiuon level
    }

    // LCEs
    list<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->executeBeforeEachReplicate(rep);
    }
}

//------------------------------------------------------------------------------
    void
Metapop::executeAfterEachReplicate(const int& rep)
{
    // traits
    map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
    for(; iter != getTraitPrototypes().end(); ++iter){
        iter->second->executeAfterEachReplicate(rep);
    }

    // LCEs
    list<LCE*>::iterator iterLCE = _theCycle.begin();
    for(; iterLCE != _theCycle.end(); ++iterLCE){
        (*iterLCE)->executeAfterEachReplicate(rep);
    }

    printEntireGeneticMap();
}

// ----------------------------------------------------------------------------------------
// densityDependend_dispRate_change
// ----------------------------------------------------------------------------------------
/** this function allows to change the dispersal rate depending on the density
 * i.e. at the specified ratio the temoral parameter is changed
 * _density_threshold[n][3] = 0: patch, 1: density, 2: change   //n is the number of changes
 */
void
Metapop::change_disp_rate_after_density(const int& gen){
    if(!_density_threshold) return;

    double density;
    Patch* curPatch;
    map<int, string>::iterator iter_tempParam;
    for(unsigned int i=0; i<_density_threshold_nbChanges; ++i){
        if(_density_threshold[i][0] == my_NAN) continue;		// if already used
        curPatch = getPatch(_density_threshold[i][0]);
        if(!curPatch->get_K()) continue;										// the patch has to have a carring capacity
        density = (double)curPatch->size(ADLTx) / curPatch->get_K();
        if(density >= _density_threshold[i][1]){            // if a change has to be made
            _density_threshold[i][0] = my_NAN;              	// set to NAN to recognise used items

            // get the x item (it is a map so scroll through it)
            iter_tempParam = _density_threshold_param[i]->get_temporal_args()->begin();
            for(int j=0; j<_density_threshold[i][3]; ++j, ++iter_tempParam){}
            string arg  = iter_tempParam->second;
            int cur_gen = iter_tempParam->first;
            int new_gen = gen+_density_threshold[i][2];
            (*_density_threshold_param[i]->get_temporal_args())[new_gen] = arg; // make the change
            _density_threshold_param[i]->get_temporal_args()->erase(iter_tempParam); // delete the old item

            // also the paramset has to be changed...
            multimap<int, map<string, Param*>* >* m = _density_threshold_param[i]->getParamSet()->getTemporalParams();
            multimap<int, map<string, Param*>* >::iterator pos = m->find(cur_gen);
            assert(pos != m->end());
            delete pos->second;
            m->erase(pos);               // delete the item
            _density_threshold_param[i]->getParamSet()->set_temporal_param(new_gen, _density_threshold_param[i]);  // add the new one
        }
    }
}

// ----------------------------------------------------------------------------------------
/** this function allows to change the dispersal rate depending on the density
 * i.e. at the specified ratio the temoral parameter is changed
 * _density_threshold[n][3] = 0: patch, 1: density, 2: delay, 3: change   //n is the number of changes
 */
    void
Metapop::set_change_disp_rate_after_density()
{
    Param* p = _paramSet->get_param("temporal_change_following_density");
    if(_density_threshold){delete[] _density_threshold; _density_threshold=NULL;}
    if(!p->isSet())return;

    // if param is specified
    if(!p->is_matrix()) fatal("Parameter 'temporal_change_following_density' must be a matrix!");
    TMatrixVar<string>* m = p->parse_matrixVarStr();

    _density_threshold_nbChanges = m->getNbRows();
    Param*    curParam;
    string paramName;
    vector<string>* curVec;
    _density_threshold  = ARRAY::new_2D<double>(_density_threshold_nbChanges, 4);
    _density_threshold_param = new Param*[_density_threshold_nbChanges];
    for(unsigned int i=0; i<_density_threshold_nbChanges; ++i){
        curVec = m->get(i);
        if(curVec->size() != 5) fatal("Parameter 'temporal_change_following_density' must be a matrix with 5 elements per row!");
        paramName = (*curVec)[0];

        // check the metapop params
        curParam = _paramSet->find_param(paramName);

        // check the LCE params
        if(!curParam){
            list<LCE*>::iterator  iterLCE = _theCycle.begin();
            for(;iterLCE != _theCycle.end(); ++iterLCE) {
                curParam = (*iterLCE)->get_paramset()->find_param(paramName);
                if(curParam) break;
            }
        }

        // check the trait params
        if(!curParam){
            map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
            for(; iter != getTraitPrototypes().end(); ++iter){
                curParam = iter->second->get_paramset()->find_param(paramName);
                if(curParam) break;
            }
        }

        // was the param found?
        if(curParam) _density_threshold_param[i] = curParam;
        else fatal("Parameter 'temporal_change_following_density' contains unknown parameter name '%s'!", paramName.c_str());

        // get the dynamic change settings
        _density_threshold[i][0] = (double)(STRING::str2int<int>((*curVec)[1])-1); // patch (here starting at 0)
        _density_threshold[i][1] = STRING::str2int<double>((*curVec)[2]);          // density
        _density_threshold[i][2] = (double)(STRING::str2int<int>((*curVec)[3]));   // delay
        _density_threshold[i][3] = (double)(STRING::str2int<int>((*curVec)[4])-1); // change (here starting at 0)

        // some tests
        if(_density_threshold[i][0]>=0 && _density_threshold[i][0]<_patchNbr) fatal("Parameter 'temporal_change_following_density' must have a valide patch id!");
        if(_density_threshold[i][2]<0) fatal("Parameter 'temporal_change_following_density' must have a postive delay!");
        if(!curParam->isTemporalParam() || curParam->get_temporal_args()->size()<_density_threshold[i][3])
            fatal("Parameter 'temporal_change_following_density' should at least contain %i temporal changes!", _density_threshold[i][3]);
    }

    delete m;
}

// ----------------------------------------------------------------------------------------
/** this function allows to reset the dynamicly changed temoral parameters to default settings
 * this fundction must be called between replicates!
 */
    void
Metapop::reset_dynamic_params()
{
    if(_currentReplicate<=1) return; 		// do it only between replicates
    if(!_density_threshold) return;    // dynamic parameters are not used

    Param* curParam;
    multimap<int, map<string, Param*>* >* pMultimap;                // pointer to multimap
    multimap<int, map<string, Param*>* >::iterator iter_multimap;   // iterator for multimap
    map<string, Param*>* pMap;                  										// pointer to the map
    map<string, Param*>::iterator iter_map;                  				// iterator for map


    for(unsigned int i=0; i<_density_threshold_nbChanges; ++i){
        curParam = _density_threshold_param[i];
        pMultimap = curParam->getParamSet()->getTemporalParams();
        iter_multimap= pMultimap->begin();

        // remove all entries of the param in the temporal multimap of the paramset
        for(; iter_multimap != pMultimap->end(); ++iter_multimap){
            pMap = iter_multimap->second;
            iter_map = pMap->find(curParam->get_name());
            if(iter_map != pMap->end()){											// item found -> delete it
                pMap->erase(iter_map);
                if(pMap->empty()){            									// check if the map is noe empty
                    delete pMap;                 									// delete the map
                    pMultimap->erase(iter_multimap);              // remocve the element from the multimap
                }
            }
        }

        // reset the param (and thus also the paramset) to the default settings
        curParam->set_temporal_args(curParam->get_tot_arg());
    }
}

// ----------------------------------------------------------------------------------------
// set_total_carrying_capacity
// ----------------------------------------------------------------------------------------
void
Metapop::set_total_carrying_capacity(){
    _total_carrying_capacity = 0;
    for(unsigned int home = 0; home < _patchNbr; ++home) {
        _total_carrying_capacity += _vPatch[home]->get_K();
    }
}

//------------------------------------------------------------------------------
/** outputs the genetic map for all the traits */
    void
Metapop::printEntireGeneticMap()
{
    if(!TTraitProto::getNbChromosome() || _service->get_genetic_map_output()) return;

    // create the file name
    string file = _service->getSimfolder()+ "genetic_map";
    if(getReplicates()>1) file += "_" + STRING::int2str_(getCurrentReplicate());
    file += ".txt";

    // open the file
    ofstream FILE(file.c_str());
    if(!FILE) error("opening file '%s' to write the genetic map!\n", file.c_str());
    else{
        // print the heading
        FILE << "#Genetic map";
        if(getReplicates()>1) FILE << " (Replicate: " << getCurrentReplicate() << ")";
        FILE << "\n#" << setfill('=') << setw(30) << "" << setfill(' ');

        // print the genetic map of each trait on a separate line
        map< string,TTraitProto* >::iterator iter = getTraitPrototypes().begin();
        for(; iter != getTraitPrototypes().end(); ++iter){
            iter->second->print_genetic_map(FILE);          // print genetic map
        }
        FILE.close();
    }
}

//------------------------------------------------------------------------------
/** Metapop::setPopulation */
    void
Metapop::setPopulation()
{
#ifdef _DEBUG
    message("Metapop::setPopulation: ");
#endif

    for(unsigned int i = 0; i < _patchNbr; ++i) {
        _vPatch[i]->reset_counters();
        _vPatch[i]->flush();
    }

    _currentAge = NONE;
    age_t requiredAge = NONE;

    //find the age class required to start the life cycle with:
    list< LCE* >::iterator IT = _theCycle.begin();
    for(; IT != _theCycle.end() && requiredAge == NONE; IT++){
        requiredAge = (*IT)->requiredAgeClass();
    }

#ifdef _DEBUG
    message("required age is: %i ",requiredAge);
#endif

    for(unsigned int i = 0; i < _patchNbr; ++i){
        _vPatch[i]->setNewGeneration(requiredAge);
    }
    _currentAge = requiredAge;

#ifdef _DEBUG
    message("loaded age is: %i (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            _currentAge, getPatchNbr(),
            size(FEM, OFFSPRG),  size(MAL, OFFSPRG),
            size(FEM, ADULTS),   size(MAL, ADULTS));
#endif

    if(_currentAge == NONE) warning("Metapop::setPopulation: generation 0 is empty!!\n");
}

//------------------------------------------------------------------------------
/** Metapop::setPopulation_FSTAT */

    void
Metapop::setPopulation_FSTAT()
{

#ifdef _DEBUG
    message("Metapop::setPopulation after FSTAT-file: ");
#endif

    for(unsigned int i = 0; i < _patchNbr; ++i) {
        _vPatch[i]->reset_counters();
        _vPatch[i]->flush();
    }

    _currentAge = NONE;
    age_t requiredAge = NONE;

    //find the age class required to start the life cycle with:
    list< LCE* >::iterator IT = _theCycle.begin();
    for(; IT != _theCycle.end() && requiredAge == NONE; IT++){
        requiredAge = (*IT)->requiredAgeClass();
    }

#ifdef _DEBUG
    message("required age is: %i ",requiredAge);
#endif

    // read the genotypes if present
    vector<unsigned int**> *vQuanti = NULL, *vNtrl = NULL, *vAny;
    string filename = TTQuantiProto::_ini_fstat_file;
    if(TTQuantiProto::_ini_fstat_file != "")   vQuanti = read_Fstat_file(TTQuantiProto::_ini_fstat_file,  FileHandler::get_base_path());
    if(TTNeutralProto::_ini_fstat_file != "")  vNtrl   = read_Fstat_file(TTNeutralProto::_ini_fstat_file, FileHandler::get_base_path());

    // if both files are given check if the characteristics of the individuals match
    if(vQuanti){
        if(vNtrl){
            if(vQuanti->size() != vNtrl->size()) fatal("inital FSTAT files: the number of individuals does not correspond between files!\n");
            for(unsigned int i=0, size=vQuanti->size(); i<size; ++i){
                if(vQuanti[i][0][0] != vNtrl[i][0][0]          // current patch
                        || vQuanti[i][1][0] != vNtrl[i][1][0]        // age
                        || vQuanti[i][2][0] != vNtrl[i][2][0]        // sex
                        || vQuanti[i][3][0] != vNtrl[i][3][0]        // id
                        || vQuanti[i][3][1] != vNtrl[i][3][1]        // natal patch
                        || vQuanti[i][4][0] != vNtrl[i][4][0]        // id mom
                        || vQuanti[i][4][1] != vNtrl[i][4][1]        // patch mom
                        || vQuanti[i][5][0] != vNtrl[i][5][0]        // id dad
                        || vQuanti[i][5][1] != vNtrl[i][5][1]){      // patch dad
                    fatal("inital FSTAT files: the charactersitics of the individuals do not match between genotype files!\n");
                }
            }
        }
        vAny = vQuanti;   // at least the quantitative genotypes are available
    }
    else vAny = vNtrl;  // only the neutral genotypes are available

    /*      - v[ind][0][0] = patch    // starting at 0 (compulsory)
     *       - v[ind][0][1] = nb_loci  // number of locus found in the file (compulsory)
     *       - v[ind][1][0] = age      // off: 0, adlt: 2
     *       - v[ind][2][0] = sex      // mal: 0, fem: 1
     *       - v[ind][3][0] = id       // id of the current patch starting at 1
     *       - v[ind][3][1] = id       // id of the natal patch starting at 1
     *       - v[ind][4][0] = id_mom   // id of mother of the current patch starting at 1
     *       - v[ind][4][1] = id_mom   // id of mother of the natal patch starting at 1
     *       - v[ind][5][0] = id_dad   // id of father of the current patch starting at 1
     *       - v[ind][5][1] = id_dad   // id of father of the natal patch starting at 1
     *       - v[ind][locus+6][allele] // allele starting at 0
     */


    Individual* cur_ind;
    unsigned int **aQuanti = NULL,   // quantitative gentyope array of a given individual
                 **aNtrl = NULL,     // neutral genotype array of a given individual
                 **array = NULL,     // any present gentoype array of a given individual used to set the characteristics
                 **aQuantiCur = NULL,// current postion in the quantitative array
                 **aNtrlCur = NULL,  // current postion in the neutral array
                 **aQuantiEnd = NULL,// after last position in the locus array
                 **aNtrlEnd = NULL;  // after last position in the locus array
    unsigned int nbLocusTotQ = vQuanti ? (*vQuanti)[0][0][1] : 0,  // total number of quantitative loci
                 nbLocusTotN = vNtrl ?   (*vNtrl)[0][0][1]   : 0,  // total number of neutral loci
                 nbInd       = vAny->size();                       // total number of individuals
    unsigned int *aMax = NULL;   // array for each patch to find the maximal index (only used when infos are present)
    ARRAY::create_1D(aMax, _patchNbr, (unsigned int)0);
    TTraitProto* pProto;
    age_idx age;
    unsigned int t;
    int p, l;
    Patch* pPatch=NULL;
    unsigned char** seq;                           // pointer to the sequence
    for(unsigned int i=0; i<nbInd; ++i){           // for each individual
        aQuanti = vQuanti ? (*vQuanti)[i] : NULL;    // get individual array of quanti
        aNtrl   = vNtrl   ? (*vNtrl)[i]   : NULL;    // get individual array of ntrl
        array   = (*vAny)[i];                        // any which is available

        // create individual
        cur_ind = _protoIndividual.clone();
        cur_ind->init();                 //allocate memory for the traits sequence

        // current patch
        pPatch = getPatch(array[0][0]);
        cur_ind->setCurrentPatch(pPatch);

        // age of the individual
        if(array[1][0] != my_NAN) age = (age_idx)array[1][0];
        else                      age = (age_idx)(requiredAge==OFFSPRG ? 0 : 1);

        // sex of the individual
        if(array[2][0] != my_NAN) cur_ind->setSex((sex_t)array[2][0]);
        else                      cur_ind->setSex(SimRunner::r.Uniform() < _sexInitRatio ? MAL : FEM);

        // ID of the individual
        if(array[3][0] != my_NAN && array[3][1] != my_NAN){
            if(array[3][1] > _patchNbr)fatal("inital FSTAT files: the patch index exceeds the number of patches!\n");
            cur_ind->setID(STRING::int2str(array[3][0]) + "_" + STRING::int2str(array[3][1]));

            // set natal patch
            cur_ind->setNatalPatch(getPatch(array[3][1]));

            // ge the maximal id of an individual
            if(aMax[array[3][1]-1] < array[3][0]) aMax[array[3][1]-1] = array[3][0];
        }
        else{
            cur_ind->setID(pPatch->get_next_IDofIndividual());
            cur_ind->setNatalPatch(pPatch); // set the current instead of the natal patch
            cur_ind->setCurrentPatch(pPatch);
        }

        // ID of the mother
        if(array[4][0] != my_NAN && array[4][1] != my_NAN){
            if(array[4][1] > _patchNbr)fatal("inital FSTAT files: the patch index exceeds the number of patches!\n");
            cur_ind->setMotherID(STRING::int2str(array[4][0]) + "_" + STRING::int2str(array[4][1]));
        }
        else cur_ind->setMotherID("");

        // ID of the father
        if(array[5][0] != my_NAN && array[5][1] != my_NAN){
            if(array[5][1] > _patchNbr)fatal("inital FSTAT files: the patch index exceeds the number of patches!\n");
            cur_ind->setFatherID(STRING::int2str(array[5][0]) + "_" + STRING::int2str(array[5][1]));
        }
        else cur_ind->setFatherID("");

        // add the indvidual to the container
        pPatch->add(cur_ind->getSex(), age, cur_ind);

        // set the genotypes
        if(aQuanti){
            aQuantiCur = aQuanti + 6;                // go to the beginning of the genotype
            aQuantiEnd = aQuantiCur + nbLocusTotQ;   // find the end of the array
        }
        if(aNtrl){
            aNtrlCur = aNtrl + 6;                    // go to the beginning of the genotype
            aNtrlEnd = aNtrlCur + nbLocusTotN;       // find the end of the array
        }
        for(t=0; t<getTraitPrototypeSize(); ++t){             // scroll through each trait // they are sorted in the same way
            pProto = &getTraitPrototype(t);                     // get the trait
            if(pProto->_type.find("quanti")!=string::npos){     // check if it is a quantitative trait
                if(aQuanti){                                      // genotypes are given
                    // check the number of loci
                    if(array[0][1] != (unsigned int)pProto->_nb_locus) fatal("initial FSTAT file (quanti): the number of loci (%i) does not correspond to the settings (%i)!\n",
                            array[0][1], pProto->_nb_locus);
                    seq = cur_ind->Traits[t]->sequence;
                    for(l = 0; l < pProto->_nb_locus && aQuantiCur!=aQuantiEnd; ++l, ++aQuantiCur) { // for each locus
                        for(p = 0; p < pProto->_ploidy; ++p) {        // for each allele
                            seq[l][p] = (*aQuantiCur)[p];
                        }
                    }
                    cur_ind->Traits[t]->set_value();
                }
                else{                                             // genotypes are not given
                    cur_ind->Traits[t]->ini_sequence(getPatch(array[0][0]));   // initialize them randomly
                }
            }
            else if(pProto->_type.find("ntrl")!=string::npos){  // check if it is a neutral trait
                if(aNtrl){                                        // genotypes are given
                    // check the number of loci
                    if(array[0][1] != (unsigned int)pProto->_nb_locus) fatal("initial FSTAT file (ntrl): the number of loci (%i) does not correspond to the settings (%i)!\n",
                            array[0][1], pProto->_nb_locus);
                    seq = cur_ind->Traits[t]->sequence;
                    for(l = 0; l < pProto->_nb_locus && aNtrlCur!=aNtrlEnd; ++l, ++aNtrlCur) {
                        for(p = 0; p < pProto->_ploidy; ++p) {
                            seq[l][p] = (*aNtrlCur)[p];
                        }
                    }
                    cur_ind->Traits[t]->set_value();
                }
                else{                                             // genotypes are not given
                    cur_ind->Traits[t]->ini_sequence(getPatch(array[0][0]));   // initialize them randomly
                }
            }
        }
        if(aQuanti) ARRAY::delete_2D(aQuanti, nbLocusTotQ+6);
        if(aNtrl)   ARRAY::delete_2D(aNtrl,   nbLocusTotN+6);
    }
    if(vQuanti) delete vQuanti;
    if(vNtrl)   delete vNtrl;

    // adjust the individual indexer of the patch to the highest present index
    for(t=0; t<_patchNbr; ++t){                                 // for each patch
        if(aMax[t]) pPatch->set_ID_individual(++aMax[t]);
    }
    delete[] aMax;

    _currentAge = requiredAge;

#ifdef _DEBUG
    message("loaded age is: %i (Patches: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            _currentAge, getPatchNbr(),
            size(FEM, OFFSPRG),  size(MAL, OFFSPRG),
            size(FEM, ADULTS),   size(MAL, ADULTS));
#endif

    if(_currentAge == NONE) warning("Metapop::setPopulation: generation 0 is empty!!\n");
}


// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
/* Metapop::reset
 * resets each Patch, all individuals are moved to the POOL
 */
void Metapop::reset()
{
    for(unsigned int i = 0; i < _patchNbr; ++i) {
        _vPatch[i]->flush();
        _vPatch[i]->reset_ID_individual();
    }
}

// ----------------------------------------------------------------------------------------
// clear
// ----------------------------------------------------------------------------------------
void Metapop::clear()
{
    unsigned int i;
    for(i = 0; i < _vPatch.size(); ++i) {
        delete _vPatch[i];
    }
    _vPatch.clear();

    for(i = 0; i < RecyclingPOOL.size(); ++i) {
        delete RecyclingPOOL[i];
    }
    RecyclingPOOL.clear();

    if(_sample_pops)       {delete[] _sample_pops;       _sample_pops=NULL;}
    if(_pSelection)        {delete _pSelection;          _pSelection=NULL;}
    if(_density_threshold_param) {delete[] _density_threshold_param; _density_threshold_param=NULL;}
    if(_density_threshold) ARRAY::delete_2D(_density_threshold, _density_threshold_nbChanges);
}
// ----------------------------------------------------------------------------------------
// getPatchPtr
// ----------------------------------------------------------------------------------------
Patch* Metapop::getPatchPtr(unsigned int patch)
{
    assert(patch < _patchNbr);
    assert(_vPatch[patch]);
    return _vPatch[patch];
}
// ----------------------------------------------------------------------------------------
// getAllIndividuals
// ----------------------------------------------------------------------------------------
deque<Individual*>* Metapop::getAllIndividuals()
{
    deque<Individual*>* all = new deque<Individual*>;
    unsigned int i, j, sizes;

    for(i = 0; i < _vPatch.size(); ++i) {
        for(j = 0, sizes = _vPatch[i]->size(MAL, ADLTx); j < sizes; ++j){
            all->push_back(_vPatch[i]->get(MAL, ADLTx, j));
        }

        for(j = 0, sizes = _vPatch[i]->size(FEM, ADLTx); j < sizes; ++j){
            all->push_back(_vPatch[i]->get(FEM, ADLTx, j));
        }

        for(j = 0, sizes = _vPatch[i]->size(FEM, OFFSx); j < sizes; ++j){
            all->push_back(_vPatch[i]->get(FEM, OFFSx, j));
        }

        for(j = 0, sizes = _vPatch[i]->size(MAL, OFFSx); j < sizes; ++j){
            all->push_back(_vPatch[i]->get(MAL, OFFSx, j));
        }
    }

    if(all->size() != size()){
        error("Metapop::getAllIndividuals:output not same size as pop (out: %i, pop: %i)\n",all->size(),size());
    }
    return all;
}

// ----------------------------------------------------------------------------------------
// getReplicateCounter
// ----------------------------------------------------------------------------------------
string
Metapop::getReplicateCounter (){
    if (_replicates == 1) return "";
    return STRING::int2str(_currentReplicate, _replicates);
}

// ----------------------------------------------------------------------------------------
string
Metapop::getReplicateCounter_ (){
    if (_replicates == 1) return "";
    return "_" + STRING::int2str(_currentReplicate, _replicates);
}

// ----------------------------------------------------------------------------------------
string
Metapop::getReplicateCounter_r (){
    if (_replicates == 1) return "";
    return "_r" + STRING::int2str(_currentReplicate, _replicates);
}

// ----------------------------------------------------------------------------------------
// getGenerationCounter
// ----------------------------------------------------------------------------------------
string
Metapop::getGenerationCounter (){
    return STRING::int2str(_currentGeneration, _generations);
}
// ----------------------------------------------------------------------------------------
string
Metapop::getGenerationCounter_ (){
    return "_" + getGenerationCounter();
}
// ----------------------------------------------------------------------------------------
string
Metapop::getGenerationCounter_g (){
    return "_g" + getGenerationCounter();
}
// ----------------------------------------------------------------------------------------
/** soft selection (at the patch level) to K */
    void
Metapop::regulate_selection_fitness_patch(age_idx AGE, unsigned int* Kmal, unsigned int* Kfem)
{
    Patch* cur_patch;
    if(Kmal && Kfem){       // the size to regulate to is given
        for(unsigned int i = 0; i < _patchNbr; ++i) {
            cur_patch = _vPatch[i];
            _pSelection->set_fitness(cur_patch, AGE); // compute the fitnesses, but don't sort or make yet the array cumulative
            cur_patch->regulate_selection_fitness(Kfem[i], _pSelection, FEM, AGE);
            cur_patch->regulate_selection_fitness(Kmal[i], _pSelection, MAL, AGE);
        }
    }
    else{                 	// the pops are regulated to carrying capacity
        for(unsigned int i = 0; i < _patchNbr; ++i) {
            cur_patch = _vPatch[i];
            _pSelection->set_fitness(cur_patch, AGE); // compute the fitnesses, but don't sort or make yet the array cumulative
            cur_patch->regulate_selection_fitness(cur_patch->get_KFem(), _pSelection, FEM, AGE);
            cur_patch->regulate_selection_fitness(cur_patch->get_KMal(), _pSelection, MAL, AGE);
        }
    }
}

// ----------------------------------------------------------------------------------------
/** soft selection (at the metapopulation level) to K */
    void
Metapop::regulate_selection_fitness_metapop(age_idx AGE)
{
    // array to buffer the fitness arrays for each patch (female and male)
    double totFitnessM = 0;                              // total fitness
    double totFitnessF = 0;                              // total fitness
    double* sumFitnessM = new double[_patchNbr];         // mean fitness of each patch
    double* sumFitnessF = new double[_patchNbr];         // mean fitness of each patch
    unsigned nbInd[2];
    unsigned nbFem, nbMal, Kmal_tot=0, Kfem_tot=0, i;
    TPatchFitness** fitnessArrays[2];
    for(i=0; i<2; ++i){
        fitnessArrays[i] = new TPatchFitness*[_patchNbr];
        nbInd[i] = 0;
    }

    // for each patch compute the fitness of the individuals
    Patch* cur_patch;
    for(i = 0; i < _patchNbr; ++i) {
        cur_patch = _vPatch[i];

        // check if the patch is populated
        nbMal = cur_patch->size(MAL, AGE);
        nbFem = cur_patch->size(FEM, AGE);
        if(!nbFem && !nbMal){
            fitnessArrays[MAL][i] = NULL;
            fitnessArrays[FEM][i] = NULL;
            continue;
        }

        // get the carrying capacities
        Kmal_tot += cur_patch->get_KMal();
        Kfem_tot += cur_patch->get_KFem();


        // compute and store the fitness of the current patch (0: male, 1: female)
        _pSelection->set_fitness(cur_patch, AGE);
        sumFitnessM[i]   = _pSelection->getSumFitness(MAL);
        totFitnessM      += sumFitnessM[i];
        sumFitnessF[i]   = _pSelection->getSumFitness(FEM);
        totFitnessF      += sumFitnessF[i];
        if(nbMal) fitnessArrays[MAL][i] = _pSelection->get_FitnessObject(MAL, true);
        if(nbFem) fitnessArrays[FEM][i] = _pSelection->get_FitnessObject(FEM, true);
    }

    // regulate each patch separately
    for(i = 0; i < _patchNbr; ++i) {
        cur_patch = _popPtr->getPatch(i);

        // regulate the males
        if(fitnessArrays[MAL][i]){
            _pSelection->set_FitnessObject(MAL, fitnessArrays[MAL][i]);    // set the fitness arrays of the current patch
            cur_patch->regulate_selection_fitness(my_round(sumFitnessM[i]/totFitnessM*Kmal_tot),
                    _pSelection, MAL, AGE);
        }

        // regulate the females
        if(fitnessArrays[FEM][i]){
            _pSelection->set_FitnessObject(FEM, fitnessArrays[FEM][i]);    // set the fitness arrays of the current patch
            cur_patch->regulate_selection_fitness(my_round(sumFitnessF[i]/totFitnessF*Kfem_tot),
                    _pSelection, FEM, AGE);
        }
    }//end_for_nbPatch

    // delete the fitness arrays
    if(sumFitnessF)      delete[] sumFitnessF;
    if(sumFitnessM)      delete[] sumFitnessM;
    for(i=0; i<2; ++i){
        if(fitnessArrays[i]) delete[] fitnessArrays[i];
    }
}

// ----------------------------------------------------------------------------------------
/** hard selection to N (note not to K!),
 * thus a indivdual with a fitness of 0.5 has a survival probability of 0.5
 */
    void
Metapop::regulate_selection_fitness_hard(age_idx AGE)
{
    Patch* cur_patch;
    for(unsigned int i = 0; i < _patchNbr; ++i) {
        cur_patch = _vPatch[i];
        _pSelection->set_fitness(cur_patch, AGE); // compute the fitnesses, but don't sort or make yet the array cumulative
        cur_patch->regulate_selection_fitness(my_round(_pSelection->getSumFitness(MAL)),	_pSelection, MAL, AGE);
        cur_patch->regulate_selection_fitness(my_round(_pSelection->getSumFitness(FEM)), _pSelection, FEM, AGE);
    }
}
// ----------------------------------------------------------------------------------------



