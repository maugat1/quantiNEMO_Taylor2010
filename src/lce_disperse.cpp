/** @file LCEdisperse.cpp
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

#include <stdlib.h>
#include <cmath>
#include "lce_disperse.h"
#include "metapop.h"
#include "individual.h"
#include "stathandler.cpp"
#include "simulation.h"


/******************************************************************************/
//                             ******** LCE_Disperse ********
/******************************************************************************/
LCE_Disperse::LCE_Disperse(int rank) : LCE("disperse","",rank),
    _disperser(NULL), _colonizer(NULL)
{
    // migration/dispersal parameters
    this->add_parameter("dispersal_model",INT2,false,true,0,3,"0");
    this->add_parameter("dispersal_border_model",INT2,false,true,0,2,"0");
    this->add_parameter("dispersal_lattice_range",INT2,false,true,0,1,"0");
    this->add_parameter("dispersal_lattice_dims",MAT, false, false,0,0,"");
    this->add_parameter("dispersal_propagule_prob",DBL,false,true,0,1,"1",true);
    this->add_parameter("dispersal_rate",DBL_MAT,false,false,0,0,"0", true);
    this->add_parameter("dispersal_rate_fem",DBL_MAT,false,false,0,0,"0",true);
    this->add_parameter("dispersal_rate_mal",DBL_MAT,false,false,0,0,"0",true);

    this->add_parameter("dispersal_k_threshold",DBL,false,true,0,1,"0");

    this->add_parameter("dispersal_k_min",DBL,false,false,0,0,"0");
    this->add_parameter("dispersal_k_max",DBL,false,false,0,0,"1");
    this->add_parameter("dispersal_k_max_growth",DBL,false,false,0,0,"-1");
    this->add_parameter("dispersal_k_growth_rate",DBL,false,false,0,0,"1e4");
    this->add_parameter("dispersal_k_symmetry",DBL,false,false,0,0,"1");

    // colonization parameters
    this->add_parameter("colonization_model",INT2,false,true,0,3,"0");
    this->add_parameter("colonization_border_model",INT2,false,true,0,2,"0");
    this->add_parameter("colonization_lattice_range",INT2,false,true,0,1,"0");
    this->add_parameter("colonization_lattice_dims",MAT, false, false,0,0,"");
    this->add_parameter("colonization_propagule_prob",DBL,false,true,0,1,"1",
            true);
    this->add_parameter("colonization_rate",DBL_MAT,false,false,0,0,"0", true);
    this->add_parameter("colonization_rate_fem",DBL_MAT,false,false,0,0,"0",
            true);
    this->add_parameter("colonization_rate_mal",DBL_MAT,false,false,0,0,"0",
            true);

    this->add_parameter("colonization_k_threshold",DBL,false,true,0,1,"0");

    this->add_parameter("colonization_k_min",DBL,false,false,0,0,"0");
    this->add_parameter("colonization_k_max",DBL,false,false,0,0,"1");
    this->add_parameter("colonization_k_max_growth",DBL,false,false,0,0,"-1");
    this->add_parameter("colonization_k_growth_rate",DBL,false,false,0,0,"1e4");
    this->add_parameter("colonization_k_symmetry",DBL,false,false,0,0,"1");
}

// -----------------------------------------------------------------------------
LCE_Disperse::~LCE_Disperse()
{
    if(_disperser) delete _disperser;
    if(_colonizer) delete _colonizer;
}

// -----------------------------------------------------------------------------
// LCE_Disperse::execute
// -----------------------------------------------------------------------------
void LCE_Disperse::execute()
{
#ifdef _DEBUG
    message("  LCE_Disperse ... ");
#endif

    preDispersal();
    if( _colonizer ) {
        _colonizer->execute();
    }
    _disperser->execute();
    postDispersal();

#ifdef _DEBUG
    unsigned int c = 0;
    for(unsigned int i = 0; i < _nb_patch; i++) {
        c += _popPtr->getPatch(i)->nbEmigrant;
    }
    message("done! (nbEmigrants: %i; Individuals(F/M): [%i/%i; %i/%i])\n",
            c,
            _popPtr->size(FEM, OFFSPRG),  _popPtr->size(MAL, OFFSPRG),
            _popPtr->size(FEM, ADULTS),   _popPtr->size(MAL, ADULTS));
#endif
}

// -----------------------------------------------------------------------------
// LCE_Disperse::preDispersal
// -----------------------------------------------------------------------------
/** prepares the dispersal, resets counters,... */
void LCE_Disperse::preDispersal()
{
    Patch* curPatch;
    for(unsigned int i = 0; i < _nb_patch; i++) {
        curPatch = _popPtr->getPatch(i);

        // remove all adult individuals
        _popPtr->getPatch(i)->flush(ADLTx);

        // reset counters
        curPatch->nbKolonisers = curPatch->get_isExtinct() ? 0 : my_NAN;
        curPatch->nbImmigrant  = 0;
        curPatch->nbEmigrant   = 0;
    }
}

// -----------------------------------------------------------------------------
// LCE_Disperse::postDispersal
// -----------------------------------------------------------------------------
/** set the demographic information about colonization and migration of the
 * current patch
 * nbEmigrant:  is not incremented here but in the function _sendEmigrants
 * nbImmigrant: is not incremented here but in the function _sendEmigrants
 */
void LCE_Disperse::postDispersal()
{
    Patch* current_patch;

    // for each patch
    for(unsigned int i = 0; i < _nb_patch; i++) {
        current_patch = _popPtr->getPatch(i);

        // change offspring to adults for each patch
        current_patch->move(OFFSx, ADLTx);

        // was the patch extinct?
        if(current_patch->get_isExtinct()) {
            // set the age since colonization
            current_patch->colnAge = 0;
            // if the patch is now inhabited we have colonizers
            assert(current_patch->nbImmigrant == current_patch->size(ADLTx));
            if(current_patch->nbImmigrant) {
                current_patch->set_isExtinct(false);
                // colonizers are also immigrants
                current_patch->nbKolonisers = current_patch->nbImmigrant;
            }
            continue;  // go to the next patch
        } else {
            // increment age since colonization for extant patches
            ++(current_patch->colnAge);
        }

        // is the patch now extinct?
        if(!current_patch->size(ADLTx)) {
            current_patch->set_isExtinct(true);
            continue;  // go to the next patch
        }

        // compute the number of residents
        current_patch->nbPhilopat = current_patch->size(ADLTx) -
            current_patch->nbImmigrant;
    }
}

// -----------------------------------------------------------------------------
// LCE_Disperse::init
// -----------------------------------------------------------------------------
bool LCE_Disperse::init (Metapop* popPtr)
{
    LCE::init(popPtr);

    // clear the Dispersers
    if(_disperser) delete _disperser;
    if(_colonizer) delete _colonizer;
    _disperser = _colonizer = NULL;
    
    // if only one patch is used there is simply no migration
    _nb_patch = popPtr->getPatchNbr();
    if( _nb_patch == 1 ) {
        _disperser = new Disperser(_paramSet, DISP_NONE);
        return true;
    }

    // create the Dispersers
    if( _paramSet->isSet("colonization_model") ||
            _paramSet->isSet("colonization_border_model") ||
            _paramSet->isSet("colonization_lattice_range") ||
            _paramSet->isSet("colonization_lattice_dims") ||
            _paramSet->isSet("colonization_propagule_prob") ||
            _paramSet->isSet("colonization_rate") ||
            _paramSet->isSet("colonization_rate_fem") ||
            _paramSet->isSet("colonization_rate_mal") ||
            _paramSet->isSet("colonization_k_threshold") ||
            _paramSet->isSet("colonization_k_min") ||
            _paramSet->isSet("colonization_k_max") ||
            _paramSet->isSet("colonization_k_max_growth") ||
            _paramSet->isSet("colonization_k_growth_rate") ||
            _paramSet->isSet("colonization_k_symmetry") )
    {
        _disperser = new Disperser(_paramSet, DISP_MIGR);
        _colonizer = new Disperser(_paramSet, DISP_COLN);
        _colonizer->init(popPtr);
    } else {
        _disperser = new Disperser(_paramSet, DISP_BOTH);
    }
    _disperser->init(popPtr);

    return true;
}

//-----------------------------------------------------------------------------
// LCE_Disperse::executeBeforeEachGeneration
//-----------------------------------------------------------------------------
void LCE_Disperse::executeBeforeEachGeneration(const int& gen)
{
    // temporal parameters
    map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);

    // if it is a temporal paramter
    if(pParam) {
        // check if a change has to be made
        if(_paramSet->updateTemporalParams(gen)) init(_popPtr);
    }
}


/******************************************************************************/
//                             ******** Disperser ********
/******************************************************************************/
Disperser::Disperser(ParamSet *paramSet, DISP_TYPE type) : _disp_model(0),
    _border_model(0), _lattice_range(0), _disp_propagule_prob(0),
    _disp_factor(0)
{
    // set the parameter set pointer
    _paramSet = paramSet;

    // set the dispersal type (none, migr, coln, both)
    _disp_type = type;

    for(int i=0; i<2; ++i) {
        _dispMatrix[i] = NULL;
        _migr_rate[i] = 0;
        _lattice_dims[i] = 0;
    }
}

// -----------------------------------------------------------------------------
Disperser::~Disperser()
{
    if(_dispMatrix[0]) delete _dispMatrix[0];
    if(_dispMatrix[1]) delete _dispMatrix[1];
    if(_disp_factor)   delete[] _disp_factor;
}

// -----------------------------------------------------------------------------
// Disperser::disperse_zeroDispersal
// -----------------------------------------------------------------------------
/** This function is set, when there is no dispersal or if there is only a
 * single population.
 */
void Disperser::disperse_zeroDispersal() {}

// -----------------------------------------------------------------------------
// Disperser::setDispMatrix
// -----------------------------------------------------------------------------
void Disperser::setDispMatrix (TMatrix* mat)
{
    checkDispMatrix(mat);
    if(_dispMatrix[0])	delete _dispMatrix[0];
    if(_dispMatrix[1])  delete _dispMatrix[1];
    _dispMatrix[0] = new TMatrix(*mat);
    _dispMatrix[1] = new TMatrix(*mat);
}

// -----------------------------------------------------------------------------
// Disperser::setDispMatrix
// -----------------------------------------------------------------------------
void Disperser::setDispMatrix (sex_t sex, TMatrix* mat)
{
    checkDispMatrix(mat);
    if(_dispMatrix[sex]) delete _dispMatrix[sex];
    _dispMatrix[sex] = new TMatrix(*mat);
}

// -----------------------------------------------------------------------------
// Disperser::checkDispMatrix
// -----------------------------------------------------------------------------
void Disperser::checkDispMatrix (TMatrix* mat)
{
    double cntr;
    unsigned int i, j;
    mat->get_dims(_lattice_dims);
    if(_lattice_dims[0] != _lattice_dims[1] || _lattice_dims[0] !=  _nb_patch) {
        fatal("The dimension of the dispersal matrix (%ix%i) does not "
                "match the number of patches (%i)!\n", _lattice_dims[0],
                _lattice_dims[1], _nb_patch);
    }
    for(i = 0; i < _lattice_dims[0]; ++i) {
        cntr = 0;
        for(j = 0; j < _lattice_dims[1]; ++j) {
            cntr += mat->get(i,j);
        }
        // (cntr != 1): floating numbers can never be the same
        if(abs(1-cntr) > 1e-4) {
            warning("The elements of row %i of the dispersal matrix do "
                    "not sum to 1!\n", i+1);
        }
    }
}

// -----------------------------------------------------------------------------
// Disperser::execute
// -----------------------------------------------------------------------------
void Disperser::execute()
{
    (this->*doDispersal)();
}

// -----------------------------------------------------------------------------
// Disperser::dispersal_matrix
// -----------------------------------------------------------------------------
/** if the dispersal rates are defined by a matrix */
void Disperser::dispersal_matrix()
{
    Patch *current_patch;
    int nbInd[2];
    double disp[2];
    unsigned int home, target;
    double factor;

    // for each patch
    for(home = 0; home < _nb_patch; ++home) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM, OFFSx);
        nbInd[MAL] = current_patch->size(MAL, OFFSx);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);

        // dispersal to each other patch
        for(target=0; target<_nb_patch; ++target) {
            // if it is the same patch, continue
            if( home == target ) continue;
            disp[FEM] = _dispMatrix[FEM]->get(home,target);
            disp[MAL] = _dispMatrix[MAL]->get(home,target);
            // send the emigrants
            _sendEmigrants(current_patch, home, target, nbInd, disp, factor);
        } // end for each target
    }//end for each home patch
}

// -----------------------------------------------------------------------------
// Disperser::_sendEmigrants
// -----------------------------------------------------------------------------
/** function to send emigrants from the home to the target patch for both sexes.
 * If the migration rate is negative, the emigrating individuals are removed
 * (absorbed)
 */
void Disperser::_sendEmigrants(Patch* curPatch, const int& home,
        const int& target, int* nbInd, double* dispRate, const double& factor)
{
    _sendEmigrants(curPatch, home, target, nbInd[FEM], dispRate[FEM], factor,
            FEM);
    _sendEmigrants(curPatch, home, target, nbInd[MAL], dispRate[MAL], factor,
            MAL);
}

// -----------------------------------------------------------------------------
// Disperser::_sendEmigrants
// -----------------------------------------------------------------------------
/** function to send emigrants form the home to the target patch for a specific
 * sex. If the migration rate is negative, the emigrating individuals are
 * removed (absorbed)
 */
void Disperser::_sendEmigrants(Patch* curPatch, const int& home,
        const int& target, const int& size, const double& dispRate,
        const double& factor, const sex_t& SEX)
{
    double d = dispRate * factor;
    int nb_emigrants;

    // check if dispersal occurs
    if(!d) return;

    // get the number of individuals remaining in the patch
    int curSize = curPatch->size(SEX,OFFSx);
    if(!curSize) return;

    // check for failed colonization/migration (colonists to nonextinct patch
    // or migrants to extinct patch)
    if( d > 0 && (_disp_type && ((_disp_type ^ DISP_MIGR) ^ DISP_COLN)) ) {
        if( (bool)(_disp_type & DISP_COLN) ^
                (bool)(_popPtr->getPatch(target)->get_isExtinct()) )
        { 
            // failed emigrants will be removed
            d = -d;
        }
    }

    // perform migration
    if(d>0) { // normal migration
        // number of emigrants defined
        if(d>=1) nb_emigrants = my_round(d);
        // emigration rate defined
        else nb_emigrants = (int)SimRunner::r.Binomial(d, size);
        // return no emigrants
        if(!nb_emigrants) return;
        // if there are not enough individuals
        if(nb_emigrants > curSize) nb_emigrants = curSize;
        // set the number of emigrants
        curPatch->nbEmigrant += nb_emigrants;
        // set the number of immigrants of the target patch
        _popPtr->getPatch(target)->nbImmigrant += nb_emigrants;
        for(; nb_emigrants; --curSize, --nb_emigrants) {
            _popPtr->move(SEX, OFFSx, home, ADLTx, target,
                    SimRunner::r.Uniform(curSize));
        }
    } else { // emigrants are removed (absorbing boundaries)
        // number of emigrants defined
        if(d<=-1) nb_emigrants = my_round((-1)*d);
        // emigration rate defined
        else nb_emigrants = (int)SimRunner::r.Binomial((-1)*d, size);
        // return no emigrants
        if(!nb_emigrants) return;
        // if there are not enough individuals
        if(nb_emigrants > curSize) nb_emigrants = curSize;
        curPatch->nbEmigrant += nb_emigrants;
        for(; nb_emigrants; --curSize, --nb_emigrants) {
            curPatch->recycle(SEX, OFFSx, SimRunner::r.Uniform(curSize));
        }
    }
}

// -----------------------------------------------------------------------------
// Disperser::disperse_island
// -----------------------------------------------------------------------------
/** island migration model with a female and a male migration rate
 * migration rate to any patch is m/(_nb_patch-1)
 */
void Disperser::disperse_island()
{
    Patch *current_patch;
    int nbInd[2];
    unsigned int home, target;
    double factor;

    // for each patch
    for(home = 0; home < _nb_patch; ++home) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);

        // dispersal to each other patch
        for(target=0; target<_nb_patch; ++target) {
            // if it is the same patch, continue
            if( home == target ) continue;
            // send emigrants
            _sendEmigrants(current_patch, home, target, nbInd, _migr_rate,
                    factor);
        } // end for each target
    }//end for each home patch
}

// -----------------------------------------------------------------------------
// Disperser::disperse_island_propagule
// -----------------------------------------------------------------------------
/** island migration model with a female and a male migration rate
 * migration rate to any patch is m/(_nb_patch-1)
 */
void Disperser::disperse_island_propagule()
{
    Patch *current_patch;
    int nbInd[2];
    unsigned int home, target, propagule;
    double factor;

    // for each patch
    for(home = 0; home < _nb_patch; ++home) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);

        // dispersal to the propagule patch
        do{ // draw the propagule patch randomly
            propagule = SimRunner::r.Uniform(_nb_patch);
        } while( home == propagule );
        _sendEmigrants(current_patch, home, propagule, nbInd,
                _migr_rate_propagule, factor);

        // migration to each other patch
        for( target=0; target<_nb_patch; ++target ) {
            // if it is the same patch, continue
            if( home == target ) continue;
            // send emigrants
            _sendEmigrants(current_patch, home, target, nbInd, _migr_rate,
                    factor);
        } // end for each target
    }//end for each home patch
}

// -----------------------------------------------------------------------------
// Disperser::disperse_1D_stepping_stone
// -----------------------------------------------------------------------------
/** island migration model with a female and a male migration rate
 * migration rate to any patch is m/2, except for the edge patches
 */
void Disperser::disperse_1D_stepping_stone()
{
    Patch *current_patch;
    int nbInd[2];
    double factor;

    // for each patch in the middle: 1-_nb_patch-2
    for(unsigned int home = 1; home < _nb_patch-1; ++home) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL]= current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home-1, nbInd, _migr_rate, factor);   // migrate left
        _sendEmigrants(current_patch, home, home+1, nbInd, _migr_rate, factor);   // migrate right
        current_patch->move(OFFSx,ADLTx);  // not migrated individuals remain in the original patch
    }//end for each home patch

    // edge patches (first patch)
    current_patch = _popPtr->getPatch(0);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL]= current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, 0, _nb_patch-1, nbInd, _migr_rateOut, factor);   // migrate left
    _sendEmigrants(current_patch, 0, 1,           nbInd, _migr_rateIn, factor);    // migrate right
    current_patch->move(OFFSx,ADLTx);  // not migrated individuals remain in the original patch

    // edge patches (last patch)
    current_patch = _popPtr->getPatch(_nb_patch-1); // last patch
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL]= current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-2, nbInd, _migr_rateIn, factor);   // migrate left
    _sendEmigrants(current_patch, _nb_patch-1, 0,           nbInd, _migr_rateOut, factor);  // migrate right
    current_patch->move(OFFSx,ADLTx);  // not migrated individuals remain in the original patch
}

// -----------------------------------------------------------------------------
// Disperser::_disperse_2D_SS_middle_4Neighbour
// -----------------------------------------------------------------------------
/** 2D stepping stone migration of the middle patches (4 neighbours)*/
void Disperser::_disperse_2D_SS_middle_4Neighbour()
{
    Patch *current_patch;
    int nbInd[2];
    unsigned int r, c, home;
    unsigned int x = _lattice_dims[0];
    unsigned int y = _lattice_dims[1];
    double factor;

    // for each patch in the middle
    for(c = 1; c < y-1; ++c) {    // for each column
        for(r = 1; r < x-1; ++r) {    // for each row
            home = c*x+r;
            current_patch = _popPtr->getPatch(home);
            nbInd[FEM]    = current_patch->size(FEM,OFFSPRG);
            nbInd[MAL]    = current_patch->size(MAL,OFFSPRG);
            factor        = (this->*get_disp_factor_funcPtr)(current_patch);

            // migrate
            _sendEmigrants(current_patch, home, home-1, nbInd, _migr_rate, factor); // left
            _sendEmigrants(current_patch, home, home+1, nbInd, _migr_rate, factor); // right
            _sendEmigrants(current_patch, home, home-x, nbInd, _migr_rate, factor); // up
            _sendEmigrants(current_patch, home, home+x, nbInd, _migr_rate, factor); // down
        }
    }//end for each home patch
}

// -----------------------------------------------------------------------------
// Disperser::_disperse_2D_SS_middle_8Neighbour
// -----------------------------------------------------------------------------
/** 2D stepping stone migration of the middle patches (8 neighbours)*/
void Disperser::_disperse_2D_SS_middle_8Neighbour()
{
    Patch *current_patch;
    int nbInd[2];
    unsigned int r, c, home;
    unsigned int x = _lattice_dims[0];
    unsigned int y = _lattice_dims[1];
    double factor;

    // for each patch in the middle
    for(c = 1; c < y-1; ++c) {    // for each column
        for(r = 1; r < x-1; ++r) {    // for each row
            home = c*x+r;
            current_patch = _popPtr->getPatch(home);
            nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
            nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
            factor     = (this->*get_disp_factor_funcPtr)(current_patch);

            // migrate
            _sendEmigrants(current_patch, home, home+1,   nbInd, _migr_rate, factor); // right
            _sendEmigrants(current_patch, home, home+1+x, nbInd, _migr_rate, factor); // right down
            _sendEmigrants(current_patch, home, home+x,   nbInd, _migr_rate, factor); // down
            _sendEmigrants(current_patch, home, home-1+x, nbInd, _migr_rate, factor); // left down
            _sendEmigrants(current_patch, home, home-1,   nbInd, _migr_rate, factor); // left
            _sendEmigrants(current_patch, home, home-1-x, nbInd, _migr_rate, factor); // left up
            _sendEmigrants(current_patch, home, home-x,   nbInd, _migr_rate, factor); // up
            _sendEmigrants(current_patch, home, home+1-x, nbInd, _migr_rate, factor); // right up
        }
    }//end for each home patch
}

// -----------------------------------------------------------------------------
// Disperser::_disperse_2D_SS_edge_4Neighbour
// -----------------------------------------------------------------------------
void Disperser::_disperse_2D_SS_edge_4Neighbour()
{
    Patch* current_patch;
    int nbInd[2];
    unsigned int x = _lattice_dims[0];
    unsigned int home;
    double factor;

    // upper edge
    for(home = 1; home <= x-2; ++home) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home+1,           nbInd, _migr_rateIn, factor);  // right
        _sendEmigrants(current_patch, home, home+x,           nbInd, _migr_rateIn, factor);  // down
        _sendEmigrants(current_patch, home, home-1,           nbInd, _migr_rateIn, factor);  // left
        _sendEmigrants(current_patch, home, _nb_patch+home-x, nbInd, _migr_rateOut, factor); // up
    }

    // right edge
    for(home = (2*x)-1; home <= _nb_patch-x-1; home+=x) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home-x+1,         nbInd, _migr_rateOut, factor); // right
        _sendEmigrants(current_patch, home, home+x,           nbInd, _migr_rateIn, factor);  // down
        _sendEmigrants(current_patch, home, home-1,           nbInd, _migr_rateIn, factor);  // left
        _sendEmigrants(current_patch, home, home-x,           nbInd, _migr_rateIn, factor);  // up
    }

    // lower edge
    for(home = _nb_patch-x+1; home <= _nb_patch-2; ++home) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home+1,           nbInd, _migr_rateIn, factor);  // right
        _sendEmigrants(current_patch, home, home-_nb_patch+x, nbInd, _migr_rateOut, factor); // down
        _sendEmigrants(current_patch, home, home-1,           nbInd, _migr_rateIn, factor);  // left
        _sendEmigrants(current_patch, home, home-x,           nbInd, _migr_rateIn, factor);  // up
    }

    // left edge
    for(home = x; home <= _nb_patch-(2*x); home+=x) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home+1,           nbInd, _migr_rateIn, factor);  // right
        _sendEmigrants(current_patch, home, home+x,           nbInd, _migr_rateIn, factor);  // down
        _sendEmigrants(current_patch, home, home+x-1,         nbInd, _migr_rateOut, factor); // left
        _sendEmigrants(current_patch, home, home-x,           nbInd, _migr_rateIn, factor);  // up
    }
}

// -----------------------------------------------------------------------------
// Disperser::_disperse_2D_SS_edge_8Neighbour
// -----------------------------------------------------------------------------
void Disperser::_disperse_2D_SS_edge_8Neighbour()
{
    Patch* current_patch;
    int nbInd[2];
    unsigned int x = _lattice_dims[0];
    unsigned int home;
    double factor;

    // upper edge
    for(home = 1; home <= x-2; ++home) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home+1,             nbInd, _migr_rateIn, factor);  // right
        _sendEmigrants(current_patch, home, home+x+1,           nbInd, _migr_rateIn, factor);  // right down
        _sendEmigrants(current_patch, home, home+x,             nbInd, _migr_rateIn, factor);  // down
        _sendEmigrants(current_patch, home, home+x-1,           nbInd, _migr_rateIn, factor);  // left down
        _sendEmigrants(current_patch, home, home-1,             nbInd, _migr_rateIn, factor);  // left
        _sendEmigrants(current_patch, home, _nb_patch+home-x-1, nbInd, _migr_rateOut, factor); // left up
        _sendEmigrants(current_patch, home, _nb_patch+home-x,   nbInd, _migr_rateOut, factor); // up
        _sendEmigrants(current_patch, home, _nb_patch+home-x+1, nbInd, _migr_rateOut, factor); // right up
    }

    // right edge
    for(home = (2*x)-1; home <= _nb_patch-x-1; home+=x) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home-x+1,           nbInd, _migr_rateOut, factor); // right
        _sendEmigrants(current_patch, home, home+1,             nbInd, _migr_rateOut, factor); // right down
        _sendEmigrants(current_patch, home, home+x,             nbInd, _migr_rateIn, factor);  // down
        _sendEmigrants(current_patch, home, home+x-1,           nbInd, _migr_rateIn, factor);  // left down
        _sendEmigrants(current_patch, home, home-1,             nbInd, _migr_rateIn, factor);  // left
        _sendEmigrants(current_patch, home, home-x-1,           nbInd, _migr_rateIn, factor);  // left up
        _sendEmigrants(current_patch, home, home-x,             nbInd, _migr_rateIn, factor);  // up
        _sendEmigrants(current_patch, home, home-(2*x)+1,       nbInd, _migr_rateOut, factor); // right up
    }

    // lower edge
    for(home = _nb_patch-x+1; home <= _nb_patch-2; ++home) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home+1,             nbInd, _migr_rateIn, factor);  // right
        _sendEmigrants(current_patch, home, home-_nb_patch+x+1, nbInd, _migr_rateOut, factor); // right down
        _sendEmigrants(current_patch, home, home-_nb_patch+x,   nbInd, _migr_rateOut, factor); // down
        _sendEmigrants(current_patch, home, home-_nb_patch+x-1, nbInd, _migr_rateOut, factor); // left down
        _sendEmigrants(current_patch, home, home-1,             nbInd, _migr_rateIn, factor);  // left
        _sendEmigrants(current_patch, home, home-x-1,           nbInd, _migr_rateIn, factor);  // left up
        _sendEmigrants(current_patch, home, home-x,             nbInd, _migr_rateIn, factor);  // up
        _sendEmigrants(current_patch, home, home-x+1,           nbInd, _migr_rateIn, factor);  // right up
    }

    // left edge
    for(home = x; home <= _nb_patch-(2*x); home+=x) {
        current_patch = _popPtr->getPatch(home);
        nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
        nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
        factor = (this->*get_disp_factor_funcPtr)(current_patch);
        _sendEmigrants(current_patch, home, home+1,             nbInd, _migr_rateIn, factor);  // right
        _sendEmigrants(current_patch, home, home+x+1,           nbInd, _migr_rateIn, factor);  // right down
        _sendEmigrants(current_patch, home, home+x,             nbInd, _migr_rateIn, factor);  // down
        _sendEmigrants(current_patch, home, home+(2*x)-1,       nbInd, _migr_rateOut, factor); // left down
        _sendEmigrants(current_patch, home, home+x-1,           nbInd, _migr_rateOut, factor); // left
        _sendEmigrants(current_patch, home, home-1,             nbInd, _migr_rateOut, factor); // left up
        _sendEmigrants(current_patch, home, home-x,             nbInd, _migr_rateIn, factor);  // up
        _sendEmigrants(current_patch, home, home-x+1,           nbInd, _migr_rateIn, factor);  // right up
    }
}

// -----------------------------------------------------------------------------
// Disperser::_disperse_2D_SS_corner_4Neighbour
// -----------------------------------------------------------------------------
void Disperser::_disperse_2D_SS_corner_4Neighbour()
{
    Patch* current_patch;
    int nbInd[2];
    unsigned int x = _lattice_dims[0];
    double factor;

    // first corner
    current_patch = _popPtr->getPatch(0);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, 0, 1,           nbInd, _migr_rateCorner, factor); // right
    _sendEmigrants(current_patch, 0, x,           nbInd, _migr_rateCorner, factor); // down
    _sendEmigrants(current_patch, 0, x-1,         nbInd, _migr_rateOut, factor);    // left
    _sendEmigrants(current_patch, 0, _nb_patch-x, nbInd, _migr_rateOut, factor);    // up

    // second corner
    current_patch = _popPtr->getPatch(x-1);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, x-1, 0,           nbInd, _migr_rateOut, factor);    // right
    _sendEmigrants(current_patch, x-1, 2*x-1,       nbInd, _migr_rateCorner, factor); // down
    _sendEmigrants(current_patch, x-1, x-2,         nbInd, _migr_rateCorner, factor); // left
    _sendEmigrants(current_patch, x-1, _nb_patch-1, nbInd, _migr_rateOut, factor);    // up

    // third corner
    current_patch = _popPtr->getPatch(_nb_patch-x);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, _nb_patch-x, _nb_patch-x+1,   nbInd, _migr_rateCorner, factor); // right
    _sendEmigrants(current_patch, _nb_patch-x, 0,               nbInd, _migr_rateOut, factor);    // down
    _sendEmigrants(current_patch, _nb_patch-x, _nb_patch-1,     nbInd, _migr_rateOut, factor);    // left
    _sendEmigrants(current_patch, _nb_patch-x, _nb_patch-(2*x), nbInd, _migr_rateCorner, factor); // up

    // forth corner
    current_patch = _popPtr->getPatch(_nb_patch-1);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-x,   nbInd, _migr_rateOut, factor);    // right
    _sendEmigrants(current_patch, _nb_patch-1, x-1,           nbInd, _migr_rateOut, factor);    // down
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-2,   nbInd, _migr_rateCorner, factor); // left
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-x-1, nbInd, _migr_rateCorner, factor); // up
}

// -----------------------------------------------------------------------------
// Disperser::_disperse_2D_SS_corner_8Neighbour
// -----------------------------------------------------------------------------
void Disperser::_disperse_2D_SS_corner_8Neighbour()
{
    Patch* current_patch;
    int nbInd[2];
    unsigned int x = _lattice_dims[0];
    double factor;

    // first corner
    current_patch = _popPtr->getPatch(0);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, 0, 1,             nbInd, _migr_rateCorner, factor); // right
    _sendEmigrants(current_patch, 0, x+1,           nbInd, _migr_rateCorner, factor); // right down
    _sendEmigrants(current_patch, 0, x,             nbInd, _migr_rateCorner, factor); // down
    _sendEmigrants(current_patch, 0, 2*x-1,         nbInd, _migr_rateOut, factor);    // left down
    _sendEmigrants(current_patch, 0, x-1,           nbInd, _migr_rateOut, factor);    // left
    _sendEmigrants(current_patch, 0, _nb_patch-1,   nbInd, _migr_rateOut, factor);    // left up
    _sendEmigrants(current_patch, 0, _nb_patch-x,   nbInd, _migr_rateOut, factor);    // up
    _sendEmigrants(current_patch, 0, _nb_patch-x+1, nbInd, _migr_rateOut, factor);    // right up

    // second corner
    current_patch = _popPtr->getPatch(x-1);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, x-1, 0,           nbInd, _migr_rateOut, factor);    // right
    _sendEmigrants(current_patch, x-1, x,           nbInd, _migr_rateOut, factor);    // right down
    _sendEmigrants(current_patch, x-1, 2*x-1,       nbInd, _migr_rateCorner, factor); // down
    _sendEmigrants(current_patch, x-1, 2*x-2,       nbInd, _migr_rateCorner, factor); // left down
    _sendEmigrants(current_patch, x-1, x-2,         nbInd, _migr_rateCorner, factor); // left
    _sendEmigrants(current_patch, x-1, _nb_patch-2, nbInd, _migr_rateOut, factor);    // left up
    _sendEmigrants(current_patch, x-1, _nb_patch-1, nbInd, _migr_rateOut, factor);    // up
    _sendEmigrants(current_patch, x-1, _nb_patch-x, nbInd, _migr_rateOut, factor);    // right up

    // third corner
    current_patch = _popPtr->getPatch(_nb_patch-x);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, _nb_patch-x, _nb_patch-x+1,     nbInd, _migr_rateCorner, factor); // right
    _sendEmigrants(current_patch, _nb_patch-x, 1,                 nbInd, _migr_rateOut, factor);    // right down
    _sendEmigrants(current_patch, _nb_patch-x, 0,                 nbInd, _migr_rateOut, factor);    // down
    _sendEmigrants(current_patch, _nb_patch-x, x-1,               nbInd, _migr_rateOut, factor);    // left down
    _sendEmigrants(current_patch, _nb_patch-x, _nb_patch-1,       nbInd, _migr_rateOut, factor);    // left
    _sendEmigrants(current_patch, _nb_patch-x, _nb_patch-x-1,     nbInd, _migr_rateOut, factor);    // left up
    _sendEmigrants(current_patch, _nb_patch-x, _nb_patch-(2*x),   nbInd, _migr_rateCorner, factor); // up
    _sendEmigrants(current_patch, _nb_patch-x, _nb_patch-(2*x)+1, nbInd, _migr_rateCorner, factor); // rigth up

    // fourth corner
    current_patch = _popPtr->getPatch(_nb_patch-1);
    nbInd[FEM] = current_patch->size(FEM,OFFSPRG);
    nbInd[MAL] = current_patch->size(MAL,OFFSPRG);
    factor = (this->*get_disp_factor_funcPtr)(current_patch);
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-x,     nbInd, _migr_rateOut, factor);    // right
    _sendEmigrants(current_patch, _nb_patch-1, 0,               nbInd, _migr_rateOut, factor);    // right down
    _sendEmigrants(current_patch, _nb_patch-1, x-1,             nbInd, _migr_rateOut, factor);    // down
    _sendEmigrants(current_patch, _nb_patch-1, x-2,             nbInd, _migr_rateOut, factor);    // left down
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-2,     nbInd, _migr_rateCorner, factor); // left
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-x-2,   nbInd, _migr_rateCorner, factor); // left up
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-x-1,   nbInd, _migr_rateCorner, factor); // up
    _sendEmigrants(current_patch, _nb_patch-1, _nb_patch-(2*x), nbInd, _migr_rateOut, factor);    // right up
}

// -----------------------------------------------------------------------------
// Disperser::disperse_2D_stepping_stone_4Neighbour
// -----------------------------------------------------------------------------
void Disperser::disperse_2D_stepping_stone_4Neighbour()
{
    // perform migration
    _disperse_2D_SS_middle_4Neighbour();
    _disperse_2D_SS_edge_4Neighbour();
    _disperse_2D_SS_corner_4Neighbour();

    // all not migrated individuals remain in the original patch
    for(unsigned int home = 0; home < _nb_patch; ++home) {
        _popPtr->getPatch(home)->move(OFFSx,ADLTx);
    }
}

// -----------------------------------------------------------------------------
// Disperser::disperse_2D_stepping_stone_8Neighbour
// -----------------------------------------------------------------------------
void Disperser::disperse_2D_stepping_stone_8Neighbour()
{
    // perform migration
    _disperse_2D_SS_middle_8Neighbour();
    _disperse_2D_SS_edge_8Neighbour();
    _disperse_2D_SS_corner_8Neighbour();

    // all not migrated individuals remain in the original patch
    for(unsigned int home = 0; home < _nb_patch; ++home) {
        _popPtr->getPatch(home)->move(OFFSx,ADLTx);
    }
}

// -----------------------------------------------------------------------------
// Disperser::init
// -----------------------------------------------------------------------------
bool Disperser::init (Metapop* popPtr)
{
    _popPtr = popPtr;
    _nb_patch = popPtr->getPatchNbr();

    // if no dispersal type
    if( _disp_type == DISP_NONE ) {
        doDispersal = &Disperser::disperse_zeroDispersal;
        return true;
    }

    //reset matrixes
    if( _dispMatrix[FEM] ) {delete _dispMatrix[FEM]; _dispMatrix[FEM]=NULL;}
    if( _dispMatrix[MAL] ) {delete _dispMatrix[MAL]; _dispMatrix[MAL]=NULL;}

    // get parameters
    if( _disp_type & DISP_MIGR ) {
        _disp_model = (int)_paramSet->getValue("dispersal_model");
        _disp_propagule_prob = _paramSet->getValue("dispersal_propagule_prob");
        _border_model = (int)_paramSet->getValue("dispersal_border_model");
        _lattice_range = (int)_paramSet->getValue("dispersal_lattice_range");
    } else {
        _disp_model = (int)_paramSet->getValue("colonization_model");
        _disp_propagule_prob =
            _paramSet->getValue("colonization_propagule_prob");
        _border_model = (int)_paramSet->getValue("colonization_border_model");
        _lattice_range = (int)_paramSet->getValue("colonization_lattice_range");
    }

    // migration is set by matrixes
    if( ((_disp_type & DISP_MIGR) && 
                ( (_paramSet->isSet("dispersal_rate")
                   || _paramSet->isSet("dispersal_rate_fem")
                   || _paramSet->isSet("dispersal_rate_mal"))
                  && (_paramSet->isMatrix("dispersal_rate")
                      || _paramSet->isMatrix("dispersal_rate_fem")
                      || _paramSet->isMatrix("dispersal_rate_mal")) ) )
            || ((_disp_type & DISP_COLN) &&
                ( (_paramSet->isSet("colonization_rate")
                   || _paramSet->isSet("colonization_rate_fem")
                   || _paramSet->isSet("colonization_rate_mal"))
                  && (_paramSet->isMatrix("colonization_rate")
                      || _paramSet->isMatrix("colonization_rate_fem")
                      || _paramSet->isMatrix("colonization_rate_mal")) ) ) )
    {
        setDispersalMatrix();
    } else { // migration is set by a single dispersal rate (or default)
        // if it is a 2D stepping stone model
        if( _disp_model == 3 ) _get_lattice_dims();
        setDispersalRate();
    }
    _setDispersalFactor();

    return true;
}

// -----------------------------------------------------------------------------
// Disperser::setDispersalFactor
// -----------------------------------------------------------------------------
/** set the genralized logistic function parameters
 * check if the function is really needed, i.e. if there is really a change
 */
void Disperser::_setDispersalFactor()
{
    if(_disp_factor) {
        delete[] _disp_factor; _disp_factor=NULL;
    }
    // get the values
    _disp_factor = new double[5];
    if( _disp_type & DISP_MIGR ) {
        _disp_factor[0] = _paramSet->getValue("dispersal_k_min");
        _disp_factor[1] = _paramSet->getValue("dispersal_k_max");
        _disp_factor[2] = _paramSet->getValue("dispersal_k_max_growth");
        _disp_factor[3] = _paramSet->getValue("dispersal_k_growth_rate");
        _disp_factor[4] = _paramSet->getValue("dispersal_k_symmetry");
        if(!_disp_factor) {
            fatal("Parameter 'dispersal_k_symmetry' cannot be zero!\n");
        }
    } else {
        _disp_factor[0] = _paramSet->getValue("colonization_k_min");
        _disp_factor[1] = _paramSet->getValue("colonization_k_max");
        _disp_factor[2] = _paramSet->getValue("colonization_k_max_growth");
        _disp_factor[3] = _paramSet->getValue("colonization_k_growth_rate");
        _disp_factor[4] = _paramSet->getValue("colonization_k_symmetry");
        if(!_disp_factor) {
            fatal("Parameter 'colonization_k_symmetry' cannot be zero!\n");
        }
    }
    // if the growth rate is too large/small, the logsitic function has a
    // size problem (exp) change to threshold function
    if( abs(_disp_factor[3]) >= STRING::str2int<double>(
                _paramSet->get_param(
                    (_disp_type & DISP_MIGR ) ?
                    "dispersal_k_growth_rate" :
                    "colonization_k_growth_rate"
                    )->get_default_arg()) )
    {
        if(_disp_factor[3]<0) { // swap min and max if slope is negative
            double temp = _disp_factor[0];
            _disp_factor[0] = _disp_factor[1];
            _disp_factor[1] = temp;
        }
        if(_disp_factor[2]<0) { // where is the threshold?
            if(abs(_disp_factor[1]-1)<1e-6) {
                get_disp_factor_funcPtr =
                    &Disperser::get_disp_factor_one;
                return;
            } else {
                get_disp_factor_funcPtr =
                    &Disperser::get_disp_factor_max;
                return;
            }
        }
        get_disp_factor_funcPtr =
            &Disperser::get_disp_factor_k_threshold;
        return;
    }
    // test if there is a change between 0 and 10 (populations may exceed
    // K...)
    double ten = generalLogisticCurve(10,_disp_factor[0],_disp_factor[1],
            _disp_factor[2],_disp_factor[3],_disp_factor[4]);
    double zero = generalLogisticCurve(0, _disp_factor[0],_disp_factor[1],
            _disp_factor[2],_disp_factor[3],_disp_factor[4]);
    if(abs(ten-zero)<1e-6) { // factor not influenced by pop density
        if(abs(ten-1)<1e-6) { // factor = 1
            get_disp_factor_funcPtr = &Disperser::get_disp_factor_one;
            return;
        }
        if(abs(ten-_disp_factor[0])<1e-6) { // factor = min
            get_disp_factor_funcPtr = &Disperser::get_disp_factor_min;
            return;
        }
        if(abs(ten-_disp_factor[1])<1e-6) { // factor = max
            get_disp_factor_funcPtr = &Disperser::get_disp_factor_max;
            return;
        }
    }
    // use the full general logistic function
    get_disp_factor_funcPtr = &Disperser::get_disp_factor_k_logistic;
}

// -----------------------------------------------------------------------------
// Disperser::_get_lattice_dims
// -----------------------------------------------------------------------------
/** get the dimensions of the 2D stepping stone matrix */
void Disperser::_get_lattice_dims()
{
    if( _disp_type & DISP_MIGR ?  _paramSet->isSet("dispersal_lattice_dims") :
            _paramSet->isSet("colonization_lattice_dims") )
    {
        TMatrix* m = (_disp_type & DISP_MIGR ?
                _paramSet->getMatrix("dispersal_lattice_dims") :
                _paramSet->getMatrix("colonization_lattice_dims"));
        // check if the dimensions of the matrix are correct
        if(m->get_dims(NULL) != 2) {
            fatal(_disp_type & DISP_MIGR ? 
                    "The parameter dispersal_lattice_dims should have a "
                    "matrix with two values!\n" :
                    "The parameter colonization_lattice_dims should have a "
                    "matrix with two values!\n");
        }
        // get the dimension of the lattice
        _lattice_dims[0] = (unsigned int)m->get(0,0);
        _lattice_dims[1] = (unsigned int)m->get(0,1);
        if(_lattice_dims[0]*_lattice_dims[1] != _nb_patch) {
            if( _disp_type & DISP_MIGR ) {
                fatal("Parameter dispersal_lattice_dims: The dimension of the "
                        "lattice (%ix%i) does not mach the number of patches "
                        "(%i)!\n", _lattice_dims[0], _lattice_dims[1],
                        _nb_patch);
            } else {
                fatal("Parameter colonization_lattice_dims: The dimension of "
                        "the lattice (%ix%i) does not mach the number of "
                        "patches (%i)!\n", _lattice_dims[0], _lattice_dims[1],
                        _nb_patch);
            }
        }
    } else {
        // if the dimensions are not set, we assume that x=y
        _lattice_dims[0] = _lattice_dims[1] = 
            (unsigned int) sqrt((double)_nb_patch);
        if(_lattice_dims[0]*_lattice_dims[1] != _nb_patch) {
            if( _disp_type & DISP_MIGR ) {
                fatal("Parameter dispersal_lattice_dims: The dimension of the "
                        "lattice (%ix%i) does not mach the number of patches "
                        "(%i)!\n", _lattice_dims[0], _lattice_dims[1],
                        _nb_patch);
            } else {
                fatal("Parameter colonization_lattice_dims: The dimension of "
                        "the lattice (%ix%i) does not mach the number of "
                        "patches (%i)!\n", _lattice_dims[0], _lattice_dims[1],
                        _nb_patch);
            }
        }
    }
    // check if it is realy a 2D and not a 1D lattice
    if(_lattice_dims[0] == 1 || _lattice_dims[1]==1) {
        _disp_model = 2; // use a 1D stepping stone model
    }
}

// -----------------------------------------------------------------------------
// Disperser::setDispersalMatrix
// -----------------------------------------------------------------------------
/** dispersal is set by a dispersal matrix
 * sex specific matrixes have precendence
 * if only one sex specific matrix is set, an error is drawn
 */
void Disperser::setDispersalMatrix()
{
    // set the dispersal matrix
    if(_paramSet->isSet(_disp_type & DISP_MIGR ? "dispersal_rate_fem" :
                "colonization_rate_fem")) {
        if(_paramSet->isSet(_disp_type & DISP_MIGR ? "dispersal_rate_fem" :
                    "colonization_rate_fem") &&
                _paramSet->isSet(_disp_type & DISP_MIGR ? "dispersal_rate_mal" :
                    "colonization_rate_mal"))
        {
            setDispMatrix(FEM, _paramSet->getMatrix(_disp_type & DISP_MIGR ?
                        "dispersal_rate_fem" : "colonization_rate_fem"));
            setDispMatrix(MAL, _paramSet->getMatrix(_disp_type & DISP_MIGR ?
                        "dispersal_rate_mal" : "colonization_rate_mal"));
        } else {
            fatal(_disp_type & DISP_MIGR ?
                    "Only one sex specific migration matrix is specified: "
                    "both are required!\n" :
                    "Only one sex specific colonization matrix is specified: "
                    "both are required!\n");
        }
    } else {
        // general dispersal matrix
        setDispMatrix(_paramSet->getMatrix(_disp_type & DISP_MIGR ?
                    "dispersal_rate" : "colonization_rate"));
    }
    // function pointer to the dispersal function
    doDispersal = &Disperser::dispersal_matrix;
}

// -----------------------------------------------------------------------------
// Disperser::setDispersalRate
// -----------------------------------------------------------------------------
/** dispersal is set by a dispersal rate
 * sex specific rates have precendence
 * if only one sex specific rate is set, an error is drawn
 */
void Disperser::setDispersalRate()
{
    if(_paramSet->isSet(_disp_type & DISP_MIGR ? "dispersal_rate_fem" :
                "colonization_rate_fem"))
    {
        if(_paramSet->isSet(_disp_type & DISP_MIGR ? "dispersal_rate_fem" :
                    "colonization_rate_fem")
                && _paramSet->isSet(_disp_type & DISP_MIGR ?
                    "dispersal_rate_mal" : "colonization_rate_mal"))
        {
            _migr_rate[FEM] = _paramSet->getValue(_disp_type & DISP_MIGR ?
                    "dispersal_rate_fem" : "colonization_rate_fem");
            _migr_rate[MAL] = _paramSet->getValue(_disp_type & DISP_MIGR ?
                    "dispersal_rate_mal" : "colonization_rate_mal");
        } else {
            fatal(_disp_type & DISP_MIGR ?
                    "Only one sex specific migration rate is specified: "
                    "both are required!\n" :
                    "Only one sex specific colonization rate is specified: "
                    "both are required!\n");
        }
    } else {
        // general dispersal rate
        _migr_rate[FEM] = _migr_rate[MAL] = 
            _paramSet->getValue(_disp_type & DISP_MIGR ? "dispersal_rate": 
                    "colonization_rate");
    }
    // if there is no migration
    if(!_migr_rate[MAL] && !_migr_rate[FEM]) {
        doDispersal = &Disperser::disperse_zeroDispersal;
        return;
    }
    switch (_disp_model) {
        case 0: // island migration model
            _migr_rate[FEM] /= _nb_patch-1;
            _migr_rate[MAL] /= _nb_patch-1;
            doDispersal = &Disperser::disperse_island;
            break;
        case 1: // island migration model with propagule pool
            // if there are only two populations use the island simple model
            if(_nb_patch == 2) {
                _disp_model = 0;
                warning("With only 2 populations it is not possible to "
                        "run a propagule dispersal model: the island "
                        "model is used!\n");
                return setDispersalRate();
            }
            _migr_rate_propagule[FEM] = 
                _migr_rate[FEM] * _disp_propagule_prob;
            _migr_rate_propagule[MAL] = 
                _migr_rate[MAL] * _disp_propagule_prob;
            _migr_rate[FEM] *= (1.0-_disp_propagule_prob)/(_nb_patch-2);
            _migr_rate[MAL] *= (1.0-_disp_propagule_prob)/(_nb_patch-2);
            doDispersal = &Disperser::disperse_island_propagule;
            break;
        case 2: // 1D steppings stone model
            int factorIn, factorOut;
            doDispersal = &Disperser::disperse_1D_stepping_stone;
            // edge effect
            switch(_border_model) {
                default:
                case 0: // circle
                    factorIn = 2; factorOut = 2;
                    break;
                case 1: // reflecting
                    factorIn = 1; factorOut = 0;
                    break;
                case 2: // absorbing
                    // negative number means removing individuals
                    factorIn = 2; factorOut = -2;
                    break;
            }
            _migr_rateIn[FEM] = _migr_rate[FEM]/factorIn;
            _migr_rateIn[MAL] = _migr_rate[MAL]/factorIn;
            _migr_rateOut[FEM] = factorOut ? _migr_rate[FEM]/factorOut : 0;
            _migr_rateOut[MAL] = factorOut ? _migr_rate[MAL]/factorOut : 0;
            _migr_rate[FEM] /= 2;
            _migr_rate[MAL] /= 2;
            break;
        case 3: // 2D stepping stone
            int factorIn4, factorOut4, factorInCorner4;
            int factorIn8, factorOut8, factorInCorner8;
            // edge effect
            switch(_border_model) {
                default:
                case 0: // torus
                    factorIn4  = 4; factorIn8  = 8;
                    factorOut4 = 4; factorOut8 = 8;
                    factorInCorner4 = 4; factorInCorner8 = 8;
                    break;
                case 1: // reflecting
                    factorIn4  = 3; factorIn8  = 5;
                    factorOut4 = 0; factorOut8 = 0;
                    factorInCorner4 = 2; factorInCorner8 = 3;
                    break;
                case 2: // absorbing
                    factorIn4  = 4; factorIn8  = 8;
                    // negative number means removing individuals
                    factorOut4 = -4; factorOut8 = -8;
                    factorInCorner4 = 4; factorInCorner8 = 8;
                    break;
            }
            // number of neighbours
            if(_lattice_range==0) { // 4 neighbours
                doDispersal = 
                    &Disperser::disperse_2D_stepping_stone_4Neighbour;
                _migr_rateIn[FEM] = _migr_rate[FEM]/factorIn4;
                _migr_rateIn[MAL] = _migr_rate[MAL]/factorIn4;
                _migr_rateOut[FEM] = 
                    factorOut4 ? _migr_rate[FEM]/factorOut4 : 0;
                _migr_rateOut[MAL] = 
                    factorOut4 ? _migr_rate[MAL]/factorOut4 : 0;
                _migr_rateCorner[FEM] = _migr_rate[FEM]/factorInCorner4;
                _migr_rateCorner[MAL] = _migr_rate[MAL]/factorInCorner4;
                _migr_rate[FEM] /= 4;
                _migr_rate[MAL] /= 4;
            } else { // 8 neighbours
                doDispersal = 
                    &Disperser::disperse_2D_stepping_stone_8Neighbour;
                _migr_rateIn[FEM] = _migr_rate[FEM]/factorIn8;
                _migr_rateIn[MAL] = _migr_rate[MAL]/factorIn8;
                _migr_rateOut[FEM] =
                    factorOut8 ? _migr_rate[FEM]/factorOut8 : 0;
                _migr_rateOut[MAL] =
                    factorOut8 ? _migr_rate[MAL]/factorOut8 : 0;
                _migr_rateCorner[FEM] = _migr_rate[FEM]/factorInCorner8;
                _migr_rateCorner[MAL] = _migr_rate[MAL]/factorInCorner8;
                _migr_rate[FEM] /= 8;
                _migr_rate[MAL] /= 8;
            }

            break;
        default:
            fatal("\nDispersal model '%i' not available!\n",_disp_model);
    }
}
