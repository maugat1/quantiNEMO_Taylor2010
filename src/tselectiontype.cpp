/** @file tselectionneutral.cpp
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


#include "tselectiontype.h"
#include "tselection.h"
#include "ttquanti.h"

// ----------------------------------------------------------------------------------------
// ini_neutral_selection
// ----------------------------------------------------------------------------------------
/** set the selection optima and intensity of the patches */
    void
TSelectionNeutral::init()
{
    _nb_selection_params = 0;

    if(_pQuantiProto->fitnessFactor_used()) _func_ptr_get_fitness = &TSelectionNeutral::get_fitnessFactor;
    else                                    _func_ptr_get_fitness = &TSelectionNeutral::get_fitness_none;
}

// ----------------------------------------------------------------------------------------
// ini_stabilizing_selection
// ----------------------------------------------------------------------------------------
/** set the selection optima and intensity of the patches */
    void
TSelectionStabilizing::init()
{
    _nb_selection_params = 2;
    _selection_pressure[FEM] = new double[_nb_selection_params];
    _selection_pressure[MAL] = new double[_nb_selection_params];
    _selection_sd            = new double[_nb_selection_params];

    // does the optima varies from generation to generation?
    if(_pSel->get_optima_sd())
        _get_selection_pressure_func_ptr[0] = &TSelectionTrait::_get_selection_pressure_var;
    else  _get_selection_pressure_func_ptr[0] = &TSelectionTrait::_get_selection_pressure_const;

    // does the intensity varies from generation to generation?
    if(_pSel->get_intensity_sd())
        _get_selection_pressure_func_ptr[1] = &TSelectionTrait::_get_selection_pressure_var;
    else  _get_selection_pressure_func_ptr[1] = &TSelectionTrait::_get_selection_pressure_const;

    // is there any fitness factor set?
    if(_pQuantiProto->fitnessFactor_used()) _func_ptr_get_fitness = &TSelectionStabilizing::get_fitnessFactor;
    else                                    _func_ptr_get_fitness = &TSelectionStabilizing::get_fitness_none;
}

// ----------------------------------------------------------------------------------------
// ini_directional_selection
// ----------------------------------------------------------------------------------------
    void
TSelectionDirectional::init()
{
    _nb_selection_params = 5;
    _selection_pressure[FEM] = new double[_nb_selection_params];
    _selection_pressure[MAL] = new double[_nb_selection_params];
    _selection_sd            = new double[_nb_selection_params];

    // does the lower asymptote varies from generation to generation?
    if(_pSel->get_min_sd())
        _get_selection_pressure_func_ptr[0] = &TSelectionTrait::_get_selection_pressure_var;
    else  _get_selection_pressure_func_ptr[0] = &TSelectionTrait::_get_selection_pressure_const;

    // does the upper asymptote varies from generation to generation?
    if(_pSel->get_max_sd())
        _get_selection_pressure_func_ptr[1] = &TSelectionTrait::_get_selection_pressure_var;
    else  _get_selection_pressure_func_ptr[1] = &TSelectionTrait::_get_selection_pressure_const;

    // does the maximal growth varies from generation to generation?
    if(_pSel->get_max_growth_sd())
        _get_selection_pressure_func_ptr[2] = &TSelectionTrait::_get_selection_pressure_var;
    else  _get_selection_pressure_func_ptr[2] = &TSelectionTrait::_get_selection_pressure_const;

    // does the growth rate varies from generation to generation?
    if(_pSel->get_growth_rate_sd())
        _get_selection_pressure_func_ptr[3] = &TSelectionTrait::_get_selection_pressure_var;
    else  _get_selection_pressure_func_ptr[3] = &TSelectionTrait::_get_selection_pressure_const;

    // does the symmetry varies from generation to generation?
    if(_pSel->get_symmetry_sd())
        _get_selection_pressure_func_ptr[4] = &TSelectionTrait::_get_selection_pressure_var;
    else  _get_selection_pressure_func_ptr[4] = &TSelectionTrait::_get_selection_pressure_const;

    // is there any fitness factor set?
    if(_pQuantiProto->fitnessFactor_used()) _func_ptr_get_fitness = &TSelectionDirectional::get_fitnessFactor;
    else                                    _func_ptr_get_fitness = &TSelectionDirectional::get_fitness_none;
}

//-----------------------------------------------------------------------------
// _get_fitness_neutral
//-----------------------------------------------------------------------------
/** if no selection acts we return a fitness of 1 (function should never be called) */
    double
TSelectionNeutral::get_fitness_none()
{
    return 1;
}

//-----------------------------------------------------------------------------
// get_fitnessFactor:
//-----------------------------------------------------------------------------
    double
TSelectionNeutral::get_fitnessFactor()
{
    return get_fitness_none() * get_fitnessFactor_individual();
}

//-----------------------------------------------------------------------------
// get_fitness_stabilizing:
//-----------------------------------------------------------------------------
/** standard stabilizing function
 * w = exp(-[(vp-optima)/omega]^2 / 2))
 * intensity is the sd
 * special is the weight when multiple tratis are used
 */
    double
TSelectionStabilizing::get_fitness_none()
{
    double value = (get_phenotype()-_selection_pressure[_curSex][0])/_selection_pressure[_curSex][1];
    value *= value;
    return exp(value/(-2));
}

//-----------------------------------------------------------------------------
// get_fitness:
//-----------------------------------------------------------------------------
    double
TSelectionStabilizing::get_fitnessFactor()
{
    return get_fitness_none() * get_fitnessFactor_individual();
}

//-----------------------------------------------------------------------------
// get_fitness_directionl
//-----------------------------------------------------------------------------
/** directional selection with a general logistic curve.
 *   A: the lower asymptote;
 *   C: the upper asymptote;
 *   M: the time of maximum growth;
 *   B: the growth rate;
 *   T: affects near which asymptote maximum growth occurs.
 */
    double
TSelectionDirectional::get_fitness_none()
{
    return generalLogisticCurve(get_phenotype(),
            _selection_pressure[_curSex][0],  // min
            _selection_pressure[_curSex][1],  // max
            _selection_pressure[_curSex][2],  // max growth
            _selection_pressure[_curSex][3],  // growth rate
            _selection_pressure[_curSex][4]); // symmetry
}

//-----------------------------------------------------------------------------
// get_fitness:
//-----------------------------------------------------------------------------
    double
TSelectionDirectional::get_fitnessFactor()
{
    return get_fitness_none() * get_fitnessFactor_individual();
}
//-----------------------------------------------------------------------------




