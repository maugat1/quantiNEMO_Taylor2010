/** @file metapopSH.h
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

#ifndef metapop_shH
#define metapop_shH

#include <list>
#include "stathandler.h"

class Individual;

/**A StatHandler for the Metapop SimComponent.*/
class MetapopSH : public StatHandler<MetapopSH> {

    double meanEmigrant,meanImmigrant,meanResidant,meanKolonisers,meanDeadDisp;
    unsigned int _compEmigrant[3], _compImmigrant[3], _compKolonisers[3], _compDeadDisp[3];
    double ObservedExtinctionRate;
    double _mean_reprod_success;
    double _var_reprod_success;

    // fitness
    double *_meanW, *_varW;

    public:

    MetapopSH( ) : _meanW(0), _varW(0){
    }

    virtual ~MetapopSH() {
        if(_meanW)    delete[] _meanW;
        if(_varW)     delete[] _varW;
    }

    virtual bool setStatRecorders(const string& token);
    virtual bool set_stat_demography(string token);
    virtual string getName() {return "MetapopSH";}

    ///@name Migration
    ///@{
    double getMeanEmigrantPerPatch           ();
    double getMeanImmigrantPerPatch           ();
    double getMeanMigrantRatio               ();
    double getMeanResidantPerPatch           ();
    double getMeanKolonisersProportion       ();
    double getMeanKolonisersPerPatch         ();
    ///@}

    ///@name Patch extinction
    ///@{
    void   setObsrvdExtinctionRate         ();
    double getObsrvdExtinctionRate         ()  {setObsrvdExtinctionRate();return ObservedExtinctionRate;}
    double get_isAlive                     ();
    double getColnAge( unsigned int i );
    ///@}

    ///@name Demography
    ///@{
    double getAdultSexRatio           ();
    double getOffsprgSexRatio         ();
    void setReproductiveStats         (bool sex);
    double getReproductiveMean        (bool sex) {setReproductiveStats(sex); return _mean_reprod_success;}
    double getReproductiveVar         (bool sex) {setReproductiveStats(sex); return _var_reprod_success;}

    double getTotNbInd                (const age_t& AGE);   // total number of individuals
    double getTotNbFem                (const age_t& AGE);   // total number of females
    double getTotNbMal                (const age_t& AGE);   // total number of males

    double getMeanNbInd               (const age_t& AGE);   // mean number of individuals across popualted pops
    double getMeanNbFem               (const age_t& AGE);   // mean number of females across popualted pops
    double getMeanNbMal               (const age_t& AGE);   // mean number of males across popualted pops

    double getNbInd                   (unsigned int i, const age_t& AGE);
    double getNbFem                   (unsigned int i, const age_t& AGE);
    double getNbMal                   (unsigned int i, const age_t& AGE);

    double getNbPops                  (const age_t& AGE);
    ///@}

    ///@name Kinship
    ///@{
    void   setKinship                 (const age_t& AGE);
    void   setKinClassCounter         (Individual *I1, Individual *I2);
    double getSibProportion           (unsigned int i, const age_t& AGE) {setKinship(AGE);return _sib_prop[i];}
    bool   set_stat_kinship(string t);
    ///@}

    // fitness
    // fitness
    double getMeanW     (unsigned int i)  {
        setMeanAndVar_W();
        return _meanW[i];
    }
    double getVarW      (unsigned int i)  {
        setMeanAndVar_W();
        return _varW[i];
    }
    double getVwB();
    double getVwW();
    void setMeanAndVar_W();
    void setMeanAndVar_W_ofPatch(Patch* crnt_patch, const unsigned int& i);

};

#endif

