/** @file lce_disperse.h
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

#ifndef lce_disperseH
#define lce_disperseH

#include <vector>
#include "lifecycleevent.h"
#include "param.h"

// enumerated type for the dispersal
typedef enum {
    DISP_NONE = 0,
    DISP_MIGR = 1,
    DISP_COLN = 2,
    DISP_BOTH = 3,
} DISP_TYPE;


/* Stores the dispersal matrices and dispersal model parameters and interface.
 * Takes care of the process of dispersal. */
class Disperser {
    protected:
        ParamSet *_paramSet;
        Metapop *_popPtr;
        // 0: None, 1: migration, 2: colonization, 3: both
        DISP_TYPE _disp_type;
        // 0: island; 1: propagule island, 2: 1D SS; 3: 2D SS
        int    _disp_model;
        // 0: torus; 1: reflecting; 2: absorbing
        int    _border_model;
        // 0: 4 neighbours; 1: 8 neighbours
        int    _lattice_range;
        double _disp_propagule_prob;
        unsigned int _nb_patch;
        // dimensions of the 2D lattice 0: nbRows; 1: nbColumns
        unsigned int _lattice_dims[2];
        // used when the migration rate differs depending on the popualtion
        // density
        //    0: min           (default: 0)
        //    1: max           (default: 1)
        //    2: max growth    (default: 0)
        //    3: growth rate   (default: 1e9)
        //    4: symmetry      (default: 1)
        double* _disp_factor;
        /**The sex-specific dispersal matrices, [0] for males, [1] for females,
         * might be used as connectivity matrix as well*/
        TMatrix* _dispMatrix[2];
        // migration rate of males [0], and females [1]
        double _migr_rate[2];
        // 2D-SS: migration from the edge to the middle (negative numbers mean
        // that the ind get lost (absorbing boundaries))
        double _migr_rateIn[2];
        // 2D-SS: migration from the edge to the outter
        double _migr_rateOut[2];
        // 2D-SS: migration from the corner to the middle
        double _migr_rateCorner[2];
        double _migr_rate_propagule[2];

        void _disperse_2D_SS_corner_4Neighbour();
        void _disperse_2D_SS_corner_8Neighbour();
        void _disperse_2D_SS_middle_4Neighbour();
        void _disperse_2D_SS_middle_8Neighbour();
        void _disperse_2D_SS_edge_4Neighbour();
        void _disperse_2D_SS_edge_8Neighbour();

        /** factor to the migration rate */
        void _setDispersalFactor();
        double (Disperser::* get_disp_factor_funcPtr)(Patch* p);
        double (Disperser::* get_coln_factor_funcPtr)(Patch* p);
        double get_disp_factor_one(Patch* p) {return 1;}
        double get_disp_factor_min(Patch* p) { return _disp_factor[0]; }
        double get_disp_factor_max(Patch* p) { return _disp_factor[1]; }
        double get_disp_factor_k_threshold(Patch* p) {
            return p->get_density(OFFSPRG)<_disp_factor[2] ?
                _disp_factor[0] : _disp_factor[1];
        }
        double get_disp_factor_k_logistic(Patch* p) {
            return generalLogisticCurve(p->get_density(OFFSPRG),
                    _disp_factor[0], _disp_factor[1], _disp_factor[2],
                    _disp_factor[3], _disp_factor[4]);
        }

        /***** General *****/
        void _sendEmigrants(Patch* curPatch,const int& home,const int& target,
                int* nbInd, double* dispRate, const double& factor);
        void _sendEmigrants(Patch* curPatch,const int& home,const int& target,
                const int& size,const double& dispRate, const double& factor,
                const sex_t& SEX);
        void _get_lattice_dims();


    public:
        // dispersal function pointer
        void (Disperser::* doDispersal) ();

        ///@name Dispersal Matrix
        ///@{
        void setDispersalMatrix();
        void setDispersalRate();
        void checkDispMatrix(TMatrix* mat);
        void setDispMatrix(TMatrix* mat);
        void setDispMatrix(sex_t sex, TMatrix* mat);
        ///@}

        int  get_disp_type() { return _disp_type; }
        int  get_disp_model() { return _disp_model; }
        string get_disp_model_str()
        {
            if(doDispersal == &Disperser::disperse_zeroDispersal)
                return "no";
            if(doDispersal == &Disperser::dispersal_matrix)
                return "matrix";
            switch(_disp_model){
                case 0: return "island";
                case 1: return "propagule island";
                case 2: return "1D stepping stone";
                case 3: return "2D stepping stone";
            }
            return"";
        }

        virtual void disperse_zeroDispersal();
        virtual void disperse_island();
        virtual void disperse_island_propagule();
        virtual void disperse_1D_stepping_stone();
        virtual void dispersal_matrix();
        virtual void disperse_2D_stepping_stone_4Neighbour();
        virtual void disperse_2D_stepping_stone_8Neighbour();

        /***** General *****/
        Disperser(ParamSet *paramSet, DISP_TYPE type = DISP_NONE);
        virtual ~Disperser();
        virtual void execute ();
        virtual bool init(Metapop* popPtr);
};


/* The base class of the dispersal LCEs, all events move offspring to the
 * post-dispersal patch containers. */
class LCE_Disperse: public LCE
{
    protected:
        Disperser *_disperser;
        Disperser *_colonizer;
        unsigned int _nb_patch;
    public:
        int isSet_disperser() { return _disperser != NULL; }
        int isSet_colonizer() { return _colonizer != NULL; }
        int  get_disp_model()
        {
            return _disperser ? _disperser->get_disp_model() : 0;
        }
        int  get_coln_model()
        {
            return _colonizer ? _colonizer->get_disp_model() : get_disp_model();
        }
        string get_disp_model_str()
        {
            if(_disperser) {
                if(_disperser->doDispersal ==
                        &Disperser::disperse_zeroDispersal)
                    return "no migration";
                if(_disperser->doDispersal == &Disperser::dispersal_matrix)
                    return "matrix";
                switch(get_disp_model()){
                    case 0: return "island";
                    case 1: return "propagule island";
                    case 2: return "1D stepping stone";
                    case 3: return "2D stepping stone";
                }
            }
            return "";
        }
        string get_coln_model_str()
        {
            if(_colonizer) {
                if(_colonizer->doDispersal ==
                        &Disperser::disperse_zeroDispersal)
                    return "no migration";
                if(_colonizer->doDispersal == &Disperser::dispersal_matrix)
                    return "matrix";
                switch(get_coln_model()){
                    case 0: return "island";
                    case 1: return "propagule island";
                    case 2: return "1D stepping stone";
                    case 3: return "2D stepping stone";
                }
            }
            return get_disp_model_str();
        }

        /***** General *****/
        LCE_Disperse(int rank = my_NAN);
        virtual ~LCE_Disperse();
        virtual void execute ();

        ///@name Implementations
        ///@{
        virtual bool init(Metapop* popPtr);
        virtual LCE_Disperse* clone () {return new LCE_Disperse();}
        virtual void loadFileServices ( FileServices* loader ) {}
        virtual void loadStatServices ( StatServices* loader ) {}
        virtual age_t removeAgeClass ( ) {return OFFSPRG;}
        virtual age_t addAgeClass ( ) {return ADULTS;}
        virtual age_t requiredAgeClass () {return OFFSPRG;}
        ///@}

        void preDispersal();
        void postDispersal();

        virtual void executeBeforeEachGeneration(const int& gen);
};


#endif //LCEDISPERSE_H

