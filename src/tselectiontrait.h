/** @file tselectiontrait.h
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

#ifndef tselectiontraitH
#define tselectiontraitH
//---------------------------------------------------------------------------
#include "simcomponent.h"


class TSelection;
class TTQuantiProto;
class Patch;
class Individual;


class TSelectionTrait{
  protected:
		TSelection* _pSel;        // pointer to the selection object (do not delete)
    Metapop*    _popPtr;      // pointer to the metapopulation

    int _traitIndex;          // index of the quantitative trait among all traits
    int _quantiIndex;         // index of the quantitative trait among all quantitative traits

    // selection pressure
    int     _nb_selection_params;
    double* _selection_pressure[5]; // if stabilizing selection:
                                    //   1. pos: optima       (default: 0)
                                    //   2. pos: intensity    (default: 1)
                                    // if directional selection:
                                    //   1. pos: min          (minimal fitness,             default: 0)
                                    //   2. pos: max          (maximal fitness,             default: 1)
																		//   3. pos: max growth   (phenotype of maximal growth, default: 1)
																		//   4. pos: growth rate  (steepness of the slope,      default: 1)
                                    //   5. pos: symmetry     (symmetry of the slope,       default: 1)
    double* _selection_sd;

		int     _selection_model;       // 0: stabilizing selection (default)
																		// 1: directional selection
																		// 2: neutral

    bool   _independent;



    // environment
    int     _environmental_model;
    double  _Ve_prop;         // Proportion of the used Ve (0: only natal patch; 1: only mating patch(default))
    double* _Ve_h2[2];        // heritability, repsectivly Ve directly  _Ve_hs[sex][patch]
    int     _Ve_model;        // 0: set VarE(0), VarE(t) constant;
                              // 1: set h2(0),  VarE(t) constant;
                              // 2: set h2(t), h2(t) constant, VarE(t) variable
                              // 3: set H2(0),  VarE(t) constant;
                              // 4: set H2(t), H2(t) constant, VarE(t) variable

    double _pheno;            // phenotype

    // current settings
    Patch*      _curPatch;
    sex_t       _curSex;
    Individual* _curInd;





	public:
		TSelectionTrait();
		virtual ~TSelectionTrait(){
			if(_selection_sd) delete[] _selection_sd;
			if(_selection_pressure[FEM]) delete[] _selection_pressure[FEM];
			if(_selection_pressure[MAL]) delete[] _selection_pressure[MAL];
		}

		TTQuantiProto* _pQuantiProto;  // pointer to the corresponding quantitative trait (do not delete)

		void init(TSelection* s, const int& trait);
    virtual void init() = 0;
		virtual void set_ve_mean_func_ptr();
		double get_fitnessFactor_individual();
		virtual double get_fitness()=0;
		virtual int get_nb_selection_params() {return _nb_selection_params;}

    // functions to set the selection pressure
		typedef double (TSelectionTrait::*_func_ptr)(double value, const int&);
		void get_selection_pressure(Patch* patch, sex_t SEX);
    _func_ptr _get_selection_pressure_func_ptr[5];
    double  _get_selection_pressure_var(double value, const int& i);
    double  _get_selection_pressure_const(double value, const int& i);

    // function pointer to get the Ve (from natal, current patch, or a mix of them)
		double (TSelectionTrait::*func_ptr_get_meanVe)();
    double (TSelectionTrait::*func_ptr_get_sdVe)();
    double  _get_null() {return 0;}
    double  _get_meanVe_current();
    double  _get_meanVe_natal();
    double  _get_meanVe_mix();
    double  _get_sdVe_current();
		double  _get_sdVe_natal();
    double  _get_sdVe_mix();
    double  get_Ve_prop()             const     {return _Ve_prop;}


    /** phenotype specific functions ********************************************/
		double   get_phenotype();
		void   (TSelectionTrait::* func_ptr_set_phenotype)(Individual* ind);
		void    set_phenotype(Individual* ind);
		void    set_phenotype_and_fitness_factor(Individual* ind);

		void    set_Ve();
		void    set_quantiHeritability();


		// getter
    int     get_Ve_model    ( )             const  {return _Ve_model;}
    bool    get_independent ( )             const  {return _independent;}
    TTQuantiProto* get_pQuantiProto()       const  {return _pQuantiProto;}

    void setPhenotype_funcPointer(bool corr);
};
#endif
