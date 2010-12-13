/** @file stathandler.h
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

#ifndef stathandlerH
#define stathandlerH

#include "stathandlerbase.h"
#include <list>


class TTraitProto;
class Individual;
class Patch;
/**A class to compute and store the summary statistics associated with a SimComponent.
 * The template type must be the type of the class that declares the methods linked into the
 * StatRecorder elements. */
template <class SH>
class StatHandler: public StatHandlerBase {

    private:
        double get_coancestry(Patch*, const sex_t&, const unsigned int&,
                Patch*, const sex_t&, const unsigned int&,
                const age_idx&);

    protected:
        /**The list of stat recorders.*/
        list<StatRecorder<SH>*> _recorders;

        TTraitProto* _SHLinkedTrait;
        int _SHLinkedTraitIndex;
        int _trait_index; // for multiple instansiations of the corresponding trait

        typedef typename list< StatRecorder<SH>* >::iterator REC_IT;

        // allele frequencies / counts
        unsigned int _nb_locus, _nb_allele, _ploidy;

        unsigned int*** _alleleCountTable;   // [patch][locus][allele]
        double***       _alleleFreqTable;    // [patch][locus][allele]
        double**        _globalAlleleFreq;   // [locus][allele]

        bool already_computed(unsigned int* array, age_t AGE=NONE);

        bool set_stat_coancestry(string t, string i, string trait);
        bool set_stat_fstat(string t, string i, string trait);
        bool set_stat_all_freq_local(string t, string i, string trait);
        bool set_stat_all_freq_global(string t, string i, string trait);
        bool set_stat_inbreeding(string t, string i, string trait);


        ///@}
        ///@name Coancestries
        ///@{
        double   Theta_FF, Theta_MM, Theta_FM;
        double   _mean_theta, _mean_alpha;
        double** _coa_matrix;

        /**Kinship classes proportions*/
        double _sib_prop[5];    // 0: non sib, 1: maternal half sibs; 2: paternal half sibs; 3: full sib; 4: selfed
        double _sib_coa[5];     // 0: non sib, 1: maternal half sibs; 2: paternal half sibs; 3: full sib; 4: selfed


        /**Gives the coancestry (probability of identity by state) of two gene sequences.
          The probability returned is the averageproportion of identical alleles per locus between the two sequences.
          @param ind1 first sequence, treated as of type (unsigned char**)
          @param ind2 second sequence, treated as of type (unsigned char**)
          */
        double Coancestry               (void** ind1, void** ind2);

        /**Computes the within and between patches coancestry coefficients.
          @param age_pos the age class index
          @param dim the dimension of the matrix to fill:
          - 1 = the diagonal (i.e. the wihtin patch coancestries or theta's)
          - 2 = the upper half (i.e. the between patch coancestries or alpha's)
          - 3 = both
          */
        void   setCoaMatrixTheta             (const age_t& AGE);
        void   setCoaMatrixAlpha             (const age_t& AGE);
        void   setSexspecific_Theta          (const age_t& AGE);

        /**Gets the given coancestry coefficient from the coancestry matrix.
          @param i combination of the row and column indexes (see setCoaMatrixRecorders()).
          \note the upper half and the diagonal of the matrix are filled, other positions are set to 0.
          */
        double getCoaTheta (unsigned int i, const age_t& AGE){setCoaMatrixTheta(AGE); return _coa_matrix[i][i];}                     // diagonal
        double getCoaAlpha (unsigned int i, const age_t& AGE){setCoaMatrixAlpha(AGE); return _coa_matrix[i%_sample_pops_size][(int)(i/_sample_pops_size)];} // pairwise
        double getMeanTheta(const age_t& AGE)                {setCoaMatrixTheta(AGE); return _mean_theta;}
        double getMeanAlpha(const age_t& AGE)                {setCoaMatrixAlpha(AGE); return _mean_alpha;}

        /**Gives the mean within females coancestry coefficient.*/
        double getTheta_FF              (const age_t& AGE)		       {setSexspecific_Theta(AGE); return Theta_FF;}

        /**Gives the mean within males coancestry coefficient.*/
        double getTheta_MM              (const age_t& AGE)		       {setSexspecific_Theta(AGE); return Theta_MM;}

        /**Gives the mean between males and females coancestry coefficient.*/
        double getTheta_FM              (const age_t& AGE)		       {setSexspecific_Theta(AGE); return Theta_FM;}

        void   setSibStats              (const age_t& AGE);
        void   setSibCoa                (Individual *I1, Individual *I2);
        double getSibProportions        (unsigned int i, const age_t& AGE) {setSibStats(AGE); return _sib_prop[i];}
        double getSibCoaMeans           (unsigned int i, const age_t& AGE) {setSibStats(AGE); return _sib_coa[i];}
        ///@}


        //fstat:
        double _ho, _hs, _ht, _hsnei, _htnei, _nb_all_global, _fst, _fis, _fit;
        double _fst_wc, _fis_wc, _fit_wc;
        double** _fst_matrix_wc, **_fst_matrix;

        double *_hsnei_locus, *_htnei_locus, *_fst_locus, *_fis_locus, *_fit_locus,  // for Nei and Chesser
        *_ho_locus, *_hs_locus, *_ht_locus;

double *_fst_WC_locus, *_fis_WC_locus, *_fit_WC_locus;                       // Weir and Cockerham

double *_nb_allele_of_pop;            // mean number of alleles (across loci) of patch i



unsigned int** _computed;    // computed[stat][option]: 0: gen; 1: rep; 2: age
unsigned int _computed_size;

public:

StatHandler( ):_SHLinkedTrait(0),  _trait_index(0),
    _alleleCountTable(0), _alleleFreqTable(0), _globalAlleleFreq(0),
    _coa_matrix(0), _fst_matrix_wc(0),_fst_matrix(0),
    _hsnei_locus(0), _htnei_locus(0), _fst_locus(0), _fis_locus(0), _fit_locus(0),
    _ho_locus(0), _hs_locus(0), _ht_locus(0), _fst_WC_locus(0), _fis_WC_locus(0), _fit_WC_locus(0),
    _nb_allele_of_pop(0){

        _computed_size = 25;
        _computed = new unsigned int*[_computed_size];
        for(unsigned int i=0; i<_computed_size; ++i){
            _computed[i] = new unsigned int[3];
            _computed[i][0] = 0;
        }
    }

void set(TTraitProto* TT);

virtual ~StatHandler ( );

virtual bool init();

virtual string getName() {return "StatHandler";}

/**Empties the _recorders list, they are destroyed in StatHandlerBase::reset().*/
virtual void clear   ( )  {_recorders.clear();}

/**Computes the stats by executing the function variables stored in the StatRecorder's.*/
virtual void execute ( );
void (StatHandler::* stat_save)();    // pointer to the correct stat function
virtual void stats_to_db( );    					// stats are stored int eh db
virtual void stats_to_file( );  					// stats are directly written to the file

/** Adds a StatRecorder to the list (for each pairwise patch combination: postfix: _pi-j) */
virtual void add_pairwisePatch (string Title, string Name, st_order Order, age_t AGE, unsigned int ARG,
        double(SH::* getStat)(void),
        double(SH::* getStatBoolArg)(bool)=0,
        double(SH::* getStatUintArg)(unsigned int)=0,
        double(SH::* getStatAGE)(const age_t&)=0,
        double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
        double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);

/** Adds a StatRecorder to the list (for each patch: postfix: _pi) */
virtual void add_perPatch (string Title, string Name, st_order Order, age_t AGE, unsigned int ARG,
        double(SH::* getStat)(void),
        double(SH::* getStatBoolArg)(bool)=0,
        double(SH::* getStatUintArg)(unsigned int)=0,
        double(SH::* getStatAGE)(const age_t&)=0,
        double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
        double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);

/** Adds a StatRecorder to the list (for each locus: postfix: _li) */
virtual void add_perLocus (string Title, string Name, st_order Order, age_t AGE, unsigned int ARG,
        double(SH::* getStat)(void),
        double(SH::* getStatBoolArg)(bool)=0,
        double(SH::* getStatUintArg)(unsigned int)=0,
        double(SH::* getStatAGE)(const age_t&)=0,
        double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
        double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);

/** Adds a StatRecorder to the list (a single stat) */
virtual void add (string Title, string Name, st_order Order, age_t AGE, unsigned int ARG,
        double(SH::* getStat)(void),
        double(SH::* getStatBoolArg)(bool)=0,
        double(SH::* getStatUintArg)(unsigned int)=0,
        double(SH::* getStatAGE)(const age_t&)=0,
        double(SH::* getStatBoolArgAGE)(bool, const age_t&)=0,
        double(SH::* getStatUintArgAGE)(unsigned int, const age_t&)=0);

// allele frequencies
void    allocateTable();
double  get_allele_freq_local(unsigned int i, const age_t& AGE);
double  get_allele_freq_global(unsigned int i, const age_t& AGE);
void    set_alleleCountTable_local(age_t AGE);
void    set_alleleCountTable_local_ofPatch(Patch*, age_t, double**&, unsigned int**&);
void    set_alleleCountTable_global(age_t AGE);
void    set_alleleCountTable_global(unsigned int nbPatch, const unsigned int& nbInd,
        unsigned int*** alleleCounts, double** globalFreqs);

double  getHarmonicMean_ofPopSize(const age_t& AGE, Patch** aPatch, unsigned int nbPatch);

///@name F-stats:
///@{

/** Computes the weighted within and between patch Fst's as well as the overall Fst (Theta) */
void   setFstat_Weir_Cockerham          (age_t AGE);      // across loci
void   setFstat_Weir_Cockerham(age_t AGE, Patch** aPatch, unsigned int nbPatch,
        double*** alleleFreqs, double** alleleFreqsGlobal,
        double& fst, double& fis, double &fit);
void   setFstat_Weir_Cockerham_perLocus (age_t AGE);      // for each locus separately
void   setFstat_Weir_Cockerham_perPatchPair(age_t AGE);
void   get_sigma_of_locus_Weir_Cockerham(double& sigma_a2,
        double& sigma_b2, double& sigma_w2, const unsigned int& nbPatch, const unsigned int& nbPops,
        double*** alleleFreqs, double** alleleFreqsGlobal, const unsigned int& l,
        const unsigned int& tot_size, unsigned int* pop_sizes, double*** ho_patch_allele, const double& nc);

/** F-statistics following Weir & Cockeram (1984) */
double getFst_WC           (const age_t& AGE) {setFstat_Weir_Cockerham(AGE); return _fst_wc;}
double getFit_WC           (const age_t& AGE) {setFstat_Weir_Cockerham(AGE); return _fit_wc;}
double getFis_WC           (const age_t& AGE) {setFstat_Weir_Cockerham(AGE); return _fis_wc;}

double getFst_WC_ij                (unsigned int i, const age_t& AGE){
    setFstat_Weir_Cockerham_perPatchPair(AGE);
    return _fst_matrix_wc[i%_sample_pops_size][(int)(i/_sample_pops_size)];
}

double getFst_WC_ofLocus  (unsigned int i, const age_t& AGE) {setFstat_Weir_Cockerham_perLocus(AGE); return _fst_WC_locus[i];}
double getFis_WC_ofLocus  (unsigned int i, const age_t& AGE) {setFstat_Weir_Cockerham_perLocus(AGE); return _fis_WC_locus[i];}
double getFit_WC_ofLocus  (unsigned int i, const age_t& AGE) {setFstat_Weir_Cockerham_perLocus(AGE); return _fit_WC_locus[i];}


/** F-statistics following Nei and Chesser 1983 */
void   setFstat_Nei_Chesser              (age_t AGE);    // F-statistics globally across patch and loci
void   setFstat_Nei_Chesser_perPatchPair (age_t AGE);    // F-statistics per pair of patches and across loci
void   setFstat_Nei_Chesser_perLocus     (age_t AGE);    // F-statistics for each locus separately across patches
void   setLociDivCounter        (age_t AGE);
unsigned int get_fixed_loci_global (age_t AGE);
double get_fixed_loci_local     (age_t AGE);
double get_fixed_loci_ofPatch(age_t AGE, unsigned int p);

double getHo                    (const age_t& AGE);
double*getHo_perLocus           (const age_t& AGE);
double*getHo_perLocus           (const age_t& AGE, const vector<Patch*>& vPatch);
double getHo_ofLocus            (unsigned int i, const age_t& AGE){return getHo_perLocus(AGE)[i];}
double getHo_ofPatch            (unsigned int i, const age_t& AGE);      
double*getHo_ofPatchperLocus    (const age_t& AGE, Patch* cur_patch, double* array=NULL);
double**getHo_ofPatchperAllele  (const age_t& AGE, Patch* cur_patch, double** array=NULL);

double getHs                    (const age_t& AGE);
double getHs_ofPatch            (unsigned int p, const age_t& AGE);
double getHsUnbiased_ofPatch    (unsigned int p, const age_t& AGE);
double getHs_ofLocus            (const age_t& AGE, const unsigned int& l);
double getHs_ofPatchandLocus    (const age_t& AGE, const unsigned int& p, const unsigned int& l);

double getHt                    (const age_t& AGE);
double getHt                    (const age_t& AGE, unsigned int* patchID, unsigned int nbPatch);   // total Ht of the given patches
double getHt_ofLocus            (const age_t& AGE, const unsigned int& l);
double getHt_ofLocus            (const age_t& AGE, const unsigned int& l, unsigned int* patchID, unsigned int nbPatch);    // total Ht of the given patches

// new getters, they should work without previously setting the stats)
double getNbAllLocal            (const age_t& AGE) {setLociDivCounter(AGE); return ARRAY::mean(_nb_allele_of_pop, _sample_pops_size);}
double getNbAll_perPatch        (unsigned int i, const age_t& AGE) {setLociDivCounter(AGE); return _nb_allele_of_pop[i];}
double getNbAllGlobal           (const age_t& AGE) {setLociDivCounter(AGE); return _nb_all_global;}
double getFixLocLocal           (const age_t& AGE) {return get_fixed_loci_local(AGE);}
double getFixLoc_perPatch       (unsigned int i, const age_t& AGE) {return get_fixed_loci_ofPatch(AGE, i);}
double getFixLocGlobal          (const age_t& AGE) {return get_fixed_loci_global(AGE);}

double getHsnei                 (const age_t& AGE) {setFstat_Nei_Chesser(AGE); return _hsnei;}
double getHtnei                 (const age_t& AGE) {setFstat_Nei_Chesser(AGE); return _htnei;}
double getFst                   (const age_t& AGE) {setFstat_Nei_Chesser(AGE); return _fst;}
double getFis                   (const age_t& AGE) {setFstat_Nei_Chesser(AGE); return _fis;}
double getFit                   (const age_t& AGE) {setFstat_Nei_Chesser(AGE); return _fit;}

double getHsnei_ofLocus         (unsigned int i, const age_t& AGE) {setFstat_Nei_Chesser_perLocus(AGE); return _hsnei_locus[i];}
double getHtnei_ofLocus         (unsigned int i, const age_t& AGE) {setFstat_Nei_Chesser_perLocus(AGE); return _htnei_locus[i];}
double getFst_ofLocus           (unsigned int i, const age_t& AGE) {setFstat_Nei_Chesser_perLocus(AGE); return _fst_locus[i];}
double getFis_ofLocus           (unsigned int i, const age_t& AGE) {setFstat_Nei_Chesser_perLocus(AGE); return _fis_locus[i];}
double getFit_ofLocus           (unsigned int i, const age_t& AGE) {setFstat_Nei_Chesser_perLocus(AGE); return _fit_locus[i];}

/**Accessor to the Fst matrix as set by setFstMatrix().*/
double getFst_ij                (unsigned int i, const age_t& AGE){
    setFstat_Nei_Chesser_perPatchPair(AGE);
    return _fst_matrix[i%_sample_pops_size][(int)(i/_sample_pops_size)];

    return my_NAN;
}

};

#endif //STATHANDLER_H

