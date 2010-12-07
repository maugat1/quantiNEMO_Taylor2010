/** @file ttquanti.h
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

#ifndef ttquantiH
#define ttquantiH

#include "ttrait.h"
#include "types.h"
#include "filehandler.h"
#include "stathandler.h"
#include "metapop.h"
#include "tree.h"
#include "simulation.h"

class TTQuantiFHvalue;
class TTQuantiProto;

/******************************************************************************/
/******************************************************************************/
/**Quantitative trait coded on discrete allelic values.
 * Alleles coded on \<char\>, max number of alleles per locus is 256.*/
class TTQuanti : public TTrait{
    friend class TTQuantiProto; // we allow to access these parameters from TTQuantiProto directly
    private:
    double _genotype;             // this is the genotype (Vg)
    double _phenotype;            // this is the phenotype (Vp = Vg + Ve)
    double _fitness_factor;       // by default 1

    TTQuantiProto* pProto;

    public:

    TTQuanti (): _phenotype(my_NAN), _fitness_factor(my_NAN), pProto(0){}

    TTQuanti(const TTQuanti& T): _phenotype(my_NAN), _fitness_factor(my_NAN), pProto(T.pProto){
        _copyTTraitParameters(T);  // copy the parameters of TTrait
    }

    virtual ~TTQuanti ();

    virtual void set_from_prototype(TTraitProto* T);

    ///@}
    ///@name Implementations
    ///@{
    virtual TTQuanti& operator= (const TTrait& T);
    virtual bool    operator== (const TTrait& T);
    virtual bool    operator!= (const TTrait& T);
    virtual void    init                 ();
    virtual void    ini_sequence         (Patch* patch);
    virtual void    reset                ();
    virtual void*   set_trait            (void* value)           {return NULL;}
    virtual void**  get_sequence         ()  const              {return (void**)sequence;}
    virtual void    set_sequence         (void** seq){
        reset();
        sequence = (unsigned char**)seq;
    }
virtual void    set_value            (); // set the genotype
virtual void    set_value            (double val)          {_phenotype = val+_genotype;}   // set the phenotype
virtual double  get_value            ()					           {return _phenotype;}
virtual double  get_genotype         ()					           {return _genotype;}
virtual double  get_phenotype        ()                    {return _phenotype;}
virtual double  get_fitnessFactor   ()                    {return _fitness_factor;}
virtual void    set_fitness_factor   () ;


virtual void    show_up              ();
virtual TTQuanti*  clone     ()                      {return new TTQuanti(*this);}

void mutate();


///@}
};



class TTQuantiSH;
class TTQuantiFH;
class TTQuantiFHvalue;
/******************************************************************************/
/******************************************************************************/
/**Prototype class for the TTQuanti trait class.
 **/
class TTQuantiProto : public TTraitProto {
    friend class TTQuanti; // we allow to access these parameters from TTQuanti directly
    private:
    int     _genetic_effect;
    int     _output;              //0: nothing, 1: allelic values
    int     _Ve_model;            // environmental model: 0: Ve directly, 1: h2; 2: h2 at every gen; 3: H2; 4: H" at every gen
    double  _Ve_prop;             // proportion of the environment from the current versus the natal patch
    int     _selection_model;     // 0: stabilizing, 1: directional 2: neutral,
    int     _Va_model;            // 0: always valable, but slow, 1: limited to random mating, but fast

    public:

    static string _ini_fstat_file; // if the inital genotypes are given   // same for all quanti traits
    inline double getAllelicValue(const int& l, const unsigned char& i){
        assert(_allelicValues[l][i] != my_NAN);
        return _allelicValues[l][i];
    }

    inline double getDominanceValue(const int& l, const unsigned char& a1, const unsigned char& a2){
        double *value;
        assert(a1<a2);
        value = &_dominanceValues[l][a1][a2];
        if(*value == my_NAN) *value = SimRunner::r.Normal(_dominance_mean, _dominance_sd);  //get the value if not yet set
        return *value;
    }

    TTQuantiSH*             _stats;
    static TTQuantiFH*      _writer;
    static TTQuantiFHvalue* _phenotyper;
    static TTQuantiFHvalue* _genotyper;

    double**             _allelicValues;      // _allelicValues[locus][allele]
    bool                 _allelic_file;       // is a alleic file passed?

    double***            _dominanceValues;    // only for dominance effect            _dominanceValues[locus][allele1][allele2]
    double               _dominance_mean;
    double               _dominance_sd;
    bool                 _dominance_file;     // is a dominance file passed?

    double*              _fitnessFactor_heterozygote;   // _fitnessFactor_heterozygote[locus]
    double*              _fitnessFactor_homozygote;     // _fitnessFactor_homozygote[locus]
    double***            _fitnessFactor;                // _fitnessFactor[locus][allele1][allele2]
    Tree<unsigned char>* _fitnessFactorTree;            //  fitness factor defined for the entire genome

    Tree<unsigned char>* _phenoTree;          // only for epistatic effect
    double               _epistatic_sd;

    // determination of the genotype
    double  (TTQuantiProto::* get_locus_genotype_func_ptr)(const int& l, const unsigned char& a1, const unsigned char& a2);
    double  get_locus_genotype_additive (const int& l, const unsigned char& a1, const unsigned char& a2);
    double  get_locus_genotype_dominance_single (const int& l, const unsigned char& a1, const unsigned char& a2);
    double  get_locus_genotype_dominance_array (const int& l, const unsigned char& a1, const unsigned char& a2);

    double  (TTQuantiProto::* get_genotype_func_ptr)(unsigned char** seq);
    double  get_genotype_none (unsigned char** seq);
    double  get_genotype_epistatic (unsigned char** seq);

    double* get_fitnessFactor_heterozygote(){return _fitnessFactor_heterozygote;}
    double* get_fitnessFactor_homozygote()  {return _fitnessFactor_homozygote;}
    double*** get_fitnessFactor_array()  {return _fitnessFactor;}
    bool    fitnessFactor_used();
    double  (TTQuantiProto::* get_fitnessFactor_func_ptr)(unsigned char** seq);
    double  get_fitnessFactor_genome(unsigned char** seq);
    double  get_fitnessFactor_explicit(unsigned char** seq);
    double  get_fitnessFactor_global(unsigned char** seq);


    void    check_allelicValues_IMM();
    void    check_mutationValues_IMM();

    void    ini_paramset             ( );   // before the normal ini()

    public:

    TTQuantiProto ( );
    TTQuantiProto (int i);
    TTQuantiProto(const TTQuantiProto& T);

    ~TTQuantiProto ( );

    //implementation of TTraitProto:
    virtual void                     init (Metapop* pMetapop);

    virtual void                      reset ();
    virtual void                      resetTotal ();

    virtual TTQuanti*          hatch ();

    virtual TTQuantiProto*      clone () {return new TTQuantiProto(*this);}

    //implementation of SimComponent:
    virtual void loadFileServices ( FileServices* loader );

    virtual void loadStatServices ( StatServices* loader );

    void read_allele_file(string name);
    void set_fitnessFactor(const string& trait, string name, double* &array);
    void set_allelicValues(const string& trait);
    void set_allelicValues(TMatrix* m, vector<map<int, string> >* pStrMatrix, double* array,
            const unsigned int& i, const unsigned int& l, const unsigned int& a, int* cols);
    void set_mutationFreq(const string& trait);
    void set_initAlleleFreq(const string& trait);
    void read_locus_file(string name);
    void read_genome_file(string name);

    void print_allelic_values (string name);
    void print_dominance_values (string name);
    void print_epistatic_values (string name);
    void print_gentoype(ostream& FILE, unsigned char** seq, const int& digit);
    bool get_next_gentoype(unsigned char** seq);

    virtual void executeBeforeEachGeneration(const int& gen);
    virtual void executeAfterEachReplicate(const int& rep);
    virtual void executeBeforeEachReplicate(const int& rep);

    int     get_selection_model()     const     {return _selection_model;}
    int     get_Va_model()            const     {return _Va_model;}

    string  get_info             ();

    void    create_regular_spaced_array(double* array, const int& size, double half_range);
    void    compute_frequencies(double* array, double* effect_array, int& size, double& sd);

};




/******************************************************************************/
/******************************************************************************/
/** A file handler to save the discrete quantitative trait genotypes in a FSTAT-like text file.
 *  The file extension is ".dat".
 */
class TTQuantiFH: public FileHandler {

    private:


    public:

        TTQuantiFH (bool rpl_per = 0, bool gen_per = 0, int rpl_occ = 0, int gen_occ = 0,TTQuantiProto* TP = NULL)
            : FileHandler (".dat"){
                FileHandler::set(rpl_per,gen_per,rpl_occ,gen_occ,0,"","",0,0,TP,1);
            }

        virtual ~TTQuantiFH ( ) { }

};


/******************************************************************************/
/******************************************************************************/
/** File handler used to save the phenotypes 0por genotypes generated from the
 * discrete quantitative trait. The file extension is
 * ".phe". !!!
 */
class TTQuantiFHvalue: public FileHandler {

    private:

        string _name;
        virtual void FHwrite (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
                Patch* current_patch, const int& patch_id);

        double (TTQuantiFHvalue::*get_value_func_ptr)(Individual* ind, const int& t);
        double get_genotype(Individual* ind, const int& t){
            return ind->getTrait(_TTidx[t])->get_genotype();
        }
        double get_phenotype(Individual* ind, const int& t){
            return ind->getTrait(_TTidx[t])->get_phenotype();
        }


    public:

        TTQuantiFHvalue (bool rpl_per = 0, bool gen_per = 0, int rpl_occ = 0, int gen_occ = 0, TTQuantiProto* TT = NULL)
            : FileHandler (".phe")
        {
            FileHandler::set(rpl_per,gen_per,rpl_occ,gen_occ,0,"","",0,0,TT,1);
        }

        virtual  void set_getter(int i){
            switch(i){
                case 0: get_value_func_ptr= &TTQuantiFHvalue::get_genotype;
                        _name = "genotypic_value";
                        set_extension(".gen");
                        break;
                case 1: get_value_func_ptr= &TTQuantiFHvalue::get_phenotype;
                        _name = "phenotypic_value";
                        set_extension(".phe");
                        break;
                default: fatal("TTQuantiFHvalue::set_getter: Not valid option for getter!\n");
            }
        }

        virtual ~TTQuantiFHvalue ( ) { }

        virtual void FHwrite  ();
};

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                               ****** TTQuantiSH ******

/*_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\*/
/**The stat handler for neutral markers. */
class TTQuantiSH: public StatHandler<TTQuantiSH> {
    private:


        double *_varA, *_varA2,
               *_meanG, *_varG,
               *_meanP, *_varP,
               _h2;

        double** _qst_matrix;

        double Theta_FF, Theta_MM, Theta_FM;
        double _mean_theta, _mean_alpha;

    public:

        TTQuantiSH (TTQuantiProto* TT) :  _varA(0), _varA2(0), _meanG(0), _varG(0),
        _meanP(0), _varP(0), _h2(0), _qst_matrix(0)
    {
        set(TT);

        // set the function to compute the additive genetic variance
        if(     TT->get_locus_genotype_func_ptr == &TTQuantiProto::get_locus_genotype_additive
                && TT->get_genotype_func_ptr       == &TTQuantiProto::get_genotype_none)
        {
            // purely additive: Va = Vg
            get_Va_ofPatch_func_ptr = &TTQuantiSH::setMeanAndVar_Vg_ofPatch;
        }
        else{
            // non additive genetic effects present
            if(TT->get_Va_model() == 0) get_Va_ofPatch_func_ptr = &TTQuantiSH::get_Va_ofPatch_regression;
            else                        get_Va_ofPatch_func_ptr = &TTQuantiSH::get_Va_ofPatch_random_mating;
        }
    }

        virtual ~TTQuantiSH ( )
        {
            if(_varA)     delete[] _varA;
            if(_varA2)    delete[] _varA2;
            if(_meanG)    delete[] _meanG;
            if(_varG)     delete[] _varG;
            if(_meanP)    delete[] _meanP;
            if(_varP)     delete[] _varP;

            ARRAY::delete_2D(_qst_matrix, _sample_pops_size);

        }

        virtual bool init ( ) ;

        virtual bool setStatRecorders (const string& token);

        // variance components
        bool   compute_alpha(double* y, int** x, int nb_ind, double* alpha,   // singular value decompostion method
                const vector<int>& availableAllele);
        void remove_private_alleles_compute_alpha(Patch* crnt_patch, const unsigned int& sizeF,
                const unsigned int& sizeM, double* alpha,
                vector<int>& availableAllele, double* arrayG,
                const age_idx& age_pos, unsigned int** allele_counts,
                const unsigned int& l);
        void (TTQuantiSH::*get_Va_ofPatch_func_ptr)(Patch*, const age_t&, double&, double&, double** freqs, unsigned int** counts);
        void get_Va_ofPatch(Patch* p, const age_t& a, double& m, double& v, double** f, unsigned int** c){(this->*get_Va_ofPatch_func_ptr)(p,a,m,v,f,c);}
        void get_Va_ofPatch_regression(Patch*, const age_t&, double&, double&, double** freqs, unsigned int** counts);    // any case, but slower
        void get_Va_ofPatch_random_mating(Patch*, const age_t&, double&, double&, double** freqs, unsigned int** counts); // only random mating

        double getVarA      (unsigned int i, const age_t& AGE)  {
            setVar_Va(AGE);
            return _varA[i];
        }
        // genotypic value
        double getMeanG     (unsigned int i, const age_t& AGE)  {
            setMeanAndVar_Vg(AGE);
            return _meanG[i];
        }
        double getVarG      (unsigned int i, const age_t& AGE)  {
            setMeanAndVar_Vg(AGE);
            return _varG[i];
        }

        // phenotypic value
        double getMeanP     (unsigned int i)  {
            setMeanAndVar_Vp();
            return _meanP[i];
        }
        double getVarP      (unsigned int i)  {
            setMeanAndVar_Vp();
            return _varP[i];
        }

        double getQst       ();

        void setQst_perPatchPair();
        double getQst_ij(unsigned int i){
            setQst_perPatchPair();
            return _qst_matrix[i%_sample_pops_size][(int)(i/_sample_pops_size)];
        }


        // between variance
        double getVgB       (const age_t& AGE);
        double getVpB       ();

        // within variance
        double getVaW       (const age_t& AGE);
        double getVgW       (const age_t& AGE);
        double getVpW       ();

        void   setVar_Va(const age_t& AGE);
        void   setMeanAndVar_Vg(const age_t& AGE);
        void   setMeanAndVar_Vg_ofPatch(Patch*, const age_t&, double&, double&, double** freqs=NULL, unsigned int** counts=NULL);
        void   setMeanAndVar_Vp();
        void   setMeanAndVar_Vp_ofPatch(Patch* crnt_patch, double& meanP, double& varP);
};

#endif //TTDISCRETEQUANTI_H

