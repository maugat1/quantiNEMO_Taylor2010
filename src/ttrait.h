/** @file ttrait.h
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

#ifndef ttraitH
#define ttraitH

#include "types.h"
#include "simcomponent.h"

#include <algorithm>
#include "random.h"

class Patch;

class TTraitProto;

/**Interface for all trait types, declares all basic trait operations.
 * This pure abstract class declares the traits interface. The precise genetic architecture
 * of the trait is not defined and is up to the designer of the trait. It uses void pointers to
 * allow users to define their own structure. Existing traits use various types for their genes
 * sequence, like \c char, \c double or \c bitset. It is up to the trait's user to know what kind
 * of structure they expect.
 * The trait objects are contained in the Individual class. Their parameters are set by their
 * TTraitProto.
 **/

class TTrait{

  /** Genetic map for the recombination:
    * Idea: if we know how long the chromosomes (in cM) are it is possible to compute
    * the mean number of recombination events. Then it is possible to compute
    * randomly the positions of the recombination events. The connection points between the
    * chromosomes are treated as recombination events with the probability of 50%.
    * Each trait knows the positions of its loci. This allows to make the inheritance within
    * each trait separately.

    * Needed elements for all the trait together (static):
    *   - vector of the cumulative lengths of the chromosomes (all chromosomes are added together)
	  *   - vector with the positions of the recombinations (f and m separately)
    *   - index for the starting chromosome (f and m separately)

	  * Needed parameters for each trait separately:
	  *   - vector of the positions ( from the start) of each loci

    * Important:
    *   - The recombination function has to be called only once for each generation of a child.
    *   - After the initial population is initialized the genetic map has to be initialized.
    *   - Before the traits are deleted the genetic map has to be deleted.
    */

public:


  unsigned char** sequence; // if the sequence is not a unsigned char it has to be redeclared

  virtual void mutate();
  virtual void inherit(TTrait* mother, TTrait* father);
  virtual void* get_allele(int loc, int all)  const;

public:
  TTraitProto* pTraitProto;

  //double get_locus_position(const int& i){return _locus_position[i];}

  // the type has to be set by the inherited class
  TTrait() : sequence(0){
  }

  TTrait(const TTrait& T) : sequence(0){
    _copyTTraitParameters(T);  // copy the parameters of TTrait
  }

  // end genetic map

protected:
  void _copyTTraitParameters(const TTrait& T); // for the copy constuctor

public:
  /**Called to allocate the trait's genotypic sequences. Called each time a new Individual is created (\c Individual::init())**/
  virtual   void            init () = 0;

  /**Called at the start of each replicate, sets the initial genotypes. Called by \c Individual::create(). **/
  virtual   void            ini_sequence (Patch* patch) = 0;

  /**Called at the end of each simulation/replicate, deallocates sequence memory. **/
  virtual   void            reset () = 0;


  /** Mutation procedure, perform mutations on the genes sequence. **/
//  virtual   void            mutate ();

  /** Called to set the phenotypic to a particular value or to give context-dependant value(s) to the trait.
    * @param value the value passed to the trait
    **/
	virtual   void*           set_trait (void* value) = 0;

	virtual   void            set_fitness_factor(){}

  virtual   void            set_from_prototype(TTraitProto* T) = 0;

  /** Called to set the sequence pointer to an existing trait
    * @param seq the existing sequence pointer
    **/
  virtual   void            set_sequence (void** seq) = 0;

  /**Tells the trait to set its phenotype from genotype, should be used instead of get_value().**/
  virtual   void            set_value () = 0;
  virtual   void            set_value (double value) = 0;

  /** Genotype to phenotype mapper.
	* @return the phenotype computed from the genotype
    **/
  virtual   double          get_value () = 0;
  virtual   double          get_phenotype (){return my_NAN;}
	virtual   double          get_fitnessFactor(){return 1;}
	virtual   double          get_genotype (){return my_NAN;}

  /** sequence accessor.
	* @return the sequence pointer
    **/
  virtual   void**          get_sequence () const = 0;

  /** Called to read one allele value at a particular locus.
   * @return the allelic value at position 'all' at locus 'loc'
   * @param loc locus position in the sequence
   * @param all which allele we want to read the value from
   **/
//  virtual   void*           get_allele   (int loc, int all) const;

  /** Writes some info to stdout. **/
  virtual   void            show_up  () = 0;

  /** Returns a copy of itself.
    \b Note: call the copy constructor of the trait which should only copy the parameters values
     not the complete state of the trait (i.e. shallow copy). The copy of the sequence data is
     made through the assignement operator!**/
  virtual   TTrait*         clone () = 0;

  ///@name Operators
  ///@{
  /**Copies the complete state of the trait from right to left side of the operator, sequence
     data included.*/
  virtual   TTrait& operator= (const TTrait&) = 0;

  /**Checks for parameters equivalence, not genetic equivalence.*/
  virtual   bool operator== (const TTrait&) = 0;
  virtual   bool operator!= (const TTrait&) = 0;

  ///@}
  virtual~TTrait ( );
};




//------------------------------------------------------------------------------
/**TTrait setter.
 * Encapsulates the methods to set the traits parameters and to generate traits. Also stores
 * the posistion of the trait in the individuals trait table.
 * This class manages the file and stat handlers through its inheritance of the SimComponent
 * interface.
*/
class TTraitProto : public SimComponent {
friend class TTrait;

   /** START of GENETIC MAP ***************************************************/
  ///@name genetic map
  ///@{
protected:
  // genetic map:
  // parameters for all the traits together (therefore static)
  // we have to know the entire length of each chromosome in centi Morgan
  // we then concatenate all the chromosomes together to a super long chromosome
  // then before a new offspring is created we have to draw new recombination positions
  // and also the starting point, i.e. selected chromosome at the start

  // no changes during simulation:
  // parameters containing the length of each chromosome (set by setGeneticMapRandom() and setGeneticMapFixed())
  static vector<double>    _chromosomeSizeVector; // vector of the length of each chromosomse in cM;
  static int               _chromosomeNb;         // number of chromosomes (_chromosomeSizeVector.size())

  // super chromosome and the mean number of recombinations for it (set by initGeneticMap())
  static double            _chromosomeLength;     // length of the super chromosome
  static double*           _chromosomeSize;       // array of the cumulative length of the chromosomses (recombination rate of 0.5);
  static double            _nbRecombinations;     // mean number of recombinations

  // the recombination positions (separate for f and m) drawn at each generation of a child (set by recombine())
  static vector<double>    _recombinationPosF;     // a list of the current recombinations for the mother
  static bool              _recombinationStartF;   // starting chromosome for the mother
  static vector<double>    _recombinationPosM;     // a list of the current recombinations for the father
  static bool              _recombinationStartM;   // starting chromosome for the father
  static void _recombine(vector<double>& vecRecombs, bool& start);

  // Vector containing the locus positions per chromosome (used to set up the super chromosome)
  vector<vector<double> > _locus_position_vector;  // set by setGeneticMapRandom() and setGeneticMapFixed()

  // super chromosome for each trait containing the locus positions
  double* _locus_position;                         // set by initTraitGeneticMap()

  // preparation to move the alleles to here
  // array containing the allele indexes
  char*     _locus_cur_allele;                     // contains the allele index, length: number of locus
  double*   _locus_mut_rate;                       // mutation rate of the locus
  double**  _locus_mut_allele_prob;                // probabilty to mutate to the allele _locu_mut_allele_prob[locus][allele]
  double*** _locus_ini_freq;                       // inital allele frequency _locus_ini_freq[locus][allele][pop]    

  TMatrix* _random_genetic_map_matrix;             // matrix of the random genetic map

  bool  set_mutation_rates(string, const string& trait);   // returns true, if the mutation rates differ among loci

public:
  /** creates the super chromosome based on the single chromosome */
  static bool initStaticGeneticMap();
  static void printGeneticMapInfo();
  static bool getNbChromosome(){return _chromosomeNb;}

  /** recombine for the female and the male separately */
  static void recombine(){
    if(_chromosomeSize){
      _recombine(_recombinationPosF, _recombinationStartF);
      _recombine(_recombinationPosM, _recombinationStartM);
    }
  }

  /** deletes the genetic map */
  static void deleteGeneticMap(){
    if(_chromosomeSize){
      delete[] _chromosomeSize;
      _chromosomeSize = NULL;
    }
  }

  /** functions to set the genetic map for each trait separately (before the call of initStaticGeneticMap()) */
  void setGeneticMapRandom(TMatrix* matrix);
  void setGeneticMapFixed(TMatrixVar<double>* matrix);
  void print_genetic_map(ofstream& FILE);

  /** function to generate the super chromosome (call after initStaticGeneticMap()) */
  void  initTraitGeneticMap();

  /** end of GENETIC MAP ******************************************************/
  ///@}
protected:

  Metapop* _popPtr;       /**The ptr to the current Metapop.*/

  double***            _initAlleleFreq;     // _initAlleleFreq[patch][locus][allele] is cumultative
  int                  _initAlleleFreqCols; // the number of columns used
  bool                 _initAlleleFreqFile;
  double**             _mutationFreq;       //  _mutationFreq[locus][allele] is cumultative
  bool                 _mutationFreqFile;
  int                  _mut_model;          // 0: RMM, 1: IMM
  int                  _ini_allele_model;   // 0: max polymorph, 1: monomorph
  vector<map<int, string> >* set_ini_allelicArray (TMatrix* mat, const int& i);
  vector<unsigned int**>* _ini_genotypes;   // vector of inital genotypes read from FSTAT file:
                                            // vector[ind], 1 value: patch, then [locus][allele]
                                            // locus and allele start already by 0
public:
  //parameters:
  int     _nb_allele;
  int     _nb_locus;
  double  _mut_rate_mean;            // 0: KAM, 1: SSM
  double* _mut_rate;
  int     _ploidy;
  string  _type;
  int     _absolute_index;      // abolute trait number across all the traits

protected:

  // parameters concerning multiple instances of the same class
  static int _nb_trait;    // total number of traits
  int     _trait_index;    // number of the instanciated class
                           // 0: if a single instanciation
                           // 1 or higher:  if several traits

  void _copyTraitPrototypeParameters(const TTraitProto& T); // for the copy constuctor

public:
  TTraitProto(): _locus_position(0), _random_genetic_map_matrix(0), _popPtr(0),  
    _initAlleleFreq(0), _initAlleleFreqFile(0), _mutationFreq(0), _mutationFreqFile(0),
    _nb_allele(0), _nb_locus(0), _mut_rate_mean(0), _mut_rate(0), _ploidy(0),
    _absolute_index(0), _trait_index(0)
    
  {
  }

  ~TTraitProto();

  virtual   TTraitProto&  operator=   (const TTraitProto& T){return setTTraitProto(T);}
  virtual   bool          operator==  (const TTraitProto& T){return isEqualTTraitProto(T);}
  virtual   bool          operator!=  (const TTraitProto& T){return !isEqualTTraitProto(T);}

  /** for the derived classes to compare directly the paramters in this base calss */
  virtual TTraitProto& setTTraitProto(const TTraitProto&);     // for operator=
  virtual bool isEqualTTraitProto(const TTraitProto&);         // for operator==
  virtual bool isUnequalTTraitProto(const TTraitProto&);       // for operator!=

  double (TTraitProto::*get_h2)(void);

  /** Inheritance procedure, creates a new trait from mother's and father's traits
   * @param mother the mother's trait
   * @param father the father's trait
   **/
  void    (TTraitProto::* _inherit_func_ptr)(TTrait* mother, TTrait* father, TTrait* child);
  void    inherit_low          (TTrait* mother, TTrait* father, TTrait* child);
  void    inherit_free         (TTrait* mother, TTrait* father, TTrait* child);

  ///@name Mutation
  ///@{
  void  (TTraitProto::* _mutate_func_ptr)(unsigned char** seq);
  void	mutate_NULL          (unsigned char** seq) { }
  void  mutate_SSM           (unsigned char** seq);
  void  mutate_SSM_single    (unsigned char** seq);
  void  mutate_KAM           (unsigned char** seq);
  void  mutate_KAM_single    (unsigned char** seq);
  void  mutate_RMM           (unsigned char** seq);
  void  mutate_RMM_single    (unsigned char** seq);
  void  mutate_IMM           (unsigned char** seq);
  void  mutate_IMM_single    (unsigned char** seq);
  ///@}


  ///@name Setters
  ///@{
  virtual void    set_type              (const string& x)   {_type=x;}
  virtual void    set_ploidy            (const int& x)      {_ploidy=x;}
  virtual void    set_nb_allele         (const int& x)      {_nb_allele=x;}
  virtual void    set_nb_locus          (const int& x)      {_nb_locus=x;}
  virtual void    set_absolute_index    (const int& x)      {_absolute_index=x;}
  virtual void    set_trait_index       (const int& x)      {_trait_index=x;}
  virtual void    set_locus_position    (double* x)         {_locus_position=x;}
  ///@}

  ///@name Getters
  ///@{
  virtual Metapop*get_popPtr()                const     {return _popPtr;}
  virtual string  get_type ()                 const     {return _type;}
  virtual int     get_ploidy()                const     {return _ploidy;}
  virtual int     get_nb_allele()             const     {return _nb_allele;}
  virtual int     get_nb_locus()              const     {return _nb_locus;}
  virtual int     get_absolute_index()        const     {return _absolute_index;}
  virtual double* get_locus_position()        const     {return _locus_position;}
  virtual double  get_locus_position(const int& i)const {return _locus_position[i];}
  ///@}

	virtual string get_info(){return "not set";}

  virtual int     get_trait_index()           const;
  virtual string  get_trait_indexStr()        const;
  virtual string  get_trait_indexStr_()       const;
  virtual string  get_trait_indexStr_t()      const;

  /**Inits the parameters, called by \c IndFactory::makePrototype().*/
  virtual   void            init (Metapop* pMetapop) = 0;

  /**Re-inits the parameters, called by \c IndFactory::makePrototype().*/
  virtual   void            reset () {return;}            // between replicates
	virtual   void            resetTotal ();

  /**Creates the trait of which it is the prototype, called by \c IndFactory::makePrototype().*/
  virtual   TTrait*         hatch () = 0;

  /**Returns a itself. \b Note: call the copy constructor and only copy the parameters state.*/
  virtual   TTraitProto* clone () = 0;

  /** These tasks are performed before/after each replicate */
  virtual void executeBeforeEachReplicate(const int& rep) {}
  virtual void executeAfterEachReplicate(const int& rep) {}

  /** These tasks are performed just before/after each generation */
  virtual void executeBeforeEachGeneration(const int& gen) {}
  virtual void executeAfterEachGeneration(const int& gen) {}

  virtual void read_allele_file (string filename);

  void    (TTraitProto::* _ini_allele_func_ptr)(unsigned char** seq, Patch* patch);
  void    ini_allele_freq_dist(unsigned char** seq, Patch* patch);
  void    ini_allele_freq_uniform(unsigned char** seq, Patch* patch);
  void    ini_allele_freq_monomorph(unsigned char** seq, Patch* patch);

  /** get any to be specified model value */
  virtual double get_model_value(){return my_NAN;}
};



#endif //TTRAIT_H

