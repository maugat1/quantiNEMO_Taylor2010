/** @file ttquanti.cpp
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

#include "ttquanti.h"
#include "stathandler.cpp"
#include "random.h"
#include "tree.cpp"
#include "output.h"
#include "ttree.h"

#include "newmatap.h"                // need matrix applications
#include "newmatio.h"                // need matrix output routines


TTQuantiFH*         TTQuantiProto::_writer;
TTQuantiFHvalue*    TTQuantiProto::_phenotyper;
TTQuantiFHvalue*    TTQuantiProto::_genotyper;


string TTQuantiProto::_ini_fstat_file;


// ----------------------------------------------------------------------------------------
// set_from_prototype (used during hatching)
// ----------------------------------------------------------------------------------------
void TTQuanti::set_from_prototype(TTraitProto* T){
  pProto = dynamic_cast<TTQuantiProto*> (T);
  pTraitProto = T;
}

// ----------------------------------------------------------------------------------------
// set_value
// ----------------------------------------------------------------------------------------
void
TTQuanti::set_value(){
	_genotype = (pProto->*(pProto->get_genotype_func_ptr))(sequence);
}   // set the genotype


// ----------------------------------------------------------------------------------------
// set_fitness_factor
// ----------------------------------------------------------------------------------------
void
TTQuanti::set_fitness_factor(){
	_fitness_factor = (pProto->*(pProto->get_fitnessFactor_func_ptr))(sequence);
}   // set the set_fitness_factor


// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTQuantiProto::TTQuantiProto( ):
  _stats(0), _allelicValues(0), _allelic_file(0),
	_dominanceValues(0), _dominance_file(0),
	_fitnessFactor_heterozygote(0), _fitnessFactor_homozygote(0), _fitnessFactor(0), _fitnessFactorTree(0),
	_phenoTree(0), get_genotype_func_ptr(0), get_fitnessFactor_func_ptr(0)
{
  _type = "quanti";
  _ploidy = 2;
  ini_paramset();
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTQuantiProto::TTQuantiProto(int i):
		_stats(0), _allelicValues(0), _allelic_file(0),
		_dominanceValues(0), _dominance_file(0),
		_fitnessFactor_heterozygote(0), _fitnessFactor_homozygote(0), _fitnessFactor(0), _fitnessFactorTree(0),
		_phenoTree(0), get_genotype_func_ptr(0), get_fitnessFactor_func_ptr(0)
{
	_trait_index = i;
	_type = "quanti"+get_trait_indexStr_();
	_ploidy = 2;
	ini_paramset();
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTQuantiProto::TTQuantiProto(const TTQuantiProto& T):
		 _stats(0), _allelicValues(T._allelicValues), _allelic_file(T._allelic_file),
		 _dominanceValues(T._dominanceValues), _dominance_file(T._dominance_file),
		 _fitnessFactor_heterozygote(T._fitnessFactor_heterozygote),
		 _fitnessFactor_homozygote(T._fitnessFactor_homozygote),
		 _fitnessFactor(T._fitnessFactor), _fitnessFactorTree(0), _phenoTree(0),
		 get_genotype_func_ptr(T.get_genotype_func_ptr),
		 get_fitnessFactor_func_ptr(T.get_fitnessFactor_func_ptr)
{
  _copyTraitPrototypeParameters(T);
	_ploidy = 2;
  ini_paramset();
}

// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TTQuantiProto::~TTQuantiProto ()
{
	if(_stats)          delete _stats;
	if(_writer)         {delete _writer;     _writer=NULL;}     //_writer is static!!!
	if(_phenotyper)     {delete _phenotyper; _phenotyper=NULL;} // _phenotyper is static!!!
	if(_genotyper)      {delete _genotyper;  _genotyper=NULL;}  // _genotyper is static!!!

	resetTotal();
}

// -----------------------------------------------------------------------------
// ini_paramset
// -----------------------------------------------------------------------------
void
TTQuantiProto::ini_paramset (){
	string trait = get_trait_indexStr_();

	set_paramset("quanti"+trait, false);

  // multiple traits
  add_parameter("quanti_nb_trait",INT2,false,false,0,0,"1");

  add_parameter("quanti_loci"+trait,INT2,true,false,0,0,"0");
  add_parameter("quanti_all"+trait,INT2,false,true,1,256,"255");

	add_parameter("quanti_allelic_var"+trait,DBL,false,false,0,0,"1");
  add_parameter("quanti_allelic_file"+trait, STR, false, false,0,0,"");
	add_parameter("quanti_ini_allele_model"+trait,INT2,false,true,0,1,"0");

  add_parameter("quanti_dominance_mean"+trait,DBL,false,false,0,0,"0");
  add_parameter("quanti_dominance_var"+trait,DBL,false,false,0,0,"0");
	add_parameter("quanti_dominance_file"+trait,STR,false,false,0,0,"");

	add_parameter("quanti_fitness_factor_heterozygote"+trait,DBL,false,false,0,0,"1");
	add_parameter("quanti_fitness_factor_homozygote"+trait,DBL,false,false,0,0,"1");

	add_parameter("quanti_epistatic_var"+trait,DBL,false,false,0,0,"0");
	add_parameter("quanti_epistatic_file"+trait,STR,false,false,0,0,"");

	add_parameter("quanti_mutation_rate"+trait,DBL,false,true,0,1,"0");
  add_parameter("quanti_mutation_shape"+trait,DBL,false,false,0,0,"0");
	add_parameter("quanti_mutation_model"+trait,INT2,false,true,0,1,"0");

  add_parameter("quanti_output"+trait,INT2,false,false,0,1,"0");

  add_parameter("quanti_save_genotype",INT2,false,true,0,2,"0");
  add_parameter("quanti_genot_sex",INT2,false,true,0,2,"0");
  add_parameter("quanti_genot_age",INT2,false,true,0,2,"0");
  add_parameter("quanti_genot_dir",STR,false,false,0,0,"");
  add_parameter("quanti_genot_logtime",INT2,false,false,0,0,"1",true);
  add_parameter("quanti_genot_script",STR,false,false,0,0,"");

  add_parameter("quanti_save_phenotype",INT2,false,true,0,2,"0");
  add_parameter("quanti_phenot_sex",INT2,false,true,0,2,"0");
	add_parameter("quanti_phenot_age",INT2,false,true,0,2,"0");
	add_parameter("quanti_phenot_dir",STR,false,false,0,0,"");
	add_parameter("quanti_phenot_logtime",INT2,false,false,0,0,"1",true);
  add_parameter("quanti_phenot_script",STR,false,false,0,0,"");

  add_parameter("quanti_save_geno_value",INT2,false,true,0,2,"0");
  add_parameter("quanti_geno_value_sex",INT2,false,true,0,2,"0");
  add_parameter("quanti_geno_value_age",INT2,false,true,0,2,"0");
  add_parameter("quanti_geno_value_dir",STR,false,false,0,0,"");
  add_parameter("quanti_geno_value_logtime",INT2,false,false,0,0,"1",true);
  add_parameter("quanti_geno_value_script",STR,false,false,0,0,"");

  add_parameter("quanti_selection_model"+trait,INT2,false,true,0,2,"0");

  // environment
  add_parameter("quanti_heritability"+trait,DBL,false,false,0,0,"1");
  add_parameter("quanti_heritability_fem"+trait,DBL,false,false,0,0,"1");
  add_parameter("quanti_heritability_mal"+trait,DBL,false,false,0,0,"1");
  add_parameter("quanti_environmental_proportion"+trait,DBL,false,true,0,1,"1");
  add_parameter("quanti_environmental_model"+trait,INT2,false,true,0,4,"0");

  // genetic map
	add_parameter("quanti_loci_positions"+trait,MAT_VAR,false,false,0,0,"");
	add_parameter("quanti_loci_positions_random"+trait,MAT,false,false,0,0,"");

  add_parameter("quanti_va_model"+trait,INT2,false,true,0,1,"0");

  add_parameter("quanti_ini_genotypes",STR,false,false,0,0,"");
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
// used between replicates
void TTQuantiProto::reset(){
	if(!_dominance_file) ARRAY::create_3D(_dominanceValues, _nb_locus, _nb_allele, _nb_allele, (double)my_NAN);
	if(_phenoTree){
		delete _phenoTree;
		_phenoTree = new Tree<unsigned char>(_nb_locus, _nb_allele);
	}
	if(_fitnessFactorTree){
		delete _fitnessFactorTree;
		_fitnessFactorTree = new Tree<unsigned char>(_nb_locus, _nb_allele);
	}
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
// used between simulations
void TTQuantiProto::resetTotal(){
	if(_phenoTree)        {delete _phenoTree;		_phenoTree=NULL;}
	if(_fitnessFactorTree){delete _fitnessFactorTree;		_fitnessFactorTree=NULL;}
	if(_dominanceValues) ARRAY::delete_3D(_dominanceValues, _nb_locus, _nb_allele);
	if(_fitnessFactor)   ARRAY::delete_3D(_fitnessFactor, _nb_locus, _nb_allele);
	if(_fitnessFactor_heterozygote) {delete[] _fitnessFactor_heterozygote; _fitnessFactor_heterozygote=NULL;}
	if(_fitnessFactor_homozygote) {delete[] _fitnessFactor_homozygote; _fitnessFactor_homozygote=NULL;}

  if(_allelicValues){
    if(_nb_locus > 1 && _allelicValues[0] != _allelicValues[1]){   // each locus has its one array
       ARRAY::delete_2D(_allelicValues, _nb_locus);
    }
    else{                                                        // all loci point to the first array
		  delete[] _allelicValues[0];
		  delete[] _allelicValues;
		  _allelicValues = NULL;
    }
	}

	TTraitProto::resetTotal();
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTQuantiProto::init (Metapop* pMetapop)
{
  _popPtr = pMetapop;
  string trait = get_trait_indexStr_();

  if(get_trait_index() <= 1){    // for sequential parameters
    if(_writer)     {delete _writer;     _writer = NULL;}
    if(_phenotyper) {delete _phenotyper; _phenotyper = NULL;}
    if(_genotyper)  {delete _genotyper;  _genotyper = NULL;}

		// FSTAT file for initial genotypes?
		_ini_fstat_file = get_parameter("quanti_ini_genotypes")->get_arg();
	}

  // get the parameters (they are validated if (and only if) they are used)
	_nb_allele          = (int)get_parameter_value("quanti_all"+trait);
  _nb_locus           = (int)get_parameter_value("quanti_loci"+trait);
  _mut_model          = (int)get_parameter_value("quanti_mutation_model"+trait);
  _epistatic_sd       = get_parameter_value("quanti_epistatic_var"+trait);
  _dominance_mean     = get_parameter_value("quanti_dominance_mean"+trait);
  _dominance_sd       = get_parameter_value("quanti_dominance_var"+trait);
  _output             = (int)get_parameter_value("quanti_output"+trait);
  _selection_model    = (int) get_parameter_value("quanti_selection_model"+trait);
  _Va_model           = (int) get_parameter_value("quanti_va_model"+trait);
  _ini_allele_model   = (int) get_parameter_value("quanti_ini_allele_model"+trait);


  //----------------------------------------------------------------------------
  // genetic map
  if(_locus_position) delete[] _locus_position;
  if(get_parameter("quanti_loci_positions"+trait)->isSet()){
    _inherit_func_ptr = &TTraitProto::inherit_low;
    _locus_position = new double[_nb_locus];
    setGeneticMapFixed(get_parameter("quanti_loci_positions"+trait)->get_matrixVar());
  }
  else if(get_parameter("quanti_loci_positions_random"+trait)->isSet()){
    _inherit_func_ptr = &TTraitProto::inherit_low;
    _locus_position = new double[_nb_locus];
    _random_genetic_map_matrix = get_parameter("quanti_loci_positions_random"+trait)->get_matrix();
    setGeneticMapRandom(_random_genetic_map_matrix);
  }
  else {
    _inherit_func_ptr = &TTraitProto::inherit_free;
    _locus_position = NULL;
  }

  //----------------------------------------------------------------------------
  // set mutation function
  if(set_mutation_rates("quanti", trait)){ // if mutation rates differ among loci
    switch(_mut_model) {
      case 0: _mutate_func_ptr = &TTraitProto::mutate_RMM_single;    break;
      case 1: _mutate_func_ptr = &TTraitProto::mutate_IMM_single;    break;
    }
  }
  else{         // all loci have the same mutation rate
    if(!_mut_rate_mean) _mutate_func_ptr = &TTraitProto::mutate_NULL;       // mutation rate = 0
    else{
      switch(_mut_model) {
        case 0: _mutate_func_ptr = &TTraitProto::mutate_RMM;    break;
        case 1: _mutate_func_ptr = &TTraitProto::mutate_IMM;    break;
      }
    }
  }

  //----------------------------------------------------------------------------
  // read the allelic file if present
  string allelicFile = get_parameter("quanti_allelic_file"+trait)->get_arg();  // if empty parameter is not set
  if(!allelicFile.empty()){
    if(!STRING::file_exists(allelicFile)){
      allelicFile = FileHandler::get_base_path()+allelicFile;
      if(!STRING::file_exists(allelicFile)) fatal("Allelic file '%s' (quanti) cannot be found!\n", allelicFile.c_str());
    }
		read_allele_file(allelicFile);
  }

  //----------------------------------------------------------------------------
  // set the allelic effects if not passed by the file
  set_allelicValues(trait);

	//----------------------------------------------------------------------------
	// set the mutation rate probabilities
  set_mutationFreq(trait);

	//----------------------------------------------------------------------------
	// set the inital allele frequencies
	set_initAlleleFreq(trait);

	//----------------------------------------------------------------------------
  // dominance effect
	string dominanceFile = get_parameter("quanti_dominance_file"+trait)->get_arg();
  if(!dominanceFile.empty()){      // if the parameter is set
    if(!STRING::file_exists(dominanceFile)){
      dominanceFile = FileHandler::get_base_path()+dominanceFile;
      if(!STRING::file_exists(dominanceFile)) fatal("Dominance file '%s' (quanti) cannot be found!\n", dominanceFile.c_str());
    }
		read_locus_file(dominanceFile);
		_dominance_file = true;
		get_locus_genotype_func_ptr = &TTQuantiProto::get_locus_genotype_dominance_array;
	}
  else {                           // no dominance effects set by file
    if(_dominance_sd){
      if(_dominance_sd < 0) fatal("Parameter 'quanti_dominance_var' must be positive!\n");
      get_locus_genotype_func_ptr = &TTQuantiProto::get_locus_genotype_dominance_array;
      ARRAY::create_3D(_dominanceValues, _nb_locus, _nb_allele, _nb_allele, (double)my_NAN);
      _dominance_sd = sqrt(_dominance_sd);
    }
    else if(_dominance_mean) get_locus_genotype_func_ptr = &TTQuantiProto::get_locus_genotype_dominance_single;
    else                     get_locus_genotype_func_ptr = &TTQuantiProto::get_locus_genotype_additive;
  }

	//----------------------------------------------------------------------------
	// epistatic effect
	string epistaticFile = get_parameter("quanti_epistatic_file"+trait)->get_arg();
	if(!epistaticFile.empty()){       // epistatic effects set by file
		if(!STRING::file_exists(epistaticFile)){
			epistaticFile = FileHandler::get_base_path()+epistaticFile;
			if(!STRING::file_exists(epistaticFile)) fatal("Epistatic file '%s' (quanti) cannot be found!\n", epistaticFile.c_str());
		}
		read_genome_file(epistaticFile);
	}

	if(_phenoTree) get_genotype_func_ptr = &TTQuantiProto::get_genotype_epistatic;
	else if(_epistatic_sd){          // epistatic effects set by its variance
		if(_epistatic_sd < 0) fatal("Parameter 'quanti_epistatic_var' must be positive!\n");
		get_genotype_func_ptr = &TTQuantiProto::get_genotype_epistatic;
		_epistatic_sd = sqrt(_epistatic_sd);
		if(!_phenoTree) _phenoTree = new Tree<unsigned char>(_nb_locus, _nb_allele);
	}
	else get_genotype_func_ptr = &TTQuantiProto::get_genotype_none;   // no epistatic effects


	//----------------------------------------------------------------------------
	// set the fitness factor
	set_fitnessFactor(trait, "quanti_fitness_factor_heterozygote", _fitnessFactor_heterozygote);
	set_fitnessFactor(trait, "quanti_fitness_factor_homozygote", _fitnessFactor_homozygote);
	if(_fitnessFactor)          get_fitnessFactor_func_ptr = &TTQuantiProto::get_fitnessFactor_explicit;
	else if(_fitnessFactorTree) get_fitnessFactor_func_ptr = &TTQuantiProto::get_fitnessFactor_genome;
	else if(_fitnessFactor_homozygote || _fitnessFactor_heterozygote){
															get_fitnessFactor_func_ptr = &TTQuantiProto::get_fitnessFactor_global;
	}
}

// ----------------------------------------------------------------------------------------
// executeBeforeEachGeneration
// ----------------------------------------------------------------------------------------
void
TTQuantiProto::executeBeforeEachGeneration(const int& gen){
  map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);

  // if it is a temporal paramter
	if(pParam){
    // check if a change has to be made
    map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
    if(pMap){
      // iterate through the map and performe the updates
      map<string, Param*>::iterator pos = pMap->begin();
      string trait = get_trait_indexStr_();

      for(; pos != pMap->end(); ++pos){
        if(pos->first == "quanti_genot_logtime" && _writer){
          _writer->set_GenerationOccurrence((unsigned int)pos->second->get_value());
        }
        else if(pos->first == "quanti_phenot_logtime" && _phenotyper){
          _phenotyper->set_GenerationOccurrence((unsigned int)pos->second->get_value());
        }
        else if(pos->first == "quanti_geno_value_logtime" && _genotyper){
          _genotyper->set_GenerationOccurrence((unsigned int)pos->second->get_value());
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void
TTQuantiProto::executeAfterEachReplicate(const int& rep){
  if(_output){
    string dir = _popPtr->get_service()->getSimfolder();
    print_allelic_values(dir+"allelic_values");
    print_dominance_values(dir+"dominance_values");
    print_epistatic_values(dir+"epistatic_values");
  }
}

// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void
TTQuantiProto::executeBeforeEachReplicate(const int& rep){
  if(rep!= 1 && _random_genetic_map_matrix) setGeneticMapRandom(_random_genetic_map_matrix);
}

// ----------------------------------------------------------------------------------------
// get_info
// ----------------------------------------------------------------------------------------
string
TTQuantiProto::get_info()
{
	string text;

	if(_trait_index) text = STRING::int2str(_trait_index) + ". quantitative trait: ";
	else             text = "Quantitative trait: ";

	text +=	STRING::int2str(_nb_locus) + " loci; "
	     + 	STRING::int2str(_nb_allele) + " alleles;";

	switch(_selection_model){
		case 0: text += " stabilizing selection";    break;
		case 1: text += " directional selection";    break;
		case 2: text += " no selection";             break;
	}

	return text;
}
// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void
TTQuantiProto::loadFileServices  (FileServices* loader)
{
	// genotype
	int choice = (int)get_parameter_value("quanti_save_genotype");
	if(choice) {
		if(!_writer) _writer = new TTQuantiFH(this);
		int logtime = (int)get_parameter_value("quanti_genot_logtime");
		if(logtime<1) fatal("Parameter 'quanti_genot_logtime' must be 1 or larger!\n");
		_writer->set(true, true, 1,
								 logtime,
								 0,
								 get_parameter("quanti_genot_dir")->get_arg(),
								 get_parameter("quanti_genot_script")->get_arg(),
								 (int)get_parameter_value("quanti_genot_sex"),
								 (int)get_parameter_value("quanti_genot_age"),
								 this, choice);

		if(get_trait_index()<=1) loader->attach(_writer);     // only the first one
	}
	else if(_writer) {
		delete _writer;
		_writer = NULL;
	}

	// phenotype
	choice = (int)get_parameter_value("quanti_save_phenotype");
	if(choice) {
		if(!_phenotyper) _phenotyper = new TTQuantiFHvalue(this);
		int logtime = (int)get_parameter_value("quanti_phenot_logtime");
		if(logtime<1) fatal("Parameter 'quanti_phenot_logtime' must be 1 or larger!\n");
			_phenotyper->set(true, true, 1,
											 logtime,
											 0,
											 get_parameter("quanti_phenot_dir")->get_arg(),
											 get_parameter("quanti_phenot_script")->get_arg(),
											 (int)get_parameter_value("quanti_phenot_sex"),
											 (int)get_parameter_value("quanti_phenot_age"),
											 this, choice);
		_phenotyper->set_getter(1); // phenotype
		if(get_trait_index()<=1) loader->attach(_phenotyper);
	}
	else {
		if(_phenotyper) {
			delete _phenotyper;
			_phenotyper = NULL;
		}
	}

	//genotypic value
	choice = (int)get_parameter_value("quanti_save_geno_value");
	if(choice) {
		if(!_genotyper) _genotyper = new TTQuantiFHvalue(this);
		int logtime = (int)get_parameter_value("quanti_geno_value_logtime");
		if(logtime<1) fatal("Parameter 'quanti_geno_value_logtime' must be 1 or larger!\n");
			_genotyper->set(true, true, 1,
											 logtime,
											 0,
											 get_parameter("quanti_geno_value_dir")->get_arg(),
											 get_parameter("quanti_geno_value_script")->get_arg(),
											 (int)get_parameter_value("quanti_geno_value_sex"),
											 (int)get_parameter_value("quanti_geno_value_age"),
											 this, choice);
			_genotyper->set_getter(0); // genotype
			if(get_trait_index()<=1) loader->attach(_genotyper);
	}
	else {
		if(_genotyper) {
      delete _genotyper;
			_genotyper = NULL;
    }
	}
}

// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void
TTQuantiProto::loadStatServices  (StatServices* loader)
{
  if(_stats)  delete _stats;
  _stats = new TTQuantiSH(this);
  loader->attach(_stats);
}

// ----------------------------------------------------------------------------------------
// clone
// ----------------------------------------------------------------------------------------
TTQuanti*
TTQuantiProto::hatch ()
{
  TTQuanti* new_trait = new TTQuanti();
  new_trait->set_from_prototype(this);
  return new_trait;
}

// ----------------------------------------------------------------------------------------
// check_allelic_array
// ----------------------------------------------------------------------------------------
/** check the mutation probabilities for the mutation model IMM */
void
TTQuantiProto::check_mutationValues_IMM()
{
  for(int l = 0; l < _nb_locus; ++l){
    if(abs(_mutationFreq[l][_nb_allele/2]) > 1e-4){
      warning("The mutation model IMM requiers that the mutation probability of the allele with effect=0 is zero (locus %i, allele %i): Automatically adjusted!\n",
        l+1, 1+_nb_allele/2);
    }
  }
}
// ----------------------------------------------------------------------------------------
// check_allelic_array
// ----------------------------------------------------------------------------------------
/** check the allelic effect array for the mutation model IMM */
void
TTQuantiProto::check_allelicValues_IMM()
{
  // the number of alleles has to be odd
  if(!(_nb_allele%2)) fatal("The mutation model IMM requires an odd number of alleles!\n");

  for(int l=0; l<_nb_locus; ++l){
    // alleles have to be regularely spaced
    double step = _allelicValues[l][1] - _allelicValues[l][0];
    for(int a = 2; a < _nb_allele; ++a){
      if(abs(_allelicValues[l][a] - _allelicValues[l][a-1] - step) > 1e-4){
        fatal("Allelic values: The effects of the alleles have to be regularly distributed if mutation model IMM is used(locus %i)!\n", l);
      }
    }

    // the effects have to be symmetrix, i.e. the middle allele has to be 0
    if(abs(_allelicValues[l][_nb_allele/2]) > 1e-4){
      fatal("Allelic values: For the IMM mutation model the allelic effects have to be symmetrically distrubted (allele %i has to be zero)!\n", 1+_nb_allele/2);
    }
  }
}

// ----------------------------------------------------------------------------------------
// create_regular_spaced_array
// ----------------------------------------------------------------------------------------
/** create a regular stepped array from -range to +range (inlcuding the edges)
  ** due to round-off errors the array is made symmetrically starting in the middle
  */
void
TTQuantiProto::create_regular_spaced_array(double* array, const int& size, double half_range)
{
  assert(array);
	double step = 2.*half_range/(size-1.);    // compute step size
  double curVal;
  double *a_up, *a_down;
  double *end = array + size;

  if(size%2){    // odd  number of alleles
    curVal = 0;
    a_up = a_down = array + (int)size/2;
  }
  else{          // event number of alleles
    curVal = step/2.;
    a_up = array + (int)size/2;
    a_down = a_up - 1;
  }

	while(1){
    *a_up   = curVal;
    *a_down = -curVal;

    // next element  (done so complicate due to code guard memory leak detection)
    ++a_up;
    if(a_up == end) break;
    --a_down;
    curVal += step;
  }
}

// ----------------------------------------------------------------------------------------
// set_allelic_valuesFreq
// ----------------------------------------------------------------------------------------
/** compute the frequencies for the normal distribution N(0, sd) for the values
	* of array effect_array and store them in the array
  * Caution: they do not yet sum up to 1!!!
	*/
void
TTQuantiProto::compute_frequencies(double* array, double* effect_array, int& size, double& sd)
{
	assert(array);
	assert(effect_array);
  double* end = array + size;
	for(; array != end; ++array, ++effect_array){
		*array = SimRunner::r.ProbDensNorm(*effect_array, 0, sd);
	}
}

// ----------------------------------------------------------------------------------------
// set_allelic_valuesFreq
// ----------------------------------------------------------------------------------------
/** read the allelic values and their frequencies from a file (used for mutation model IMM and RMM)
  */
void
TTQuantiProto::read_allele_file (string filename)
{
	#ifdef _DEBUG
	 message("  TTQuantiProto::read_allele_file (%s) ...",filename.c_str());
	#endif

	// read the file
  TMatrix mat;
  map<string, int> fileInfo;
  bool hasLocus = true;
  int i, l, a, nbValues=0;
  vector<map<int, string> >* pStrMatrix = NULL;
  double* array = NULL;

  // read the file
  fileInfo = mat.read_matrix(filename);
  int nbLines = mat.getNbRows();
  int nbCols  = mat.getNbCols();

  const int colsSize = 5;
  int cols[colsSize]        = {my_NAN,my_NAN,my_NAN,my_NAN,my_NAN};  // 0: locus; 1: allele; 2: allelic_value; 3: mut_freq; 4: ini_freq;
  string colsName[colsSize] = {"col_locus","col_allele","col_allelic_value","col_mut_freq","col_ini_freq"};

	// get col values form file info
	map<string, int>::iterator pos;
	for(i=0; i<colsSize; ++i){
		pos = fileInfo.find(colsName[i]);
		if(pos != fileInfo.end()){
			cols[i] = pos->second - 1;
			fileInfo.erase(pos);
		}
	}

	// ouput all not used columns
	for(pos = fileInfo.begin(); pos != fileInfo.end(); ++pos){
		warning("Allelic file for trait '%s' (%s): Column '%s' (%i) is not a valid column names!\n",
						_type.c_str(), filename.c_str(), pos->first.c_str(), pos->second);
	}

	// check if all cols are available and also dimensions of the matrix are met adn create the arrays
	if(cols[0] == my_NAN) hasLocus = false; // fatal("Allelic file (%s): Locus column is missing!", filename.c_str());
	if(cols[1] == my_NAN) fatal("Allelic file (%s): Allele column is missing!", filename.c_str());
	if(cols[2] != my_NAN) {
		ARRAY::create_2D(_allelicValues, _nb_locus, _nb_allele, (double)my_NAN);
		_allelic_file = true;
	}
	if(cols[3] != my_NAN){    // mutation frequency
		ARRAY::create_2D(_mutationFreq, _nb_locus, _nb_allele, (double)my_NAN);
		_mutationFreqFile = true;
	}

	if(cols[4] != my_NAN){     // inital frequencies
		pStrMatrix = set_ini_allelicArray(&mat, cols[4]);  // _initAlleleFreqCols is set
		if(_initAlleleFreqCols) array = new double[_initAlleleFreqCols];
	}
	else _initAlleleFreqCols = 0;

	for (i=0; i<colsSize; ++i){
		if(cols[i]!=my_NAN && cols[i]>nbCols) fatal("Allelic file (%s): File info is wrong: matrix has only %i columns and not %i!\n", filename.c_str(), nbCols, cols[i]);
	}

	// read the table: copy the allelic effects and their frequencies
	for(i=0; i<nbLines; ++i){
		if(hasLocus) l = (int)mat.get(i,cols[0]);    // if locus specific values
		else         l = 1;                          // if all loci have the same settings
		a = (int)mat.get(i,cols[1]);    // allele

		// test if out of range
		if(hasLocus && l>_nb_locus){
			fatal("Allelic values: Locus %i (allele %i) is out of range (only %i loci specified)!\n",
						l, a, _nb_locus);
		}
		if(a>_nb_allele){
			fatal("Allelic values: Allele %i of locus %i is out of range (only %i alleles specified)!\n",
						a, l, _nb_allele);
		}

		// get the array position (-1 as arrays start at 0)
		--l; --a;

		// check if the combination was already read
		if(   (_allelicValues  && _allelicValues[l][a]     != my_NAN)
			 || (_mutationFreq   && _mutationFreq[l][a]      != my_NAN)
			 || (_initAlleleFreq && _initAlleleFreq[0][l][a] != my_NAN)){

			fatal("Allelic values: Allelic value was already specified for locus %i and allele %i!\n", l+1, a+1);
		}

		// set the values
		if(hasLocus){         // locus specific settings
			set_allelicValues(&mat, pStrMatrix, array, i, l, a, cols);
			++nbValues;
		}
		else {                // all loci have the same settings
			for(l = 0; l < _nb_locus; ++l){
				set_allelicValues(&mat, pStrMatrix, array, i, l, a, cols);
				++nbValues;
			}
		}
	}

	// check if we have enough values
	int totValues = _nb_locus*_nb_allele;
  if(nbValues != totValues){
    fatal("Allelic file: %i of %i allelic values are not set by the file '%s'!\n",
          totValues-nbValues, totValues, filename.c_str());
  }

  if(_initAlleleFreqCols) delete[] array;

  // control the arrays and make them cumulative
  if(_allelic_file && _mutate_func_ptr == &TTQuantiProto::mutate_IMM){
    check_allelicValues_IMM();
  }

  if(_mutationFreqFile){
    if(_mutate_func_ptr == &TTQuantiProto::mutate_IMM){
      check_mutationValues_IMM();
    }
    for(int l = 0; l < _nb_locus; ++l){
      ARRAY::make_frequency(_mutationFreq[l], _nb_allele);  // adjust the sum to 1
      ARRAY::cumulative(_mutationFreq[l], _nb_allele);      // make it cumulative
    }
  }

  if(_initAlleleFreqCols){
    for(int c = 0; c < _initAlleleFreqCols; ++c) {
      for(int l = 0; l < _nb_locus; ++l){
        ARRAY::make_frequency(_initAlleleFreq[c][l], _nb_allele);  // adjust the sum to 1
        ARRAY::cumulative(_initAlleleFreq[c][l], _nb_allele);      // make it cumulative
      }
    }
	}

	#ifdef _DEBUG
	 message("done!\n");
	#endif
}

// ----------------------------------------------------------------------------------------
// set_mutationFreq
// ----------------------------------------------------------------------------------------
/** set the probabilities to mutate to a certain allele */
void
TTQuantiProto::set_mutationFreq(const string& trait){
	if(_mutationFreqFile || _mutate_func_ptr == &TTraitProto::mutate_NULL) return;

  // the allelic effects have to be fixed automatically
  Param* p = get_parameter("quanti_allelic_var"+trait);

  // if the variance of the alleic effects varies between loci
  if(p->is_matrix()){     // each locus has an individual mutation freq array
    // get the values
    int nb = p->get_matrix()->get_dims(NULL);
    double* values = p->get_matrix()->get();

    // check the number of variances with the number of loci
    if(nb>_nb_locus) warning("There are more allelic effect variances than loci defined! Only a part of the allelic effect variances is considered!\n");
    else if(_nb_locus % nb) warning("The number of allelic effect variances is not a entire subset of the number of loci!\n");


    // create the mutation freq matrix: each locus has other settings
		ARRAY::create_2D(_mutationFreq, _nb_locus, _nb_allele);

    // get the variances for each locus
    double sd;
    for(int l=0; l<_nb_locus; ++l){
      sd = values[l % nb];
		  if(sd<0) fatal("Parameter 'quanti_allelic_var%s' must be positive!\n", trait.c_str());
		  sd = sqrt(sd);

      // compute the density
		  compute_frequencies(_mutationFreq[l], _allelicValues[l], _nb_allele, sd);

      //if IMM the "middle allele" has to have a freq of zero!
	  	if(_mutate_func_ptr == &TTraitProto::mutate_IMM){
		    _mutationFreq[l][_nb_allele/2] = 0;
		  }

      // check that the frequencies sum up to 1 and make the frequencies cumulative
      ARRAY::make_frequency(_mutationFreq[l], _nb_allele);
      ARRAY::cumulative(_mutationFreq[l], _nb_allele);
    }
    return;
  }

  // all loci have the same settings -> all loci point to the same allele array (caution, when deleting)
  double sd  = get_parameter_value("quanti_allelic_var"+trait);
	if(sd < 0) fatal("Parameter 'quanti_allelic_var%s' must be positive!\n", trait.c_str());
	sd = sqrt(sd);

	ARRAY::create_1D(_mutationFreq, _nb_locus);
	_mutationFreq[0] = new double[_nb_allele];         // create the array for the first locus
	for(int l = 0; l < _nb_locus; ++l){
	  _mutationFreq[l] = _mutationFreq[0];             // all loci point to this first locus array
	}

  // fill the array with values
	compute_frequencies(_mutationFreq[0], _allelicValues[0], _nb_allele, sd);

  //if IMM the "middle allele" has to have a freq of zero!
	if(_mutate_func_ptr == &TTraitProto::mutate_IMM){
	  _mutationFreq[0][_nb_allele/2] = 0;
	}

  // check that the frequencies sum up to 1 and make the frequencies cumulative
  ARRAY::make_frequency(_mutationFreq[0], _nb_allele);
  ARRAY::cumulative(_mutationFreq[0], _nb_allele);
}


// ----------------------------------------------------------------------------------------
// set_initAlleleFreq
// ----------------------------------------------------------------------------------------
/** set the inital allelic values for a maximal "polymorphism" */
void
TTQuantiProto::set_initAlleleFreq(const string& trait){
  // if the initial frequencies are given explicitly by a file
	if(_initAlleleFreqFile){
		_ini_allele_func_ptr = &TTraitProto::ini_allele_freq_dist;
    return;
  }

  // if all populations are monomorph
  if(_ini_allele_model){
		_ini_allele_func_ptr = &TTraitProto::ini_allele_freq_monomorph;
    return;
  }

  // if the populations are maximal polymorph
  _ini_allele_func_ptr = &TTraitProto::ini_allele_freq_dist;

  // the allelic effects have to be fixed automatically
  Param* p = get_parameter("quanti_allelic_var"+trait);

  // if all loci have the same polymorphism
  if(p->is_matrix()){     // each locus has another setting
    // get the values
    int nb = p->get_matrix()->get_dims(NULL);
    double* values = p->get_matrix()->get();

    // check the number of variances with the number of loci
    if(nb>_nb_locus) warning("There are more allelic effect variances than loci defined! Only a part of the allelic effect variances is considered!\n");
    else if(_nb_locus % nb) warning("The number of allelic effect variances is not a entire subset of the number of loci!\n");


    // create the initial allele frequency matrix: each locus has other settings
    _initAlleleFreqCols = 1;      // all patches have the same initial allele frequencies
		ARRAY::create_3D(_initAlleleFreq, _initAlleleFreqCols, _nb_locus, _nb_allele); //_initAlleleFreq[patch][locus][allele]

    // get the variances for each locus
    double sd;
    for(int l=0; l<_nb_locus; ++l){
      sd = values[l % nb];
		  if(sd<0) fatal("Parameter 'quanti_allelic_var%s' must be positive!\n", trait.c_str());
		  sd = sqrt(sd);

      // compute the density
		  compute_frequencies(_initAlleleFreq[0][l], _allelicValues[l], _nb_allele, sd);

      // check that the frequencies sum up to 1 and make the frequencies cumulative
      ARRAY::make_frequency(_initAlleleFreq[0][l], _nb_allele);
      ARRAY::cumulative(_initAlleleFreq[0][l], _nb_allele);
    }
    return;
  }

  // all loci have the same settings -> all loci point to the same allele array (caution, when deleting)
  double sd  = get_parameter_value("quanti_allelic_var"+trait);
	if(sd < 0) fatal("Parameter 'quanti_allelic_var%s' must be positive!\n", trait.c_str());
	sd = sqrt(sd);

  // the inital freuqncies are the same as the mutation propablities for the RMM
  _initAlleleFreqCols = 1;      // all patches have the same initial allele frequencies
	ARRAY::create_2D(_initAlleleFreq, _initAlleleFreqCols, _nb_locus);  // all are the same -> use only one locus array (caution, when deleting)
	_initAlleleFreq[0][0] = new double[_nb_allele];   // create the array for the first locus
	for(int l = 1; l < _nb_locus; ++l){
	  _initAlleleFreq[0][l] = _initAlleleFreq[0][0];       // all loci point to this first locus array
	}

  // fill the array with values
  compute_frequencies(_initAlleleFreq[0][0], _allelicValues[0], _nb_allele, sd);

  // check that the frequencies sum up to 1 and make the frequencies cumulative
  ARRAY::make_frequency(_initAlleleFreq[0][0], _nb_allele);
  ARRAY::cumulative(_initAlleleFreq[0][0], _nb_allele);
}

// ----------------------------------------------------------------------------------------
// set_fitnessFactor
// ----------------------------------------------------------------------------------------
void
TTQuantiProto::set_fitnessFactor(const string& trait, string name, double* &array)
{
	Param* p = get_parameter(name+trait);
	if(array) {delete[] array; array=NULL;}
	if(p->is_matrix()){      // locus specific settings
		// get the values
		int nb = p->get_matrix()->get_dims(NULL);
		double* values = p->get_matrix()->get();

		// check the number of values with the number of loci
		if(nb>_nb_locus) warning("There are more fitness factors ('%s') than loci defined! Only a part of the fitness factors is considered!\n", name.c_str());
		else if(_nb_locus % nb) warning("The number of fitness factors ('%s') is not a entire subset of the number of loci!\n", name.c_str());

		// create the array
		array = new double[_nb_locus];

		// get the values for each locus
		for(int l=0; l<_nb_locus; ++l){
			array[l] = values[l % nb];
		}
	}
	else{              	// common settings for all loci
		double value = p->get_value();
		if(value != 1){   // create the array only if it is not the default value
			array = new double[_nb_locus];
			for(int l=0; l<_nb_locus; ++l){
				array[l] = value;
			}
		}
	}
}

// ----------------------------------------------------------------------------------------
// set_allelicValues
// ----------------------------------------------------------------------------------------
void
TTQuantiProto::set_allelicValues(const string& trait){
	if(_allelic_file) return;       // the allelic effects are set by the alleic file

	// controls
	if(_nb_allele < 2) fatal("%i alleles makes no sense to the simulation!\n", _nb_allele);
	if(_mut_model == 1){		// IMM
		if(!(_nb_allele%2)) fatal("The IMM model requires an odd number of alleles!\n");
		if(_nb_allele < 51) fatal("The IMM model requires at least 51 alleles per locus!\n");
	}

	// the allelic effects have to be fixed automatically
	Param* p = get_parameter("quanti_allelic_var"+trait);

	// if the variance of the alleic effects varies between loci
	if(p->is_matrix()){
		// get the values
		int nb = p->get_matrix()->get_dims(NULL);
		double* values = p->get_matrix()->get();

		// check the number of variances with the number of loci
		if(nb>_nb_locus) warning("There are more allelic effect variances than loci defined! Only a part of the allelic effect variances is considered!\n");
		else if(_nb_locus % nb) warning("The number of allelic effect variances is not a entire subset of the number of loci!\n");

		// create the allelic file: each locus has other settings
		ARRAY::create_2D(_allelicValues, _nb_locus, _nb_allele);

		// get the variances for each locus
		double range, sd;
		for(int l=0; l<_nb_locus; ++l){
			sd = values[l % nb];
			if(sd<0) fatal("Parameter 'quanti_allelic_var%s' must be positive!\n", trait.c_str());
			sd = sqrt(sd);

			// create the regularly spaced allelic effect array
			if(_nb_allele < 6) range = sd * (_nb_allele-1);
			else range = sd * ((_mut_model == 1) ? 20 : 6);   //IMM?
			create_regular_spaced_array(_allelicValues[l], _nb_allele, range);
		}
		return;
	}

	// all loci have the same allelic variance
	double sd  = get_parameter_value("quanti_allelic_var"+trait);
	if(sd < 0) fatal("Parameter 'quanti_allelic_var%s' must be positive!\n", trait.c_str());
	sd = sqrt(sd);

	// controls specific for the IMM model:
	if(_mut_model == 1){    // IMM?
		if(!(_nb_allele%2)) fatal("The IMM model requires an odd number of alleles!\n");
		if(_nb_allele < 51) fatal("The IMM model requires at least 51 alleles per locus!\n");
	}

	// all loci have the same settings -> all loci point to the same allele array (caution, when deleting)
	ARRAY::create_1D(_allelicValues, _nb_locus);
	_allelicValues[0] = new double[_nb_allele];    // create the allele array for the first locus
	for(int l = 0; l < _nb_locus; ++l){
		_allelicValues[l] = _allelicValues[0];      // all loci point to this first allele array
	}

	// create the regularly spaced allelic effect array
	double range;
	if(_nb_allele < 2) fatal("%i alleles makes no sense to the simulation!\n", _nb_allele);
	if(_nb_allele < 6) range = sd * (_nb_allele-1);
	else range = sd * ((_mut_model == 1) ? 20 : 6);  // IMM?
	create_regular_spaced_array(_allelicValues[0], _nb_allele, range);
}

// ----------------------------------------------------------------------------------------
// set_allelicValues
// ----------------------------------------------------------------------------------------
/** set the values for a single row */
void
TTQuantiProto::set_allelicValues(TMatrix* mat, vector<map<int, string> >* pStrMatrix, double* array,
						const unsigned int& i, const unsigned int& l, const unsigned int& a, int* cols){
	if(_allelicValues)  _allelicValues[l][a]  = mat->get(i,cols[2]); // allelic value
  if(_mutationFreq)  _mutationFreq[l][a]  = mat->get(i,cols[3]);   // mut  frequency
  if(_initAlleleFreq){
    if(pStrMatrix){
      STRING::strMatrix2array((*pStrMatrix)[i].find(cols[4])->second, array, _initAlleleFreqCols);
      for(int p=0; p<_initAlleleFreqCols; ++p){
		    _initAlleleFreq[p][l][a] = array[p]; // init frequencies per population
      }
    }
    else _initAlleleFreq[0][l][a] = mat->get(i,cols[4]); // init frequencies
  }
}

// ----------------------------------------------------------------------------------------
// read_locus_file
// ----------------------------------------------------------------------------------------
/** read dominance values and fitness factors from a file
	* not set combinations will be rezessive
	*/
void
TTQuantiProto::read_locus_file(string filename)
{
	#ifdef _DEBUG
	 message("  TTQuantiProto::read_locus_file (%s) ...",filename.c_str());
	#endif
	int i, l, a1, a2;
	bool hasLocus = true;

	// read the file
	TMatrix mat;
	map<string, int> fileInfo = mat.read_matrix(filename);      // read the file
	int nbLines = mat.getNbRows();
	int nbCols  = mat.getNbCols();

	int colsSize=5;
	int cols[5]        = {my_NAN,my_NAN,my_NAN,my_NAN};  // 0: locus; 1: allele1; 2: allele2; 3: value; 4: fitness_factor
	string colsName[5] = {"col_locus","col_allele1","col_allele2","col_dominance","col_fitness_factor"};

	// get col values form file info
	map<string, int>::iterator pos;
	for(i=0; i<colsSize; ++i){
		pos = fileInfo.find(colsName[i]);
		if(pos != fileInfo.end()){
			cols[i] = pos->second - 1;
			fileInfo.erase(pos);
		}
	}

	// ouput all not used columns
	for(pos = fileInfo.begin(); pos != fileInfo.end(); ++pos){
		warning("Dominance file for trait '%s' (%s): Column '%s' (%i) is not a valid column names!\n",
						_type.c_str(), filename.c_str(), pos->first.c_str(), pos->second);
	}

	// check if all cols are available and also dimensions of the matrix are met adn create the arrays
	if(cols[0] == my_NAN) hasLocus = false; // fatal("Dominance file (%s): Locus column is missing!", filename.c_str());
	if(cols[1] == my_NAN) fatal("Dominace file (%s): 1. allele column is missing!", filename.c_str());
	if(cols[2] == my_NAN) fatal("Dominace file (%s): 2. allele column is missing!", filename.c_str());
	if(cols[3] != my_NAN) {  // dominance values
		ARRAY::create_3D(_dominanceValues, _nb_locus, _nb_allele, _nb_allele, (double) my_NAN);
		_dominance_file = true;
	}
	if(cols[4] != my_NAN) {  // fitness factor
		ARRAY::create_3D(_fitnessFactor, _nb_locus, _nb_allele, _nb_allele, (double) my_NAN);
	}

	for (i=0; i<colsSize; ++i){
		if(cols[i]!=my_NAN && cols[i]>nbCols) fatal("Dominance file (%s): File info is wrong: matrix has only %i columns and not %i!\n", filename.c_str(), nbCols, cols[i]);
	}

	// read the table: copy the dominance values
  for(i=0; i<nbLines; ++i){
    if(hasLocus) l  = (int)mat.get(i,cols[0]);    // if locus specific values
		else         l = 1;                           // if all loci have the same settings
    a1 = (int)mat.get(i,cols[1]);    // allele 1
    a2 = (int)mat.get(i,cols[2]);    // allele 2

    // test if out of range
    if(l>_nb_locus){
      fatal("Dominance values: Locus %i is out of range (only %i loci specified)!\n",
            l, _nb_locus);
    }
    if(a1>_nb_allele){
      fatal("Dominance values: Allele %i of locus %i is out of range (only %i alleles specified)!\n",
            a1, l, _nb_allele);
    }
    if(a2>_nb_allele){
      fatal("Dominance values: Allele %i of locus %i is out of range (only %i alleles specified)!\n",
            a2, l, _nb_allele);
    }
		if(a1 > a2){ // wrong order of the alleles: -> swap
      int temp = a1;
      a1 = a2;
      a2 = temp;
    }

    // get the array position (-1 as arrays start at 0)
    --l; --a1; --a2;

		// check if the combination was already read
		if(   (_dominanceValues && _dominanceValues[l][a1][a2] != my_NAN)
			 || (_fitnessFactor   && _fitnessFactor[l][a1][a2]   != my_NAN)){

			fatal("Dominance values: Value was already specified for locus %i and alleles %i and %i!\n", l+1, a1+1, a2+1);
		}

    // set the values
    if(hasLocus){       // if all loci have the same settings
			if(_dominanceValues) _dominanceValues[l][a1][a2] = mat.get(i,cols[3]);
			if(_fitnessFactor)   _fitnessFactor[l][a1][a2]   = mat.get(i,cols[4]);
		}
		else{               // all loci have the same settings
			for(l=0; l<_nb_locus; ++l){
				if(_dominanceValues) _dominanceValues[l][a1][a2] = mat.get(i,cols[3]);
				if(_fitnessFactor)   _fitnessFactor[l][a1][a2]   = mat.get(i,cols[4]);
			}
    }
	}

	#ifdef _DEBUG
	 message("done!\n");
	#endif
}

// ----------------------------------------------------------------------------------------
// read_genome_file
// ----------------------------------------------------------------------------------------
/** read epistatic values from a file */
void
TTQuantiProto::read_genome_file(string filename)
{
	#ifdef _DEBUG
	 message("  TTQuantiProto::read_genome_file (%s) ...",filename.c_str());
	#endif

	// read the file
  TMatrix mat;
  map<string, int> fileInfo;
  int i, l, nbValues=0;

  // read the file
  fileInfo = mat.read_matrix(filename);
  int nbLines = mat.getNbRows();
  int nbCols  = mat.getNbCols();
  vector<map<int, string> >* pStrMatrix = mat.get_strMatrix();

	const int nbParam = 4;
	int colsSize=nbParam;
	int cols[nbParam]        = {my_NAN,my_NAN,my_NAN,my_NAN};  // 0: genotype; 1: epistaticVal; 2: genotypicVal; 3: fitnessFactor
	string colsName[nbParam] = {"col_genotype","col_epistatic_value","col_genotypic_value", "col_fitness_factor"};

  // get col values form file info
  map<string, int>::iterator pos;
	for(i=0; i<colsSize; ++i){
    pos = fileInfo.find(colsName[i]);
		if(pos != fileInfo.end()){
			cols[i] = pos->second - 1;
			fileInfo.erase(pos);        // the element is correct remove it from the list
    }
  }

  // ouput all not used columns
	for(pos = fileInfo.begin(); pos != fileInfo.end(); ++pos){
    warning("Epistatic file for trait '%s' (%s): Column '%s' (%i) will not be considered!\n",
            _type.c_str(), filename.c_str(), pos->first.c_str(), pos->second);
  }

	// check if all cols are available and also dimensions of the matrix are met and create the arrays
  if(cols[0] == my_NAN) fatal("Epistatic file (%s): Genotype column is missing!", filename.c_str());

	if(cols[1] != my_NAN || cols[2] != my_NAN){      // epistatic or genotypic value
		if(_phenoTree)  delete _phenoTree;
		_phenoTree = new Tree<unsigned char>(_nb_locus, _nb_allele);
	}
	else if(_phenoTree){delete _phenoTree; _phenoTree=NULL;}

	if(cols[1] != my_NAN && cols[2] != my_NAN){
		warning("Epistatic file: Both genotypic and epistatic values are present: only the genotypic values are considered!", filename.c_str());
		cols[1] = my_NAN;
	}

	if(cols[2] != my_NAN){
		// allelic and dominance values are not used: remove them
		_allelic_file = _dominance_file = false;
		ARRAY::delete_2D(_allelicValues, _nb_locus);
		ARRAY::delete_3D(_dominanceValues, _nb_locus, _nb_allele);
	}

	if(cols[3] != my_NAN){        // fitness factor
		if(_fitnessFactorTree)  delete _fitnessFactorTree;
		_fitnessFactorTree = new Tree<unsigned char>(_nb_locus, _nb_allele);
	}
	else if(_fitnessFactorTree){delete _fitnessFactorTree; _fitnessFactorTree=NULL;}

	for (i=0; i<colsSize; ++i){
		if(cols[i] != my_NAN && cols[i] > nbCols) fatal("Epistatic file (%s): File info is wrong: matrix has only %i columns and not %i!\n", filename.c_str(), nbCols, cols[i]);
	}


	unsigned int length, digit=-1;       // number of characters representing an allele
	unsigned char** genotype = ARRAY::new_2D<unsigned char>(_nb_locus, _ploidy);
	double value;
  string text, t;
  int locus;

	// get the values
	for(i=0; i<nbLines; ++i){
    // get the genotype
    istringstream LINE;       // make a new stream
		if(mat.get(i,cols[0]) != my_NAN) fatal("Could not read epistatic file '%s': Genotype must be written within brackets!\n", filename.c_str());
		assert((*pStrMatrix)[i].find(cols[0]) != (*pStrMatrix)[i].end());
    text = (*pStrMatrix)[i].find(cols[0])->second;
    text = text.substr(1, text.length()-2); // remove the start and ending bracket
    LINE.str(text);                         // allocate the matrix to the new stream

    // for each genotype
		locus = 0;
    do{
      LINE >> text;
      length = text.length();
			if(length != 2*digit){
				digit = length/2;
				if(length != 2*digit) fatal("Could not read epistatic file '%s': Genotype could not be read at line %i!\n", filename.c_str(), i+1);
      }

      // get first allele of locus i
      t = text.substr(0, digit);
      if(!STRING::is_integer(t)) fatal("File '%s' could not be read at line %i!\n", filename.c_str(), i+1);
      genotype[locus][0] = (unsigned char) STRING::str2int<unsigned int>(t)-1;
      if((unsigned int)genotype[locus][0] >= (unsigned int)_nb_allele){
        fatal("Epistatic values: Allele 1 of locus %i (line %i) is out of range (only %i alleles specified)!\n", locus+1, i+1, _nb_allele);
      }

      // get second allele of locus i
      t = text.substr(digit, digit);
      if(!STRING::is_integer(t)) fatal("File '%s' could not be read at line %i: locus %i is not a number!\n", filename.c_str(), i+1, locus+1);
      genotype[locus][1] = (unsigned char) STRING::str2int<unsigned int>(t)-1;
      if((unsigned int)genotype[locus][1] >= (unsigned int)_nb_allele){
        fatal("Epistatic values: Allele 2 of locus %i (line %i) is out of range (only %i alleles specified)!\n", locus+1, i+1, _nb_allele);
      }
      ++locus;
    }while(!LINE.eof());

    // control if the genotype has the correct number of loci
    if(locus!=_nb_locus) fatal("File '%s' could not be read at line %i: wrong number of loci specified (%i instead of %i!\n", filename.c_str(), i+1, locus+1, _nb_locus);

		// control if the genotype was already set

		// set the genotype
		if(cols[1] != my_NAN){ // set the epistatic value
			if(_phenoTree->get_value(genotype) != my_NAN){
				fatal("Epistatic file %s: value for genotype of line %i was already specified!\n", filename.c_str(), i+1);
			}
			value = mat.get(i,cols[1]);
			if(value == my_NAN) fatal("Epistatic file '%s': epistatic values are not available!\n", filename.c_str());
			_phenoTree->set_value(genotype, get_genotype_none(genotype)+value);
		}
		if(cols[2] != my_NAN){ // set directly the genotypic value
			if(_phenoTree->get_value(genotype) != my_NAN){
				fatal("Epistatic file %s: value for genotype of line %i was already specified!\n", filename.c_str(), i+1);
			}
			value = mat.get(i,cols[2]);
			if(value == my_NAN) fatal("Epistatic file '%s': genotypic values are not available!\n", filename.c_str());
			_phenoTree->set_value(genotype, value);
		}
		if(cols[3] != my_NAN){ // set the fitness factor
			if(_fitnessFactorTree->get_value(genotype) != my_NAN){
				fatal("Epistatic file %s: value for genotype of line %i was already specified!\n", filename.c_str(), i+1);
			}
			value = mat.get(i,cols[3]);
			if(value != my_NAN) _fitnessFactorTree->set_value(genotype, value);
		}
		++nbValues;
	}

	// check if we have enough values
	int totValues = (int)pow((double)(_nb_allele*(_nb_allele+1) / 2), _nb_locus);
	if(nbValues != totValues && (cols[1] != my_NAN || cols[2] != my_NAN)){
		fatal("Epistatic values: %i of %i epistatic values are not set by the file '%s'!\n",
					totValues-nbValues, totValues, filename.c_str());
	}

	ARRAY::delete_2D(genotype, _nb_locus);

	#ifdef _DEBUG
	 message("done!\n");
	#endif
}

// ----------------------------------------------------------------------------------------
// print_allelic_values
// ----------------------------------------------------------------------------------------
/** print only the allelic values which were used */
void
TTQuantiProto::print_allelic_values(string name){

	if(!_allelicValues) return;  // if not used do not print it

	string rpl = _popPtr->getReplicateCounter();
	int p;

	// create the filename
	string filename = name
									+ get_trait_indexStr_t()
                  + _popPtr->getReplicateCounter_r()
                  + ".txt";
  ofstream FILE(filename.c_str());

  FILE.width(12);
  FILE.setf(ios::left,ios::adjustfield);
  FILE << setprecision(4);

  // write title
  FILE << "#Allelic values ";
  if(get_trait_index() || !rpl.empty()){
    FILE << "(";
    if(get_trait_index())                 FILE << "Trait: " << get_trait_indexStr();
    if(get_trait_index() && !rpl.empty()) FILE << ", ";
    if(!rpl.empty())                      FILE << "Replicate: " << rpl;
    FILE << ")";
  }
  FILE << "\n################################################################\n";

	// write file info
	int col = 1;
	FILE << "\n[FILE_INFO]{";
	FILE << "\n    col_locus "         << col++;
	FILE << "\n    col_allele "        << col++;
	FILE << "\n    col_allelic_value " << col++;
	FILE << "\n    col_mut_freq "      << col++;
	FILE << "\n    col_ini_freq "      << col;
  FILE << "\n}";

  FILE << "\n\n# " << _nb_locus  << " loci"
       << "\n# "   << _nb_allele << " alleles per locus" ;


  // write heading
	FILE << "\n\n#locus" << "\tallele" << "\tvalue" << "\tmut_freq" << "\tini_freq";

  // write the values line by line
  for(int l=0; l<_nb_locus; ++l){
    for(int a=0; a<_nb_allele; ++a){
      if(_allelicValues[l][a] != my_NAN){
        FILE << "\n" << (l+1)        // locus
             << "\t" << (a+1)        // allele
						 << "\t" << _allelicValues[l][a]; // allelic value

        if(_mutationFreq){
					if(a) FILE << "\t" << (_mutationFreq[l][a]-_mutationFreq[l][a-1]);  // the array is cumulative!
					else  FILE << "\t" << _mutationFreq[l][a];
				}
				else    FILE << "\t" << "0"; 		// mutation rate is zero

				if(_initAlleleFreq){
          if(_initAlleleFreqCols == 1){
            if(a) FILE << "\t" << (_initAlleleFreq[0][l][a]-_initAlleleFreq[0][l][a-1]);  // the array is cumulative!
            else  FILE << "\t" << _initAlleleFreq[0][l][a];
          }
          else{
            FILE << "\t{";
            for(p=0; p<_initAlleleFreqCols; ++p){
              if(p) FILE << "\t";
              if(a) FILE << (_initAlleleFreq[p][l][a]-_initAlleleFreq[p][l][a-1]);  // the array is cumulative!
              else  FILE << _initAlleleFreq[p][l][a];
            }
            FILE << "}";
          }
				}
				else FILE << "\t" << ((_nb_allele/2)== a ? "1" : "0"); // monomorph inital populations
			}
		}
  }
  FILE.close();
}

// ----------------------------------------------------------------------------------------
// print_dominance_values
// ----------------------------------------------------------------------------------------
/** print only the dominance values which were used */
void
TTQuantiProto::print_dominance_values(string name){

	if(!(_dominanceValues || _fitnessFactor)) return;   // if not used don't make the output

  string rpl = _popPtr->getReplicateCounter();

  // create the filename
  string filename = name
                  + get_trait_indexStr_t()
                  + _popPtr->getReplicateCounter_r()
                  + ".txt";
  ofstream FILE(filename.c_str());

  // write title
  FILE << "#Dominance values";
  if(get_trait_index() || !rpl.empty()){
    FILE << " (";
    if(get_trait_index())                 FILE << "Trait: " << get_trait_indexStr();
    if(get_trait_index() && !rpl.empty()) FILE << ", ";
    if(!rpl.empty())                      FILE << "Replicate: " << rpl;
    FILE << ")";
  }
  FILE << "\n################################################################\n";

	// write file info
	int col = 1;
	FILE << "\n[FILE_INFO]{";
	FILE << "\n    col_locus " << col++;
	FILE << "\n    col_allele1 " << col++;
	FILE << "\n    col_allele2 " << col++;
	if(_dominanceValues) FILE << "\n    col_dominance " << col++;
	if(_fitnessFactor)   FILE << "\n    col_fitness_factor " << col;
	FILE << "\n}";

	FILE << "\n\n# " << _nb_locus  << " loci"
			 << "\n# "   << _nb_allele << " alleles per locus"
			 << "\n# vg(l) = a1 + a2 + k(a2-a1)   with a1 < a2";

	// write heading
	FILE << "\n\n#locus"
			 << "\tall_1"
			 << "\tall_2";
	if(_dominanceValues) FILE << "\tvalue";
	if(_fitnessFactor)   FILE << "\tfitness";

	// write the values line by line
  int l, a1, a2;
	for(l=0; l<_nb_locus; ++l){
		for(a1=0; a1<_nb_allele; ++a1){
			for(a2=a1; a2<_nb_allele; ++a2){
				// are there values present?
				if((_dominanceValues && _dominanceValues[l][a1][a2]!=my_NAN)
					 || (_fitnessFactor && _fitnessFactor[l][a1][a2]!=my_NAN)) continue;

				// write the allele info
				FILE << "\n" << (l+1)                             // locus
						 << "\t" << (a1+1)                            // allele 1
						 << "\t" << (a2+1);                           // allele 2

				// write dominance value if needed
				if(_dominanceValues){
					if(_dominanceValues[l][a1][a2]!=my_NAN) FILE << "\t" << _dominanceValues[l][a1][a2];
					else                                    FILE << "\tNaN";
				}

				// write fitness factor if needed
				if(_fitnessFactor){
					if(_fitnessFactor[l][a1][a2]!=my_NAN)   FILE << "\t" << _fitnessFactor[l][a1][a2];
					else                                    FILE << "\tNAN";
        }
			}
		}
	}

	FILE << "\n";
	FILE.close();
}

// ----------------------------------------------------------------------------------------
// print_epistatic_values
// ----------------------------------------------------------------------------------------
/** print only the epistatic values which were used */
void
TTQuantiProto::print_epistatic_values(string name){

	if(!(_phenoTree || _fitnessFactorTree)) return;   // if not used don't make the output

  string rpl = _popPtr->getReplicateCounter();

	// create the filename
	string filename = name
                  + get_trait_indexStr_t()
                  + _popPtr->getReplicateCounter_r()
                  + ".txt";
  ofstream FILE(filename.c_str());

  int digit = STRING::getNbDigits(_nb_allele);

	// write title
	FILE << "#Epistatic values";
	if(get_trait_index() || !rpl.empty()){
    FILE << " (";
    if(get_trait_index())                 FILE << "Trait: " << get_trait_indexStr();
    if(get_trait_index() && !rpl.empty()) FILE << ", ";
    if(!rpl.empty())                      FILE << "Replicate: " << rpl;
    FILE << ")";
  }
  FILE << "\n################################################################\n";

  // write file info
	int col = 1;
	FILE << "\n[FILE_INFO]{";
	FILE << "\n    col_genotype " << col++;
	if(_phenoTree){
		if(_allelicValues){
			FILE << "\n    col_epistatic_value " << col++;
			FILE << "\n    # col_genotypic_value " << col++;
		}
		else FILE << "\n    col_genotypic_value " << col++;
	}
	if(_fitnessFactorTree)   FILE << "\n    col_fitness_factor " << col;
	FILE << "\n}";

	FILE << "\n\n# " << _nb_locus  << " loci"
			 << "\n# "   << _nb_allele << " alleles per locus";

	// write heading
	FILE << "\n\n#genotype";
	if(_phenoTree){
		if(_allelicValues) FILE << "\tepistaticVal";
		FILE << "\tgenotypicVal";
	}
	if(_fitnessFactorTree) FILE <<"\tfitnessFactor";

	// create the first possible genotype to explore all possible genotypes
	unsigned char** seq = ARRAY::new_2D<unsigned char>(_nb_locus, _ploidy, (unsigned char) 0);

	double val_pheno, val_fitness;
	do{
		// get the values, i.e. check if the genotype is set
		if(_phenoTree)         val_pheno   = _phenoTree->get_value(seq);
		else                   val_pheno   = my_NAN;
		if(_fitnessFactorTree) val_fitness = _fitnessFactorTree->get_value(seq);
		else                   val_fitness = my_NAN;

		// if genotype is not set do not plot it
		if(val_pheno != my_NAN || val_fitness != my_NAN){
			// print the genotype
			FILE << "\n{";
			print_gentoype(FILE, seq, digit);
			FILE << "}";

			// print the values
			if(_phenoTree){
				if(_allelicValues) FILE << "\t" << val_pheno - get_genotype_none(seq);
				FILE << "\t" << val_pheno;
			}
			if(_fitnessFactorTree) FILE << "\t" << val_fitness;
		}
	}while(get_next_gentoype(seq));      // get the next genotype

	FILE.close();

  ARRAY::delete_2D(seq, _nb_locus);
}

//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
void
TTQuantiProto::print_gentoype(ostream& FILE, unsigned char** seq, const int& digit)
{
	int l, a;
	for(l=0; l<_nb_locus; ++l){
		if(l) FILE << " ";
		for(a=0; a<_ploidy; ++a){
			FILE.fill('0');
			FILE.width(digit);
			FILE<<(unsigned int)(seq[l][a]+1);
		}
	}
}

//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
/** generatets the next possible genotype (seq is altered)
	* returns true if this was possible and false if the last possible genotype is reached
	*/
bool
TTQuantiProto::get_next_gentoype(unsigned char** seq)
{
	int l, p, a;

	for(l=_nb_locus-1; l>=0; --l){
			// increment the last allele and check if it is valable
			++seq[l][1];
			if((unsigned int)seq[l][1] < (unsigned int)_nb_allele) return true;

			// increment the second allele and reset the least allele to the same one (note 0102 = 0201!!!)
			++seq[l][0];
			if((unsigned int)seq[l][0] < (unsigned int)_nb_allele){
				seq[l][1] = seq[l][0];
				return true;
			}

			// the locus has to be switched: reset both alleles to zero
			seq[l][0] = seq[l][1] = 0;
	}
	return false;
}
//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTQuanti&
TTQuanti::operator= (const TTrait& T)
{
  pTraitProto->setTTraitProto(*T.pTraitProto);
  const TTQuanti& TN = dynamic_cast<const TTQuanti&> (T);
  if(this != &TN) {
    pProto = TN.pProto;

    reset();
    init();

    int i, nbLocus  = pProto->_nb_locus;
    int j, nbPloidy = pProto->_ploidy;
    for(i = 0; i < nbLocus; ++i) {
      for(j = 0; j < nbPloidy; ++j) {
        sequence[i][j] = TN.sequence[i][j];
      }
    }
  }

  return *this;
}

//----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool
TTQuanti::operator== (const TTrait& T)
{
  if(*pTraitProto != *T.pTraitProto) return false;
  const TTQuanti& TN = dynamic_cast<const TTQuanti&> (T);
  if(this != &TN || *pProto != *TN.pProto) return false;
  return true;
}

//----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool
TTQuanti::operator!= (const TTrait& T)
{
  if(!((*this) == T)) return true;
  else                return false;
}

// -----------------------------------------------------------------------------
// mutate
// -----------------------------------------------------------------------------
void
TTQuanti::mutate(){
  (pProto->*(pProto->_mutate_func_ptr))(sequence);
}

// ----------------------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------------------
TTQuanti::~TTQuanti()
{
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void
TTQuanti::init ()
{
  if(sequence) fatal("TTQuanti::init::sequence is not NULL !\n");

  int size = pProto->_nb_locus;
  sequence = new unsigned char* [size];
  for(int i=0; i < size; i++){
    sequence[i] = new unsigned char[2];
  }
}

// ----------------------------------------------------------------------------------------
// ini_sequence
// ----------------------------------------------------------------------------------------
/** initialization of the sequence */
void
TTQuanti::ini_sequence (Patch* patch)
{
	(pProto->*(pProto->_ini_allele_func_ptr))(sequence, patch);
} // ini_sequence

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void
TTQuanti::reset()
{
	_phenotype      = my_NAN;
	_genotype       = my_NAN;
	_fitness_factor = my_NAN;

}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** Prototype ********

// ----------------------------------------------------------------------------------------
// get_locus_genotype
// ----------------------------------------------------------------------------------------
double
TTQuantiProto::get_locus_genotype_additive(const int& l, const unsigned char& a1, const unsigned char& a2){
	return (getAllelicValue(l, a1) + getAllelicValue(l, a2));
}

 /* G = 2a1 + (1+k)(a2-a1) = a1 + a2 + k(a2-a1)       // a1 < a2
  * k < -1:  underdominant
  * k = -1:  the smaller allele (a1) is dominant
  * k = 0:   purely additive
  * k = 1:   the bigger allele (a2) is dominant
  * k > 1:   overdominance
  */
double
TTQuantiProto::get_locus_genotype_dominance_array(const int& l, const unsigned char& a1, const unsigned char& a2){
	double a_1 = getAllelicValue(l, a1);    // get effect of allele 1
	double a_2 = getAllelicValue(l, a2);    // get effect of allele 2
	if(a_1<a_2) return a_1 + a_2 + getDominanceValue(l, a1, a2)*(a_2 - a_1);
	if(a_1>a_2) return a_2 + a_1 + getDominanceValue(l, a2, a1)*(a_1 - a_2);
	return a_1 + a_2;                       // homozygote: 2a1
}

double
TTQuantiProto::get_locus_genotype_dominance_single(const int& l, const unsigned char& a1, const unsigned char& a2){
  double a_1 = getAllelicValue(l, a1);    // get effect of allele 1
  double a_2 = getAllelicValue(l, a2);    // get effect of allele 2
  if(a_1<a_2) return a_1 + a_2 + _dominance_mean*(a_2 - a_1);
  if(a_1>a_2) return a_2 + a_1 + _dominance_mean*(a_1 - a_2);
  return a_1 + a_2;                       // homozygote: 2a1
}

// ----------------------------------------------------------------------------------------
// get_genotype
// ----------------------------------------------------------------------------------------
double
TTQuantiProto::get_genotype_none(unsigned char** seq){
	double sum=0;
	for(int l=0; l<_nb_locus; ++l){
		sum += (this->*get_locus_genotype_func_ptr)(l, seq[l][0], seq[l][1]);
	}
	return sum;
}

double
TTQuantiProto::get_genotype_epistatic(unsigned char** seq){
	double val  = _phenoTree->get_value(seq);
	if(val != my_NAN) return val;

	// if not yet computed compute it
	val = get_genotype_none(seq);
	val += SimRunner::r.Normal(0, _epistatic_sd);
	_phenoTree->set_value(seq, val);
	return val;
}

// ----------------------------------------------------------------------------------------
// get_selectionCoef_explcit
// ----------------------------------------------------------------------------------------
double
TTQuantiProto::get_fitnessFactor_explicit(unsigned char** seq){
	assert(_fitnessFactor);

	double product=1, value;
	unsigned char a1, a2;
	for(int l=0; l<_nb_locus; ++l){
		a1 = seq[l][0];     // get allele 1
		a2 = seq[l][1];     // get allele 2
		if(a1 > a2) value = _fitnessFactor[l][a2][a1];
		else        value = _fitnessFactor[l][a1][a2];

		if(value != my_NAN){ // value is explicitly defined
			product *= value;
		}
		else{                // generally defined
			if(a1==a2) product *= _fitnessFactor_homozygote   ? _fitnessFactor_homozygote[l]   : 1;
			else       product *= _fitnessFactor_heterozygote ? _fitnessFactor_heterozygote[l] : 1;
		}
	}
	return product;
}

// ----------------------------------------------------------------------------------------
// get_fitnessFactor_global
// ----------------------------------------------------------------------------------------
double
TTQuantiProto::get_fitnessFactor_global(unsigned char** seq){
	double product=1;
	for(int l=0; l<_nb_locus; ++l){
		if(seq[l][0]==seq[l][1]) product *= _fitnessFactor_homozygote   ? _fitnessFactor_homozygote[l]   : 1;
		else                     product *= _fitnessFactor_heterozygote ? _fitnessFactor_heterozygote[l] : 1;
	}
	return product;
}

// ----------------------------------------------------------------------------------------
// get_fitnessFactor_genome
// ----------------------------------------------------------------------------------------
double
TTQuantiProto::get_fitnessFactor_genome(unsigned char** seq){
	double val  = _fitnessFactorTree->get_value(seq);
	if(val != my_NAN) return val;

	// if not set for the genome get the explicit fitness factor
	if(_fitnessFactor) return get_fitnessFactor_explicit(seq);
	return get_fitnessFactor_global(seq);
}

// ----------------------------------------------------------------------------------------
// fitnessFactor_used
// ----------------------------------------------------------------------------------------
/** is a selection coefficent set? */
bool
TTQuantiProto::fitnessFactor_used(){
	return (get_fitnessFactor_func_ptr != NULL);
}


// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void
TTQuanti::show_up ()
{
  message("\n  Trait's type: discretequanti\n\
       locus: %i\n\
     alleles: %i\n\
    sequence:",pProto->_nb_locus,pProto->_nb_allele);

  for(int i = 0; (i < pProto->_nb_locus && i < 10); i++){
    message("\n              %i %i",(int)sequence[i][0],(int)sequence[i][1]);
  }
  message("\n");
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                 ******** Phenotyper / Genotyper ********

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void
TTQuantiFHvalue::FHwrite ()
{
  int patchNbr = get_pop_ptr()->getPatchNbr();
  int t;
  Patch* current_patch;

  // get the number of occupied patches
  int nb_patches_alive = 0;
  for (int i=0; i<patchNbr; ++i) {
		if(!get_pop_ptr()->getPatch(i)->get_isExtinct()) ++nb_patches_alive;
	}
	if(!nb_patches_alive) return;

	string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + get_extension();

  #ifdef _DEBUG
   message("TTQuantiFHvalue::FHwrite (%s)\n",filename.c_str());
  #endif

	ofstream FILE (filename.c_str(), ios::out);

  if(!FILE) fatal("could not open Phenotypes output file!!\n");

  // write the heading line
  FILE << nb_patches_alive << " " << _nb_trait << "\n";

  // write the names of the traits
  for(t=0; t<_nb_trait; ++t){
    FILE << _name << "_trait-" << (t+1) << "\n";
  }

  //FILE.width(12);
  FILE.setf(ios::left,ios::adjustfield);
  FILE.precision(4);

  unsigned int a, mask;
  for (int i=0; i<patchNbr; ++i) {
    current_patch = get_pop_ptr()->getPatch(i);
		for(a = 0, mask=1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
			if(mask & _age){
        if(_sex != 2) FHwrite(static_cast<age_idx>(a), FEM, FILE, current_patch, i);
        if(_sex != 1) FHwrite(static_cast<age_idx>(a), MAL, FILE, current_patch, i);
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void TTQuantiFHvalue::FHwrite (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
      Patch* current_patch, const int& patch_id)
{
  Individual* ind;
  unsigned int nbInd = current_patch->size(cur_sex, cur_age);
  if(!nbInd) return;

  double value;

  for (unsigned int j = 0; j < nbInd; ++j){
	  ind = current_patch->get(cur_sex,cur_age, j);
    FILE << (patch_id+1);
    for (int t=0; t<_nb_trait; ++t){
      value = (this->*get_value_func_ptr)(ind, t);
      if(value == my_NAN) FILE << "\t" << "NaN";
			else                FILE << "\t" << value;
    }
	  if(get_fstat_choice()==2){
			FILE << "\t";
			write_individual_info_to_stream(FILE, ind, cur_age, cur_sex, '\t');
		}
    FILE << "\n";
  }
}


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** StatHandler ********

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool TTQuantiSH::init()
{
  StatHandler<TTQuantiSH>::init();
  allocateTable();
  return true;
}

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTQuantiSH::setStatRecorders(const string& t)
{
  // group quanti
			 if(t=="q.varA_p")  add_perPatch("Additive genetic variance","q.varA",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getVarA);
	else if(t=="q.meanG_p") add_perPatch("Genetic mean","q.meanG",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getMeanG);
	else if(t=="q.varG_p")  add_perPatch("Genetic variance","q.varG",FLAT,ADULTS,0,0,0,0,0,0,&TTQuantiSH::getVarG);
	else if(t=="q.meanP_p") add_perPatch("Phenotypic mean","q.meanP",FLAT,ADULTS,0,0,0,&TTQuantiSH::getMeanP);
	else if(t=="q.varP_p")  add_perPatch("Phenotypic variance","q.varP",FLAT,ADULTS,0,0,0,&TTQuantiSH::getVarP);

  else if(t=="q.VaW") add("Additive genetic variance within patches","q.VaW",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVaW);
  else if(t=="q.VgB") add("Genetic variance between patches","q.VgB",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVgB);
  else if(t=="q.VgW") add("Genetic variance within patches","q.VgW",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVgW);
  else if(t=="q.VpB") add("Phenotypic variance between patches","q.VpB",FLAT,ADULTS,0,&TTQuantiSH::getVpB);
	else if(t=="q.VpW") add("Phenotypic variance within patches","q.VpW",FLAT,ADULTS,0,&TTQuantiSH::getVpW);
	else if(t=="q.qst") add("Qst","q.qst",FLAT,ADULTS,0,&TTQuantiSH::getQst);
	else if(t=="q.qst_pair") add_pairwisePatch("Qst","q.qst_pair",FLAT,ADULTS,0,0,0,&TTQuantiSH::getQst_ij);

	else if(t=="quanti"){
		add("Genetic variance between patches","q.VgB",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVgB);
   	add("Genetic variance within patches","q.VgW",FLAT,ADULTS,0,0,0,0,&TTQuantiSH::getVgW);
   	add("Phenotypic variance between patches","q.VpB",FLAT,ADULTS,0,&TTQuantiSH::getVpB);
   	add("Phenotypic variance within patches","q.VpW",FLAT,ADULTS,0,&TTQuantiSH::getVpW);
	}

  else if(set_stat_coancestry     (t, "q", "quanti")){}     // adults and offspring
  else if(set_stat_fstat          (t, "q", "quanti")){}     // adults and offspring
  else if(set_stat_all_freq_local (t, "q", "quanti")){}     // adults and offspring
  else if(set_stat_all_freq_global(t, "q", "quanti")){}     // adults and offspring
  else return false;

	return true;
}

// ----------------------------------------------------------------------------------------
// setVarAdditiveGenetic
// ----------------------------------------------------------------------------------------
/** computation of the additive genetic variance (for any case) for each patch */
void
TTQuantiSH::setVar_Va(const age_t& AGE){
  // check if the table has already been computed
  if(already_computed(_computed[20], AGE)) return;

  if(!_varA)   _varA   = new double[_sample_pops_size];
  set_alleleCountTable_local(AGE);        // compute allele frequencies
  setMeanAndVar_Vg(AGE);                  // compute genetic mean and variance
  double mean; // not really used...

  // for each population
  for(unsigned int i = 0; i < _sample_pops_size; ++i) {
    (this->*get_Va_ofPatch_func_ptr)(_sample_pops[i], AGE, mean, _varA[i], _alleleFreqTable[i],_alleleCountTable[i]);
  }
}

// ----------------------------------------------------------------------------------------
// setVarAdditiveGenetic
// ----------------------------------------------------------------------------------------
/** computation of the additive genetic variance (only for random matings)
  * much quicker than the version for non-random matings
  * meanA is nto computed. But needed for function pointer
  */
void
TTQuantiSH::get_Va_ofPatch_random_mating(Patch* crnt_patch, const age_t& AGE,
                                      double& meanA, double& varA, double** freqs, unsigned int** counts)
{
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  unsigned int sizeF = crnt_patch->size(FEM, age_pos),
               sizeM = crnt_patch->size(MAL, age_pos),
               size  = sizeF + sizeM;
  if(!size) {varA = my_NAN; return;}

  double G, meanG=0;
  double *curElem;
  TTree<double>       condMeanGG(_nb_locus, _nb_allele, 0.0);
  TTree<unsigned int> condMeanGGsize(_nb_locus, _nb_allele, 0);
  unsigned int i, l, a, a1, a2;
  double* aGeno = new double[size];
  unsigned char** genes;
  Individual* ind;

  // females
  for(i = 0; i < sizeF; ++i) {
    ind   = crnt_patch->get(FEM, age_pos, i);
    meanG += G = ind->getTraitGenotype(_SHLinkedTraitIndex);          // genotype
    genes = (unsigned char**)ind->getTrait(_SHLinkedTraitIndex)->get_sequence();  // sequence
    for (l = 0; l < _nb_locus; ++l){
      a1 = genes[l][0];          // 1. allele
      a2 = genes[l][1];          // 2. allele

      if(a1<a2){
        condMeanGG.get(l,a1,a2) += G;
        condMeanGGsize.get(l,a1,a2) += 1;
      }
      else{
        condMeanGG.get(l,a2,a1) += G;
        condMeanGGsize.get(l,a2,a1) += 1;
      }
    }
  } // for each female

  // males
  for(i = 0; i < sizeM; ++i) {
    ind   = crnt_patch->get(MAL, age_pos, i);
    meanG += G = ind->getTraitGenotype(_SHLinkedTraitIndex);          // genotype
    genes = (unsigned char**)ind->getTrait(_SHLinkedTraitIndex)->get_sequence();  // sequence
    for (l = 0; l < _nb_locus; ++l){
      a1 = genes[l][0];          // 1. allele
      a2 = genes[l][1];          // 2. allele

      if(a1<a2){
        condMeanGG.get(l,a1,a2) += G;
        condMeanGGsize.get(l,a1,a2) += 1;
      }
      else{
        condMeanGG.get(l,a2,a1) += G;
        condMeanGGsize.get(l,a2,a1) += 1;
      }
    }
  } // for each male

  // for each present genotype
  meanG /= size;
  double** alphaStar = ARRAY::new_2D(_nb_locus, _nb_allele, -meanG);
  for(curElem = condMeanGG.first(l,a1,a2); curElem; curElem = condMeanGG.next(l,a1,a2)){
    // correct the dominance effects for the sample number
    (*curElem) /= condMeanGGsize.get(l,a1,a2);

    // compute the average excess (alphaStar) which is in this case (random mating) identical to the additive effects
    alphaStar[l][a1]   += (*curElem) * freqs[l][a2];
    if(a1 != a2){                 // if the two alles are not identical ...
      alphaStar[l][a2] += (*curElem) * freqs[l][a1];
    }
  }

  // compute the additive variance
  varA = 0;
  for(l=0; l<_nb_locus; ++l){
    for(a = 0; a < _nb_allele; ++a){    // for each allele
      if(freqs[l][a]){
        varA += alphaStar[l][a]*alphaStar[l][a]*freqs[l][a];        // only valable when random mating
      }
    }
  }
  varA *= 2;    // we have diploid individuals ...

  // delete arrays
  ARRAY::delete_2D(alphaStar, _nb_locus);
  delete[] aGeno;
}

// ----------------------------------------------------------------------------------------
// setVarAdditiveGenetic
// ----------------------------------------------------------------------------------------
/** computation of the additive genetic variance (for any case)  */
void
TTQuantiSH::get_Va_ofPatch_regression(Patch* crnt_patch, const age_t& AGE,
                                      double& meanA, double& varA, double** freqs, unsigned int** counts)
{
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  unsigned int sizeF = crnt_patch->size(FEM, age_pos),
               sizeM = crnt_patch->size(MAL, age_pos),
               size  = sizeF + sizeM;
  if(!size) {varA = my_NAN; return;}

  unsigned int i, m, l;
  double *curElem;
  double* arrayG = new double[size];
  unsigned char a1, a2;
  unsigned char** genes;
  unsigned int a, b;
  double G, meanG=0;

  TTree<double>       condMeanGG(_nb_locus, _nb_allele, 0.0);
  TTree<unsigned int> condMeanGGsize(_nb_locus, _nb_allele, 0);

  // make matrix of predictor values (each column for an allele, each row for an individual)
  int*** predMatrix = ARRAY::new_3D(_nb_locus, _nb_allele, size, (int) 0); // array[locus][allele][individual]

  Individual* ind;

  // females
  for(i = 0; i < sizeF; ++i) {
    ind   = crnt_patch->get(FEM, age_pos, i);
    meanG += G = ind->getTraitGenotype(_SHLinkedTraitIndex);                               // genotype
    genes = (unsigned char**)ind->getTrait(_SHLinkedTraitIndex)->get_sequence();  // sequence
    for (l = 0; l < _nb_locus; ++l){
      a1 = genes[l][0];          // 1. allele
      a2 = genes[l][1];          // 2. allele

      if(a1<a2){
        condMeanGG.get(l,a1,a2) += G;
        condMeanGGsize.get(l,a1,a2) += 1;
      }
      else{
        condMeanGG.get(l,a2,a1) += G;
        condMeanGGsize.get(l,a2,a1) += 1;
      }

      // create the predictor matrix
      ++predMatrix[l][a1][i];        // [locus][allele][individual]
      ++predMatrix[l][a2][i];
    }
  } // for each female

  // males
  for(m = 0; m < sizeM; ++i, ++m) {
    ind   = crnt_patch->get(MAL, age_pos, m);
    meanG += G = ind->getTraitGenotype(_SHLinkedTraitIndex);                               // genotype
    genes = (unsigned char**)ind->getTrait(_SHLinkedTraitIndex)->get_sequence();  // sequence
    for (l = 0; l < _nb_locus; ++l){
      a1 = genes[l][0];          // 1. allele
      a2 = genes[l][1];          // 2. allele

      if(a1<a2){
        condMeanGG.get(l,a1,a2) += G;
        condMeanGGsize.get(l,a1,a2) += 1;
      }
      else{
        condMeanGG.get(l,a2,a1) += G;
        condMeanGGsize.get(l,a2,a1) += 1;
      }

      // create the predictor matrix
      ++predMatrix[l][a1][i];    // [locus][allele][individual]
      ++predMatrix[l][a2][i];
    }
  } // for each male

  // correct the genotype by the mean genotype
  meanG /= size;
  for(i = 0; i < size; ++i) {
    arrayG[i] -= meanG;
  }

  // compute allele frequencies from the number of alleles
  vector<int>* availableAllele = new vector<int>[_nb_locus];
  for(l=0; l<_nb_locus; ++l){
    for(a = 0; a < _nb_allele; ++a){         // compute number of alleles
      if(freqs[l][a]){
        availableAllele[l].push_back(a);    // get the present alleles
      }
    }
  }

  // compute the additive effects (alpha, time consuming!)
  double** alpha = new double*[_nb_locus];
  for(l=0; l<_nb_locus; ++l){
    alpha[l] = new double[availableAllele[l].size()];   // only the present alleles have to be computed!
    if(!compute_alpha(arrayG, predMatrix[l], size, alpha[l], availableAllele[l])){
      remove_private_alleles_compute_alpha(crnt_patch,sizeF,sizeM, alpha[l],availableAllele[l],arrayG,age_pos,counts,l);
    }
  }

  // compute the average excess (alphaStar)
  double** alphaStar = ARRAY::new_2D(_nb_locus, _nb_allele, -meanG);
  for(curElem = condMeanGG.first(l,a,b); curElem; curElem = condMeanGG.next(l,a,b)){   // for each genotype
    // correct the dominance effects for the sample number
    (*curElem) /= condMeanGGsize.get(l,a,b);

    // compute the average excess (alphaStar)
    alphaStar[l][a]     += (*curElem) * freqs[l][b];
    if(a != b){                 // if the two alles are not identical ...
      alphaStar[l][b]     += (*curElem) * freqs[l][a];
    }
  }

  // compute the additive variance  var = 2*sum(p a a*)
  varA = 0;
  unsigned int a_size, nbFailure = 0;
  for(l=0; l<_nb_locus; ++l){
    if(alpha[l][0] == my_NAN) ++nbFailure;                    // sum up to failures in computing alpha
    for(a = 0, a_size = availableAllele[l].size(); a < a_size; ++a){    // for each allele
      a1 = availableAllele[l][a];
      varA += alpha[l][a]*alphaStar[l][a1]*freqs[l][a1];          // any case
      // var += alphaStar[l][a]*alphaStar[l][a]*freqs[l][a1];     // limited to random mating
    }
  }
  varA *= 2;    // we have diploid individuals ...

  // if more than 10% of the loci were not able to be computed set Va to NAN
  if(nbFailure > 0.1*_nb_locus) varA = my_NAN;

  // clean up
  delete[] availableAllele;
  delete[] arrayG;
  ARRAY::delete_2D(alpha, _nb_locus);
  ARRAY::delete_2D(alphaStar, _nb_locus);
  ARRAY::delete_3D(predMatrix, _nb_locus, _nb_allele);
}

// ----------------------------------------------------------------------------------------
// compute_alpha
// ----------------------------------------------------------------------------------------
   // traditional sum of squares and products method of calculation
   // compute the alpha
   // returns false if the regression was not able to compute due to
   // individuals which have two private alleles
bool
TTQuantiSH::compute_alpha(double* y, int** x, int nb_ind,
                    double* alpha, const vector<int>& availableAllele){

  int nbAllele = availableAllele.size();  // number of alleles present in the population

  // test if the population is fixed for one allele
  if(nbAllele <=1){
    for (int i = 0; i<nbAllele; ++i){
      alpha[i] = 0;  // if singular matrix
    }
    if(nbAllele) return true;
    else         return false;    // the computation was not successful
  }

  // test if the population is fixed for two alleles
  if(nbAllele ==2){
    int i;
    for (i = 0; i<nb_ind; ++i){
      if(x[availableAllele[0]][i] != 1) break;   // if not all individuals have both alleles
    }
    if(i==nb_ind){      // if singular matrix
      for (i = 0; i<nbAllele; ++i){
        alpha[i] = 0;
      }
      return true;
    }
  }

  // create the matrixes
  Matrix X(nb_ind, nbAllele);          // make matrix of predictor values
  ColumnVector Y(nb_ind);
  try{
    // load the predictive matrix
    for(int i=0; i<nbAllele; ++i){
      X.Column(i+1) << x[availableAllele[i]];
    }

    // load Y values
    Y << y;

    // calculate estimate
    ColumnVector A = (X.t() * X).i() * (X.t() * Y);  // .i(): inverse; .t(): transpose

    for (int i = 0; i<nbAllele; ++i){
      alpha[i] = A(i+1);  // all the other values
    }
  }catch(...){
    // warning("%s!\n", BaseException::what());
    // cout << "\nmatrix y - x:\n\n" << (Y|X) << endl;

    // the problem is that "X.t() * X" can result in a singular matrix
    // thus there is no solution to the regression
    for (int i = 0; i<nbAllele; ++i){
      alpha[i] = 0;  // no additive effect at this locus
    }
    return false;
  }
  return true;
}

// ----------------------------------------------------------------------------------------
// remove_private_alleles_compute_alpha
// ----------------------------------------------------------------------------------------
/** this function is called when the regression cannot be resolved
  * (i.e. when t(X)*X results in a singular matrix det(t(X)*X) = 0)
  * this function remnoves all individuals havin two private alleles, i.e.
  * alleles which are unique to this individual.
  * Then the addtitive effects (alpha) are recomputed.
  * If more than 10 percent of the individuals have to be removed or the matrix is
  * again singular the computation is assumed to be unsuccessful
  */
void
TTQuantiSH::remove_private_alleles_compute_alpha(Patch* crnt_patch, const unsigned int& sizeF,
                                                 const unsigned int& sizeM, double* alpha,
                                                 vector<int>& availableAllele, double* arrayG,
                                                 const age_idx& age_pos, unsigned int** allele_counts,
                                                 const unsigned int& l)
{
  // the regression was not possible to compute: remove all individuals with two private alleles
  unsigned int i;
  unsigned char* g;
  vector<unsigned char*> geno;    // vector with all "good" locus genotypes
  vector<double>   vectorG1;      // vector of the corrected genotypic values
  vector<int> rem_all;            // vector of alleles to remove
  bool stop = false;

  // check each individual if it has two private alleles
  for(i = 0; i<sizeF; ++i){       // for each female
    g = (unsigned char*)crnt_patch->get(FEM, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence()[l];
    if(allele_counts[l][g[0]]==1 && allele_counts[l][g[1]]==1){ // if both alleles are private          {
      rem_all.push_back(g[0]);
      rem_all.push_back(g[1]);
    }
    else {
      geno.push_back(g);              // add the locus genotype
      vectorG1.push_back(arrayG[i]);
    }
  }
  for(i = 0; i<sizeM; ++i){       // for each male
    g = (unsigned char*)crnt_patch->get(MAL, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence()[l];
    if(allele_counts[l][g[0]]==1 && allele_counts[l][g[1]]==1){ // if both alleles are private          {
      rem_all.push_back(g[0]);
      rem_all.push_back(g[1]);
    }
    else{
      geno.push_back(g);              // add the locus genotype
      vectorG1.push_back(arrayG[sizeF+i]);
    }
  }

  // if more then 10% of the individuls shas to be removed consider it as not computable
  if(geno.size() > 0.1*(sizeF+sizeM)) stop = true;
  else{
    // create the new available allele vector
    vector<int> availableAllele1;
    unsigned int r=0;
    unsigned int s=availableAllele.size();
    sort(rem_all.begin(), rem_all.end()); // sort the vector
    for(i=0; i<s; ++i){
      if(r < rem_all.size() && availableAllele[i] == rem_all[r]) ++r;
      else availableAllele1.push_back(availableAllele[i]);
    }

    // create the new predMatrix
    unsigned int size1 = geno.size();
    double* arrayG1 = new double[size1];
    int** predMatrix1 = ARRAY::new_2D(_nb_allele, size1, (int) 0);        // array[allele][individual]

    for(i=0; i<size1; ++i){                   // for each individual
      ++predMatrix1[geno[i][0]][i];           // allele 1
      ++predMatrix1[geno[i][1]][i];           // allele 2
      arrayG1[i] = vectorG1[i];
    }

    // compute the alpha again
    if(!compute_alpha(arrayG1, predMatrix1, size1, alpha, availableAllele1)) stop = true;
  }

  if(stop){// if the least-square regression had again no solution, set the first additive effect to NAN
    alpha[0]=my_NAN;

    // generate a warning the first 10 times the problem occures and then every 100 time
    static unsigned int passes;      // is initialised with 0
    static unsigned int replicate;   // needed to reset the counter when a new replicate is computed
    if(replicate != _current_replicate){ passes=0; replicate = _current_replicate;}
    ++passes;
    if(passes < 10 || !(passes%100)){
      warning("Va could not be correctly estimated (%i. time at generation %i, see manual parameter 'quanti_va_model')!\n", passes, _current_generation);
    }
  }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarGenotyp
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vg(const age_t& AGE){
  // check if the table has already been computed
  if(already_computed(_computed[21], AGE)) return;

  if(!_meanG) _meanG = new double[_sample_pops_size];
  if(!_varG)  _varG  = new double[_sample_pops_size];

  // for each population
  for(unsigned int i = 0; i < _sample_pops_size; ++i) {
    setMeanAndVar_Vg_ofPatch(_sample_pops[i], AGE, _meanG[i], _varG[i]);
  }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarGenotyp
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vg_ofPatch(Patch* crnt_patch, const age_t& AGE,
                                    double& meanG, double& varG, double** freqs, unsigned int** counts){
  // create a temporar array with all individuals
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  unsigned int sizeF = crnt_patch->size(FEM, age_pos),
               sizeM = crnt_patch->size(MAL, age_pos),
               size  = sizeF + sizeM;
  if(!size){
    meanG = varG = my_NAN;
    return;
  }

  double* array = new double[size];

  unsigned int f, m;
  for(f = 0; f < sizeF; ++f) {
    array[f] = crnt_patch->get(FEM, age_pos, f)->getTraitGenotype(_SHLinkedTraitIndex);
  }
  for(m = 0; m < sizeM; ++m, ++f) {
    array[f] = crnt_patch->get(MAL, age_pos, m)->getTraitGenotype(_SHLinkedTraitIndex);
  }

  // compute mean and var
  meanG = ARRAY::mean(array, size);
  varG  = ARRAY::var(array, size, meanG);
  delete[] array;
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarPhenotype
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vp(){
  // check if the table has already been computed
  if(already_computed(_computed[22])) return;

  if(!_meanP) _meanP = new double[_sample_pops_size];
  if(!_varP)  _varP  = new double[_sample_pops_size];

  // for each population
  for(unsigned int i = 0; i < _sample_pops_size; ++i) {
    setMeanAndVar_Vp_ofPatch(_sample_pops[i], _meanP[i], _varP[i]);
  }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarPhenotype
// ----------------------------------------------------------------------------------------
void
TTQuantiSH::setMeanAndVar_Vp_ofPatch(Patch* crnt_patch, double& meanP, double& varP){
  unsigned int sizeF = crnt_patch->size(FEM, ADLTx),
               sizeM = crnt_patch->size(MAL, ADLTx),
               size  = sizeF + sizeM;

  // if the patch is empty or the phneotype is not yet computed -> stop
  if(  !((sizeF && crnt_patch->get(FEM, ADLTx, 0)->getTraitPhenotype(_SHLinkedTraitIndex) != my_NAN)
    || (sizeM && crnt_patch->get(MAL, ADLTx, 0)->getTraitPhenotype(_SHLinkedTraitIndex) != my_NAN))){
    meanP = varP = my_NAN;
    return;
  }

  double* array = new double[size];

  unsigned int f, m;
  for(f = 0; f < sizeF; ++f) {
    array[f] = crnt_patch->get(FEM, ADLTx, f)->getTraitPhenotype(_SHLinkedTraitIndex);
  }
  for(m = 0; m < sizeM; ++f, ++m) {
    array[f] = crnt_patch->get(MAL, ADLTx, m)->getTraitPhenotype(_SHLinkedTraitIndex);
  }

  // compute mean and var
  meanP = ARRAY::mean(array, size);
  varP  = ARRAY::var(array, size, meanP);
  delete[] array;
}

// ----------------------------------------------------------------------------------------
// getQst
// ----------------------------------------------------------------------------------------
/** calculation of QST:
  * QST = Vb/(Vb+2*h2*Vp) = Vb/(Vb+2*Va)
  * Vb  = the variance of the phenotypic population means
  * Vp  = the mean of the phenotypic population variance
  * Va  = the mean of the additive genotypic population variance
  */
double
TTQuantiSH::getQst(){
  if(get_current_nbPops(ADLTx)<2) return my_NAN;        // at least two populations are needed

  setVar_Va(ADULTS);
  setMeanAndVar_Vp();

  double Vb = ARRAY::var(_meanP, _sample_pops_size);    // round of problems
  if(Vb == my_NAN) return my_NAN;

  double Va = ARRAY::mean(_varA, _sample_pops_size);
  if(Va == my_NAN) return my_NAN;

  return  (2*Va+Vb) ? Vb/(2*Va+Vb) : my_NAN;
}

// ----------------------------------------------------------------------------------------
// setQst_perPatchPair
// ----------------------------------------------------------------------------------------
/** calculation of QST for all pairwise combinations:
  * QST = Vb/(Vb+2*h2*Vp) = Vb/(Vb+2*Va)
	* Vb  = the variance of the phenotypic population means
  * Vp  = the mean of the phenotypic population variance
  * Va  = the mean of the additive genotypic population variance
	*/
void
TTQuantiSH::setQst_perPatchPair()
{
	if(!_qst_matrix) ARRAY::create_2D(_qst_matrix, _sample_pops_size, _sample_pops_size);

	setVar_Va(ADULTS);
	setMeanAndVar_Vp();

	unsigned int i, j;
	double Vb, Va;
	double array[2];		// temporare array


  for(i=0; i<_sample_pops_size-1; ++i){
		if(_meanP[i] == my_NAN){                   // both pops have to be populated
			for(j=i+1; j<_sample_pops_size; ++j){
				_qst_matrix[i][j] = my_NAN;
			}
			continue;
		}
		array[0] = _meanP[i];

		for(j=i+1; j<_sample_pops_size; ++j){
			if(_meanP[j] == my_NAN){                   // both pops have to be populated
				_qst_matrix[i][j] = my_NAN;
				continue;
			}
			array[1] = _meanP[j];
			Vb = ARRAY::var(array, 2);         		// var of means
			Va = (_varA[i] +_varA[j])/2.0;    // mean of vars

			_qst_matrix[i][j] = (2*Va+Vb) ? Vb/(2*Va+Vb) : my_NAN;
		}
	}
}

// ----------------------------------------------------------------------------------------
// getVgB
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVgB(const age_t& AGE){
  if(get_current_nbPops(AGE)<2) return my_NAN;        // at least two populations are needed

  setMeanAndVar_Vg(AGE);

  return ARRAY::var(_meanG, _sample_pops_size);
}

// ----------------------------------------------------------------------------------------
// getVpB
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVpB(){
  if(get_current_nbPops(ADLTx)<2) return my_NAN;        // at least two populations are needed

  setMeanAndVar_Vp();

  return ARRAY::var(_meanP, _sample_pops_size);
}
// ----------------------------------------------------------------------------------------
// getVaW
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVaW(const age_t& AGE){
  setVar_Va(AGE);

  return ARRAY::mean(_varA, _sample_pops_size);
}

// ----------------------------------------------------------------------------------------
// getVgW
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVgW(const age_t& AGE){
  setMeanAndVar_Vg(AGE);

  return ARRAY::mean(_varG, _sample_pops_size);
}

// ----------------------------------------------------------------------------------------
// getVpW
// ----------------------------------------------------------------------------------------
double
TTQuantiSH::getVpW(){
  setMeanAndVar_Vp();

  return ARRAY::mean(_varP, _sample_pops_size);
}

