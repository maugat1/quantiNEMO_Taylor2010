/** @file ttneutral.cpp
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
#include "ttneutral.h"
#include "random.h"
#include "stathandler.cpp"

string TTNeutralProto::_ini_fstat_file;

TTNeutralFH*         TTNeutralProto::_writer;


// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTNeutralProto::TTNeutralProto( ) : _stats(0)
{
  _type = "ntrl";
  ini_paramset();
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTNeutralProto::TTNeutralProto(int i) : _stats(0)
{
	_trait_index = i;
  _type = "ntrl"+get_trait_indexStr_();
  ini_paramset();
}

// ----------------------------------------------------------------------------------------
// cstor
// ----------------------------------------------------------------------------------------
TTNeutralProto::TTNeutralProto(const TTNeutralProto& T): _stats(0)
{
  _copyTraitPrototypeParameters(T);
  ini_paramset();
}

// ----------------------------------------------------------------------------------------
// ini_paramset
// ----------------------------------------------------------------------------------------
void
TTNeutralProto::ini_paramset(){
  string trait = get_trait_indexStr_();

  set_paramset("ntrl"+trait, false);

  add_parameter("ntrl_nb_trait",INT2,false,false,0,0,"1");

  add_parameter("ntrl_loci"+trait,INT2,true,false,0,0,"0");
  add_parameter("ntrl_all"+trait,INT2,false,true,1,256,"255");

  add_parameter("ntrl_allelic_file"+trait, STR, false, false,0,0,"");
	add_parameter("ntrl_ini_allele_model"+trait,INT2,false,true,0,1,"0");

  add_parameter("ntrl_mutation_rate"+trait,DBL,false,true,0,1,"0");
  add_parameter("ntrl_mutation_shape"+trait,DBL,false,false,0,0,"0");
  add_parameter("ntrl_mutation_model"+trait,INT2,false,true,0,1,"0");

	add_parameter("ntrl_save_genotype",INT2,false,true,0,2,"0");
	add_parameter("ntrl_genot_sex",INT2,false,true,0,2,"0");
	add_parameter("ntrl_genot_age",INT2,false,true,0,2,"0");
	add_parameter("ntrl_genot_dir",STR,false,false,0,0,"");
	add_parameter("ntrl_genot_logtime",INT2,false,false,0,0 ,"1",true);
	add_parameter("ntrl_genot_script",STR,false,false,0,0,"");

  // genetic map
  add_parameter("ntrl_loci_positions"+trait,MAT,false,false,0,0,"");
  add_parameter("ntrl_loci_positions_random"+trait,MAT,false,false,0,0,"");

  // genotype file as input?
  add_parameter("ntrl_ini_genotypes",STR,false,false,0,0,"");
}

// ----------------------------------------------------------------------------------------
// dstor
// ----------------------------------------------------------------------------------------
TTNeutralProto::~TTNeutralProto ()
{
	if(_stats)          delete _stats;
	if(_writer)         {delete _writer; _writer=NULL;}         //_writer is static!!!

	resetTotal();
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTNeutralProto::init (Metapop* pMetapop)
{
  _popPtr = pMetapop;
  string trait = get_trait_indexStr_();

  _nb_allele = (int)get_parameter_value("ntrl_all"+trait);
  _nb_locus =  (int)get_parameter_value("ntrl_loci"+trait);
  _ploidy = 2;

  // only done for the first trait if several traits are simulated
  if(get_trait_index() <= 1){
    if(_writer)     {delete _writer;     _writer = NULL;}

    // FSTAT file for initial genotypes?
    _ini_fstat_file = get_parameter("ntrl_ini_genotypes")->get_arg();
  }

  // recombination
  if(get_parameter("ntrl_loci_positions"+trait)->isSet()){
    _inherit_func_ptr = &TTraitProto::inherit_low;
    _locus_position = new double[_nb_locus];
    setGeneticMapFixed(get_parameter("ntrl_loci_positions"+trait)->get_matrixVar());
  }
  else if(get_parameter("ntrl_loci_positions_random"+trait)->isSet()){
    _inherit_func_ptr = &TTraitProto::inherit_low;
    _locus_position = new double[_nb_locus];
    _random_genetic_map_matrix = get_parameter("ntrl_loci_positions_random"+trait)->get_matrix();
     setGeneticMapRandom(_random_genetic_map_matrix); 
  }
  else {
    _inherit_func_ptr = &TTraitProto::inherit_free;
    _locus_position = NULL;
  }


  // read the allelic file if available
  string allelicFile = get_parameter("ntrl_allelic_file"+trait)->get_arg();  // if empty parameter is not set
  if(!allelicFile.empty()) {
    if(!STRING::file_exists(allelicFile)){
      allelicFile = FileHandler::get_base_path()+allelicFile;
      if(!STRING::file_exists(allelicFile)) fatal("Allelic file '%s' (ntrl) cannot be found!\n", allelicFile.c_str());
    }
    read_allele_file(allelicFile);
  }

  // mutation
  if(set_mutation_rates("ntrl", trait)){ // if mutation rates differ among loci
    if(_mutationFreqFile){
      switch(_mut_model) {
        case 0: _mutate_func_ptr = &TTraitProto::mutate_RMM_single;    break;
        case 1: _mutate_func_ptr = &TTraitProto::mutate_IMM_single;    break;
      }
    }
    else{
      switch(_mut_model) {
        case 0: _mutate_func_ptr = &TTraitProto::mutate_KAM_single;    break;
        case 1: _mutate_func_ptr = &TTraitProto::mutate_SSM_single;    break;
      }
		}
  }
  else{           // same mutation rates for all loci
    if(!_mut_rate_mean) _mutate_func_ptr = &TTraitProto::mutate_NULL;
    else{
      if(_mutationFreqFile){
        switch(_mut_model) {
          case 0: _mutate_func_ptr = &TTraitProto::mutate_RMM;    break;
          case 1: _mutate_func_ptr = &TTraitProto::mutate_IMM;    break;
        }
      }
      else{
        switch(_mut_model) {
          case 0: _mutate_func_ptr = &TTraitProto::mutate_KAM;    break;
          case 1: _mutate_func_ptr = &TTraitProto::mutate_SSM;    break;
        }
		  }
    }
  }

  // initial alleles?
  if(_initAlleleFreqFile) _ini_allele_func_ptr = &TTraitProto::ini_allele_freq_dist;
  else{
    _ini_allele_model  = (int) get_parameter_value("ntrl_ini_allele_model"+trait);
    if(_ini_allele_model) _ini_allele_func_ptr = &TTraitProto::ini_allele_freq_monomorph;
    else                  _ini_allele_func_ptr = &TTraitProto::ini_allele_freq_uniform;
  }
}

// ----------------------------------------------------------------------------------------
// executeBeforeEachGeneration
// ----------------------------------------------------------------------------------------
void
TTNeutralProto::executeBeforeEachGeneration(const int& gen){
  map<string, Param*>* pParam = _paramSet->getTemporalParams(gen);

  // if it is a temporal parameter
  if(pParam){
    // check if a change has to be made
    map<string, Param*>* pMap = _paramSet->updateTemporalParams(gen);
    if(pMap){
      // iterate through the map and performe the updates
      map<string, Param*>::iterator pos = pMap->begin();
      string trait = get_trait_indexStr_();

      for(; pos != pMap->end(); ++pos){
       if(pos->first == "ntrl_genot_logtime"+trait && _writer){
          _writer->set_GenerationOccurrence((unsigned int)pos->second->get_value());
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
// get_info
// ----------------------------------------------------------------------------------------
string
TTNeutralProto::get_info()
{
	string text;

	if(_trait_index) text = STRING::int2str(_trait_index) + ". neutral marker type: ";
	else             text = "Neutral marker type: ";

	text +=	STRING::int2str(_nb_locus) + " loci; "
	     + 	STRING::int2str(_nb_allele) + " alleles";

	return text;
}
// ----------------------------------------------------------------------------------------
// loadFileServices
// ----------------------------------------------------------------------------------------
void TTNeutralProto::loadFileServices  (FileServices* loader)
{
	//writer
	int choice = (int)get_parameter_value("ntrl_save_genotype");
	if(choice) {
		if(!_writer) _writer = new TTNeutralFH(this);
		int logtime = (int)get_parameter_value("ntrl_genot_logtime");
		if(logtime<1) fatal("Parameter 'ntrl_genot_logtime' must be 1 or larger!\n");
		_writer->set(true, true, 1,
								 logtime,
								 0,
								 get_parameter("ntrl_genot_dir")->get_arg(),
								 get_parameter("ntrl_genot_script")->get_arg(),
								 (int)get_parameter_value("ntrl_genot_sex"),
								 (int)get_parameter_value("ntrl_genot_age"),
								 this, choice);

		if(get_trait_index()<=1) loader->attach(_writer);
	}
  else if(_writer) {
		delete _writer;
		_writer = NULL;
	}
}

// ----------------------------------------------------------------------------------------
// loadStatServices
// ----------------------------------------------------------------------------------------
void TTNeutralProto::loadStatServices  (StatServices* loader)
{
  //allocate the stat handler
  if(_stats)  delete _stats;
  _stats = new TTNeutralSH(this);
  loader->attach(_stats);
}

// ----------------------------------------------------------------------------------------
// clone
// ----------------------------------------------------------------------------------------
TTNeutral* TTNeutralProto::hatch ()
{
  TTNeutral* new_trait = new TTNeutral();
  new_trait->set_from_prototype(this);
  return new_trait;
}
//----------------------------------------------------------------------------------------
// operator=
// ----------------------------------------------------------------------------------------
TTNeutral&
TTNeutral::operator= (const TTrait& T)
{
  pTraitProto->setTTraitProto(*T.pTraitProto);
  const TTNeutral& TN = dynamic_cast<const TTNeutral&> (T);

  if(this != &TN) {
    pProto = TN.pProto;

    reset();
    init();

    for(int i = 0; i < pProto->_nb_locus; ++i) {
      for(int j = 0; j < pProto->_ploidy; ++j) {
        sequence[i][j] = TN.sequence[i][j];
      }
    }
  }

  return *this;
}

//----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
void*
TTNeutral::get_allele(int loc, int all)  const{
  return ( !(loc<pProto->_nb_locus)||!(all<pProto->_ploidy) ? 0 : (void*)&sequence[loc][all]);
}

//----------------------------------------------------------------------------------------
// operator==
// ----------------------------------------------------------------------------------------
bool TTNeutral::operator== (const TTrait& T)
{
  if(*pTraitProto != *T.pTraitProto) return false;

  const TTNeutral& TN = dynamic_cast<const TTNeutral&> (T);

  if(this != &TN || *pProto != *TN.pProto) return false;

  return true;
}

//----------------------------------------------------------------------------------------
// operator!=
// ----------------------------------------------------------------------------------------
bool TTNeutral::operator!= (const TTrait& T)
{
  if(!((*this) == T))
    return true;
  else
    return false;
}

// ----------------------------------------------------------------------------------------
// Destructor
// ----------------------------------------------------------------------------------------
TTNeutral::~TTNeutral()
{
}

// ----------------------------------------------------------------------------------------
// set_from_prototype (used during hatching)
// ----------------------------------------------------------------------------------------
void TTNeutral::set_from_prototype(TTraitProto* T){
  pProto = dynamic_cast<TTNeutralProto*> (T);
  pTraitProto = T;
}

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
void TTNeutral::init ()
{
  if(sequence != NULL) fatal("TTNeutral::init::sequence is not NULL !\n");

  int nbLocus = pProto->_nb_locus;
  int ploidy = pProto->_ploidy;
  sequence = new unsigned char* [nbLocus];
  for(int i=0; i < nbLocus; i++){
    sequence[i] = new unsigned char[ploidy];
  }
}
// ----------------------------------------------------------------------------------------
// ini_sequence
// ----------------------------------------------------------------------------------------
void TTNeutral::ini_sequence (Patch* patch)
{
  (pProto->*(pProto->_ini_allele_func_ptr))(sequence, patch);
}
// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void TTNeutral::reset()
{
}


// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void TTNeutral::show_up ()
{
  message("\n  Trait's type: ntrl\n\
       locus: %i\n\
     alleles: %i\n\
    sequence:",pProto->_nb_locus,pProto->_nb_allele);
  for(int i = 0; (i < pProto->_nb_locus && i < 10); i++){
    message("\n              %i %i",(int)sequence[i][0],(int)sequence[i][1]);
  }
  message("\n");
}


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** StatHandler ********/

// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool TTNeutralSH::init()
{
	StatHandler<TTNeutralSH>::init();

	allocateTable();
  return true;
}

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool TTNeutralSH::setStatRecorders(const string& t)
{
       if(set_stat_coancestry     (t, "n", "ntrl")){}     // adults and offspring
  else if(set_stat_fstat          (t, "n", "ntrl")){}     // adults and offspring
  else if(set_stat_all_freq_local (t, "n", "ntrl")){}     // adults and offspring
  else if(set_stat_all_freq_global(t, "n", "ntrl")){}     // adults and offspring
  else return false;

	return true;
}
