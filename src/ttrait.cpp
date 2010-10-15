/** @file ttrait.cpp
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

#include <cmath>
#include "ttrait.h"
#include "simulation.h"
#include "param.h"
#include <iomanip>
#include "patch.h"
#include "metapop.h"


#include <algorithm>
using namespace std;

//------------------------------------------------------------------------------
/** initialization of static class members */
double*         TTraitProto::_chromosomeSize = NULL;
vector<double>  TTraitProto::_recombinationPosF;
vector<double>  TTraitProto::_recombinationPosM;
bool            TTraitProto::_recombinationStartF=0;
bool            TTraitProto::_recombinationStartM=0;
int             TTraitProto::_chromosomeNb=0;
double          TTraitProto::_chromosomeLength=0;
double          TTraitProto::_nbRecombinations=0;
vector<double>  TTraitProto::_chromosomeSizeVector;

//------------------------------------------------------------------------------
TTraitProto::~TTraitProto(){
	resetTotal();
}

//------------------------------------------------------------------------------
/** same as the copy constructor */
// ----------------------------------------------------------------------------------------
void
TTraitProto::_copyTraitPrototypeParameters(const TTraitProto& T){
	fatal("Not implemented!\n");
}

//------------------------------------------------------------------------------
/** same as the copy constructor */
// ----------------------------------------------------------------------------------------
void
TTraitProto::resetTotal(){
	if(_locus_position) {delete[] _locus_position; _locus_position=NULL;}

  if(_initAlleleFreq){
    if(_nb_locus == 1 || _initAlleleFreq[0][0] != _initAlleleFreq[0][1]){   // each locus has its one array
       ARRAY::delete_3D(_initAlleleFreq, _initAlleleFreqCols, _nb_locus);
    }
    else{                                                        // all loci point to the first array
		  delete[] _initAlleleFreq[0][0];
		  ARRAY::delete_2D(_initAlleleFreq, _initAlleleFreqCols);
    }
	}

  if(_mutationFreq){
    if(_nb_locus == 1 || _mutationFreq[0] != _mutationFreq[1]){   // each locus has its one array
       ARRAY::delete_2D(_mutationFreq, _nb_locus);
    }
    else{                                                        // all loci point to the first array
		  delete[] _mutationFreq[0];
		  delete[] _mutationFreq;
		  _mutationFreq = NULL;
    }
	}

  if(_mut_rate) {delete[] _mut_rate; _mut_rate=NULL;}
}

// ----------------------------------------------------------------------------------------
/** get_trait_index */
// ----------------------------------------------------------------------------------------
int
TTraitProto::get_trait_index() const {
  return _trait_index;
}

string
TTraitProto::get_trait_indexStr() const {
  if(!_trait_index) return "";
  return STRING::int2str(_trait_index);
}

string
TTraitProto::get_trait_indexStr_() const {
  if(!_trait_index) return "";
  return STRING::int2str_(_trait_index);
}

string
TTraitProto::get_trait_indexStr_t() const {
  if(!_trait_index) return "";
  return STRING::int2str_(_trait_index, (int)0, "t");
}

// ----------------------------------------------------------------------------------------
// sequence initialization
// ----------------------------------------------------------------------------------------
/** initialization of the sequences based on the _initAlleleFreq array (settings file) */
void
TTraitProto::ini_allele_freq_dist(unsigned char** seq, Patch* patch){
  int p, l;
  int cur_patch = patch->get_ID()%_initAlleleFreqCols;
  for(l = 0; l < _nb_locus; ++l) {
    for(p = 0; p < _ploidy; ++p) {
			seq[l][p] = (unsigned char)SimRunner::r.AfterDistribution(_initAlleleFreq[cur_patch][l], _nb_allele);
    }
  }
}

/** initialization of the sequences: alleles are randomly distributed -> maximal variance */
void
TTraitProto::ini_allele_freq_uniform(unsigned char** seq, Patch* patch){
  int p, l;
  for(p = 0; p < _ploidy; ++p) {
    for(l = 0; l < _nb_locus; ++l) {
       seq[l][p] = (unsigned char)SimRunner::r.Uniform(_nb_allele);
    }
  }
}

/** initialization of the sequences: populations are fixed for the middle allele -> minimal variance */
void
TTraitProto::ini_allele_freq_monomorph(unsigned char** seq, Patch* patch){
	int p, l;
	unsigned char allele = (unsigned char)_nb_allele/2;
	for(p = 0; p < _ploidy; ++p) {
		for(l = 0; l < _nb_locus; ++l) {
			 seq[l][p] = allele;
    }
	}
}


// ----------------------------------------------------------------------------------------
// set_ini_allelicArray
// ----------------------------------------------------------------------------------------
/** set up the inital frequency array/ get the dimensions (nothing else!)
  */
vector<map<int, string> >*
TTraitProto::set_ini_allelicArray (TMatrix* mat, const int& i){
  // check if the frequencies are identical for all patches
  vector<map<int, string> >* pStrMatrix = NULL;
  if(mat->get(0, i) != my_NAN) _initAlleleFreqCols = 1;  // no matrix: all pops have the same init freqs
  else{
    pStrMatrix = mat->get_strMatrix();
    assert(pStrMatrix);
    map<int, string>::iterator pos = (*pStrMatrix).begin()->find(i);  // get the corresponding map element from the first row
    assert(pos != (*pStrMatrix).begin()->end());
    //find the size of the matrix
		_initAlleleFreqCols = STRING::strMatrix2vector<double>(pos->second).size();

    // check the number of carrying capacities with the number of patches
		int nbPatch = _popPtr->getPatchNbr();
    if(_initAlleleFreqCols>nbPatch) warning("There are more inital allele frequencies defined than patches! Only a part of the inital allele frequencies is considered!\n");
    else if(nbPatch % _initAlleleFreqCols) warning("The number of inital allele frequencies is not an entire subset of the number of patches!\n");
  }
  ARRAY::create_3D(_initAlleleFreq, _initAlleleFreqCols, _nb_locus, _nb_allele, (double)my_NAN);
  _initAlleleFreqFile = true;
  return pStrMatrix;
}

// ----------------------------------------------------------------------------------------
// set_allelic_valuesFreq
// ----------------------------------------------------------------------------------------
/** read the allelic values and their frequencies from a file (used for mutation model IMM and RMM)
  */
void
TTraitProto::read_allele_file (string filename){
  // read the file
  TMatrix mat;
  map<string, int> fileInfo;
  bool hasLocus = true;
  int i, l, a, p, nbValues=0;
  vector<map<int, string> >* pStrMatrix = NULL;
  double* array = NULL;

  // read the file
	fileInfo = mat.read_matrix(filename);
  int nbLines = mat.getNbRows();
	int nbCols  = mat.getNbCols();

  const int colsSize = 4;
  int cols[colsSize]        = {my_NAN,my_NAN,my_NAN,my_NAN};  // 0: locus; 1: allele; 2: mut_freq; 3: ini_freq;
  string colsName[colsSize] = {"col_locus","col_allele","col_mut_freq", "col_ini_freq"};

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
    warning("Allelic file for trait '%s' (%s): Column '%s' (%i) will not be considered!\n",
            _type.c_str(), filename.c_str(), pos->first.c_str(), pos->second);
  }

  // check if all cols are available and also dimensions of the matrix are met and create the arrays
  if(cols[0] == my_NAN) hasLocus = false; // fatal("Allelic file (%s): Locus column is missing!", filename.c_str());
  if(cols[1] == my_NAN) fatal("Allelic file (%s): Allele column is missing!", filename.c_str());
  if(cols[2] != my_NAN){    // mutation frequency
    ARRAY::create_2D(_mutationFreq, _nb_locus, _nb_allele, (double)my_NAN);
    _mutationFreqFile = true;
  }

  if(cols[3] != my_NAN){     // inital frequencies
    pStrMatrix = set_ini_allelicArray (&mat, cols[2]);
    if(_initAlleleFreqCols) array = new double[_initAlleleFreqCols];
  }
  else fatal("Allelic file (%s): No information is passed (initial frequencies are not specified)!\n", filename.c_str());

  for (i=0; i<colsSize; ++i){
    if(cols[i] > nbCols) fatal("Allelic file (%s): File info is wrong: matrix has only %i columns and not %i!\n", filename.c_str(), nbCols, cols[i]);
  }

  // copy the allelic effects and their frequencies
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
    if(   (_mutationFreq  && _mutationFreq[l][a] != my_NAN)
       || (_initAlleleFreq && _initAlleleFreq[0][l][a] != my_NAN)){

      fatal("Allelic values: Allelic value was already specified for locus %i and allele %i!\n", l+1, a+1);
    }

    // set the values
    if(hasLocus){         // locus specific settings
      if(_mutationFreq)  _mutationFreq[l][a]  = mat.get(i,cols[2]); // mut  frequency
      if(_initAlleleFreq){
        if(pStrMatrix){
          STRING::strMatrix2array((*pStrMatrix)[i].find(cols[3])->second, array, _initAlleleFreqCols);
          for(p=0; p<_initAlleleFreqCols; ++p){
            _initAlleleFreq[p][l][a] = array[p]; // init frequencies per population
          }
        }
        else _initAlleleFreq[0][l][a] = mat.get(i,cols[3]); // init frequencies
      }
      ++nbValues;
    }
    else {                // all loci have the same settings
      for(l = 0; l < _nb_locus; ++l){
        if(_mutationFreq)  _mutationFreq[l][a]  = mat.get(i,cols[2]); // mut  frequency
        if(_initAlleleFreq){
          if(pStrMatrix){
            STRING::strMatrix2array((*pStrMatrix)[i].find(cols[3])->second, array, _initAlleleFreqCols);
            for(p=0; p<_initAlleleFreqCols; ++p){
              _initAlleleFreq[p][l][a] = array[p]; // init frequencies per population
            }
          }
          else _initAlleleFreq[0][l][a] = mat.get(i,cols[3]); // init frequencies
        }
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
	if(_mutationFreqFile){
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
}

// ----------------------------------------------------------------------------------------
// set_mutation_rates
// ----------------------------------------------------------------------------------------
/** this function sets the mutation rates dependig on the parameters. If the
  * mutation rates differ among loci, a true is returned.
  */
bool
TTraitProto::set_mutation_rates(string type, const string& trait){
  _mut_model     = (int)get_parameter_value(type+"_mutation_model"+trait);
  if(_mut_rate) {delete[] _mut_rate; _mut_rate=NULL;}

  // if the single mutation rates are given as a matrix
  if(get_parameter(type+"_mutation_rate"+trait)->is_matrix()){
    TMatrix* m = get_parameter(type+"_mutation_rate"+trait)->get_matrix();
    int nb = m->get_dims(NULL);
    double* v = m->get();
    // check the number of mutation rates with the number of patches
    if(nb>_nb_locus) warning("There are more mutation rates than loci defined! Only a part of the mutation rates is considered!\n");
    else if(_nb_locus % nb) warning("The number of defined mutation rates is not a entire subset of the number of loci!\n");
    _mut_rate = new double[_nb_locus];
    for(int l=0; l<_nb_locus; ++l){
      _mut_rate[l] = v[l % nb];
    }
    return true;
  }

  // if the single mutation rates are given by the variance
  int mut_shape  = (int) get_parameter_value(type+"_mutation_shape"+trait);
  _mut_rate_mean = get_parameter_value(type+"_mutation_rate"+trait);
  if(mut_shape){   // mutation rates vary between loci
    if(mut_shape<0) fatal("The parameter '%s_mutation_shape%s' has to be positive!\n", type.c_str(), trait.c_str());
    _mut_rate = new double[_nb_locus];
    for(int l=0; l<_nb_locus; ++l){
      _mut_rate[l] = SimRunner::r.Gamma(mut_shape)/mut_shape*_mut_rate_mean;
    }
    return true;
  }

  // the mutation rates do not differ among loci
  return false;
}

// ----------------------------------------------------------------------------------------
// mutate_SSM
// ----------------------------------------------------------------------------------------
void
TTraitProto::mutate_SSM(unsigned char** seq)
{
  unsigned int NbMut;
  unsigned char * curPos;

	for(NbMut = SimRunner::r.Poisson(_ploidy*_nb_locus*_mut_rate_mean) ; NbMut != 0; --NbMut) {
		curPos = &seq[SimRunner::r.Uniform( _nb_locus)][SimRunner::r.Uniform( _ploidy)];

    //alleles values are from 0 to NtrlAll - 1 !!!
		if(SimRunner::r.Bool() && *curPos < _nb_allele-1){   // Direction
      *curPos += 1; //one step to the right
    }
    else if(*curPos > 0){ // !direction || all==_nb_allele
      *curPos -= 1; //one step to the left
    }
    else{ //!direction && all == 0
      *curPos += 1; // on step to the right if the curfent allele is 0
    }
  }
}

// ----------------------------------------------------------------------------------------
// mutate_SSM
// ----------------------------------------------------------------------------------------
void
TTraitProto::mutate_SSM_single(unsigned char** seq)
{
  int l, a;
  unsigned char * curPos;

  for(l=0; l<_nb_locus; ++l){ // for each locus
    for(a=0; a<_ploidy; ++a){ // for each allele
      if(SimRunner::r.Uniform()<_mut_rate[l]){
		    curPos = &seq[l][a];
        //alleles values are from 0 to NtrlAll - 1 !!!
				if(SimRunner::r.Bool() && *curPos < _nb_allele-1){   // Direction
          *curPos += 1; //one step to the right
        }
        else if(*curPos > 0){ // !direction || all==_nb_allele
          *curPos -= 1; //one step to the left
        }
        else{ //!direction && all == 0
          *curPos += 1; // on step to the right if the curfent allele is 0
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
// mutate_KAM
// ----------------------------------------------------------------------------------------
void
TTraitProto::mutate_KAM(unsigned char** seq)
{
  unsigned int NbMut;
  unsigned char mut;
  unsigned char * curPos;

	for(NbMut = SimRunner::r.Poisson(_ploidy*_nb_locus*_mut_rate_mean) ; NbMut != 0; NbMut--) {
		curPos = &seq[SimRunner::r.Uniform( _nb_locus)][SimRunner::r.Uniform( _ploidy)];

    //assign an arbitrary allele value:
    do{
			mut = (unsigned char) (SimRunner::r.Uniform(_nb_allele));
    } while (mut == *curPos);  // it has to change
    *curPos = mut;
  }
}

// ----------------------------------------------------------------------------------------
// mutate_KAM
// ----------------------------------------------------------------------------------------
void
TTraitProto::mutate_KAM_single(unsigned char** seq)
{
  int l, a;
  unsigned char mut;
  unsigned char * curPos;

  for(l=0; l<_nb_locus; ++l){ // for each locus
    for(a=0; a<_ploidy; ++a){ // for each allele
      if(SimRunner::r.Uniform()<_mut_rate[l]){
        curPos = &seq[l][a];
        //assign an arbitrary allele value:
        do{
		  	  mut = (unsigned char) (SimRunner::r.Uniform(_nb_allele));
        } while (mut == *curPos);  // it has to change
        *curPos = mut;
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
// mutate_IMM
// ----------------------------------------------------------------------------------------
/** Incremental mutation: at each mutation an allelic effect is drawn proportional to their frequencies
  * This drawn allelic affect is added to the present one.
  * If the new allele is out of range, nothing is made (following fl)
  */
void
TTraitProto::mutate_IMM(unsigned char** seq)
{
  unsigned int NbMut,mutLocus;
  int mutStep, newPos;
  unsigned char * curPos;

	for(NbMut = SimRunner::r.Poisson(_ploidy*_nb_locus*_mut_rate_mean) ; NbMut != 0; NbMut--) {
		mutLocus  = SimRunner::r.Uniform( _nb_locus);
		curPos = &seq[mutLocus][SimRunner::r.Uniform( _ploidy)];

    // perform the mutation
		mutStep = SimRunner::r.AfterDistribution(_mutationFreq[mutLocus], _nb_allele) - _nb_allele/2;  // _nb_allele/2 is the middle index
    newPos = mutStep+(*curPos);
    if (newPos>=0 && newPos<_nb_allele){   // do nothing if the new allele is out of range
      *curPos += mutStep;
    }
  }
}

// ----------------------------------------------------------------------------------------
// mutate_IMM
// ----------------------------------------------------------------------------------------
/** Incremental mutation: at each mutation an allelic effect is drawn proportional to their frequencies
  * This drawn allelic affect is added to the present one.
  * If the new allele is out of range, nothing is made (following fl)
  */
void
TTraitProto::mutate_IMM_single(unsigned char** seq)
{
  int l, a;
  int mutStep, newPos;
  unsigned char * curPos;

  for(l=0; l<_nb_locus; ++l){ // for each locus
    for(a=0; a<_ploidy; ++a){ // for each allele
      if(SimRunner::r.Uniform()<_mut_rate[l]){
        curPos = &seq[l][a];
        // perform the mutation
		    mutStep = SimRunner::r.AfterDistribution(_mutationFreq[l], _nb_allele) - _nb_allele/2;  // _nb_allele/2 is the middle index
        newPos = mutStep+(*curPos);
        if (newPos>=0 && newPos<_nb_allele){   // do nothing if the new allele is out of range
          *curPos += mutStep;
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
// mutate_RMM
// ----------------------------------------------------------------------------------------
/** At each mutation a new allelic effect is drawn proportional to their frequencies
  * RMM: rekursive mutation model
  */
void
TTraitProto::mutate_RMM(unsigned char** seq)
{
  unsigned int NbMut,mutLocus,i;
  unsigned char mutPos;
  unsigned char * curPos;

	for(NbMut = SimRunner::r.Poisson(_ploidy*_nb_locus*_mut_rate_mean) ; NbMut != 0; NbMut--) {
		mutLocus = SimRunner::r.Uniform( _nb_locus);
		curPos = &seq[mutLocus][SimRunner::r.Uniform( _ploidy)];

    //assign an arbitrary allele value:
    i=0;
    do{
      mutPos = (unsigned char) (SimRunner::r.AfterDistribution(_mutationFreq[mutLocus], _nb_allele));
      ++i;                                 // used to asure that the loop ends
   		} while (mutPos == *curPos && i<1e4);  // it has to change (except if after 1e4 trials it has not yet changed)
    *curPos = mutPos;
  }
}

// ----------------------------------------------------------------------------------------
// mutate_RMM
// ----------------------------------------------------------------------------------------
/** At each mutation a new allelic effect is drawn proportional to their frequencies
  * RMM: rekursive mutation model
  */
void
TTraitProto::mutate_RMM_single(unsigned char** seq)
{
  int l, a, i;
  unsigned char mutPos;
  unsigned char * curPos;

  for(l=0; l<_nb_locus; ++l){ // for each locus
    for(a=0; a<_ploidy; ++a){ // for each allele
      if(SimRunner::r.Uniform()<_mut_rate[l]){
        curPos = &seq[l][a];
        //assign an arbitrary allele value:
        i=0;
        do{
          mutPos = (unsigned char) (SimRunner::r.AfterDistribution(_mutationFreq[l], _nb_allele));
          ++i;                                 // used to asure that the loop ends
        } while (mutPos == *curPos && i<1e4);  // it has to change (except if after 1e4 trials it has not yet changed)
        *curPos = mutPos;
      }
    }
  }
}

//------------------------------------------------------------------------------
/** for the derived classes to compare directly the paramters in this base calss */
TTraitProto&
TTraitProto::setTTraitProto(const TTraitProto& T){             // for operator=
  _copyTraitPrototypeParameters(T);
  return *this;
}

//------------------------------------------------------------------------------
bool
TTraitProto::isEqualTTraitProto(const TTraitProto& T){         // for operator==
  if(_type              != T._type              )  return false;
  if(_ploidy            != T._ploidy            )  return false;
  if(_nb_allele         != T._nb_allele         )  return false;
  if(_nb_locus          != T._nb_locus          )  return false;
  if(_mut_rate          != T._mut_rate          )  return false;
  if(_absolute_index    != T._absolute_index    )  return false;
  if(_trait_index       != T._trait_index       )  return false;
  if(_locus_position    != T._locus_position    )  return false;

  return true;
}

//------------------------------------------------------------------------------
bool
TTraitProto::isUnequalTTraitProto(const TTraitProto& T){       // for operator!=
  return !(*this == T);
}

//------------------------------------------------------------------------------
/** print the genetci map of the trait to the stream */
void
TTraitProto::print_genetic_map(ofstream& FILE){
  unsigned int c, size, l;

  // for each chromosome
  FILE << "\n" << setw(10) << left << _type  << "{";
  for(c=0, size=_locus_position_vector.size(); c<size; ++c){
    FILE << "{";
    // print each locus position
    for(l=0; l<_locus_position_vector[c].size(); ++l){
      if(l) FILE << " ";
      FILE << _locus_position_vector[c][l];
    }
    FILE << "}";
  }
  FILE << "}";
}

//------------------------------------------------------------------------------
/** sets the locus positions according to the input parameters.
  * The input matrix consists of a list of chromosomes defined by there lenght in cM.
  * The locus positions are randomly assigned to positions on these chromosomes.
  * The _chromosomeSizeVector (static) and _locus_position_vector (per trait) is filled.
  */
void
TTraitProto::setGeneticMapRandom(TMatrix* matrix){
  // check the maximum size of each chromosome
  unsigned int nbChromosome = matrix->get_dims(NULL); // nbr of chromosome
  unsigned int c;
  int l;
  double upperLimit, lowerLimit;

  double size, tot_size=0;
  for(c=0; c<nbChromosome; ++c){
    size=matrix->get(0, c);
    if(_chromosomeSizeVector.size()<=c) _chromosomeSizeVector.push_back(0);  // create the new chromosome
    if(size == my_NAN ) continue;               // this chromsome contains no loci: jump to the next one
    if(_chromosomeSizeVector[c]<size) _chromosomeSizeVector[c] = size;       // resize it if it is bigger
    tot_size += size;
  }

  // set randomly the locus positions on a temp vector and sort it
  double* vecTemp = new double[_nb_locus];
  for(l=0; l<_nb_locus; ++l){
		vecTemp[l] = tot_size*SimRunner::r.Uniform();  //position between 0 and tot_size
  }
  // sort(vecTemp, vecTemp+nbLoci);    // no idea why this does not work!!!!
  ARRAY::quicksort(vecTemp, _nb_locus);

  // copy the positions to the single chromosomes
  // this has to be done as another trait could have defined a bigger chromosome
 
  _locus_position_vector.clear();
  l = 0;
  upperLimit = lowerLimit = 0;
  for(c=0; c<_chromosomeSizeVector.size() && l<_nb_locus; ++c){   // for each chromosome
    _locus_position_vector.push_back(vector<double>());           // create the new trait chromosome
    lowerLimit = upperLimit;
    upperLimit += matrix->get(0, c);                              // get the size of the current chromosome
    if(upperLimit-lowerLimit == my_NAN ){
      upperLimit = lowerLimit;
      continue;               // this chromsome contains no loci: jump to the next one
    }

    // get all the loci assigned to this chromosome
    while(l<_nb_locus && upperLimit>=vecTemp[l]){
      _locus_position_vector[c].push_back(vecTemp[l]-lowerLimit);
      ++l;
    }
  }

  delete[] vecTemp;
}

//------------------------------------------------------------------------------
/** sets the locus positions according to the input parameters.
  * The input matrix consists of the exact positions of the loci defined in cM from the BEGINNING on.
  * The _chromosomeSizeVector (static) and _locus_position_vector (per trait) is filled.
  */
void
TTraitProto::setGeneticMapFixed(TMatrixVar<double>* matrix){
  unsigned int c, l, size;
  unsigned int nbChromosome = matrix->getNbRows();
  double chromosomeSize;
  int nbLoci=0;

  _locus_position_vector.clear();

  // for each chromosome
  for(c=0; c<nbChromosome; ++c){
    // check if the static chromosome alredy exists, otherwise create it
    if(_chromosomeSizeVector.size()<=c) _chromosomeSizeVector.push_back(0);

    // create a chromosome for the trait vector
    _locus_position_vector.push_back(vector<double>());

    // add the postions of the loci to this chromosome
    size = matrix->getNbCols(c);
    if(!size) continue;               // if the chromosome is empty
    nbLoci += size;
    for(l=0; l<size; ++l){
      _locus_position_vector[c].push_back(matrix->get(c, l));
    }
    sort(_locus_position_vector[c].begin(), _locus_position_vector[c].end());
    chromosomeSize = _locus_position_vector[c].back();

    // if the current cromosome is bigger than the static chromosome resize it
    if(_chromosomeSizeVector[c] < chromosomeSize) _chromosomeSizeVector[c] = chromosomeSize;
  }

  // control if the number of loci is correct
  if(nbLoci != _nb_locus){
    fatal("Parameter '%s' contains %i instead of %i loci!",
          _type.replace(_type.find("_"),0,"_loci_positions").c_str(), nbLoci, _nb_locus);
  }
}

// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
/** inheritance with recombination */
void
TTraitProto::inherit_low (TTrait* mother, TTrait* father, TTrait* child)
{
  int i, curChrom;
  vector<double>::iterator nextRecomb;
  unsigned char** mother_seq = (unsigned char**)mother->get_sequence();
  unsigned char** father_seq = (unsigned char**)father->get_sequence();
  unsigned char** child_seq  = (unsigned char**)child->get_sequence();

  // get gamete of mother
  curChrom = _recombinationStartF;          // get the starting chromosome
  nextRecomb = _recombinationPosF.begin();  // get the positon of the first recombination event
  for(i = 0; i < _nb_locus; ++i) {
    while(_locus_position[i] > *nextRecomb){             // the last position of recomb is always bigger or equal
      curChrom = curChrom ? 0 : 1;
      ++nextRecomb;
    }
    child_seq[i][0] = mother_seq[i][curChrom];
  }

  // get gamete of father
  curChrom = _recombinationStartM;          // get the starting chromosome
  nextRecomb = _recombinationPosM.begin();  // get the positon of the first recombination event
  for(i = 0; i < _nb_locus; ++i) {
    while(_locus_position[i] > *nextRecomb){            // the last position of recomb is always bigger or equal
      curChrom = curChrom ? 0 : 1;
      ++nextRecomb;
    }
    child_seq[i][1] = father_seq[i][curChrom];
  }
}

// ----------------------------------------------------------------------------------------
// inherit
// ----------------------------------------------------------------------------------------
/** inheritance for unlinked loci */
void
TTraitProto::inherit_free(TTrait* mother, TTrait* father, TTrait* child)
{
  unsigned char** mother_seq = (unsigned char**)mother->get_sequence();
  unsigned char** father_seq = (unsigned char**)father->get_sequence();
  unsigned char** child_seq  = (unsigned char**)child->get_sequence();

  for(int i = 0; i < _nb_locus; ++i) {
		child_seq[i][0] = mother_seq[i][SimRunner::r.Bool()];
		child_seq[i][1] = father_seq[i][SimRunner::r.Bool()];
  }
}

//------------------------------------------------------------------------------
/** destrucor */
TTrait::~TTrait ( ){
  if(sequence){
    for(int i = 0; i<pTraitProto->_nb_locus; ++i){
      if(sequence[i]) delete[] sequence[i];
    }
    delete[] sequence;
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void
TTrait::mutate(){
  (pTraitProto->*(pTraitProto->_mutate_func_ptr))(sequence);
}

void
TTrait::inherit(TTrait* mother, TTrait* father) {
  (pTraitProto->*(pTraitProto->_inherit_func_ptr)) (mother, father, this);
}

void*
TTrait::get_allele(int loc, int all)  const  {
  if(loc<pTraitProto->_nb_locus && all<pTraitProto->_nb_allele) return (void*)&sequence[loc][all];
  return 0;
}

//------------------------------------------------------------------------------
/** The super chromosome is created:
  *   - _chromosomeNb (int):        number of chromosomes
  *   - _chromosomeLength (double): length of super chromosome
  *   - _chromosomeSize (double*):  positions of the change from one to the other chromosome
  *   - _nbRecombinations (double): mean number of recombinations for the total super chromosome
  * If no genetic map is inizilised false is returned, othersie true
  */
bool
TTraitProto::initStaticGeneticMap(){
  if(_chromosomeSize) delete[] _chromosomeSize;

  // get the dimensions of the super chromosome (a cumulative size vector of the chromosomes)
  _chromosomeNb = _chromosomeSizeVector.size();
  if(!_chromosomeNb){
    return false;
  }
  _chromosomeSize = new double[_chromosomeNb];

  // set the positions of the changes of one to the other chromosomes
  double cur_size;   // size of the current chromosome
  for(int i=0; i<_chromosomeNb; ++i){
    cur_size = _chromosomeSizeVector[i];
    if(!i) _chromosomeSize[i] = cur_size;                         // first chromosome
    else   _chromosomeSize[i] = cur_size + _chromosomeSize[i-1];  // cumulative size
  }
  _chromosomeLength=_chromosomeSize[_chromosomeNb-1];             // total length
  _nbRecombinations=_chromosomeLength/100;
  return true;
}

//------------------------------------------------------------------------------
/** prints to the console the genetic map */
void
TTraitProto::printGeneticMapInfo(){
  if(!_chromosomeNb){
    message("\n    not used: all loci are independent");
    return;
  }

  message("\n    %i chromosoms of", _chromosomeNb);
  for(int i=0; i<_chromosomeNb; ++i){
    message(" %g", _chromosomeSizeVector[i]);
  }
  message(" cM (in total %g cM)\n", _chromosomeLength);
}

//------------------------------------------------------------------------------
/** create the super chromosome (containing the specific postions of the loci)
  * for each trait
  */
void
TTraitProto::initTraitGeneticMap(){
  double start=0;
  // for each chromosome
  for (unsigned int c=0, l=0; c<_locus_position_vector.size(); ++c){
    // for each locus
    for (unsigned int i=0; i<_locus_position_vector[c].size(); ++i, ++l){
      _locus_position[l] = start + _locus_position_vector[c][i];
    }
    start += _chromosomeSizeVector[c];
  }

}

//------------------------------------------------------------------------------
/** compute the recombination positions */
void
TTraitProto::_recombine(vector<double>& vecRecombs, bool& start){
  // clear first the old recombination events
  vecRecombs.clear();

  // get the new recombination events
	for(int i=SimRunner::r.Poisson(_nbRecombinations); i>0; --i){
     vecRecombs.push_back(_chromosomeLength*SimRunner::r.Uniform());
  }

  // add the choice of the starting chromosome as a recombination of the super chromosome
  // first one:
	start = SimRunner::r.Bool();

  // for all postions where a chromosome ends we have a recombination rate of 0.5
  for(int i=0; i<_chromosomeNb-1; ++i){
    if(SimRunner::r.Bool()) vecRecombs.push_back(_chromosomeSize[i]);
  }

  // add the end of the chomosome as a recombination event (easier for the inheritance function)
  vecRecombs.push_back(_chromosomeLength);

  // sort the vector
  sort(vecRecombs.begin(), vecRecombs.end());
}


//------------------------------------------------------------------------------
/** same as the copy constructor */
void
TTrait::_copyTTraitParameters(const TTrait& T){
  if(sequence) delete[] sequence; sequence = NULL;
  pTraitProto   = T.pTraitProto;
}

//------------------------------------------------------------------------------


