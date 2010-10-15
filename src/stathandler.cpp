/** @file stathandler.cpp
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


#include "stathandler.h"
#include "ttrait.h"
#include "metapop.h"
#include "metapop_sh.h"
#include "output.h"

//------------------------------------------------------------------------------
/**Adds a StatRecorder to the list, it is also added to the StatHandlerBase::_stats list.
		Two types of function variables are passed to this function. The "getter" and "setter".
		A "getter" returns a double value that will be stored in the StatRecorder structure. It
		may or may not take an argument. Only one getter should be passed to the stat recorder.
		The setter is used to set variables in the SH class which are then read by the getter.
		A setter and a getter may be given together but a setter alone will issue an error at
		runtime as no getter is present in the stat recorder.
   * @param Title 						the stat title
   * @param Name 							the stat name (headers in the output file)
   * @param Order 						stat table ordering flag
   * @param AGE 							age on which the stat should be processed
   * @param ARG 							the argument to pass to the S function
   * @param getSt 						function ptr to a S getter
   * @param getStBoolArg 			function ptr to a S getter with boolean argument
   * @param getStUintArg 			function ptr to a S getter with unsigned int argument
   * @param getStAGE 					function ptr to a S getter with age argument
   * @param getStBoolArgAGE 	function ptr to a S getter with boolean and age argument
   * @param getStUintArgAGE		function ptr to a S getter with unsigned int and age argument
   */
template <class SH>
void StatHandler<SH>::add_pairwisePatch (string Title, string Name, st_order Order, age_t AGE, unsigned int ARG,
					double(SH::* getStat)(void),
          double(SH::* getStatBoolArg)(bool),
					double(SH::* getStatUintArg)(unsigned int),
					double(SH::* getStatAGE)(const age_t&),
          double(SH::* getStatBoolArgAGE)(bool, const age_t&),
					double(SH::* getStatUintArgAGE)(unsigned int, const age_t&))
{
	assert(!ARG && (getStatUintArg || getStatUintArgAGE));
	if(_sample_pops_size==1) return add(Title,Name,Order,AGE,ARG,getStat,getStatBoolArg,getStatUintArg,
                                      getStatAGE,getStatBoolArgAGE,getStatUintArgAGE);

	// several patches are available
	string t1, t2, n1, n2;
	unsigned int i, j, id1, id2;
	for(i=0; i<_sample_pops_size; ++i){
    id1 = _sample_pops[i]->get_ID()+1;
    n1 = STRING::int2str_(id1, _pop->getPatchNbr(), "p") + "-";      		 // "_p1-2"   first patch (name)
    t1 = " between patch " + STRING::int2str(id1, (unsigned int)1) + " and ";   // first patch (description)
		for(j=i+1; j<_sample_pops_size; ++j){
      id2 = _sample_pops[j]->get_ID()+1;
			n2  = STRING::int2str(id2, _pop->getPatchNbr());                   //           second patch (name)
			t2  = STRING::int2str(id2);                              //           second patch (description)
			add(Title+t1+t2,Name+n1+n2,Order,AGE,i+j*_sample_pops_size,getStat,getStatBoolArg,getStatUintArg,
          getStatAGE,getStatBoolArgAGE,getStatUintArgAGE);
  	}
	}
}

//------------------------------------------------------------------------------
template <class SH>
void StatHandler<SH>::add_perPatch (string Title, string Name, st_order Order, age_t AGE, unsigned int ARG,
					double(SH::* getStat)(void),
          double(SH::* getStatBoolArg)(bool),
					double(SH::* getStatUintArg)(unsigned int),
					double(SH::* getStatAGE)(const age_t&),
          double(SH::* getStatBoolArgAGE)(bool, const age_t&),
					double(SH::* getStatUintArgAGE)(unsigned int, const age_t&))
{
	assert(!ARG && (getStatUintArg || getStatUintArgAGE));
	if(_sample_pops_size==1) return add(Title,Name,Order,AGE,ARG,getStat,getStatBoolArg,getStatUintArg,
                                      getStatAGE,getStatBoolArgAGE,getStatUintArgAGE);

	// several patches are available
	string t, n;
  unsigned int id;
	for(unsigned int i=0; i<_sample_pops_size; ++i){
    id = _sample_pops[i]->get_ID()+1;
		n = STRING::int2str_(id, _pop->getPatchNbr(), "p");     // name: "_p1"
		t = " of population " + STRING::int2str(id);  // description
		add(Title+t,Name+n,Order,AGE,i,getStat,getStatBoolArg,getStatUintArg,
        getStatAGE,getStatBoolArgAGE,getStatUintArgAGE);
	}
}

//------------------------------------------------------------------------------
template <class SH>
void StatHandler<SH>::add_perLocus (string Title, string Name, st_order Order, age_t AGE, unsigned int ARG,
					double(SH::* getStat)(void),
          double(SH::* getStatBoolArg)(bool),
					double(SH::* getStatUintArg)(unsigned int),
					double(SH::* getStatAGE)(const age_t&),
          double(SH::* getStatBoolArgAGE)(bool, const age_t&),
					double(SH::* getStatUintArgAGE)(unsigned int, const age_t&))
{
	assert(!ARG && (getStatUintArg || getStatUintArgAGE));
	if(_nb_locus==1) return add(Title,Name,Order,AGE,ARG,getStat,getStatBoolArg,getStatUintArg,
                              getStatAGE,getStatBoolArgAGE,getStatUintArgAGE);

	// several patches are available
	string t, n;
	for(unsigned int i=0; i<_nb_locus; ++i){
		n = STRING::int2str_(i+1, _nb_locus, "l");     // name: "_l1"
		t = " of locus " + STRING::int2str(i+1);  // description
		add(Title+t,Name+n,Order,AGE,i,getStat,getStatBoolArg,getStatUintArg,
        getStatAGE,getStatBoolArgAGE,getStatUintArgAGE);
	}
}

//------------------------------------------------------------------------------
template <class SH>
void StatHandler<SH>::add (string Title, string Name, st_order Order, age_t AGE, unsigned int ARG,
					double(SH::* getStat)(void),
          double(SH::* getStatBoolArg)(bool),
					double(SH::* getStatUintArg)(unsigned int),
					double(SH::* getStatAGE)(const age_t&),
          double(SH::* getStatBoolArgAGE)(bool, const age_t&),
					double(SH::* getStatUintArgAGE)(unsigned int, const age_t&))
{
	StatRecorder<SH>* new_rec = new StatRecorder<SH>(this->get_nrows(),this->get_ncols(), this->get_index());
	if(_trait_index) Name += STRING::int2str_(_trait_index, (int)0, "t");   // "_t2"
	new_rec->set(Title,Name,Order,AGE,ARG,getStat,getStatBoolArg,getStatUintArg,
               getStatAGE,getStatBoolArgAGE,getStatUintArgAGE);
	_recorders.push_back(new_rec);
	StatHandlerBase::add(new_rec);

	#ifdef _DEBUG
	 if(_recorders.size()>1) message(" ");
	 message("%s", Name.c_str());
	#endif
}

//------------------------------------------------------------------------------
template <class SH>
StatHandler<SH>::~StatHandler ( ) {
  if(_hsnei_locus)   delete[] _hsnei_locus;
  if(_htnei_locus)   delete[] _htnei_locus;
  if(_fst_locus)     delete[] _fst_locus;
  if(_fis_locus)     delete[] _fis_locus;
  if(_fit_locus)     delete[] _fit_locus;
  if(_fst_WC_locus)  delete[] _fst_WC_locus;
  if(_fis_WC_locus)  delete[] _fis_WC_locus;
  if(_fit_WC_locus)  delete[] _fit_WC_locus;
  if(_ho_locus)      delete[] _ho_locus;
  if(_hs_locus)      delete[] _hs_locus;
  if(_ht_locus)      delete[] _ht_locus;
  if(_nb_allele_of_pop)delete[] _nb_allele_of_pop;

  ARRAY::delete_3D(_alleleCountTable, _sample_pops_size, _nb_locus);

  ARRAY::delete_2D(_fst_matrix_wc, _sample_pops_size);
  ARRAY::delete_2D(_fst_matrix, _sample_pops_size);
  ARRAY::delete_2D(_coa_matrix, _sample_pops_size);
  ARRAY::delete_3D(_alleleFreqTable, _sample_pops_size, _nb_locus);
  ARRAY::delete_2D(_globalAlleleFreq, _nb_locus);

  REC_IT pos = _recorders.begin();
  for(;pos != _recorders.end(); ++pos){
    delete *pos;
  }
  _recorders.clear();

  ARRAY::delete_2D(_computed, _computed_size);
}

//------------------------------------------------------------------------------
template <class SH>
void StatHandler<SH>::execute ( ) {
	(this->*stat_save)();
}

//------------------------------------------------------------------------------
/** stats are stored in the db */
template <class SH>
void StatHandler<SH>::stats_to_db ( )
{
	for(REC_IT rec = _recorders.begin(); rec != _recorders.end(); rec++) {
		if( !((*rec)->setVal(_pop->getCurrentAge(),_pop->getCurrentGeneration(),
												 _pop->getCurrentReplicate(),dynamic_cast<SH*>(this))) ) {
			fatal("StatHandler::stats_to_db returns false\n");
		}
	}
}

//------------------------------------------------------------------------------
/** stats are directly writen to the file */
template <class SH>
void StatHandler<SH>::stats_to_file ( )
{
	// open the file each time
	ofstream FH(get_service()->get_file_name_stats().c_str(), ios::app);
	if(!FH) fatal("Could not open stat output file \"%s\"!\n",get_service()->get_file_name_stats().c_str());

	FH.width(12);
	FH.setf(ios::left,ios::adjustfield);
	FH << setprecision(4);
	FH << get_service()->get_pop_ptr()->getCurrentReplicate() << "\t";
	FH.width(12);
	FH << get_service()->get_pop_ptr()->getCurrentGeneration() << "\t";

	for(REC_IT rec = _recorders.begin(); rec != _recorders.end(); rec++) {
		if( !((*rec)->setValDirectly(_pop->getCurrentAge(),_pop->getCurrentGeneration(),
												 _pop->getCurrentReplicate(),dynamic_cast<SH*>(this), FH)) ) {
			fatal("StatHandler::stats_to_file returns false\n");
		}
	}

	FH << "\n";
	FH.close();    // close the file
}

//------------------------------------------------------------------------------
template <class SH>
void StatHandler<SH>::set(TTraitProto* TT){
	_trait_index = TT->get_trait_index();
	_SHLinkedTrait = TT;
  _SHLinkedTraitIndex = TT->get_absolute_index();
	_pop = TT->get_popPtr();
  _nb_locus  = _SHLinkedTrait->get_nb_locus();
	_nb_allele = _SHLinkedTrait->get_nb_allele();
  _ploidy    = _SHLinkedTrait->get_ploidy();
}

//------------------------------------------------------------------------------
template <class SH>
bool StatHandler<SH>::init(){
	StatHandlerBase::init();

	// only if each stat is outputted solely the file is directly written
	if(get_service()->get_save_choice() == 6) stat_save = &StatHandler::stats_to_file;
	else																			stat_save = &StatHandler::stats_to_db;

	return true;
}

//------------------------------------------------------------------------------
	/** check if it has already been computed: return true if it was already computed
  * the array has to have a length of three
  * if the age is not used it is set to NONE (0)
  */
template <class SH>
bool StatHandler<SH>::already_computed(unsigned int* array, age_t AGE){
	if(   array[0] == _current_generation
     && array[1] == _current_replicate
   	 && array[2] == AGE){
  	return true;
	}

	array[0] = _current_generation;
	array[1] = _current_replicate;
	array[2] = AGE;
	return false;
}

// ----------------------------------------------------------------------------------------
//////////////////////// FSTAT functions //////////////////////////////////////////////////
// ----------------------------------------------------------------------------------------
// allocateTable
// ----------------------------------------------------------------------------------------
/** tables are only available for sampled pops */
template <class SH>
void StatHandler<SH>::allocateTable ()
{
  ARRAY::create_3D(_alleleCountTable, _sample_pops_size, _nb_locus, _nb_allele);
  ARRAY::create_3D(_alleleFreqTable,  _sample_pops_size, _nb_locus, _nb_allele);
  ARRAY::create_2D(_globalAlleleFreq, _nb_locus, _nb_allele);
}

// ----------------------------------------------------------------------------------------
// get_allele_freq
// ----------------------------------------------------------------------------------------
/** function to retrieve the local allele frequencies with only a single index
  * Caution: only the sampled pops are available !
  */
template <class SH>
double StatHandler<SH>::get_allele_freq_local(unsigned int i, const age_t& AGE){

  set_alleleCountTable_local(AGE);

  unsigned int size_p = _nb_locus * _nb_allele;

  unsigned int p = i/size_p;
  unsigned int l = (i-(p*size_p))/_nb_allele;
  unsigned int a = i-l*_nb_allele-p*size_p;

  return _alleleFreqTable[p][l][a];

}

// ----------------------------------------------------------------------------------------
// get_allele_freq_global
// ----------------------------------------------------------------------------------------
/** function to retrieve the global allele frequencies with only a single index
  * Caution: only the sampled pops are available
  */
template <class SH>
double StatHandler<SH>::get_allele_freq_global(unsigned int i, const age_t& AGE){

  set_alleleCountTable_global(AGE);

  unsigned int l = i/_nb_allele;
  unsigned int a = i-l*_nb_allele;

  return _globalAlleleFreq[l][a];

}

// ----------------------------------------------------------------------------------------
// set_alleleCountTable_local (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the allele frequencies for each patch separately
  * Caution: only the sampled pops are available
  */
template <class SH>
void StatHandler<SH>::set_alleleCountTable_local(age_t AGE)
{
  // check if the table has already been  created
  if(already_computed(_computed[4], AGE)) return;

  // get the local frequencies and counts for each patch
  for(unsigned int k = 0; k < _sample_pops_size; ++k) {
    set_alleleCountTable_local_ofPatch(_sample_pops[k], AGE, _alleleFreqTable[k], _alleleCountTable[k]);
  }
}

// ----------------------------------------------------------------------------------------
// set_alleleCountTable_local (per patch)
// ----------------------------------------------------------------------------------------
/** computation of the allele frequencies for each patch separately
  * Caution: only the sampled pops are available
  */
template <class SH>
void StatHandler<SH>::set_alleleCountTable_local_ofPatch(Patch* crnt_patch, age_t AGE,
                                              double**& freqs, unsigned int**& counts)
{
  if(!counts) ARRAY::create_2D(counts, _nb_locus, _nb_allele, (unsigned int) 0);
  else        ARRAY::reset_2D (counts, _nb_locus, _nb_allele, (unsigned int) 0);
  if(!freqs)  ARRAY::create_2D(freqs,  _nb_locus, _nb_allele);

  unsigned int i, l, a, p, sizeF, sizeM, size;
  unsigned char** genes;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

  sizeF = crnt_patch->size(FEM, age_pos);    // number of females
  sizeM = crnt_patch->size(MAL, age_pos);    // number of males
  size = (sizeF + sizeM) * _ploidy;          // total number of alleles

 // if patch is empty all frequencies are zero
 if(!size){
   for (l = 0; l < _nb_locus; ++l){
     for (a = 0; a < _nb_allele; ++a){
       freqs[l][a] = 0;        // counts are already set to zero
     }
   }
   return;
 }

 // get female allele counts
 for(i = 0; i < sizeF; ++i) {
   genes = (unsigned char**)crnt_patch->get(FEM, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
   for (l = 0; l < _nb_locus; ++l){
     for (p = 0; p < _ploidy; ++p){
       ++counts[l][genes[l][p]];
     }
   }
 }

 // get male allele counts
 for(i = 0; i < sizeM; ++i) {
   genes = (unsigned char**)crnt_patch->get(MAL, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
   for (l = 0; l < _nb_locus; ++l){
     for (p = 0; p < _ploidy; ++p){
       ++counts[l][genes[l][p]];
     }
   }
 }

 // compute allele frequencies
 for (l = 0; l < _nb_locus; ++l){
   for (a = 0; a < _nb_allele; ++a){
     freqs[l][a] = (double)counts[l][a] / size;
   }
 }
}

// ----------------------------------------------------------------------------------------
// set_alleleCountTable_global (globally)
// ----------------------------------------------------------------------------------------
/** computation of the allele frequencies in the entire metapopulation
  * Caution: only the sampled pops are available
  */
template <class SH>
void StatHandler<SH>::set_alleleCountTable_global(age_t AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[5], AGE)) return;

  // get the local allele frequencies if necessary
  set_alleleCountTable_local(AGE);

  set_alleleCountTable_global(_sample_pops_size, _current_nbInds[AGE==ADULTS ? ADLTx : OFFSx],  _alleleCountTable, _globalAlleleFreq);
}

// ----------------------------------------------------------------------------------------
// set_alleleCountTable_global (globally)
// ----------------------------------------------------------------------------------------
/** computation of the allele frequencies in the entire metapopulation
  * Caution: only the sampled pops are available
  */
template <class SH>
void StatHandler<SH>::set_alleleCountTable_global(unsigned int nbPatch, const unsigned int& nbInd,
                                                  unsigned int*** alleleCounts, double** globalFreqs)
{
  unsigned int p, l, a, count;

  for (l = 0; l < _nb_locus; ++l){
    for (a = 0; a < _nb_allele; ++a){
      // sum up the allele counts of a patch
      count = 0;
      for(p = 0; p < nbPatch; ++p){
        count += alleleCounts[p][l][a];
      }
      globalFreqs[l][a] = (double)count/(2*nbInd);
    }
  }
}

// ----------------------------------------------------------------------------------------
// setFstat
// ----------------------------------------------------------------------------------------
/** compute global F-satistics averaged across loci and pops following Nei & Chesser (1983)
  * Hs  = Mean_patch(Mean_loci(1-SUM_allele(p^2))  // p: allele freq of locus l and patch p
        NaN if metapop empty
  * Ht  = Mean_loci(1-SUM_allele(mean(p)^2))       // p: allele freq of locus l
        NaN if metapop empty
  * Fst = 1-Hs/Ht
				NaN if metapop empty, 0 if Ht=0 (all individuals are monomorph for a SINGLE allele)
  * Fis = Ho/Hs
				NaN if metapop empty, 0 if Hs=0 (all individual of a patch are monomorph for an allele (may differ among patches))
	* Fit = Ho/Ht
        NaN if metapop empty, 0 if Ht=0 (all individuals are monomorph for a SINGLE allele)
  */
template <class SH>
void StatHandler<SH>::setFstat_Nei_Chesser(age_t AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[6], AGE)) return;

	age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
	unsigned int nbpatch = _current_nbPops[age_pos];

  // if metapop is empty
	if(!nbpatch){
    _hsnei = _htnei = _fis = _fit = _fst = my_NAN;
		return;
  }

	// compute harmonic mean of patch sizes:
	double harmonic = getHarmonicMean_ofPopSize(AGE, _sample_pops, _sample_pops_size);

	getHo(AGE);       // mean observed heterozygostiy
	getHs(AGE);       // mean expected heterozyogsity

	// Fis
	_hsnei = harmonic!=1 ? harmonic/(harmonic-1.0)*(_hs-(_ho/(2.0*harmonic))) : _hs;   //Nei's corrections:
	_fis = _hsnei ? 1.0-(_ho/_hsnei)    : 0;   // monomorphic = ht=0: after definition 0 and not NaN!

	if(nbpatch>1){                              // more than one colonized patch
		getHt(AGE);
		_htnei = _ht + (_hsnei/(harmonic*nbpatch))-(_ho/(2.0*harmonic*nbpatch));  //Nei's corrections:
		_fit = _htnei ? 1.0-(_ho/_htnei)    : 0;   // monomorphic = ht=0: after definition 0 and not NaN!
		_fst = _htnei ? 1.0-(_hsnei/_htnei) : 0;   // monomorphic = ht=0: after definition 0 and not NaN!
	}
	else _ht =_htnei = _fit = _fst = my_NAN;
}

// ----------------------------------------------------------------------------------------
// setFstatMatrix
// ----------------------------------------------------------------------------------------
/** computes the pairwise Fst values and stores them in a matrix
  * Fst computed following Nei & Chesser (1983)
  */
template <class SH>
void StatHandler<SH>::setFstat_Nei_Chesser_perPatchPair(age_t AGE)
{

  // check if the table has already been computed
  if(already_computed(_computed[7], AGE)) return;

  double harmonic, hsnei, htnei, ho, ht, hs;
  unsigned int nbpatch, i, j;
  unsigned int aPatchID[2];                         // array for the two patch id's
  Patch* aPatch[2];                                 // array for the two patches

  if(!_fst_matrix) ARRAY::create_2D(_fst_matrix, _sample_pops_size, _sample_pops_size);

  double * hs_pop = new double[_sample_pops_size];
  double * ho_pop = new double[_sample_pops_size];
  for(i = 0; i< _sample_pops_size; ++i){
    ho_pop[i] = getHo_ofPatch (i, AGE);
    hs_pop[i] = getHs_ofPatch (i, AGE);

  }
  // for each pair of pops
  // we have to check all pops, since the Fst after Nei is also defined for a
  // single populated pop, but not if both pops are empty
  for(i = 0; i< _sample_pops_size-1; ++i){
    aPatchID[0] = i;                                  // get the corresponding patch id
    aPatch[0] = _sample_pops[i];                      // get the corresponding patch
    for(j = i+1; j< _sample_pops_size; ++j){
      aPatchID[1] = j;                                // get the corresponding patch id
      aPatch[1] = _sample_pops[j];                    // get the corresponding patch

      // if neither of the patches is populated
      if(hs_pop[i] == my_NAN && hs_pop[j] == my_NAN){
        _fst_matrix[i][j] =  my_NAN;
        continue;
      }

      harmonic = getHarmonicMean_ofPopSize(AGE, aPatch, 2);

      // compute hs, ho and ht
      if(hs_pop[i] != my_NAN){
        if(hs_pop[j] != my_NAN){         // both patches are populated
          hs = (hs_pop[i] + hs_pop[j])/2.0;
          ho = (ho_pop[i] + ho_pop[j])/2.0;
          ht = getHt(AGE, aPatchID, 2);
          nbpatch = 2;
        }
        else{                            // only patch 1 is populated
          hs = hs_pop[i];
          ho = ho_pop[i];
          ht = hs;
          nbpatch = 1;
        }
      }
      else{                              // only patch 2 is populated
        hs = hs_pop[j];
        ho = ho_pop[j];
        ht = hs;
        nbpatch = 1;
      }

      //Nei's corrections:
      hsnei = harmonic!=1 ? harmonic/(harmonic-1.0)*(hs-(ho/(2.0*harmonic))) : hs;
      htnei = ht + (hsnei/(harmonic*nbpatch))-(ho/(2.0*harmonic*nbpatch));

			_fst_matrix[i][j] = htnei ? 1.0-(hsnei/htnei) : 0;  // monomorphic = ht=0: after definition 0 and not NaN!
		}
  }
  delete[] hs_pop;
  delete[] ho_pop;
}

// ----------------------------------------------------------------------------------------
// setFstat
// ----------------------------------------------------------------------------------------
// following Nei & Chesser (1983)
template <class SH>
void StatHandler<SH>::setFstat_Nei_Chesser_perLocus(age_t AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[8], AGE)) return;

  unsigned int i, nbpatch = get_current_nbPops(AGE);

  if(!_hsnei_locus){
    _hsnei_locus = new double[_nb_locus];
    _htnei_locus = new double[_nb_locus];
    _fis_locus   = new double[_nb_locus];
    _fit_locus   = new double[_nb_locus];
    _fst_locus   = new double[_nb_locus];
    _hs_locus    = new double[_nb_locus];
    _ht_locus    = new double[_nb_locus];
  }

  // if the patches are empty
	if(!nbpatch){
    for (i = 0; i < _nb_locus; ++i){
      _hsnei_locus[i] = _htnei_locus[i] = _fis_locus[i] = _fit_locus[i] = _fst_locus[i] = my_NAN;
    }
    return;
  }

  // compute harmonic mean of patch sizes:
  double harmonic = getHarmonicMean_ofPopSize(AGE, _sample_pops, _sample_pops_size);

  getHo_perLocus(AGE);  // _ho: computed for all loci at a time

  for (i = 0; i < _nb_locus; ++i){
		_hs_locus[i] = getHs_ofLocus(AGE, i);  // compute Hs
		_hsnei_locus[i] = harmonic != 1 ? harmonic/(harmonic-1.0)*(_hs_locus[i]-(_ho_locus[i]/(2.0*harmonic))) : _hs_locus[i];
		_fis_locus[i] = _hsnei_locus[i] ? 1.0-(_ho_locus[i]/_hsnei_locus[i]) : 0; // monomorphic = ht=0: after definition 0 and not NaN!

		if(nbpatch>1){       // more than one colonized patch
			_ht_locus[i] = getHt_ofLocus(AGE, i);  // compute Ht
			_htnei_locus[i] = _ht_locus[i] + (_hsnei_locus[i]/(harmonic*nbpatch))-(_ho_locus[i]/(2.0*harmonic*nbpatch));
			_fit_locus[i] = _htnei_locus[i] ? 1.0-(_ho_locus[i]/_htnei_locus[i])    : 0; // monomorphic = ht=0: after definition 0 and not NaN!
			_fst_locus[i] = _htnei_locus[i] ? 1.0-(_hsnei_locus[i]/_htnei_locus[i]) : 0; // monomorphic = ht=0: after definition 0 and not NaN!
		}
		else _ht_locus[i]=_htnei_locus[i]=_fit_locus[i]=_fst_locus[i] = my_NAN;
	}
}

// ----------------------------------------------------------------------------------------
// setLociDivCounter
// ----------------------------------------------------------------------------------------
/** compute the mean number of alleles per locus of the:
  *   - entire metapop: _nb_all_global
  *   - array of the mean number of alleles across loci of pop i: _nb_allele_of_pop[pop]
  */
template <class SH>
void StatHandler<SH>::setLociDivCounter (age_t AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[9], AGE)) return;

  if(!_nb_allele_of_pop) _nb_allele_of_pop = new double[_sample_pops_size];

  unsigned int nbpatch = get_current_nbPops(AGE);
  if(!nbpatch){
    _nb_all_global = my_NAN;
    ARRAY::reset_1D(_nb_allele_of_pop, _sample_pops_size, (double)my_NAN);
    return;
  }

  unsigned int i, j, k, val;
  double patch_mean;
  bool **pop_div = ARRAY::new_2D(_nb_locus, _nb_allele, (bool) 0);

  // get the mean number of alleles per locus and patch
  for(k = 0; k < _sample_pops_size; ++k) {
    if(!_sample_pops[k]->size(AGE)) continue;   // if the pop is empty skip it
    patch_mean = 0;
    for(i = 0; i < _nb_locus; ++i){
	    for(j = 0; j < _nb_allele; ++j) {
        val = _alleleCountTable[k][i][j];
		    patch_mean += (val != 0);
        pop_div[i][j] |= (val != 0);
      }
    }
    //add mean nb of alleles per locus for Patch k to the pop mean
    _nb_allele_of_pop[k] = patch_mean/_nb_locus;
  }

  // get the mean number of alleles per locus in the entire metapopulation
  _nb_all_global = 0;
  for(i = 0; i < _nb_locus; ++i){
    for(j = 0; j < _nb_allele; ++j){
      _nb_all_global += pop_div[i][j];
    }
  }
  _nb_all_global /= _nb_locus;

  ARRAY::delete_2D(pop_div, _nb_locus);
}

// ----------------------------------------------------------------------------------------
// get_fixed_loci_local
// ----------------------------------------------------------------------------------------
/** get the mean number of fixed loci per population */
template <class SH>
double StatHandler<SH>::get_fixed_loci_local(age_t AGE)
{
  // get the local allele frequencies if necessary
  set_alleleCountTable_local(AGE);

  //number of fixed loci, local and global counters:
   double fix_loc_local = 0;

  for (unsigned int p = 0; p < _sample_pops_size; ++p){
    fix_loc_local += get_fixed_loci_ofPatch(AGE, p);
  }
  return (double)fix_loc_local/get_current_nbPops(AGE);
}

// ----------------------------------------------------------------------------------------
// get_fixed_loci_local
// ----------------------------------------------------------------------------------------
/** get the mean number of fixed loci per population */
template <class SH>
double StatHandler<SH>::get_fixed_loci_ofPatch(age_t AGE, unsigned int p)
{
  // get the local allele frequencies if necessary
  set_alleleCountTable_local(AGE);

  //number of fixed loci, local and global counters:
  unsigned int l, a, fix_loc_local = 0;

  for (l = 0; l < _nb_locus; ++l){
    for (a = 0; a < _nb_allele; ++a){
      if(_alleleFreqTable[p][l][a] == 1) ++fix_loc_local;
    }
  }
  return (double)fix_loc_local;
}

// ----------------------------------------------------------------------------------------
// get_fixed_loci_global
// ----------------------------------------------------------------------------------------
/** get the number of fixed loci in the entire metapopulation */
template <class SH>
unsigned int StatHandler<SH>::get_fixed_loci_global(age_t AGE)
{
  // get the global allele frequencies if necessary
  set_alleleCountTable_global(AGE);

  unsigned int l, a, fix_loc_global = 0;
  for (l = 0; l < _nb_locus; ++l){
    for (a = 0; a < _nb_allele; ++a){
      if(_globalAlleleFreq[l][a] == 1) ++fix_loc_global;
    }
  }
  return fix_loc_global;
}

// ----------------------------------------------------------------------------------------
// getHarmonicMean_ofPopSize
// ----------------------------------------------------------------------------------------
// compute harmonic mean of patch sizes:
template <class SH>
double StatHandler<SH>::getHarmonicMean_ofPopSize(const age_t& AGE, Patch** aPatch, unsigned int nbPatch)
{
  double harmonic = 0;
  unsigned int nbind, p, nbPopFull=0;
  for(p=0; p<nbPatch; ++p){
    nbind = aPatch[p]->size(AGE);
	  if(nbind){
      harmonic += 1.0/nbind;
      ++nbPopFull;
    }
  }
  return (double) nbPopFull / harmonic;
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get mean observed heterozygosity */
template <class SH>
double StatHandler<SH>::getHo (const age_t& AGE)
{
  // compute the observed heterozygosity for each locus
  getHo_perLocus (AGE);

  //  if the metapop is empty
  if(_ho_locus[0] == my_NAN){
    _ho = my_NAN;
    return _ho;
  }

  // compute the mean observed heterozygosity
  _ho = 0;
	for(unsigned int k = 0; k < _nb_locus; ++k) {
    _ho += _ho_locus[k];
  }

  _ho /= _nb_locus;
  return _ho;
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get observed heterozygoxity for each loci individually */
template <class SH>
double* StatHandler<SH>::getHo_perLocus (const age_t& AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[10], AGE)) return _ho_locus;

  unsigned int l, nbPatch = get_current_nbPops(AGE);
  double* ho = new double[_nb_locus];                 // temp array
  if(_ho_locus) ARRAY::reset_1D(_ho_locus, _nb_locus, (double)0);
  else          ARRAY::create_1D(_ho_locus, _nb_locus, (double)0);

  if(!nbPatch){      // if metapop is empty
    ARRAY::reset_1D(_ho_locus, _nb_locus, (double)my_NAN);
    return _ho_locus;
  }

  // compute Ho per patch
  Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
  for(; cur != end; ++cur){
    getHo_ofPatchperLocus (AGE, *cur, ho);     // compute ho for each locus of the patch i
    if(ho[0]!=my_NAN){                         // check if patch is empty
      for (l = 0; l < _nb_locus; ++l) {
        _ho_locus[l] += ho[l];                 // sum up the ho of each patch
      }
    }
  }

  for (l = 0; l < _nb_locus; ++l) {
    _ho_locus[l] /= nbPatch;            // compute mean across patches
  }

  delete[] ho;

  return _ho_locus;
}
// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get observed heterozygosity of patch i */
template <class SH>
double StatHandler<SH>::getHo_ofPatch (unsigned int i, const age_t& AGE)
{
  Patch* cur_patch = _sample_pops[i];

  // compute the observed heterozygosity for each locus of that patch
  double* ho_locus = getHo_ofPatchperLocus (AGE, cur_patch);

  //  if the metapop is empty
  if(ho_locus[0] == my_NAN){
    delete[] ho_locus;
    return my_NAN;
  }
  
  // compute the mean observed heterozygosity
  double ho = 0;
	for(unsigned int l = 0; l < _nb_locus; ++l) {
    ho += ho_locus[l];
  }

  delete[] ho_locus;

  return ho /= _nb_locus;
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get observed heterozygoxity for each locus of the selected patch individually
  * if no array is passed then the returned array has to be deleted !
  */
template <class SH>
double* StatHandler<SH>::getHo_ofPatchperLocus (const age_t& AGE, Patch* cur_patch, double* array)
{
  if(!array) ARRAY::create_1D(array, _nb_locus, (double)0); // create the array if not passed
  else       ARRAY::reset_1D( array, _nb_locus, (double)0);

  unsigned int i, l;
  unsigned char** genes;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  unsigned int sizeF = cur_patch->size(FEM, age_pos),
               sizeM = cur_patch->size(MAL, age_pos);

  // if patch is empty
  if(!sizeF && !sizeM){                     // if patch is empty
    ARRAY::reset_1D(array, _nb_locus, (double)my_NAN);
    return array;
  }


  // count heterozygote females
	for(i = 0; i < sizeF; ++i) {
	  genes = (unsigned char**)cur_patch->get(FEM, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
	  for (l = 0; l < _nb_locus; ++l){
      array[l] += (genes[l][0] != genes[l][1]);
    }
	}

  // count heterozygote males
	for(i = 0; i < sizeM; ++i) {
	  genes = (unsigned char**)cur_patch->get(MAL, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
	  for (l = 0; l < _nb_locus; ++l){
      array[l] += (genes[l][0] != genes[l][1]);
    }
	}

  // compute the mean heterozygosity for each locus
  for (l = 0; l < _nb_locus; ++l) {
    array[l] /= sizeM+sizeF;
  }

  return array;
}

// ----------------------------------------------------------------------------------------
/** get observed heterozygosity of the patch for each allele separately array[locus][allele]
  * if no array is passed then the returned array has to be deleted !
  */
template <class SH>
double** StatHandler<SH>::getHo_ofPatchperAllele (const age_t& AGE, Patch* cur_patch, double** array)
{
  if(!array) ARRAY::create_2D(array, _nb_locus, _nb_allele, (double)0); // create the array if not passed
  else       ARRAY::reset_2D( array, _nb_locus, _nb_allele, (double)0);

  unsigned int i, l;
  unsigned char** genes;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);
  unsigned int sizeF = cur_patch->size(FEM, age_pos),
               sizeM = cur_patch->size(MAL, age_pos);

  // if patch is empty
  if(!sizeF && !sizeM){                     // if patch is empty
    ARRAY::reset_2D(array, _nb_locus, _nb_allele, (double)my_NAN);
    return array;
  }

  // count heterozygote females
	for(i = 0; i < sizeF; ++i) {
	  genes = (unsigned char**)cur_patch->get(FEM, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
	  for (l = 0; l < _nb_locus; ++l){
      if(genes[l][0] != genes[l][1]){
        ++array[l][genes[l][0]];
        ++array[l][genes[l][1]];
      }
    }
	}

  // count heterozygote males
	for(i = 0; i < sizeM; ++i) {
	  genes = (unsigned char**)cur_patch->get(MAL, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence();
	  for (l = 0; l < _nb_locus; ++l){
      if(genes[l][0] != genes[l][1]){
        ++array[l][genes[l][0]];
        ++array[l][genes[l][1]];
      }
    }
	}

  // compute the mean heterozygosity for each locus
  for (l = 0; l < _nb_locus; ++l) {
    for (i = 0; i < _nb_allele; ++i) {
      array[l][i] /= (sizeM+sizeF);
    }
  }

  return array;
}

// ----------------------------------------------------------------------------------------
// setHo
// ----------------------------------------------------------------------------------------
/** get observed heterozygoxity for each loci individually
  * Ho = mean(sum_patch(a1 != a2)/nb_ind_of_patch)
  * the returned array must be deleted!
  */
template <class SH>
double* StatHandler<SH>::getHo_perLocus (const age_t& AGE, const vector<Patch*>& vPatch)
{
  // check if the table has already been computed
  if(already_computed(_computed[11], AGE)) return _ho_locus;

  unsigned int size= vPatch.size(), sizeM, sizeF, j, i, k;
  int* ho = new int[_nb_locus];                 // temp parameter
  double* ho_locus = new double[_nb_locus];     // is returned and must be deleted in the calling function
  unsigned char** genes;
  Patch* current_patch;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

  // compute Ho per patch
  for(vector<Patch*>::iterator pos = vPatch.begin(); pos != vPatch.end(); ++pos) {

    // reset all ho to zero
    ARRAY::reset_1D(ho, _nb_locus, (int)0);

    // females
	  sizeF = (*pos)->size(FEM, age_pos);
	  for(j = 0; j < sizeF; ++j) {
	    genes = (unsigned char**)(*pos)->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
	    for (k = 0; k < _nb_locus; ++k){
        ho[k] += (genes[k][0] != genes[k][1]);
      }
	  }

    // males
    sizeM = (*pos)->size(MAL, age_pos);
	  for(j = 0; j < sizeM; ++j) {
	    genes = (unsigned char**)(*pos)->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence();
	    for (k = 0; k < _nb_locus; ++k){
        ho[k] += (genes[k][0] != genes[k][1]);
      }
	  }

    // sum up the mean ho per locus for each patch
    for (k = 0; k < _nb_locus; ++k) {
      if(sizeM+sizeF) _ho_locus[k] += (double)ho[k]/(sizeM+sizeF);
    }
  }

  // compute mean across patches
  for (k = 0; k < _nb_locus; ++k) {
    _ho_locus[k] /= vPatch.size();
  }

  delete[] ho;

  return ho_locus;
}

// ----------------------------------------------------------------------------------------
// getHs
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity (across loci and patches) hs = 1 - sum(p^2) */
template <class SH>
double StatHandler<SH>::getHs (const age_t& AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[12], AGE)) return _hs;

  // check if the metapop is empty
  if(!get_current_nbPops(AGE)){
    _hs = my_NAN;
    return _hs;
  }

  // compute mean expected heterozygosity
  _hs = 0;
	for (unsigned int l = 0; l < _nb_locus; ++l) {
    _hs += getHs_ofLocus(AGE, l);
  }
  _hs /= _nb_locus;

  return _hs;
}

// ----------------------------------------------------------------------------------------
// getHs_perPatch_perLocus
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity of locus l across patches: hs = 1 - sum(p^2) */
template <class SH>
double StatHandler<SH>::getHs_ofLocus (const age_t& AGE, const unsigned int& l)
{
  unsigned int nbPatch = get_current_nbPops(AGE);
  if(!nbPatch) return my_NAN;        // if metapop is empty

  double val, hs = 0;

	for (unsigned int p = 0; p < _sample_pops_size; ++p) {
    val = getHs_ofPatchandLocus(AGE, p, l);
    if(val != my_NAN) hs += val;
  }
  return hs/nbPatch;
}

// ----------------------------------------------------------------------------------------
// getHs_perPatch
// ----------------------------------------------------------------------------------------
/** return mean UNBIASED expected heterozygosity of patch p across loci (Nei 1987, eq. 7.39, p. 164) */
template <class SH>
double StatHandler<SH>::getHsUnbiased_ofPatch (unsigned int p, const age_t& AGE)
{
  // check if this pop is inhabited
  unsigned int size = _sample_pops[p]->size(AGE);
  if(size<2) return my_NAN;

  double hs = 0;
	for (unsigned int l = 0; l < _nb_locus; ++l) {
    hs += getHs_ofPatchandLocus(AGE, p, l);
  }
  return (double)size/(size-1)*(hs/_nb_locus - (getHo_ofPatch(p, AGE)/(2*size)));
}

// ----------------------------------------------------------------------------------------
// getHs_perPatch_perPatch
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity of patch p across loci: hs = 1 - sum(p^2) */
template <class SH>
double StatHandler<SH>::getHs_ofPatch (unsigned int p, const age_t& AGE)
{
  // check if this pop is inhabited
  if(!_sample_pops[p]->size(AGE)) return my_NAN;

  double hs = 0;
	for (unsigned int l = 0; l < _nb_locus; ++l) {
    hs += getHs_ofPatchandLocus(AGE, p, l);
  }
  return hs/_nb_locus;
}

// ----------------------------------------------------------------------------------------
// getHs_perPatchandLocus
// ----------------------------------------------------------------------------------------
/** return expected heterozygosity of patch p and locus l: hs = 1 - sum(p^2)
  * returns NaN if the pop is empty
  */
template <class SH>
double StatHandler<SH>::getHs_ofPatchandLocus (const age_t& AGE, const unsigned int& p, const unsigned int& l)
{
  // get the local allele frequencies if necessary
  set_alleleCountTable_local(AGE);

  double freq, hs = 1;
	for (unsigned int a = 0; a < _nb_allele; ++a) {
    freq = _alleleFreqTable[p][l][a];
	  if(freq) hs -= freq*freq;
	}
  if(hs==1) return my_NAN; // no alleles were present, respectively the pop is empty
  return hs;
}

// ----------------------------------------------------------------------------------------
// getHt
// ----------------------------------------------------------------------------------------
/** return mean total expected heterozygosity (across loci) ht = 1 - sum(p^2) */
template <class SH>
double StatHandler<SH>::getHt (const age_t& AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[13], AGE)) return _ht;

  _ht = 0;
	for (unsigned int l = 0; l < _nb_locus; ++l) {
    _ht += getHt_ofLocus(AGE, l);
  }
  _ht /= _nb_locus;
  return _ht;
}

// ----------------------------------------------------------------------------------------
// getHt
// ----------------------------------------------------------------------------------------
/** return mean total expected heterozygosity (across loci) ht = 1 - sum(mean(p)^2)
 *  across the given patches
 *  patchID is not the id of the patch but the id in the smaple_pops container
 */
template <class SH>
double StatHandler<SH>::getHt (const age_t& AGE, unsigned int* patchID, unsigned int nbPatch)
{
  double ht = 0;
	for (unsigned int l = 0; l < _nb_locus; ++l) {
    ht += getHt_ofLocus(AGE, l, patchID, nbPatch);
  }
  return ht /= _nb_locus;
}

// ----------------------------------------------------------------------------------------
// getHt_perPatch_perLocus
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity of locus l across all patches: hs = 1 - sum(mean(p)^2)
  * used by Nei
  */
template <class SH>
double StatHandler<SH>::getHt_ofLocus (const age_t& AGE, const unsigned int& l)
{
  // get the local allele frequencies if necessary
	set_alleleCountTable_global(AGE);

  double hs = 1, mean_freq;
  unsigned int a, p;
	for(a = 0; a < _nb_allele; ++a) {
		mean_freq = _globalAlleleFreq[l][a];
		hs  -= mean_freq*mean_freq;
	}
  return hs;
}

// ----------------------------------------------------------------------------------------
// getHt_perPatch_perLocus
// ----------------------------------------------------------------------------------------
/** return mean expected heterozygosity of locus l across patches: hs = 1 - sum(mean(p)^2)
  * used by Nei for the given patches
  * the array patchID containa the id of the patch in relation to sample_pops and not the id of the patch!
  */
template <class SH>
double StatHandler<SH>::getHt_ofLocus (const age_t& AGE, const unsigned int& l, unsigned int* patchID, unsigned int nbPatch)
{
  // get the local allele frequencies if necessary
  set_alleleCountTable_local(AGE);

  double hs = 1, mean_freq;
  unsigned int a;
  unsigned int* pos, *end = patchID + nbPatch;
	for(a = 0; a < _nb_allele; ++a) {
    mean_freq=0;
	  for(pos = patchID; pos != end; ++pos) {
      mean_freq += _alleleFreqTable[*pos][l][a];
    }
    mean_freq /= nbPatch;
	  hs  -= mean_freq*mean_freq;
	}
  return hs;
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham
// ----------------------------------------------------------------------------------------
/** computes the global F-statistics following Weir and Cockerham (1984) */
template <class SH>
void StatHandler<SH>::setFstat_Weir_Cockerham(age_t AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[14], AGE)) return;

  if(get_current_nbPops(AGE) < 2){
    _fst_wc = _fis_wc = _fit_wc = my_NAN;
  }

  set_alleleCountTable_global(AGE);

  setFstat_Weir_Cockerham(AGE, _sample_pops, _sample_pops_size, _alleleFreqTable, _globalAlleleFreq,
                             _fst_wc, _fis_wc, _fit_wc);
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham_perLocus
// ----------------------------------------------------------------------------------------
/** computes the global F-statistics following Weir and Cockerham (1984) */
template <class SH>
void StatHandler<SH>::setFstat_Weir_Cockerham_perLocus(age_t AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[15], AGE)) return;

  unsigned int l, nbPatch = get_current_nbPops(AGE);

  if(!_fis_WC_locus){
    _fis_WC_locus   = new double[_nb_locus];
    _fit_WC_locus   = new double[_nb_locus];
    _fst_WC_locus   = new double[_nb_locus];
   }

  // if the patches are empty
  if(!nbPatch){
    for (l = 0; l < _nb_locus; ++l){
      _fis_WC_locus[l] = _fit_WC_locus[l] = _fst_WC_locus[l] = my_NAN;
    }
    return;
  }

  set_alleleCountTable_global(AGE);

  unsigned int* pop_sizes   = new unsigned int[_sample_pops_size];
  double*** ho_patch_allele = new double**[_sample_pops_size];
  unsigned int tot_size=0, cur_size, p, nbPop = 0;
  double a2, b2, w2, tot_square_size=0, nc;

  // get pop sizes
  for(p = 0; p < _sample_pops_size; ++p) {
    pop_sizes[p] = cur_size = _sample_pops[p]->size(AGE);
    if(!cur_size){
      ho_patch_allele[p] = NULL;
      continue;       // if pop is empty
    }
    ++nbPop;                      // get the number of populated pops
    tot_size += cur_size;
    ho_patch_allele[p] = getHo_ofPatchperAllele (AGE, _sample_pops[p]);
    tot_square_size += cur_size*cur_size;
  }

  if(nbPop < 2){
    for(l=0; l<_nb_locus; ++l){
      _fis_WC_locus[l] = _fit_WC_locus[l] = _fis_WC_locus[l] = my_NAN;
    }
    delete[] pop_sizes;
    ARRAY::delete_3D(ho_patch_allele, _sample_pops_size, _nb_locus);
    return;
  }

  nc = (tot_size - (tot_square_size/tot_size))/(nbPop-1);


  // for each locus
  for(l = 0; l < _nb_locus; ++l) {                      // for each locus
    a2 = b2 = w2 = 0;

    get_sigma_of_locus_Weir_Cockerham(a2, b2, w2, _sample_pops_size,nbPop,
              _alleleFreqTable, _globalAlleleFreq, l, tot_size, pop_sizes, ho_patch_allele, nc);

    if(a2 || b2 || w2){
      _fst_WC_locus[l] = a2/(a2 + b2 + w2);
      _fit_WC_locus[l] = (a2 + b2)/(a2 + b2 + w2);
    }
    else _fst_WC_locus[l] = _fit_WC_locus[l] = my_NAN;

    if(b2 || w2) _fis_WC_locus[l] = b2/(b2 + w2);
    else         _fis_WC_locus[l] = my_NAN;
  }


  delete[] pop_sizes;
  ARRAY::delete_3D(ho_patch_allele, _sample_pops_size, _nb_locus);
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham_perLocus
// ----------------------------------------------------------------------------------------
/** computes the pairwise Fst following Weir and Cockerham (1984) */
template <class SH>
void StatHandler<SH>::setFstat_Weir_Cockerham_perPatchPair(age_t AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[23], AGE)) return;

  unsigned int l;

  if(!_fst_matrix_wc) ARRAY::create_2D(_fst_matrix_wc, _sample_pops_size, _sample_pops_size);

  set_alleleCountTable_local(AGE);

  unsigned int* pop_sizes   = new unsigned int[_sample_pops_size];
  double*** ho_patch_allele = new double**[_sample_pops_size];
  unsigned int tot_size, nbPop = 2, i, j;
  double a2, b2, w2, tot_square_size, nc;
  unsigned int    cur_pop_sizes[2];
  double**        cur_ho_patch_allele[2];
  double**        cur_allele_freq_local[2];
  unsigned int* * cur_count_allele[2];
  double**        cur_allele_freq_global = ARRAY::new_2D<double>(_nb_locus, _nb_allele);

  // get pop sizes and Ho
  for(i = 0; i < _sample_pops_size; ++i) {
    pop_sizes[i] = _sample_pops[i]->size(AGE);
    if(pop_sizes[i]) ho_patch_allele[i] = getHo_ofPatchperAllele (AGE, _sample_pops[i]);
    else             ho_patch_allele[i] = NULL;
  }

  for(i=0; i<_sample_pops_size-1; ++i){
    if(!pop_sizes[i]){                   // both pops have to be populated
      for(j=i+1; j<_sample_pops_size; ++j){
        _fst_matrix_wc[i][j] = my_NAN;
      }
      continue;
    }
    cur_ho_patch_allele[0]   = ho_patch_allele[i];
    cur_pop_sizes[0]         = pop_sizes[i];
    cur_allele_freq_local[0] = _alleleFreqTable[i];
    cur_count_allele[0]      = _alleleCountTable[i];

    for(j=i+1; j<_sample_pops_size; ++j){
      if(!pop_sizes[j]){                 // both pops have to be populated!
        _fst_matrix_wc[i][j] = my_NAN;
        continue;
      }
      cur_ho_patch_allele[1]   = ho_patch_allele[j];
      cur_pop_sizes[1]         = pop_sizes[j];
      cur_allele_freq_local[1] = _alleleFreqTable[j];
      cur_count_allele[1]      = _alleleCountTable[j];

      tot_size = pop_sizes[i] + pop_sizes[j];
      tot_square_size = pop_sizes[i]*pop_sizes[i] + pop_sizes[j]*pop_sizes[j];

      // compute the global allele frequencies
      set_alleleCountTable_global(2, tot_size,  cur_count_allele, cur_allele_freq_global);

      nc = (tot_size - (tot_square_size/tot_size))/(nbPop-1);


      // loop over each lcous
      a2 = b2 = w2 = 0;
      for(l = 0; l < _nb_locus; ++l) {                      // for each locus
        get_sigma_of_locus_Weir_Cockerham(a2, b2, w2, 2, nbPop,
              cur_allele_freq_local, cur_allele_freq_global, l, tot_size, cur_pop_sizes, cur_ho_patch_allele, nc);
      }

      if(a2 || b2 || w2) _fst_matrix_wc[i][j] = a2/(a2 + b2 + w2);
      else               _fst_matrix_wc[i][j] = my_NAN;
    }
  }

  delete[] pop_sizes;
  ARRAY::delete_3D(ho_patch_allele, _sample_pops_size, _nb_locus);
  ARRAY::delete_2D(cur_allele_freq_global, _nb_locus);

}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham
// ----------------------------------------------------------------------------------------
/** computes the global F-statistics following Weir and Cockerham (1984) for the given patches */
template <class SH>
void StatHandler<SH>::setFstat_Weir_Cockerham(age_t AGE, Patch** aPatch, unsigned int nbPatch,
                                                double*** alleleFreqs, double** alleleFreqsGlobal,
                                                double& fst, double& fis, double &fit)
{
  unsigned int* pop_sizes   = new unsigned int[nbPatch];
  double*** ho_patch_allele = new double**[nbPatch];
  unsigned int tot_size=0, cur_size, p, l, nbPop = 0;
  double a2, b2, w2, tot_square_size=0, nc;

  // get pop sizes
  for(p = 0; p < nbPatch; ++p) {
    pop_sizes[p] = cur_size = aPatch[p]->size(AGE);
    if(!cur_size){
      ho_patch_allele[p] = NULL;
      continue;       // if pop is empty
    }
    ++nbPop;                      // get the number of populated pops
    tot_size += cur_size;
    ho_patch_allele[p] = getHo_ofPatchperAllele (AGE, aPatch[p]);
    tot_square_size += cur_size*cur_size;
  }

  if(nbPop<=1){
    fst = fis = fit = my_NAN;
    delete[] pop_sizes;
    ARRAY::delete_3D(ho_patch_allele, nbPatch, _nb_locus);
    return;
  }

  nc = (tot_size - (tot_square_size/tot_size))/(nbPop-1);


  // compute sigma(a2), sigma(b2) and sigma(w2) for each locus and sum them up
  a2 = b2 = w2 = 0;
  for(l = 0; l < _nb_locus; ++l) {
    get_sigma_of_locus_Weir_Cockerham(a2,b2,w2,nbPatch,nbPop,
                 alleleFreqs,alleleFreqsGlobal,l,tot_size,pop_sizes,ho_patch_allele, nc);
  }

  if(a2 || b2 || w2){
    fst = a2/(a2 + b2 + w2);
    fit = (a2 + b2)/(a2 + b2 + w2);
  }
  else fst = _fit_wc = my_NAN;

  if(b2 || w2) fis = b2/(b2 + w2);
	else         fis = my_NAN;

	delete[] pop_sizes;
  ARRAY::delete_3D(ho_patch_allele, nbPatch, _nb_locus);
}

// ----------------------------------------------------------------------------------------
// setFstat_Weir_Cockerham
// ----------------------------------------------------------------------------------------
/** computes the trhee sigma for each locus separately (F-statistics following Weir and Cockerham (1984)) */
template <class SH>
void StatHandler<SH>::get_sigma_of_locus_Weir_Cockerham(double& sigma_a2,
          double& sigma_b2, double& sigma_w2, const unsigned int& nbPatch, const unsigned int& nbPop,
          double*** alleleFreqs, double** alleleFreqsGlobal, const unsigned int& l,
          const unsigned int& tot_size, unsigned int* pop_sizes, double*** ho_patch_allele, const double& nc)
{
  double p_freq, sum1, sum2, nom, bloc, val;
  unsigned int a, p;
  for(a = 0; a < _nb_allele; ++a) {                   // for each allele
    p_freq = alleleFreqsGlobal[l][a];                 // the global allele frequency of allele a (or weighted frequency)

    sum1 = tot_size*p_freq*(1-p_freq);

    sum2 = 0;
    nom  = 0;
    for(p=0; p<nbPatch; ++p){
      if(!pop_sizes[p]) continue;
      val = alleleFreqs[p][l][a] - p_freq;       // the allele frequency of allele a in patch p
      sum2 += pop_sizes[p]*val*val;
      nom +=ho_patch_allele[p][l][a]*pop_sizes[p];
    }

    bloc = (sum1-sum2-(nom/4))/(tot_size-nbPop);

    sigma_a2 += (sum2/(nbPop-1) - bloc)/nc;                 // sigma(a)^2
    sigma_b2 += bloc - (nom/(4*tot_size));                  // sigma(b)^2
    sigma_w2 += nom/(2*tot_size);                           // sigma(w)^2
  }
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                          ******** coancestries analysis ********
// ----------------------------------------------------------------------------------------
// coancestry
// ----------------------------------------------------------------------------------------
template <class SH>
double StatHandler<SH>::Coancestry(void** ind1, void** ind2)
{
  unsigned int p = 0;
  unsigned char **seq1 = (unsigned char**)ind1;
  unsigned char **seq2 = (unsigned char**)ind2;

  for (unsigned int k = 0; k < _nb_locus; ++k){
	  p += !(seq1[k][0]^seq2[k][0])
      + !(seq1[k][0]^seq2[k][1])
      + !(seq1[k][1]^seq2[k][0])
      + !(seq1[k][1]^seq2[k][1]);
  }

  return (double)p/(4.*_nb_locus);
}
// ----------------------------------------------------------------------------------------
// setCoaMatrixTheta
// ----------------------------------------------------------------------------------------
/** computes the pairwise coancestry matrix: diagonal (theta)*/
template <class SH>
void StatHandler<SH>::setCoaMatrixTheta (const age_t& AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[16], AGE)) return;

  unsigned int Fsize, Msize, tot_size, wt, i, j, k;
  Patch *P1;
  double coa;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

  // create the matrix if not present
  if(!_coa_matrix) ARRAY::create_2D(_coa_matrix, _sample_pops_size, _sample_pops_size);

  //first fill the diagonale: within deme coancestry (theta)
  _mean_theta = 0;
  wt = 0;

  for(i = 0; i < _sample_pops_size; ++i) {
    P1 = _sample_pops[i];
    Fsize = P1->size(FEM, age_pos);
    Msize = P1->size(MAL, age_pos);
    tot_size = Fsize + Msize;
    coa = 0;
    if(tot_size) {
      //fem-fem coa
      for(j = 0; j < Fsize-1; ++j){
        for(k = j+1; k < Fsize; ++k){
          coa += Coancestry(P1->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                            P1->get(FEM, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence());
        }
      }
      //mal-mal coa
      for(j = 0; j < Msize-1; ++j){
        for(k = j+1; k < Msize; ++k){
          coa += Coancestry(P1->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                            P1->get(MAL, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence());
        }
      }
      //fem-mal coa
      coa += get_coancestry(P1, FEM, Fsize, P1, MAL, Msize, age_pos);

      coa /= tot_size*(tot_size -1)/2.0;
    }//end if

    _coa_matrix[i][i] = coa;
    _mean_theta += tot_size * coa;
    wt += tot_size;
  }//end for patchNbr

  _mean_theta /= wt;        //weighted average
}


// ----------------------------------------------------------------------------------------
// setCoaMatrixAlpha
// ----------------------------------------------------------------------------------------
/** computes the pairwise coancestry matrix: upper half part (alpha)*/
template <class SH>
void StatHandler<SH>::setCoaMatrixAlpha (const age_t& AGE)
{
  //fill the first upper half of the matrix: between deme coancestry (alpha)
  // check if the table has already been computed
  if(already_computed(_computed[17], AGE)) return;

  unsigned int Fsize1, Msize1, Fsize2, Msize2, tot_size, wt, i, l;
  double coa;
  Patch *P1, *P2;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

  // create the matrix if not present
  if(!_coa_matrix) ARRAY::create_2D(_coa_matrix, _sample_pops_size, _sample_pops_size);

  _mean_alpha = 0;
  wt = 0;

  for(i = 0; i < _sample_pops_size-1; ++i) {   // first patch
    P1 = _sample_pops[i];
    Fsize1 = P1->size(FEM, age_pos);
    Msize1 = P1->size(MAL, age_pos);
    tot_size = Fsize1 + Msize1;  		// size of patch 1
  	if(!tot_size){      // if patch is empty continue with the next patch
    	for(l = i+1; l < _sample_pops_size; ++l) {   // set the corrresponding elements to zero
      	_coa_matrix[i][l] = 0;
      }
    	continue;
    }
    for(l = i+1; l < _sample_pops_size; ++l) {   // second patch
      P2 = _sample_pops[l];
    	Fsize2 = P2->size(FEM, age_pos);
   		Msize2 = P2->size(MAL, age_pos);
      tot_size += Fsize2 + Msize2;  // size of patch 1 and patch 2
      coa = 0;
      if(Fsize2 + Msize2) {                // if patch is empty continue with the next patch
        coa += get_coancestry(P1, FEM, Fsize1, P2, FEM, Fsize2, age_pos); // fem-fem coa
        coa += get_coancestry(P1, FEM, Fsize1, P2, MAL, Msize2, age_pos); // fem-mal coa
        coa += get_coancestry(P1, MAL, Msize1, P2, FEM, Fsize2, age_pos); // mal-fem coa
        coa += get_coancestry(P1, MAL, Msize1, P2, MAL, Msize2, age_pos); // mal-mal coa
        coa /= (Fsize1*Fsize2) + (Fsize1*Msize2) + (Msize1*Fsize2) + (Msize1*Msize2);
      }//endif

      _coa_matrix[i][l] = coa;
      _mean_alpha += tot_size * coa;
      wt += tot_size;
    }//end for P2
  }//end for P1
  _mean_alpha /= wt;           //weighted average
}

// ----------------------------------------------------------------------------------------
// setSexspecific_Theta
// ----------------------------------------------------------------------------------------
template <class SH>
double StatHandler<SH>::get_coancestry(Patch* P1, const sex_t& SEX1, const unsigned int& size1,
 																						 Patch* P2, const sex_t& SEX2, const unsigned int& size2,
                                             const age_idx& age_pos){
  unsigned int i, j;
  double sum = 0;
	for(i = 0; i < size1; ++i){           // for each individual of patch 1
		for(j = 0; j < size2; ++j){         // for each individual of patch 2
      sum += Coancestry(P1->get(SEX1, age_pos, i)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                                         P2->get(SEX2, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence());
		}
	}
  return sum;
}

// ----------------------------------------------------------------------------------------
// setSexspecific_Theta
// ----------------------------------------------------------------------------------------
template <class SH>
void StatHandler<SH>::setSexspecific_Theta(const age_t& AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[18], AGE)) return;

  unsigned int Fsize,Msize,i,j,k, FFsize, MMsize, FMsize;
  Patch *P1;
  double mean, grand_mean;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

  Theta_FF = 0;
  Theta_MM = 0;
  Theta_FM = 0;

  _mean_theta = 0;

  for(i = 0; i < _sample_pops_size; ++i) {
    P1 = _sample_pops[i];
	  Fsize = P1->size(FEM, age_pos);
	  Msize = P1->size(MAL, age_pos);

    FFsize = Fsize*(Fsize-1)/2;
    MMsize = Msize*(Msize-1)/2;
    FMsize = Fsize*Msize;

    grand_mean = 0;

	  if(Fsize || Msize) {
      if(Fsize) {
        mean = 0;
        for(j = 0; j < Fsize-1;++j){
          for(k = j+1; k < Fsize;++k){
            mean += Coancestry(P1->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
            									 P1->get(FEM, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence());
          }
        }
        Theta_FF += mean/FFsize;
      }
      grand_mean += mean;
      if(Msize) {
        mean = 0;
        for(j = 0; j < Msize-1;++j){
          for(k = j+1; k < Msize;++k){
            mean += Coancestry(P1->get(MAL, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                               P1->get(MAL, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence());
          }
        }
        Theta_MM += mean/MMsize;
      }
      grand_mean += mean;
      if(Fsize && Msize) {
        mean = 0;
        for(j = 0; j < Fsize;++j){
          for(k = 0; k < Msize;++k){
            mean += Coancestry(P1->get(FEM, age_pos, j)->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                               P1->get(MAL, age_pos, k)->getTrait(_SHLinkedTraitIndex)->get_sequence());
          }
        }
        Theta_FM += mean/FMsize;
      }
      grand_mean += mean;
      _mean_theta += grand_mean / (FFsize + MMsize + FMsize);
    }
  }

  _mean_theta /= _sample_pops_size;
  Theta_FF    /= _sample_pops_size;
  Theta_MM    /= _sample_pops_size;
  Theta_FM    /= _sample_pops_size;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                          ******** kinship analysis ********
// ----------------------------------------------------------------------------------------
// setSibsStats
// ----------------------------------------------------------------------------------------
/** sets the sib coancestry and at the same time the kinship as it does n ot cost a lot of extra effort
  * (the kinship can separately been set by the function setKinship()
  */
template <class SH>
void StatHandler<SH>::setSibStats(const age_t& AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[19], AGE)) return;

  unsigned int i, j, k, Fsize, Msize;
  Individual *I1, *I2;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

  // if it is the first generation stop here
  if(_pop->getCurrentGeneration() == 1){
    for(i = 0; i < 5; ++i) {
	    _sib_prop[i] =  _sib_coa[i]  = my_NAN;
    }
    return;
  }

  //counters initialization
  for(i = 0; i < 5; ++i) {
    _sib_prop[i] =  _sib_coa[i]  = 0.0;
  }

  Patch **cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
  for(;cur != end; ++cur){
		if((Fsize = (*cur)->size(FEM, age_pos)) != 0) {
	    for(j = 0; j < Fsize -1; ++j) {
      	I1 = (*cur)->get(FEM, age_pos, j);
			  for(k = j+1; k < Fsize; ++k) {
          I2 = (*cur)->get(FEM, age_pos, k);
          setSibCoa(I1, I2);
			  }
        //selfed offspring counter:
        if(I1->getIsSelfed()) _sib_prop[4]++;
	    }
    }

		if((Msize = (*cur)->size(MAL, age_pos)) != 0) {
      for(j = 0; j < Msize -1; ++j) {
        I1 = (*cur)->get(MAL, age_pos, j);
        for(k = j+1; k < Msize; ++k) {
          I2 = (*cur)->get(MAL, age_pos, k);
          setSibCoa(I1, I2);
        }
        //selfed offspring counter:
        if(I1->getIsSelfed()) _sib_prop[4]++;
      }
    }

    //male-female
    for(j = 0; j < Msize; ++j) {
      I1 = (*cur)->get(MAL, age_pos, j);
      for(k = 0; k < Fsize; ++k) {
        I2 = (*cur)->get(FEM, age_pos, k);
        setSibCoa(I1, I2);
      }
    }
  }

  double tot = _sib_prop[0] + _sib_prop[1] + _sib_prop[2] + _sib_prop[3];

  for(i = 0 ; i < 4; ++i) {
	  _sib_coa[i] = ( (_sib_prop[i]) ? _sib_coa[i]/_sib_prop[i] : my_NAN);
	  _sib_prop[i] /= tot;
  }
  _sib_prop[4] /= _pop->size(AGE);
}

// ----------------------------------------------------------------------------------------
// setSibCoa
// ----------------------------------------------------------------------------------------
/** sets the sib coancestry and at the same time the kinship as it does n ot cost a lot of extra effort
  * (the kinship can separately been set by the function setKinship()
  */
template <class SH>
void StatHandler<SH>::setSibCoa(Individual *I1, Individual *I2)
{
	double coa = Coancestry(I1->getTrait(_SHLinkedTraitIndex)->get_sequence(),
                          I2->getTrait(_SHLinkedTraitIndex)->get_sequence());
	if(I1->getMotherID() == I2->getMotherID()){
		if(I1->getFatherID() == I2->getFatherID()) {_sib_prop[3]++; _sib_coa[3]+=coa;}  // full sibs
    else 																	     {_sib_prop[1]++; _sib_coa[1]+=coa;}  // maternal half sibs
	}
	else{
		if(I1->getFatherID() == I2->getFatherID()) {_sib_prop[2]++; _sib_coa[2]+=coa;}  // paternal half sibs
		else                                       {_sib_prop[0]++; _sib_coa[0]+=coa;}  // non sibs
	}
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                          ******** stat options ********
// ----------------------------------------------------------------------------------------
// set_stat_coancestry
// ----------------------------------------------------------------------------------------
/** coancestry stat options */
template <class SH>
bool StatHandler<SH>:: set_stat_coancestry(string t, string i, string trait){

  // check the prefix
  assert(i.length() == 1);
  if(t[0] != i[0]) return false; // first character

  // get the age
  int pos = t.find('.', 2);  // find the second '.'
  string ageStr, token = t.substr(2, pos-2);
  age_t AGE;
  if(token == "adlt")     {AGE = ADULTS;  ageStr = "adult";}
  else if(token == "off") {AGE = OFFSPRG; ageStr = "offsprg";}
  else return false;

  t = t.substr(pos+1);  // get the search token

	     if(t == "theta") 	  add("Within patch coancestry ("+ageStr+", "+trait+")",i+"."+token+".theta",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getMeanTheta);
	else if(t == "alpha")		  add("Between patch coancestry ("+ageStr+", "+trait+")",i+"."+token+".alpha",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getMeanAlpha);
	else if(t == "thetaFF")   add("Mean within patch, within females coancestry ("+ageStr+", "+trait+")",i+"."+token+".thetaFF",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getTheta_FF);
	else if(t == "thetaMM")   add("Mean within patch, within males coancestry ("+ageStr+", "+trait+")",i+"."+token+".thetaMM",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getTheta_MM);
	else if(t == "thetaFM")   add("Mean within patch, between sexes coancestry ("+ageStr+", "+trait+")",i+"."+token+".thetaFM",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getTheta_FM);
  else if(t == "coa.fsib")  add("Coancestry of full-sib ("+ageStr+", "+trait+")",i+"."+token+".coa.fsib",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getSibCoaMeans);
  else if(t == "coa.phsib") add("Coancestry of paternal half-sib ("+ageStr+", "+trait+")",i+"."+token+".coa.phsib",FLAT,AGE,1,0,0,0,0,0,&StatHandler<SH>::getSibCoaMeans);
	else if(t == "coa.mhsib")	add("Coancestry of maternal half-sib ("+ageStr+", "+trait+")",i+"."+token+".coa.mhsib",FLAT,AGE,2,0,0,0,0,0,&StatHandler<SH>::getSibCoaMeans);
  else if(t == "coa.nsib")  add("Coancestry of non-sib ("+ageStr+", "+trait+")",i+"."+token+".coa.nsib",FLAT,AGE,3,0,0,0,0,0,&StatHandler<SH>::getSibCoaMeans);
	else if(t == "theta_p")   add_perPatch("Within coancestry ("+ageStr+", "+trait+")",i+"."+token+".theta.p", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getCoaTheta);
	else if(t == "alpha_pair")add_pairwisePatch("Mean coancestry ("+ageStr+", "+trait+")",i+"."+token+".alpha.pair", FLAT, AGE, 0, 0, 0, 0, 0, 0, &StatHandler<SH>::getCoaAlpha);

  else if(t == "coa"){
    add("Within patch coancestry ("+ageStr+", "+trait+")",i+"."+token+".theta",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getMeanTheta);
    add("Between patch coancestry ("+ageStr+", "+trait+")",i+"."+token+".alpha",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getMeanAlpha);
    add("Mean within patch, within females coancestry ("+ageStr+", "+trait+")",i+"."+token+".thetaFF",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getTheta_FF);
    add("Mean within patch, within males coancestry ("+ageStr+", "+trait+")",i+"."+token+".thetaMM",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getTheta_MM);
    add("Mean within patch, between sexes coancestry ("+ageStr+", "+trait+")",i+"."+token+".thetaFM",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getTheta_FM);
    add("Coancestry of full-sib ("+ageStr+", "+trait+")",i+"."+token+".coa.fsib",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getSibCoaMeans);
    add("Coancestry of paternal half-sib ("+ageStr+", "+trait+")",i+"."+token+".coa.phsib",FLAT,AGE,1,0,0,0,0,0,&StatHandler<SH>::getSibCoaMeans);
    add("Coancestry of maternal half-sib ("+ageStr+", "+trait+")",i+"."+token+".coa.mhsib",FLAT,AGE,2,0,0,0,0,0,&StatHandler<SH>::getSibCoaMeans);
    add("Coancestry of non-sib ("+ageStr+", "+trait+")",i+"."+token+".coa.nsib",FLAT,AGE,3,0,0,0,0,0,&StatHandler<SH>::getSibCoaMeans);
  }
  else return false;
  return true;
}


// ----------------------------------------------------------------------------------------
// set_stat_fstat
// ----------------------------------------------------------------------------------------
/** fstat stat options */
template <class SH>
bool StatHandler<SH>:: set_stat_fstat(string t, string i, string trait){

  // check the prefix
  assert(i.length() == 1);
  if(t[0] != i[0]) return false; // first character

  // get the age
  int pos = t.find('.', 2);  // find the second '.'
  string ageStr, token = t.substr(2, pos-2);
  age_t AGE;
  if(token == "adlt")     {AGE = ADULTS;  ageStr = "adult";}
  else if(token == "off") {AGE = OFFSPRG; ageStr = "offsprg";}
  else return false;

  t = t.substr(pos+1);  // get the search token

  // genetic diversity
	     if(t == "nbAll")  	add("Mean number of alleles/locus ("+ageStr+", "+trait+")",i+"."+token+".nbAll",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getNbAllGlobal);
	else if(t == "meanAll") add("Mean number of alleles/locus/patch ("+ageStr+", "+trait+")",i+"."+token+".meanAll",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getNbAllLocal);
	else if(t == "nbFixLoc") 	 add("Total number of fixed loci ("+ageStr+", "+trait+")",i+"."+token+".nbFixLoc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFixLocGlobal);
	else if(t == "meanFixLoc") add("Mean number of fixes loci/patch ("+ageStr+", "+trait+")",i+"."+token+".meanFixLoc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFixLocLocal);
	else if(t == "ho") 			add("Ho (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".ho",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getHo);
	else if(t == "hs") 	    add("Hs (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".hs",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getHsnei);
	else if(t == "ht")    	add("Ht (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".ht",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getHtnei);
	else if(t == "nbAll_p") add_perPatch("Mean number of alleles/locus ("+ageStr+", "+trait+")",i+"."+token+".nbAll",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getNbAll_perPatch);
	else if(t == "nbFixLoc_p") add_perPatch("Total number of fixed loci ("+ageStr+", "+trait+")",i+"."+token+".nbFixLoc",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getFixLoc_perPatch);
	else if(t == "ho_p") 		add_perPatch("Ho (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".ho",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getHo_ofPatch);
	else if(t == "hs_p") 	  add_perPatch("Hs (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".hs",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getHsUnbiased_ofPatch);

  else if(t == "gendiv") {
		add("Mean number of alleles/locus ("+ageStr+", "+trait+")",i+"."+token+".nbAll",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getNbAllGlobal);
		add("Mean number of alleles/locus/patch ("+ageStr+", "+trait+")",i+"."+token+".meanAll",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getNbAllLocal);
		add("Total number of fixed loci ("+ageStr+", "+trait+")",i+"."+token+".nbFixLoc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFixLocGlobal);
		add("Mean number of fixes loci/patch ("+ageStr+", "+trait+")",i+"."+token+".meanFixLoc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFixLocLocal);
		add("Ho (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".ho",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getHo);
	  add("Hs (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".hs",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getHsnei);
		add("Ht (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".ht",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getHtnei);
  }

  // F-stats following Nei and Chesser
	else if(t == "fst") 		add("Fst (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fst",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFst);
	else if(t == "fis") 		add("Fis (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fis",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFis);
	else if(t == "fit")     add("Fit (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fit",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFit);
	else if(t == "fst_pair")   add_pairwisePatch("Fst (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fst.pair",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getFst_ij);

  else if(t == "fstat") {
		add("Fst (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fst",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFst);
		add("Fis (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fis",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFis);
	  add("Fit (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fit",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFit);
  }

  // F-stats following Weir and Cockerham
  else if(t == "fst.wc")	add("Fst (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fst.wc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFst_WC);
  else if(t == "fis.wc")	add("Fis (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fis.wc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFis_WC);
  else if(t == "fit.wc")	add("Fit (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fit.wc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFit_WC);
  else if(t == "fst.wc_pair")	add_pairwisePatch("Fst(Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fst.wc.pair",FLAT,AGE, 0,0,0,0,0,0,&StatHandler<SH>::getFst_WC_ij);

	else if(t == "fstat.wc") {
		add("Fst (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fst.wc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFst_WC);
		add("Fis (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fis.wc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFis_WC);
		add("Fit (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fit.wc",FLAT,AGE,0,0,0,0,&StatHandler<SH>::getFit_WC);
	}

	// statistics per locus
	else if(t == "ho_l") 		add_perLocus("Ho (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".ho",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getHo_ofLocus);
	else if(t == "hs_l")    add_perLocus("Hs (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".hs",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getHsnei_ofLocus);
	else if(t == "ht_l")    add_perLocus("Ht (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".ht",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getHtnei_ofLocus);
	else if(t == "fst_l") 	add_perLocus("Fst (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fst",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getFst_ofLocus);
	else if(t == "fis_l") 	add_perLocus("Fis (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fis",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getFis_ofLocus);
	else if(t == "fit_l")   add_perLocus("Fit (Nei & Chesser, 1983)("+ageStr+", "+trait+")",i+"."+token+".fit",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getFit_ofLocus);
	else if(t == "fst.wc_l")add_perLocus("Fst (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fst.wc",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getFst_WC_ofLocus);
	else if(t == "fis.wc_l")add_perLocus("Fis (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fis.wc",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getFis_WC_ofLocus);
	else if(t == "fit.wc_l")add_perLocus("Fit (Weir & Cockerham, 1984)("+ageStr+", "+trait+")",i+"."+token+".fit.wc",FLAT,AGE,0,0,0,0,0,0,&StatHandler<SH>::getFit_WC_ofLocus);

	else return false;
  return true;
}

// ----------------------------------------------------------------------------------------
// set_stat_all_freq_local
// ----------------------------------------------------------------------------------------
/** local allele freuqencies for each patch stat options */
template <class SH>
bool StatHandler<SH>:: set_stat_all_freq_local(string t, string i, string trait){

  // check the prefix
  assert(i.length() == 1);
  if(t[0] != i[0]) return false; // first character

  // get the age
  int pos = t.find('.', 2);  // find the second '.'
  string ageStr, token = t.substr(2, pos-2);
  age_t AGE;
  if(token == "adlt")     {AGE = ADULTS;  ageStr = "adult";}
  else if(token == "off") {AGE = OFFSPRG; ageStr = "offsprg";}
  else return false;

  t = t.substr(pos+1);  // get the search token

  // group allele frequencies
  if(t=="a.freq") {
    unsigned int a, l, p;
    string _i, n;
    Patch* cur;

    for(p=0; p<_sample_pops_size; ++p){
      cur = _sample_pops[p];
      for(l=0; l<_nb_locus; ++l){
        for(a=0; a<_nb_allele; ++a){
          _i=i+".a_freq";
          n = "Local allele frequency";
          if(_sample_pops_size>1){
            _i += STRING::int2str_(cur->get_ID()+1, _pop->getPatchNbr(), "p");    // "_p1"
            n  += " of population " + STRING::int2str(cur->get_ID()+1);
          }
          if(_nb_locus>1){
            _i += STRING::int2str_(l+1, _nb_locus, "l");    // "_l1"
            if(!n.empty()) n += ", ";
            n  += " locus " + STRING::int2str(l+1);
          }
          if(_nb_allele>1){
            _i += STRING::int2str_(a+1, _nb_allele, "a");    // "_a1"
            if(!n.empty()) n += ", and ";
            n  += " allele " + STRING::int2str(a+1);
          }
          n += " ("+ageStr+", "+trait+")";
          add(n,_i,FLAT,AGE, p*(_nb_allele*_nb_locus) + l*_nb_allele + a,0,0,0,0,0,&StatHandler<SH>::get_allele_freq_local);
        }
      }
    }
  }
  else return false;
  return true;
}

// ----------------------------------------------------------------------------------------
// set_stat_all_freq_global
// ----------------------------------------------------------------------------------------
/** globel allele freuqencies for each locus and allele stat options */
template <class SH>
bool StatHandler<SH>:: set_stat_all_freq_global(string t, string i, string trait){

  // check the prefix
  assert(i.length() == 1);
  if(t[0] != i[0]) return false; // first character

  // get the age
  int pos = t.find('.', 2);  // find the second '.'
  string ageStr, token = t.substr(2, pos-2);
  age_t AGE;
  if(token == "adlt")     {AGE = ADULTS;  ageStr = "adult";}
  else if(token == "off") {AGE = OFFSPRG; ageStr = "offsprg";}
  else return false;

  t = t.substr(pos+1);  // get the search token

  // group allele frequencies
  if(t=="a.freq.global") {
    unsigned int a, l;
    string _i, n;

    for(l=0; l<_nb_locus; ++l){
      for(a=0; a<_nb_allele; ++a){
        _i=i+".a_freq";
        n = "Global allele frequency of ";
        if(_nb_locus>1){
          _i += STRING::int2str_(l+1, _nb_locus, "l");    // "_l1"
          if(!n.empty()) n += ", ";
          n  += " locus " + STRING::int2str(l+1);
        }
        if(_nb_allele>1){
          _i += STRING::int2str_(a+1, _nb_allele, "a");    // "_a1"
          if(!n.empty()) n += ", and ";
          n  += " allele " + STRING::int2str(a+1);
        }
        n += " ("+ageStr+", "+trait+")";
        add(n,_i,FLAT,AGE, l*_nb_allele + a,0,0,0,0,0,&StatHandler<SH>::get_allele_freq_global);
      }
    }
  }
  else return false;
  return true;
}
