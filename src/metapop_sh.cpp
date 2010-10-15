/** @file metapop_sh.cpp
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


#include <sstream>
#include "metapop_sh.h"
#include "metapop.h"
#include "stathandler.cpp"


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** StatHandler ********

// ----------------------------------------------------------------------------------------
// setStatsRecorders
// ----------------------------------------------------------------------------------------
bool MetapopSH::setStatRecorders(const string& t)
{
  // patch extinction
       if(t == "ext.rate") add("Extinction rate","ext.rate",FLAT,ALL,0,&MetapopSH::getObsrvdExtinctionRate);

  // Fecundity
	else if(t == "fem.meanFec") add("Mean fecundity of females","fem.meanFec",FLAT,ADULTS,1,0,&MetapopSH::getReproductiveMean);
	else if(t == "fem.varFec")  add("Variance of the fecundity of females","fem.varFec",FLAT,ADULTS,1,0,&MetapopSH::getReproductiveVar);
	else if(t == "mal.meanFec") add("Mean fecundity of males","mal.meanFec",FLAT,ADULTS,0,0,&MetapopSH::getReproductiveMean);
	else if(t == "mal.varFec")  add("Variance of the fecundity of males","mal.varFec",FLAT,ADULTS,0,0,&MetapopSH::getReproductiveVar);

	else if(t == "fecundity"){
		add("Mean realized female fecundity","fem.real.fec",FLAT,ADULTS,1,0,&MetapopSH::getReproductiveMean);
		add("Variance of the mean realized female fecundity","fem.var.fec",FLAT,ADULTS,1,0,&MetapopSH::getReproductiveVar);
		add("Mean realized male fecundity  ","mal.real.fec",FLAT,ADULTS,0,0,&MetapopSH::getReproductiveMean);
		add("Variance of the mean realized male fecundity","mal.var.fec",FLAT,ADULTS,0,0,&MetapopSH::getReproductiveVar);
  }

  // migration
  else if(t == "emigrants")  add("Mean number of emigrants/patch","emigrants",FLAT,ADULTS,0,&MetapopSH::getMeanEmigrantPerPatch);
	else if(t == "immigrants") add("Mean number of immigrants/patch","immigrants",FLAT,ADULTS,0,&MetapopSH::getMeanImmigrantPerPatch);
	else if(t == "residents")  add("Mean number of residents/patch","residents",FLAT,ADULTS,0,&MetapopSH::getMeanResidantPerPatch);
	else if(t == "immigrate")  add("Mean effective immigration rate","immigrate",FLAT,ADULTS,0,&MetapopSH::getMeanMigrantRatio);
	else if(t == "colonisers") add("Mean number of colonizers per extinct patch","colonisers",FLAT,ADULTS,0,&MetapopSH::getMeanKolonisersPerPatch);
	else if(t == "colon.rate") add("Mean effective colonization rate of extinct patches","colon.rate",FLAT,ADULTS,0,&MetapopSH::getMeanKolonisersProportion);

  else if(t == "migration"){
  	add("Mean number of emigrants/patch","emigrants",FLAT,ADULTS,0,&MetapopSH::getMeanEmigrantPerPatch);
    add("Mean number of immigrants/patch","immigrants",FLAT,ADULTS,0,&MetapopSH::getMeanImmigrantPerPatch);
    add("Mean number of residents/patch","residents",FLAT,ADULTS,0,&MetapopSH::getMeanResidantPerPatch);
    add("Mean effective immigration rate","immigrate",FLAT,ADULTS,0,&MetapopSH::getMeanMigrantRatio);
    add("Mean number of colonizers per extinct patch","colonisers",FLAT,ADULTS,0,&MetapopSH::getMeanKolonisersPerPatch);
    add("Mean effective colonization rate of extinct patches","colon.rate",FLAT,ADULTS,0,&MetapopSH::getMeanKolonisersProportion);
  }

	// fitness
  else if(t=="VwW")     add("Fitness variance of adults within patches ","VwW",FLAT,ADULTS,0,&MetapopSH::getVwW);
  else if(t=="VwB")     add("Fitness variance of adults between patches ","VwB",FLAT,ADULTS,0,&MetapopSH::getVwB);
	else if(t=="meanW_p")	add_perPatch("Fitness mean of adults","meanW",FLAT,ADULTS,0,0,0,&MetapopSH::getMeanW);
	else if(t=="varW_p")  add_perPatch("Fitness variance of adults","varW",FLAT,ADULTS,0,0,0,&MetapopSH::getVarW);

	else if(t=="fitness") {
		add("Fitness variance of adults within patches ","VwW",FLAT,ADULTS,0,&MetapopSH::getVwW);
		add("Fitness variance of adults between patches ","VwB",FLAT,ADULTS,0,&MetapopSH::getVwB);
	}

  else if(set_stat_kinship(t)){}     // adults and offspring
  else if(set_stat_demography(t)){}

  else return false;

  return true;
}

// ----------------------------------------------------------------------------------------
// set_stat_demography
// ----------------------------------------------------------------------------------------
/** fstat stat options */
bool MetapopSH::set_stat_demography(string t)
{
  // get the age
  int pos = t.find('.', 0);  // find the second '.'
  string ageStr, token = t.substr(0, pos);
  age_t AGE;
  if(token == "adlt")     {AGE = ADULTS;  ageStr = "adults";}
  else if(token == "off") {AGE = OFFSPRG; ageStr = "offsprings";}
  else return false;

  t = t.substr(pos+1);  // get the search token

	// Demography
       if(t == "nbInd")    add("Total number of "+ageStr,token+".nbInd",FLAT,AGE,0,0,0,0,&MetapopSH::getTotNbInd);
  else if(t == "nbFem")    add("Total number of female "+ageStr,token+".nbFem",FLAT,AGE,0,0,0,0,&MetapopSH::getTotNbFem);
  else if(t == "nbMal")    add("Total number of male "+ageStr,token+".nbMal",FLAT,AGE,0,0,0,0,&MetapopSH::getTotNbMal);
  else if(t == "meanInd")  add("Mean number of "+ageStr,token+".meanInd",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbInd);
  else if(t == "meanFem")  add("Mean number of female "+ageStr,token+".meanFem",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbFem);
  else if(t == "meanMal")  add("Mean number of male "+ageStr,token+".meanMal",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbMal);
  else if(t == "sexRatio") add("Sex ratio of "+ageStr,token+".sexRatio",FLAT,AGE,0,&MetapopSH::getAdultSexRatio);
	else if(t == "nbPops")   add("Number of inhabited patches by "+ageStr,token+".nbPops",FLAT,AGE,0,0,0,0,&MetapopSH::getNbPops);

  else if(t == "nbInd_p")  add_perPatch("Number of "+ageStr,token+".nbInd",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbInd);
  else if(t == "nbFem_p")  add_perPatch("Number of female "+ageStr,token+".nbFem",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbFem);
  else if(t == "nbMal_p")  add_perPatch("Number of male "+ageStr,token+".nbMal",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getNbMal);

  else if(t == "demo"){
    add("Total number of "+ageStr,token+".nbInd",FLAT,AGE,0,0,0,0,&MetapopSH::getTotNbInd);
    add("Total number of female "+ageStr,token+".nbFem",FLAT,AGE,0,0,0,0,&MetapopSH::getTotNbFem);
    add("Total number of male "+ageStr,token+".nbMal",FLAT,AGE,0,0,0,0,&MetapopSH::getTotNbMal);
    add("Mean number of "+ageStr,token+".meanInd",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbInd);
    add("Mean number of female "+ageStr,token+".meanFem",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbFem);
    add("Mean number of male "+ageStr,token+".meanMal",FLAT,AGE,0,0,0,0,&MetapopSH::getMeanNbMal);
    add("Sex ratio of "+ageStr,token+".sexRatio",FLAT,AGE,0,&MetapopSH::getAdultSexRatio);
	  add("Number of inhabited patches by "+ageStr,token+".nbPops",FLAT,AGE,0,0,0,0,&MetapopSH::getNbPops);
  }
  else return false;
  return true;
}

// ----------------------------------------------------------------------------------------
// set_stat_kinship
// ----------------------------------------------------------------------------------------
/** kinship stat options */
bool MetapopSH:: set_stat_kinship(string t){

  // get the age
  int pos = t.find('.');  // find the first '.'
  string ageStr, token = t.substr(0, pos);
  age_t AGE;
  if(token == "adlt")     {AGE = ADULTS;  ageStr = "adult";}
  else if(token == "off") {AGE = OFFSPRG; ageStr = "offsprg";}
  else return false;

  t = t.substr(pos+1);  // get the search token

       if(t == "fsib")  add("Mean proportion of full-sib ("+ageStr+")",token+".fsib",FLAT,AGE,3,0,0,0,0,0,&MetapopSH::getSibProportion);
  else if(t == "phsib") add("Mean proportion of paternal half-sib ("+ageStr+")",token+".phsib",FLAT,AGE,2,0,0,0,0,0,&MetapopSH::getSibProportion);
  else if(t == "mhsib") add("Mean proportion of maternal half-sib ("+ageStr+")",token+".mhsib",FLAT,AGE,1,0,0,0,0,0, &MetapopSH::getSibProportion);
  else if(t == "nsib")	add("Mean proportion of non-sib ("+ageStr+")",token+".nsib",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getSibProportion);
  else if(t == "self")  add("Mean proportion of selfed offspring ("+ageStr+")",token+".self",FLAT,AGE,4,0,0,0,0,0,&MetapopSH::getSibProportion);

  else if(t == "kinship"){
    add("Mean proportion of full-sib ("+ageStr+")",token+".fsib",FLAT,AGE,3,0,0,0,0,0,&MetapopSH::getSibProportion);
    add("Mean proportion of paternal half-sib ("+ageStr+")",token+".phsib",FLAT,AGE,2,0,0,0,0,0,&MetapopSH::getSibProportion);
    add("Mean proportion of maternal half-sib ("+ageStr+")",token+".mhsib",FLAT,AGE,1,0,0,0,0,0, &MetapopSH::getSibProportion);
    add("Mean proportion of non-sib ("+ageStr+")",token+".nsib",FLAT,AGE,0,0,0,0,0,0,&MetapopSH::getSibProportion);
    add("Mean proportion of selfed offspring ("+ageStr+")",token+".self",FLAT,AGE,4,0,0,0,0,0,&MetapopSH::getSibProportion);
  }
  else return false;
  return true;
}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                         ******** Migrants analysis ********
// ----------------------------------------------------------------------------------------
// getMeanEmigrantPerPatch
// ----------------------------------------------------------------------------------------
/** compute the mean number of emigrants per patch */
double MetapopSH::getMeanEmigrantPerPatch()
{
  // before any migration
  if(_current_generation == 1){
    meanEmigrant = my_NAN;
    return meanEmigrant;
  }

  meanEmigrant = 0;
  Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
  for(; cur != end; ++cur){
    meanEmigrant += (*cur)->nbEmigrant;
  }
  meanEmigrant /= _sample_pops_size;
  return meanEmigrant;
}
// ----------------------------------------------------------------------------------------
// getMeanImmigrantPerPatch
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanImmigrantPerPatch()
{
  // check if the table has already been computed
  if(already_computed(_computed[0])) return meanImmigrant;

  if(_current_nbPops[ADLTx] && _current_generation > 1){
    meanImmigrant = 0;
    Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
    for(; cur != end; ++cur){
      if(!(*cur)->get_isExtinct()) meanImmigrant += (*cur)->nbImmigrant;
    }
    meanImmigrant /= _current_nbPops[ADLTx];
  }
  else meanImmigrant = my_NAN;
  return meanImmigrant;
}
// ----------------------------------------------------------------------------------------
// getMeanResidantPerPatch
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanResidantPerPatch()
{
  // check if the table has already been computed
  if(already_computed(_computed[1])) return meanResidant;

  if(_current_nbPops[ADLTx] && _current_generation > 1){
    meanResidant = 0;
    Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
    for(; cur != end; ++cur){
      if(!(*cur)->get_isExtinct()) meanResidant += (*cur)->nbPhilopat;
    }
    meanResidant /= _current_nbPops[ADLTx];
  }
  else meanResidant = my_NAN;
  return meanResidant;
}
// ----------------------------------------------------------------------------------------
// getMeanMigrantRatio
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanMigrantRatio()
{
  // before any migration
  if(!_current_nbPops[ADLTx] || _current_generation == 1) return my_NAN;

  getMeanImmigrantPerPatch(); // set the variable meanImmigrant
  getMeanResidantPerPatch();  // set the variable meanResidant
  return ((meanImmigrant+meanResidant) ? meanImmigrant/(meanImmigrant+meanResidant) : 0);
}
// ----------------------------------------------------------------------------------------
// getMeanMigrantProportion
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanKolonisersProportion()
{
  double mean = 0;

  if(!_current_nbPops[ADLTx] || _current_generation == 1) return my_NAN;

  Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
  for(; cur != end; ++cur){
    if((*cur)->nbKolonisers != my_NAN) {
      mean += (double)(*cur)->nbKolonisers / (*cur)->get_K();
    }
  }
  return  mean/_current_nbPops[ADLTx];
}
// ----------------------------------------------------------------------------------------
// getMeanKolonisersPerPatch
// ----------------------------------------------------------------------------------------
double MetapopSH::getMeanKolonisersPerPatch()
{
  if(_current_nbPops[ADLTx] && _current_generation > 1){
    double colonizers;
    meanKolonisers = 0;
    Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
    for(; cur != end; ++cur){
      colonizers = (*cur)->nbKolonisers;
      if(colonizers != my_NAN) meanKolonisers += colonizers;
    }
    meanKolonisers /= _current_nbPops[ADLTx];
  }
  else meanKolonisers = my_NAN;
  return meanKolonisers;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                     ****** Patch extinction analysis ********

// ----------------------------------------------------------------------------------------
// setObsrvdExtinctionRate
// ----------------------------------------------------------------------------------------
void MetapopSH::setObsrvdExtinctionRate ()
{
  ObservedExtinctionRate = 0;
  Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
  for(; cur != end; ++cur){
	  ObservedExtinctionRate += (*cur)->isEmpty();
  }

  ObservedExtinctionRate /= _sample_pops_size;
}


/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                   ****** Fecundity and matings analysis ********

// ----------------------------------------------------------------------------------------
// getMeanRealizedFecundity
// ----------------------------------------------------------------------------------------
void MetapopSH::setReproductiveStats(bool sex)  // 0: MAL, 1: FEM
{
  // check if the table has already been computed
	if(already_computed(_computed[24], sex)) return;

	unsigned int j, size, tot_size = _current_nbIndsS[sex][ADLTx], v = 0;
	double *stat;

  if(!tot_size) {
		_mean_reprod_success = my_NAN;
    _var_reprod_success = my_NAN;
		return ;
  }
  stat = new double [tot_size];

  Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
  for(; cur != end; ++cur){
    for(j = 0, size = (*cur)->size((sex_t)sex, ADLTx); j < size; ++j){
      stat[v] = (*cur)->get((sex_t)sex, ADLTx, j)->getTotRealizedFecundity();
      v++;
    }
  }

	_mean_reprod_success = ARRAY::mean(stat, tot_size);
	_var_reprod_success  = ARRAY::var(stat, tot_size, _mean_reprod_success);
  delete [] stat;
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ****** kinship analysis ********
// ----------------------------------------------------------------------------------------
// setKinship
// ----------------------------------------------------------------------------------------
/** compute the mean kinship */
void MetapopSH::setKinship(const age_t& AGE)
{
  // check if the table has already been computed
  if(already_computed(_computed[2], AGE)) return;

	unsigned int j,k;
	unsigned int Msize, Fsize;
   Individual *I1,*I2;
  age_idx age_pos = (AGE == ADULTS ? ADLTx : OFFSx);

  // if it is the first generation stop here
  if(_current_generation == 1){
    for(j = 0; j < 5; ++j) {
	    _sib_prop[j] = my_NAN;
    }
    return;
  }

  //reset counters
  for(j = 0; j < 5; ++j){
    _sib_prop[j] = 0.0;
  }

  Patch** cur = _sample_pops, **end = _sample_pops + _sample_pops_size;
  for(; cur != end; ++cur){
    //male-male
		if ((Msize = (*cur)->size(MAL, age_pos)) != 0) {
		  for(j = 0; j < Msize -1; ++j) {
        I1 = (*cur)->get(MAL, age_pos, j);
			  for(k = j+1; k < Msize; ++k) {
          I2 = (*cur)->get(MAL, age_pos, k);
          setKinClassCounter(I1, I2);
        } //end for k < size

        //selfed offspring counter:
        if(I1->getIsSelfed()) _sib_prop[4]++;
		  } //end for j < size-1

      //don't forget the last one!
      if((*cur)->get(MAL, age_pos, Msize -1)->getIsSelfed()) _sib_prop[4]++;
		}//endif

      //female-female
		if ((Fsize = (*cur)->size(FEM, age_pos)) != 0) {
		  for(j = 0; j < Fsize -1; ++j) {
        I1 = (*cur)->get(FEM, age_pos, j);
        for(k = j+1; k < Fsize; ++k) {
          I2 = (*cur)->get(FEM, age_pos, k);
          setKinClassCounter(I1, I2);
        } //end for k < size

        //selfed offspring counter:
        if(I1->getIsSelfed()) _sib_prop[4]++;
		  } //end for j < size-1

      //don't forget the last one!
      if((*cur)->get(FEM, age_pos, Fsize -1)->getIsSelfed()) _sib_prop[4]++;
		}//endif

    //male-female
    for(j = 0; j < Msize; ++j) {
      I1 = (*cur)->get(MAL, age_pos, j);
      for(k = 0; k < Fsize; ++k) {
        I2 = (*cur)->get(FEM, age_pos, k);
        setKinClassCounter(I1, I2);
      } //end for k
    } //end for j
  } //end for i < patchNbr

  //total number of pairwise comparisons:
  double tot = _sib_prop[0] + _sib_prop[1] + _sib_prop[2] + _sib_prop[3];

  for(j = 0 ; j < 4; ++j){
    _sib_prop[j] /= tot;
  }

  _sib_prop[4] /= _current_nbInds[age_pos];
}

// ----------------------------------------------------------------------------------------
// setKinClassCounter
// ----------------------------------------------------------------------------------------
/** sets the kinship */
void MetapopSH::setKinClassCounter(Individual *I1, Individual *I2)
{
	if(I1->getMotherID() == I2->getMotherID()){
		if(I1->getFatherID() == I2->getFatherID()) _sib_prop[3]++;   // full sibs
    else 																	     _sib_prop[1]++;   // maternal half sibs
	}
	else{
		if(I1->getFatherID() == I2->getFatherID()) _sib_prop[2]++;   // paternal half sibs
		else                                       _sib_prop[0]++;   // non sibs
	}
}
// ----------------------------------------------------------------------------------------


double MetapopSH::getOffsprgSexRatio () {return _current_nbIndsS[FEM][OFFSx] ? (double)_current_nbIndsS[MAL][OFFSx]/_current_nbIndsS[FEM][OFFSx] : my_NAN;}
double MetapopSH::getAdultSexRatio   () {return _current_nbIndsS[FEM][ADLTx] ? (double)_current_nbIndsS[MAL][ADLTx]/_current_nbIndsS[FEM][ADLTx] : my_NAN;}

double MetapopSH::getTotNbInd        (const age_t& AGE){return get_current_nbInds(AGE);}
double MetapopSH::getTotNbFem        (const age_t& AGE){return get_current_nbIndsS(AGE,FEM);}
double MetapopSH::getTotNbMal        (const age_t& AGE){return get_current_nbIndsS(AGE,MAL);}

double MetapopSH::getMeanNbInd       (const age_t& AGE){return get_current_nbPops(AGE) ? (double)get_current_nbInds(AGE)/get_current_nbPops(AGE) : my_NAN;}
double MetapopSH::getMeanNbFem       (const age_t& AGE){return get_current_nbPops(AGE) ? (double)get_current_nbIndsS(AGE,FEM)/get_current_nbPops(AGE) : my_NAN;}
double MetapopSH::getMeanNbMal       (const age_t& AGE){return get_current_nbPops(AGE) ? (double)get_current_nbIndsS(AGE,MAL)/get_current_nbPops(AGE) : my_NAN;}

double MetapopSH::getNbInd           (unsigned int i, const age_t& AGE){return _sample_pops[i]->size(AGE);}
double MetapopSH::getNbFem           (unsigned int i, const age_t& AGE){return _sample_pops[i]->size(FEM, AGE);}
double MetapopSH::getNbMal           (unsigned int i, const age_t& AGE){return _sample_pops[i]->size(MAL, AGE);}

double MetapopSH::getNbPops          (const age_t& AGE){return get_current_nbPops(AGE);}
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** Fitness ********

// ----------------------------------------------------------------------------------------
// setMeanAndVarFitness
// ----------------------------------------------------------------------------------------
void
MetapopSH::setMeanAndVar_W(){
  // check if the table has already been computed
  if(already_computed(_computed[3])) return;

  if(!_meanW) _meanW = new double[_sample_pops_size];
  if(!_varW)  _varW  = new double[_sample_pops_size];

  // for each population
  for(unsigned int i = 0; i < _sample_pops_size; ++i) {
    setMeanAndVar_W_ofPatch(_sample_pops[i], i);
  }
}

// ----------------------------------------------------------------------------------------
// setMeanAndVarFitness
// ----------------------------------------------------------------------------------------
void
MetapopSH::setMeanAndVar_W_ofPatch(Patch* crnt_patch, const unsigned int& i){
  unsigned int sizeF = crnt_patch->size(FEM, ADLTx),
               sizeM = crnt_patch->size(MAL, ADLTx),
               size  = sizeF + sizeM;

  // if the patch is empty or the fitness is not yet computed -> stop
  if(!((sizeF && crnt_patch->get(FEM, ADLTx, 0)->getFitness() != my_NAN)
    ||(sizeM && crnt_patch->get(MAL, ADLTx, 0)->getFitness() != my_NAN))){
    _meanW[i] = _varW[i] = my_NAN;
    return;
  }

  unsigned int f, m;
  double* array = new double[size];

  for(f = 0; f < sizeF; ++f) {
		array[f] = crnt_patch->get(FEM, ADLTx, f)->getFitness();
  }
  for(m = 0; m < sizeM; ++m, ++f) {
    array[f] = crnt_patch->get(MAL, ADLTx, m)->getFitness();
  }

  // compute mean and var
  _meanW[i] = ARRAY::mean(array, size);
  _varW[i]  = ARRAY::var(array, size, _meanW[i]);
  delete[] array;
}

// ----------------------------------------------------------------------------------------
// getVwB    fitness variance between patches
// ----------------------------------------------------------------------------------------
double
MetapopSH::getVwB(){
  setMeanAndVar_W();
  if(_sample_pops_size<2) return my_NAN;        // at least two populations are needed
  return ARRAY::var(_meanW, _sample_pops_size);
}

// ----------------------------------------------------------------------------------------
// getVwB   fitness variance within patches
// ----------------------------------------------------------------------------------------
double
MetapopSH::getVwW(){
  setMeanAndVar_W();
  return ARRAY::mean(_varW, _sample_pops_size);
}


