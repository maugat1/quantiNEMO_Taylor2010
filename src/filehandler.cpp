/** @file filehandler.cpp
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


#ifdef _BORLAND_
	#include <dir.h>
#else
   #include <sys/types.h>
   #include <sys/stat.h>
   #include <errno.h>
#endif

#include <cstdlib>
#include <cmath>
#include <dirent.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <dirent.h>
#include "output.h"
#include "filehandler.h"
#include "lce_misc.h"
#include "metapop.h"

using namespace std;



string FileHandler::_base_path;
// ----------------------------------------------------------------------------------------
// init
// ----------------------------------------------------------------------------------------
bool FileHandler::init ()
{
	Metapop *pop = get_pop_ptr();

  //check if the occurrences exceed the pop parameters:
  if(_GenerationOccurrence > pop->getGenerations()){
    _GenerationOccurrence = pop->getGenerations();
  }

  if(_ReplicateOccurrence > pop->getReplicates()){
    _ReplicateOccurrence = pop->getReplicates();
  }

	// open the file
	string filename = get_path() + get_service()->getGenerationReplicateFileName() + get_extension();

	//check if the basefilename is already used on disk:
  ifstream ifExist;
  ifExist.setstate(ios::failbit);

	ifExist.open(filename.c_str(),ios::in);

  if(ifExist.is_open() && !_service->getOverwriteFiles()) {
		warning("a file named \"%s\" exists already in folder \"%s\" !!\n",(get_service()->getGenerationReplicateFileName() + get_extension()).c_str(),get_path().c_str());
		ifExist.close();
    return false;
  }
	return true;
}

// ----------------------------------------------------------------------------------------
// set_path if it does not exist
// ----------------------------------------------------------------------------------------
void FileHandler::check_path (string& path)
{
	if(path.size()){
		if(path[path.length()-1] != SEP) path += SEP;

		// for each folder
		string::size_type cur, next;
		string curFolder;
		cur=0;
		next=path.find(SEP, cur);
		while(next != string::npos){
			curFolder = path.substr(0, next+1);
			DIR *dirname=opendir(curFolder.c_str());
			if(!dirname) {      // if folder does not exist
				#ifdef _BORLAND_
				 if((mkdir(curFolder.c_str())) == -1){
				#else
				 string cmd = "mkdir -p \"" + curFolder + "\"";      // "" are needed if a folder name contains spaces
				 if ( system( cmd.c_str() ) < 0 ){
				#endif

					error("could not create directory \"%s\", saving in wd.\n",path.c_str());
          path = "";
          return;
				}
			}
      else closedir(dirname);

      // go to next folder
      cur = next+1;
      next=path.find(SEP, cur);
		}
  }
}

// ----------------------------------------------------------------------------------------
// get_filename
// ----------------------------------------------------------------------------------------
string& FileHandler::get_filename (){
  _current_filename = _path + _service->getReplicateFileName() + _extension;
  return _current_filename;
}

// ----------------------------------------------------------------------------------------
// update
// ----------------------------------------------------------------------------------------
void FileHandler::update()
{
	Metapop* _pop = get_pop_ptr();
	_current_replicate = _pop->getCurrentReplicate();
	_current_generation = _pop->getCurrentGeneration();

	if(((_isReplicatePeriodic && !(_current_replicate % _ReplicateOccurrence))
														|| (_current_replicate == _ReplicateOccurrence))
		&&((_isGenerationPeriodic && !(_current_generation % _GenerationOccurrence))
														|| (_current_generation == _GenerationOccurrence))){
		FHwrite();
	}
}

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
/** Writes genotype in a FSTAT-like format, the age class and the sex are added after the pop id.*/
void FileHandler::FHwrite ()
{
  int patchNbr = get_pop_ptr()->getPatchNbr();
  Patch* current_patch;
  int t;

  // get the number of occupied patches
  int nb_patches_alive = 0;
  for (int i=0; i<patchNbr; ++i) {
    if(!get_pop_ptr()->getPatch(i)->get_isExtinct()) ++nb_patches_alive;
	}
	if(!nb_patches_alive) return;  // if all patches are extinct

  // get the number of digits to visualize
  int max_allele=0;    // find the max number of alleles across all traits
  for(t = 0; t<_nb_trait; ++t){
    if(max_allele < _trait[t]->get_nb_allele()) max_allele = _trait[t]->get_nb_allele();
  }


  // open the file
  string filename = get_path() + this->get_service()->getGenerationReplicateFileName() + get_extension();
  #ifdef _DEBUG
   message("FileHandler::FHwrite (%s)\n",filename.c_str());
  #endif
  ofstream FILE (filename.c_str(), ios::out);
  if(!FILE) fatal("Could not open FSTAT output file '%s'!\n", filename.c_str());

  // write the heading line
  int position = STRING::getNbDigits(max_allele);
  int tot_nb_loci=0;
  for(t=0; t<_nb_trait; ++t){
    tot_nb_loci += _trait[t]->get_nb_locus();
  }
  FILE << nb_patches_alive << " " << tot_nb_loci << " " << max_allele << " " << position << "\n";

  // write the names of the loci
  for(t=0; t<_nb_trait; ++t){
    for (int l = 0; l < _trait[t]->get_nb_locus(); ++l){
      FILE << "trait-" << (t+1) << "_locus-" << (l+1) << "\n";
    }
  }

  // write all individuals of all patches
  int nbPatchDigits = STRING::getNbDigits(patchNbr);
  unsigned int a, mask;
  for (int i=0; i<patchNbr; ++i) {
    current_patch = get_pop_ptr()->getPatch(i);
    for(a = 0, mask=1; a < NB_AGE_CLASSES; ++a, mask <<= 1) {
      if(mask & _age){
        if(_sex != 2) FHwrite(static_cast<age_idx>(a), FEM, FILE, current_patch, i, nbPatchDigits, position);
        if(_sex != 1) FHwrite(static_cast<age_idx>(a), MAL, FILE, current_patch, i, nbPatchDigits, position);
      }
    }
  }

  FILE.close();

  // call an external program and pass the name of this file
  if(!_script.empty()) execute_script(_script, filename);
}

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
void FileHandler::FHwrite (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
      Patch* current_patch, const int& patch_id, const int& nbPatchDigit, const int& position){

  unsigned char** seq;
  Individual *ind;
  int ploidy, nb_locus;

  for (unsigned int j = 0, nbInd = current_patch->size(cur_sex, cur_age); j < nbInd; ++j) {
    FILE << setfill('0') << setw(nbPatchDigit) << (patch_id+1) << setfill(' ') << " ";
	  ind = current_patch->get(cur_sex, cur_age, j);
    for(int t=0; t<_nb_trait; ++t){   // for multiple instanciations of a trait
      ploidy   = _trait[t]->get_ploidy();
      nb_locus = _trait[t]->get_nb_locus();
	    seq = (unsigned char**)ind->getTrait(_TTidx[t])->get_sequence();

      for(int k = 0; k < nb_locus; ++k) {
        for (int l = 0; l < ploidy; ++l) {
          FILE.fill('0');
          FILE.width(position);
          FILE<<(unsigned int)(seq[k][l]+1);
        }
        FILE<<" ";
      }
    }
	  if(_fstat_choice==2){
			write_individual_info_to_stream(FILE, ind, cur_age, cur_sex, ' ');
		}
		FILE << "\n";
	}
}

// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
/** writes the genotype of the individual to the stream (for multiple instances of the same trait)*/
void FileHandler::FHwriteIndividual2Stream(ofstream& FILE, Individual* ind, const int& ploidy, const int& position){
	int nb_locus;
	unsigned char** seq;
  for(int t=0; t<_nb_trait; ++t){
    nb_locus = _trait[t]->get_nb_locus();
		seq = (unsigned char**)ind->getTrait(_TTidx[t])->get_sequence();
		for(int k = 0; k < nb_locus; ++k) {
			for (int l = 0; l < ploidy; ++l) {
				FILE.fill('0');
				FILE.width(position);
				FILE<<(unsigned int)(seq[k][l]+1);
			}
			FILE<<" ";
		}
	}
}


// ----------------------------------------------------------------------------------------
// FHwrite
// ----------------------------------------------------------------------------------------
/** writesw the information of the individual to the stream */
void FileHandler::write_individual_info_to_stream(ostream& FILE, Individual* ind,
																									const age_idx& cur_age, const sex_t& cur_sex, char sep){
	FILE << cur_age+1 << sep
			 << cur_sex   << sep
			 << ind->getID() << sep;
	if(ind->getMotherID() != "")    FILE << ind->getMotherID() << sep;
	else                            FILE << "NaN" << sep;
	if(ind->getFatherID() != "")    FILE << ind->getFatherID();
	else                            FILE << "NaN";
	if(ind->getFitness() == my_NAN) FILE << sep << "NaN";
	else                            FILE << sep << ind->getFitness();
}

// ----------------------------------------------------------------------------------------
/** launch the script at the end of all simulations */
void FileHandler::execute_script(string scriptname, string filename)
{
  ifstream script(scriptname.c_str(),ios::in);
  string cmd;

  if(!script.is_open()) {
	  error("could not open script %s!\n", scriptname.c_str());
	  return;
  }


  #ifdef _BORLAND_
    cmd = scriptname + " " + filename;
  #else
    cmd = "./" + scriptname + " " + filename;
  #endif

  #ifdef _DEBUG
   message("Executing shell script \"%s\" ",cmd.c_str());
  #endif
  fflush(stdout);

  if(system(cmd.c_str()) < 0){
	  error("Execution of `sh %s' failed: %s\n",cmd.c_str(),strerror(errno));
	  return;
  }

  #ifdef _DEBUG
   message("...done\n");
  #endif
}
// ----------------------------------------------------------------------------------------


