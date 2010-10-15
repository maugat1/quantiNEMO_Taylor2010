/** @file output.cpp
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
	#include <stdlib.h>
  #include <process.h>
#include <stdio.h>
#include <errno.h>

#else
  #include <unistd.h>
  #include <sys/types.h>
 // #include <sys/wait.h>
  #include <errno.h>
#endif


#include <stdarg.h>
#include "output.h"
#include <fstream>
#include "simulation.h"
using namespace std;

void message (const char* message, ...)
{
	va_list ap;
  va_start(ap, message);
	vprintf(message,ap);
	va_end(ap);
}

void warning (const char* message, ...)
{
	va_list ap;
  va_start(ap, message);
	printf("\n***WARNING*** ");
	vprintf(message,ap);
	va_end(ap);
}

void error (const char* message, ...)
{
	va_list ap;
  va_start(ap, message);
	fprintf(stderr,"\n***ERROR*** ");
  vfprintf(stderr,message,ap);
	va_end(ap);
}

void fatal (const char* message, ...)
{
	va_list ap;
  va_start(ap, message);
	printf("\n***FATAL-ERROR*** ");
  vfprintf(stderr,message,ap);
	va_end(ap);

	throw message;
}

//------------------------------------------------------------------------------
string getElapsedTime(clock_t time)
{
  int e_time = time / CLOCKS_PER_SEC;
  int hour = e_time / 3600;
  int min  = ((e_time % 3600) / 60);
  int sec  = (e_time % 3600) % 60;

  ostringstream o;
  if(hour<10) o << setfill('0') << setw(2) << hour;
  else o << hour;
  o << ":" << setfill('0') << setw(2) << min;
  o << ":" << setfill('0') << setw(2) << sec;

  return o.str();
}

//------------------------------------------------------------------------------
/** similar to sample in R:
  * returns nb UNIQUE integer numbers between 0 and max-1
  * the returned array must be deleted
  * the returned array is sorted!!!
  */
unsigned int* sample(const unsigned int& nb, const unsigned int& max)
{
  assert(nb<=max);

  // create a temp vector (0, 1, 2, ... max-1)
  vector<unsigned int> vec;
  for(unsigned int i = 0; i< max; ++i){
    vec.push_back(i);
  }

  unsigned int* array = new unsigned int[nb];     // array to return
  if(nb<max/2){                                   // draw the "good" ones
    unsigned int pos;
    for(unsigned int i = 0; i< nb; ++i){
      pos = SimRunner::r.Uniform(max-i);          // get the pos
      array[i] = vec[pos];                        // store the value to the array
      vec.erase(vec.begin()+pos);                 // remove the element
    }
  }
  else{                                           // remove the bad ones
    for(unsigned int i = 0; i< max-nb; ++i){
      vec.erase(vec.begin()+SimRunner::r.Uniform(max-i));
    }
    for(unsigned int i = 0; i< nb; ++i){          // copy the vector to the array
      array[i] = vec[i];
    }
  }
  return array;
}

// ----------------------------------------------------------------------------------------
// read_Fstat_file
// ----------------------------------------------------------------------------------------
/** this function reads any FSTAT file, a vector is returned which contains each individual.
  * Each individual is determined by a double unsigned int array:
  *  - first the suplement info is stored if available:
  *       - v[ind][0][0] = patch       // starting at 0 (compulsory)
  *       - v[ind][0][1] = nb_loci     // number of locus found int eh file (compulsory)
  *       - v[ind][1][0] = age         // off: 0, adlt: 2
  *       - v[ind][2][0] = sex         // mal: 0, fem: 1
  *       - v[ind][3][0] = id          // id of the current patch starting at 1
  *       - v[ind][3][1] = natal patch // id of the natal patch starting at 1
  *       - v[ind][4][0] = id_mom      // id of mother of the current patch starting at 1
  *       - v[ind][4][1] = moms patch  // id of mother of the natal patch starting at 1
  *       - v[ind][5][0] = id_dad      // id of father of the current patch starting at 1
  *       - v[ind][5][1] = dads patch  // id of father of the natal patch starting at 1
  *  - then the loci are stored: v[ind][loc+][all] (compulsory
  * Note, that the arrays of the vector has to be deleted
  */
vector<unsigned int**>*
read_Fstat_file(string filename, string dir){
  vector<unsigned int**>* v = new vector<unsigned int**>();     // vector containing the individuals
  try{
    if(!STRING::file_exists(filename)){
      filename = dir + filename;
			if(!STRING::file_exists(filename))fatal("FSTAT file '%s' could not be opened!\n", filename.c_str());
    }
    ifstream FILE(filename.c_str());

    // read the heading line and test if the settings correspond to the ini file
    int nbPop, nbLoci, maxAllele, nbDigit, i;
    string text;
    FILE >> nbPop >> nbLoci >> maxAllele >> nbDigit >> ws;
    if(maxAllele > maxAllele) fatal("FSTAT file: the number of alleles does not correspond to the settings!\n");


    // skip all locus names
    for(i = 0; i < nbLoci; ++i) {
      FILE.ignore(2048, '\n');  // remove line by line
    }
    // read individual by individual
    unsigned int** curInd;
    char c;
    unsigned int val;
    string text1;
    while(FILE.good() && !FILE.eof()){      // for each individual
      curInd = ARRAY::new_2D(nbLoci+6, 2, (unsigned int) my_NAN);

      FILE >> val;                               // get the pop
      curInd[0][0] = val-1;
      curInd[0][1] = nbLoci;

      // get the alleles
      for(i=0; i<nbLoci; ++i){
        FILE >> text;
        if((int)text.length() != 2*nbDigit) fatal("FSTAT file: The number of digits is wrong (line %i, locus %i)!\n", v->size()+1, i+1);
        curInd[i+6][0] = STRING::str2int<unsigned int>(text.substr(0,nbDigit))-1;
        curInd[i+6][1] = STRING::str2int<unsigned int>(text.substr(nbDigit,nbDigit))-1;
      }

      // check if the suplement columns are present
      do{
        FILE.get(c);
      }while(isspace(c) && c != EOL && c != FILE.eof());

      if(c != EOL && c != FILE.eof()){    // the suplement columns are present
        FILE.putback(c);
				FILE >> val;  curInd[1][0] = val-1;   // age
				FILE >>       curInd[2][0];                      // sex
        FILE >> text; curInd[3][0] = STRING::str2int<unsigned int>(text.substr(0,text.find('_')));               // id
                      curInd[3][1] = STRING::str2int<unsigned int>(text.substr(text.find('_')+1, text.length()));
        FILE >> text; curInd[4][0] = STRING::str2int<unsigned int>(text.substr(0,text.find('_')));               // id_mother
                      curInd[4][1] = STRING::str2int<unsigned int>(text.substr(text.find('_')+1, text.length()));
        FILE >> text; curInd[5][0] = STRING::str2int<unsigned int>(text.substr(0,text.find('_')));               // id_father
                      curInd[5][1] = STRING::str2int<unsigned int>(text.substr(text.find('_')+1, text.length()));
				FILE >> text;		// remove the fitness: fitness is not considered
				FILE >> ws;
      }
      v->push_back(curInd);      // add the individual to the vector
      FILE >> ws;
    }

    FILE.close();
	}catch(...) {throw("Problems with reading the FSTAT file!");}

	return v;
}






