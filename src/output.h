/** @file output.h
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

#ifndef outputH
#define outputH

#include <cstdlib>

/** end of line character */
#ifdef _MAC_
  #define EOL '\r'      // \r (only macs)
#else
  #define EOL '\n'      // \r\n or \n  (PC and UNIX)
#endif


#include <assert.h>
#include <time.h>

#include "tstring.h"
#include <cmath>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <map>
#include "tarray.h"
#include "types.h"

using namespace std;

extern void message (const char* message, ...);
extern void warning (const char* message, ...);
extern void error   (const char* message, ...);
extern void fatal   (const char* message, ...);

extern string getElapsedTime(clock_t time);



//------------------------------------------------------------------------------
/** returns the rounded value with the specified number of digits */
inline int my_round(double value)
{
    return (int) (value<0 ? value-0.5 : value+0.5);
}

inline double my_round(double value, int digits)
{
    if(!digits) return (int) (value<0 ? value-0.5 : value+0.5);
    return floor(value * pow((double)10, digits) + 0.5) * pow((double)10, -digits);
}

//------------------------------------------------------------------------------
/** similar to sample in R:
	* returns nb UNIQUE integer numbers between 0 and max-1
	* the returned array must be deleted
	* the returned array is sorted!!!
	*/
unsigned int* sample(const unsigned int& nb, const unsigned int& max);

/** similar to sample in R WITHOUT replacement:
	* input:  -array to sample from (array is changed)
	*         -size: total size of the array
	*         -nb: points to sample
	*         -order: does the order of the samples matters? (true=yes)
	* output: -nothing, but the array is modified
						 the sampled items are the "nb" first elements of the array
	*
template <typename T>
void sample(T* array, const unsigned int& size, const unsigned int& nb, bool order=true){
	if(size<nb) nb=size;
	if(!nb) return;                                 // as it could cause loop problems

	if(order || nb<size/2){                         // draw the "good" ones (start left)
		unsigned int pos, i;
		for(i = 0; i< nb; ++i){
			pos = SimRunner::r.Uniform(size-i);         // get the pos
			if(pos+i != i) ARRAY::swap(array, i, pos+i);
		}
	}
	else{                                           // draw the bad ones (start right)
		unsigned int pos, i;
		for(i = size-1; i>=nb; --i){
			pos = SimRunner::r.Uniform(i+1);         		// get the pos
			if(pos != i) ARRAY::swap(array, i, pos);
		}
	}
}
*/
//------------------------------------------------------------------------------
/** returns the new population size assuming logistic growth limited by the carrying capacity */
inline unsigned int logisticGrowth(const double& r, const unsigned int& K, const unsigned int& N)
{
   if(!K) return 0;
   return my_round(N + r*N*(1.0 - (double)N/K));
}

//------------------------------------------------------------------------------
/** returns the new population size assuming logistic growth limited by the carrying capacity
  * The Beverton-Hold model is implemented:  newN[i]<-N*(1+r)*K/((1+r)*N-N+K)
  * Beverton RJH and Holt SJ (1957). On the Dynamics of Exploited Fish Populations.
  * The model (repectivel r) was modified such as r=0 results in a stable population size
  */
inline double beverton_hold(const double& r, const unsigned int& K, const unsigned int& N)
{
  if(!K) return 0;
	return N*K*(1+r)/(N*(1+r)-N+K);
}


//------------------------------------------------------------------------------
/** that is an implementation of the general logistic curve
  * (Richards, F.J. 1959 A flexible growth function for empirical use. J. Exp. Bot. 10: 290--300.)
  * The curve is defined by 5 parameters.
  * Do to size problems a slope of more than (+/-)1e4 is regarded as instantenious change from min to max
  */
inline double generalLogisticCurve(const double& x,     // time
                                   const double& min,   // the lower asymptote
                                   const double& max,   // the upper asymptote
                                   const double& max_r, // the time of maximum growth
                                   const double& r,     // the growth rate
                                   const double& s)     // affects near which asymptote maximum growth occurs
{

  // check if the slope exceeds the limits
  if(r>=1e4){         // slope is too big            => use instantenious change
    if(x<max_r)      return min;
    else if(x>max_r) return max;
    else             return min+((max-min)/pow(1+s,1.0/s));   // if x == max_r
	}

	if(r<=-1e4){       // slope is too big (negative) => use instantenious change
		if(x<max_r)      return max;
    else if(x>max_r) return min;
    else             return min+((max-min)/pow(1+s,1.0/s));   // if x == max_r
  }

  // the slope is small enough to use the generalized logistic function
  assert(s);
  return min+((max-min)/pow(1+(s*(double)exp((long double)(max_r-x)*r)),1.0/s));
}

// ----------------------------------------------------------------------------------------
// read_Fstat_file
// ----------------------------------------------------------------------------------------
/** this function reads any FSTAT file, a vector is returned which contains each individual.
  * Each individual is determined by a double unsigned int array:
  *  - first the loci are stored: v[ind][loc][all]
  *  - Then the suplement info is stored if available:
  *       - v[ind][_nb_locus]  [0] = patch  // startign at 0
  *       - v[ind][_nb_locus+1][0] = age    // off: 0, adlt: 2
  *       - v[ind][_nb_locus+2][0] = sex    // mal: 0, fem: 1
  *       - v[ind][_nb_locus+3][0] = id     // id of the current patch starting at 1
  *       - v[ind][_nb_locus+3][1] = id     // id of the natal patch starting at 1
  *       - v[ind][_nb_locus+4][0] = id_mom // id of mother of the current patch starting at 1
  *       - v[ind][_nb_locus+4][1] = id_mom // id of mother of the natal patch starting at 1
  *       - v[ind][_nb_locus+5][0] = id_dad // id of father of the current patch starting at 1
  *       - v[ind][_nb_locus+5][1] = id_dad // id of father of the natal patch starting at 1
  * Note, that the arrays of the vector has to be deleted
  */
vector<unsigned int**>* read_Fstat_file(string filename, string dir);
//------------------------------------------------------------------------------
#endif
