/** @file types.h
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

#ifndef typesH
#define typesH

// use the assert macro only when debugging
#ifndef _DEBUG
  #define NDEBUG  // deactivate assert
#endif

//#define _DEBUG

#include <new>
#include <string.h>
#include <fstream>
#include <limits>
using namespace std;

/**Sex types, males are always 0 and females 1!!**/
typedef enum {
  MAL=0, FEM=1
} sex_t;

/**Array index of the age classes in the patch sizes and containers arrays.**/
typedef enum {
  OFFSx=0, ADLTx=1
}age_idx;

/**Age class flags.*/  // binary representation
typedef unsigned int age_t;
#define NONE 0        // 00000000   /** No age flag.*/
#define OFFSPRG 1     // 00000001   /** Offspring age class flag.*/
#define ADULTS 2      // 00000010   /** Adults age class flag (breeders).*/
#define ALL 3         // 00000011   /** Adults and offspring flag */


#define my_NAN 99999

#define NB_AGE_CLASSES 2

extern char SEP;

/**Ordering type used to record statistics in the StatRecorders.**/
typedef enum {
  FLAT = 2,     // values are stored for each generation and replicate
  GEN  = 4,     // values are stored for each generation: replicate values are added
  RPL  = 6      // values are stored for each replicate: generation values are added
}st_order;

/**Trait types**/
typedef string trait_t;
/**Max number of characters in the trait's type descriptor.*/
#define TRAIT_T_MAX 5
#define DELE "delet"
#define DISP "disp"
#define FDISP "fdisp"
#define MDISP "mdisp"
#define NTRL "ntrl"
#define DQUANT "quanti"

/**Param's types**/
typedef enum {
	DBL,INT2,STR,MAT,DIST,MAT_VAR,
	INT_MAT,DBL_MAT,STR_MAT
}param_t;

#endif

