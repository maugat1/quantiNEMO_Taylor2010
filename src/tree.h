/** @file tree.h
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

#ifndef treeH
#define treeH

#include <iostream>
using namespace std;

#include "node.h"
#include "output.h"


/**A class to compute phenotypes with epitasis, is an aggregate of Node.*/
template <class T>
class Tree {

private:
  /**The depth of the tree as determined by the number of locus of the trait*/
  unsigned int _nb_locus;
  /**The number of branches per node, determined by the number of possible genotypes at a locus.*/
  unsigned int _nb_branches;
  /**The number of allelic states of the trait.*/
  unsigned int _nb_all;
  /**The first Node, root of the tree.*/
  Node _root;
  /**A  nb_all x nb_all matrix used to convert a locus genotype into a unique value (the coordinate of that locus).*/
  unsigned int** _mapper;
  /**a _nb_branches*ploidy matrix to get back the allelic combination from the _coord.*/
  T** _un_mapper;
  /**A table of length = number of locus, the coordinate of the genotype in the tree after its mapping.*/
  unsigned int * _coord;

public:

  Tree (unsigned int nbloc, unsigned int nball);
  ~Tree ();
  /**Gives the phenotype of the genotype given in argument.*/
  double get_value(T** genotype);
  double get_first(T** genotype);
  double get_next (T** genotype);
  void   set_value(T** genotype, double value);

};

#endif //TREE_H

