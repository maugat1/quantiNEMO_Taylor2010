/** @file ttree.h
*
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
//---------------------------------------------------------------------------

#ifndef ttreeH
#define ttreeH

#include "output.h"
//---------------------------------------------------------------------------
template <class T>
class TNode2{
private:
  T**          _array;
  T**          _end_pos;
  T**          _cur_pos;
  unsigned int _size;
  T            _default;

public:
  TNode2(unsigned int s, T d){
    _array = new T*[s];
    _size = s;
    _default = d;
    _end_pos = _array + _size;
    _cur_pos = _array;
    for(; _cur_pos != _end_pos; ++_cur_pos){
      *_cur_pos = NULL;
    }
  }

  ~TNode2(){
    _cur_pos = _array;
    for(; _cur_pos != _end_pos; ++_cur_pos){
      if(*_cur_pos) delete (*_cur_pos);
    }
    delete[] _array;
  }

  void add(unsigned int a2, T val){
    if(!_array[a2]){
      _array[a2] = new T;
      *_array[a2] = _default;
    }
    *_array[a2] += val;
  }

  void set(unsigned int a2, T val){
    if(!_array[a2]) new T;
    *_array[a2] = val;
  }

  T& get(unsigned int a2){
    if(!_array[a2]){
      _array[a2] = new T;
      *_array[a2] = _default;
    }
    return *_array[a2];
  }

  // returns a pointer to the next element and changes a2 to the correct value.
  // returns NULL if there are no more elements
  T* next(unsigned int& a2){
    if(a2>=_size){a2=0; return NULL;}
    while(!_array[a2]){
      ++a2;
      if(a2>=_size){a2 = 0; return NULL;}
    }
    return _array[a2];
  }
};

//---------------------------------------------------------------------------
template <class T>
class TNode1{
private:
  TNode2<T>**  _array;
  TNode2<T>**  _end_pos;
  TNode2<T>**  _cur_pos;
  unsigned int _size;
  T            _default;

public:
  TNode1(unsigned int s, T d){
    _array = new TNode2<T>*[s];
    _size = s;
    _default = d;
    _end_pos = _array + _size;
    _cur_pos = _array;
    for(; _cur_pos != _end_pos; ++_cur_pos){
      *_cur_pos = NULL;
    }
  }

  ~TNode1(){
    _cur_pos = _array;
    for(; _cur_pos != _end_pos; ++_cur_pos){
      if(*_cur_pos) delete (*_cur_pos);
    }
    delete[] _array;
  }

  void add(unsigned int a1, unsigned int a2, T val){
    if(!_array[a1]) _array[a1] = new TNode2<T>(_size, _default);
    _array[a1]->add(a2,val);
  }

  void set(unsigned int a1, unsigned int a2, T val){
    if(!_array[a1]) _array[a1] = new TNode2<T>(_size, _default);
    _array[a1]->set(a2,val);
  }

  T& get(unsigned int a1, unsigned int a2){
    if(!_array[a1]) _array[a1] = new TNode2<T>(_size, _default);
    return _array[a1]->get(a2);
  }

  // returns a pointer to the next element and changes a1 and a2 to the correct value.
  // returns NULL if there are no more elements
  T* next(unsigned int& a1, unsigned int& a2){
    do{
      if(_array[a1]){
        T* elem = _array[a1]->next(a2);
        if(elem) return elem;
      }
      do{
        ++a1;
        if(a1>=_size){a2=0; a1=0; return NULL;}
      }while(!_array[a1]);
    }while(1);
    return NULL;
  }
};

template <class T>
class TTree{
private:
  TNode1<T>**  _array;
  TNode1<T>**  _end_pos;
  TNode1<T>**  _cur_pos;
  unsigned int _nb_loci;
  unsigned int _nb_allele;
  T            _default;

  T* _find_next(unsigned int& l, unsigned int& a1, unsigned int& a2){
    do{
      if(_array[l]){
        T* elem = _array[l]->next(a1, a2);
        if(elem) return elem;
      }
      do{
        ++l;
        if(l>=_nb_loci){l=0; return NULL;}
      }while(!_array[l]);
    }while(1);
    return NULL;
  }

public:
  TTree(unsigned int l, unsigned int a, T d=0){
    _array = new TNode1<T>*[l];
    _nb_loci = l;
    _nb_allele = a;
    _default = d;
    _end_pos = _array + l;
    _cur_pos = _array;
    for(; _cur_pos != _end_pos; ++_cur_pos){
      *_cur_pos = NULL;
    }
  }

  ~TTree(){
    _cur_pos = _array;
    for(; _cur_pos != _end_pos; ++_cur_pos){
      if(*_cur_pos) delete (*_cur_pos);
    }
    delete[] _array;
  }

  void add(unsigned int l, unsigned int a1, unsigned int a2, T val){
    if(!_array[l]) _array[l] = new TNode1<T>(_nb_allele, _default);
    _array[l]->add(a1,a2,val);
  }

  void set(unsigned int l, unsigned int a1, unsigned int a2, T val){
    if(!_array[l]) _array[l] = new TNode1<T>(_nb_allele, _default);
    _array[l]->set(a1,a2,val);
  }

	// get the value of locus l and alleles a1 and a2
  T& get(unsigned int l, unsigned int a1, unsigned int a2){
    if(!_array[l]) _array[l] = new TNode1<T>(_nb_allele, _default);
    return _array[l]->get(a1,a2);
  }

  // returns a pointer to the next element and changes l, a1 and a2 to the correct value.
  // returns NULL if there are no more elements
  T* next(unsigned int& l, unsigned int& a1, unsigned int& a2){
    ++a2; // increment the last allele
    return _find_next(l, a1, a2);
  }

  // returns a pointer to the first element and changes l, a1 and a2 to the correct value.
  // returns NULL if the tree is empty
  T* first(unsigned int& l, unsigned int& a1, unsigned int& a2){
    l = a1 = a2 = 0;
    return _find_next(l, a1, a2);
  }
};

#endif
