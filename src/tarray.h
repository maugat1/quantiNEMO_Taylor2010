/** @file tarray.h
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

#ifndef tarrayH
#define tarrayH
//---------------------------------------------------------------------------

#include "output.h"
#include "types.h"


class ARRAY{
private:

public:
//------------------------------------------------------------------------------
// array functions
template <typename T>
static void delete_3D(T*** &array, int l, int a){
  try{
    if(array){
      int i, j;
      for(i=0; i<l; ++i){
        if(array[i]){
          for(j=0; j<a; ++j){
            if(array[i][j]) delete[] array[i][j];
          }
          delete[] array[i];
        }
      }
      delete[] array;
      array = NULL;
    }
  }catch(...) {throw ("Out of memory!\n");}
}

template <typename T>
static void create_3D(T*** &array, int l, int a, int t){    // with initialization
  try{
    if(array) delete_3D<T>(array, l, a);
    int i, j;
    array = new T** [l];
    for(i=0; i < l; ++i){
      array[i] = new T*[a];
      for(j=0; j < a; ++j){
        array[i][j] = new T[t];
      }
    }
  }catch(...) {throw ("Out of memory!\n");}
}

template <typename T>
static void create_3D(T*** &array, int l, int a, int t, T init){    // with initialization
  try{
    if(array) delete_3D<T>(array, l, a);
    int i, j, k;
    array = new T** [l];
    for(i=0; i < l; ++i){
      array[i] = new T*[a];
      for(j=0; j < a; ++j){
        array[i][j] = new T[t];
        for(k = 0; k < t; ++k){
          array[i][j][k] = init;       //init values
        }
      }
    }
  }catch(...) {throw ("Out of memory!\n");}
}

template <typename T>
static T*** new_3D(int l, int a, int t){
  T*** array = NULL;
  create_3D(array, l, a, t);
  return array;
}

template <typename T>
static T*** new_3D(int l, int a, int t, T init){
  T*** array = NULL;
  create_3D(array, l, a, t, init);
  return array;
}

template <typename T>
static void reset_3D(T*** &array, int l, int a, int t, T init){
  assert(array);
  int i, j, k;
  for(i=0; i < l; ++i){
    assert(array[i]);
    for(j=0; j < a; ++j){
      assert(array[i][j]);
      for(k = 0; k < t; ++k){
        array[i][j][k] = init;       //init values
      }
    }
  }
}

// 2Darray functions
template <typename T>
static void delete_2D(T** &array, int l){
  if(array){
    int i;
    for(i=0; i<l; ++i){
      if(array[i]) delete[] array[i];
    }
    delete[] array;
    array = NULL;
  }
}

template <typename T>
static void create_2D(T** &array, int l, int a){
  try{
    if(array) delete_2D<T>(array, l);
    int i;
    array = new T* [l];
    for(i=0; i < l; ++i){
      array[i] = new T[a];
    }
  }catch(...) {throw ("Out of memory!\n");}
}

template <typename T>
static void create_2D(T** &array, int l, int a, T init){   // with initialization
  try{
    if(array) delete_2D<T>(array, l);
    int i, j;
    array = new T* [l];
    for(i=0; i < l; ++i){
      array[i] = new T[a];
      for(j=0; j < a; ++j){
        array[i][j] = init;       //init values
      }
    }
  }catch(...) {throw ("Out of memory!\n");}
}

template <typename T>
static T** new_2D(int l, int a){
  T** array = NULL;
  create_2D(array, l, a);
  return array;
}
template <typename T>
static T** new_2D(int l, int a, T init){
  T** array = NULL;
  create_2D(array, l, a, init);
  return array;
}


template <typename T>
static void reset_2D(T** &array, int l, int a, T init){
  assert(array);
  int i, j;
  for(i=0; i < l; ++i){
    assert(array[i]);
    for(j=0; j < a; ++j){
      array[i][j] = init;       //init values
    }
  }
}

// 1Darray functions
template <typename T>
static void delete_1D(T* &array){
  if(array){
    delete[] array;
    array = NULL;
  }
}

template <typename T>
static void create_1D(T* &array, int l){
  try{
    if(array) delete_1D<T>(array);
    array = new T [l];
  }catch(...) {throw ("Out of memory!\n");}
}

template <typename T>
static void create_1D(T* &array, int l, T init){    // with initialization
  try{
    if(array) delete_1D<T>(array);
    int i;
    if(l) array = new T [l];
    for(i=0; i < l; ++i){
      array[i] = init;       //init values
    }
  }catch(...) {throw ("Out of memory!\n");}
}

template <typename T>
static T* new_1D(int l){
  T* array = NULL;
  create_3D(array, l);
  return array;
}

template <typename T>
static T* new_1D(int l, T init){
  T* array = NULL;
  create_1D(array, l, init);
  return array;
}

template <typename T>
static void reset_1D(T* &array, int l, T init){
  assert(array);
  int i;
  for(i=0; i < l; ++i){
    array[i] = init;       //init values
  }
}


// swap two elements in an array
template <typename T>
static void swap(T* a, const int& i, const int& j){
  T temp = a[j];
  a[j] = a[i];
  a[i] = temp;
}

private:
template <typename T>
static int increment(const T& i1, const T& i2){
  if( i1 < i2) return 1;
  else return -1;
}
template <typename T>
static int decrement(const T& i1, const T& i2){
  if( i1 < i2) return -1;
  else return 1;
}
public:
//******************************************************************************
// quicksort for one dimensional array
//******************************************************************************
private:
template <typename T>
static void quicksort( T* a, int left, int right, int (*comp)(const T&, const T&)){
  if (left < right) {
    int lp, rp, i = left + 1, j = left + 1;
    T x = a[left];
    while (j <= right) {
      if((*comp)(a[j], x)<0){
      //if (a[j] < x) {
        swap(a, i, j);
        i++;
      }
      j++;
    }
    swap(a, left, i-1);
    lp = i - 2;
    rp = i;

    quicksort(a, left, lp, comp);
    quicksort(a, rp, right, comp);
  }
}
public:

//------------------------------------------------------------------------------
template <typename T>
static void quicksort( T* a, int n, bool incremental=true){
  if(incremental) quicksort<T>(a, 0, n-1, &decrement<T>);
  else            quicksort<T>(a, 0, n-1, &increment<T>);
}

//******************************************************************************
// quicksortIndex for one dimensional array where the index of the sorted array is returned
//******************************************************************************
/** get the sorted index
  * it is possible to pass an array for the output (aIndex), otherwise the output
  * array is created and has to be deleted afterwards.
  * The original unsorted array wil not be affected by the sort
  */
private:
template <typename T>
static void quicksortIndex( T* a, int left, int right, int* index, int (*comp)(const T&, const T&)){
  if (left < right) {
    int lp, rp, i = left + 1, j = left + 1;
    T x = a[index[left]];
    while (j <= right) {
      if((*comp)(a[index[j]], x)<0){
        swap(index, i, j);
        i++;
      }
      j++;
    }
    swap(index, left, i-1);
    lp = i - 2;
    rp = i;

    quicksortIndex(a, left, lp, index, comp);
    quicksortIndex(a, rp, right, index, comp);
  }
}
public:

//------------------------------------------------------------------------------
template <typename T>
static int* quicksortIndex( T* a, int n, int* index, bool incremental=true){
  // set the indexes (1,2,3,4,...size-1)
  if(!index) index = new int[n];
  for(int i = 0; i<n; ++i){
    index[i] = i;
  }

  if(incremental) quicksortIndex(a, 0, n-1, index, &decrement<T>);
  else            quicksortIndex(a, 0, n-1, index, &increment<T>);
  return index;
}

//******************************************************************************
// quicksortBoth for one dimensional array where the index AND the array are sorted
//******************************************************************************
/** sort the array a and change the elements of array index at the same time
  * Both arrays have to have the same size
  */
private:
template <typename T, typename S>
static void quicksortBoth( T* a, int left, int right, S* index, int (*comp)(const T&, const T&)){
  if (left < right) {
    int lp, rp, i = left + 1, j = left + 1;
    T x = a[left];
    while (j <= right) {
      if((*comp)(a[j], x)<0){
        swap(a, i, j);
        swap(index, i, j);
        i++;
      }
      j++;
    }
    swap(a, left, i-1);
    swap(index, left, i-1);
    lp = i - 2;
    rp = i;

    quicksortBoth(a, left, lp, index, comp);
    quicksortBoth(a, rp, right, index, comp);
  }
}
public:

//------------------------------------------------------------------------------
template <typename T, typename S>
static void quicksortBoth( T* a, int n, S* index, bool incremental=true){
  if(incremental) quicksortBoth<T,S>(a, 0, n-1, index, &decrement<T>);
  else            quicksortBoth<T,S>(a, 0, n-1, index, &increment<T>);
}

//------------------------------------------------------------------------------
/** makes the array a cumulative array */
template <typename T>
static void cumulative(T* a, int n){
  for(int i = 1; i<n; ++i){
    a[i] += a[i-1];
  }
}
//------------------------------------------------------------------------------
/** return a cumulative array to a non cumulative array */
template <typename T>
static void undo_cumulative(T* a, int n){
	for(int i = n-1; i>0; --i){
		a[i] -= a[i-1];
	}
}
//------------------------------------------------------------------------------
/** the smallest values become the biggest ones
	* a value of zero is approximated by 1e-9
	*/
template <typename T>
static void reverse_cumulative(T* a, int n){
	a[0] = a[0] ? (1.0/a[0]) : 1e9; // first element
  for(int i = 1; i<n; ++i){
		a[i] = a[i] ? (1.0/a[i] + a[i-1]) : (1e9 + a[i-1]);
  }
}

//------------------------------------------------------------------------------
/** return a cumulative array to a non cumulative array
  * it is not possible to return to the absolut values, but the relation between the values is maintained
  */
template <typename T>
static void undo_reverse_cumulative(T* a, int n){
	for(int i = n-1; i>0; --i){
		a[i] = 1.0/(a[i]-a[i-1]);
	}
	a[0] = 1.0/a[0]; // last element
}

//******************************************************************************
// compute mean, variance, standard variance
//******************************************************************************
/** compute the sum of an array */
template <typename T>
static T sum(T* array, const int& size){
  if(!size) return my_NAN;

  T sum=0;
  T* end = &array[size];          // get the end
  for(; array != end; ++array){
    if(*array == my_NAN) continue;   // continue if the value is missing
    sum += *array;
  }
  return  sum;
}

//------------------------------------------------------------------------------
/** compute the mean of an array */
template <typename T>
static double mean(T* array, const int& size){
  if(!size) return my_NAN;

  double sum=0;
  unsigned int nb=0;
  T* end = array + size;          // get the end

  for(; array != end; ++array){   // for each element
    if(*array == my_NAN) continue;// if it is not a missing value
    sum += *array;
    ++nb;
  }
  return  nb ? sum/nb : my_NAN;
}

//------------------------------------------------------------------------------
/** compute the unbiased estimated variance of an array */
template <typename T>
static double estimateVariance(T* array, const int& size){
  if(size<=1) return my_NAN;

  double sum=0.0, sq_sum=0.0;
  unsigned int nb = 0;
  T* end = array + size;          // get the end

  for(; array != end; ++array) {
    if(*array == my_NAN) continue;
    sum    += (*array);
    sq_sum += (*array) * (*array);
    ++nb;
  }

  return (nb>1) ? (double)(sq_sum - sum*sum/nb) / (nb-1) : my_NAN;
}

//------------------------------------------------------------------------------
/** compute the biased variance of an array */
template <typename T>
static double var(T* array, const int& size, double mu=my_NAN){
  if(mu == my_NAN) mu = mean(array, size);  // if the mean has not be computed
  if(mu == my_NAN) return my_NAN;                    // the array has no valid values

  double diff, sum=0;
  unsigned int nb=0;
  T* end = array + size;  // get the end

  for(; array != end; ++array) {                // for each element
    if(*array == my_NAN) continue;              // if it is not a valid element
    diff = (*array)-mu;
    sum  += diff*diff;
    ++nb;
  }

  sum /=nb;
  return (abs(sum) > 1e-30) ? sum : 0;                                 // no control as it was alredy checked
}

//------------------------------------------------------------------------------
/** compute the standard deviation (sd) of an array */
template <typename T>
static double sd(T* array, const int& size){

  // special cases
  if(size<=1) return my_NAN;

  return sqrt(var(array, size));
}

//------------------------------------------------------------------------------
/** transforms the array that the sum of the array is 1
  * returns true, if this was already the case, and false if the array had to be adjusted
	*/
template <typename T>
static bool make_frequency(T* array, const int& size){

  double Sum = (double) sum(array, size);
  if(abs(1-Sum) < 1e-4) return true;

  // if all values are zero
  if(!Sum) Sum = 1.0/size;

	T* end = array + size;          	// get the end

  for(; array != end; ++array) {
		(*array) = (*array)/Sum;
		if(*array<1e-13) *array = 0;		// due to accuracy problems
	}
	return false;
}

};


#endif
