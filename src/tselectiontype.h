/** @file tselectiontype.h
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

#ifndef tselectiontypeH
#define tselectiontypeH
//---------------------------------------------------------------------------
#include "tselectiontrait.h"
class TSelection;

class TSelectionNeutral: public TSelectionTrait{
  private:
		double (TSelectionNeutral::* _func_ptr_get_fitness)();
		double get_fitness_none();
		double get_fitnessFactor();

	public:
		TSelectionNeutral(TSelection* s, const int& t){
			TSelectionTrait::init(s, t);
			init();
		}
		~TSelectionNeutral(){}
		void init();

		double get_fitness(){return (this->*_func_ptr_get_fitness)();}
};

class TSelectionStabilizing: public TSelectionTrait{
	private:
		double (TSelectionStabilizing::* _func_ptr_get_fitness)();
		double get_fitness_none();
		double get_fitnessFactor();

	public:
		TSelectionStabilizing(TSelection* s, const int& t){
			TSelectionTrait::init(s, t);
			init();
		}
		~TSelectionStabilizing(){}
		void init();

		double get_fitness(){return (this->*_func_ptr_get_fitness)();}
};

class TSelectionDirectional: public TSelectionTrait{
	private:
		double (TSelectionDirectional::* _func_ptr_get_fitness)();
		double get_fitness_none();
		double get_fitnessFactor();

	public:
		TSelectionDirectional(TSelection* s, const int& t){
			TSelectionTrait::init(s, t);
			init();
		}
		~TSelectionDirectional(){}
		void init();

		double get_fitness(){return (this->*_func_ptr_get_fitness)();}
};
#endif
