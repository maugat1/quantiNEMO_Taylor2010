/** @file service.h
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

#ifndef serviceH
#define serviceH

#include <list>
#include "handler.h"
#include "output.h"
using namespace std;


class SimComponent;

/**Interface for the simulation services (files and stats).
 * Implements the observer pattern. Notify the observers (Handler) to update their state.
 * Contains the observer list. Provides interface to attach the observers to the service.
**/
class Service {
private:
	/**The list of observers. */
  list<Handler*> _observers;

public:

  Service( ) { }

  virtual ~Service( ) { }

  /**Inits internals. */
  virtual bool init ( ) = 0;

  /**Notifies all observers to update their state. */
	virtual void notify ( ) {
		#ifdef _DEBUG
      message(" %i ... ",_observers.size());
		#endif

		list< Handler* >::iterator HIT = _observers.begin();
		while(HIT != _observers.end()) {
			(*HIT)->update();
	    HIT++;
	  }
	}

	virtual list<Handler*>* get_observers( ){return &_observers;}
	
  /**Interface to used by a simulation component to load its obervers onto a service provider.*/
  virtual void load ( SimComponent* sc ) = 0;

	/**Adds an observer to the list. */
  virtual void attach ( Handler* h ) {
    _observers.push_back(h);
  }

  /**Clears the observers list. */
  virtual void reset ( ) {
    _observers.clear();
  }

};

#endif //SERVICE_H






