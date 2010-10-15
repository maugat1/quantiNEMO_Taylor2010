/** @file main.cpp
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


#include "version.h"
#include "simulation.h"
#include "stathandler.cpp"
#include "filehandler.h"

#include <iostream>
using namespace std;

char SEP;

//--------------------------------------------------------------------
int main (int argc, char **argv)
{
	SimRunner *theSim = NULL;
	try{
		message("\n**************************************************\n"   \
				"\n                q u a n t i N E M O"                     \
				"\n**************************************************"     \
				"\n*     Release: %i.%i.%i%s [%s; %s]     *"                \
				"\n*    Copyright (C) 2008 Samuel Neuenschwander    *"     \
				"\n* http://www.unil.ch/popgen/softwares/quantinemo *"     \
				"\n**************************************************\n",
				RELEASE,REVISION,MINOR_VERSION,TEMP_VERSION,VERSION_DATE,VERSION_TIME);

#ifdef _DEBUG
		message("\n***** DEBUG MODE *****\n\n");
#endif

		// get the working directory (remove the exe name)
		string dir = argv[0];
		SEP = '\\';
		string::size_type  pos = dir.find_last_of(SEP);
		if(pos == string::npos){
			SEP = '/';
			pos = dir.find_last_of(SEP);
			if(pos == string::npos) argv[0][0]     = '\0';
			else                    argv[0][pos+1] = '\0';
		}
		else                      argv[0][pos+1] = '\0';

		// read quantiNemo.ini
		theSim = new SimRunner();
		theSim->readInputFile(argc, argv);

		FileHandler::set_base_path(argv[0]);

		if( !theSim->run() ) fatal(" could not run the sim!\n");
		delete theSim;


		message("\nquantiNEMO terminated successfully!\n");
	}
	catch(const bad_alloc& x){
		if(theSim) delete theSim;
		error("Out of memory (%s)!\n", x.what());
		error("\nquantiNEMO was not able to run successfully!\n");
	}
	catch(const char* text){
		if(theSim) delete theSim;
		error("%s!\n", text);
		error("\nquantiNEMO was not able to run successfully!\n");
	}
	catch(...){
		if(theSim) delete theSim;
		error("\nquantiNEMO was not able to run successfully!\n");
	}

#ifdef _DEBUG
	// this allows to stop the program before the console window is closed
	int w;
	cin >> w;
#endif
	return 0;
}




