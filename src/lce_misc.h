/** @file lce_misc.h
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

#ifndef lce_miscH
#define lce_miscH


#include "types.h"
#include "lifecycleevent.h"
#include "fileservices.h"
#include "statservices.h"
#include "filehandler.h"
#include "stathandler.h"


// LCE_FileServicesNotifier
//
/**Event used to notify all file handlers to update their state through the FileServices::notify() interface.*/
class LCE_FileServicesNotifier: public LCE{

    private:
        FileServices* _service;


    public:
        LCE_FileServicesNotifier(int rank = my_NAN)
            : LCE("save_files","",rank),_service(0){
            }

        virtual ~LCE_FileServicesNotifier( ) { }

        virtual void execute ();

        virtual LCE_FileServicesNotifier* clone ( ) {return new LCE_FileServicesNotifier();}

        //SimComponent overrides:
        virtual void  loadFileServices ( FileServices* loader ) {
            _service = loader;
        }
        virtual void  loadStatServices (StatServices* loader ) {}
        virtual age_t removeAgeClass ( ) {return 0;}
        virtual age_t addAgeClass ( ) {return 0;}
        virtual age_t requiredAgeClass () {return 0;}
};

// LCE_StatFH
//
/**FileHandler of the LCE_StatServiceNotifier class, writes the recorded stats to txt files.*/
class LCE_StatFH : public FileHandler {
    public:
        StatServices* _statService;

    public:

        LCE_StatFH () : FileHandler(".txt"), _statService(0) {
            FileHandler:: set(false,false,0,0,0,"","",0,0,NULL, 1);
        }

        ~LCE_StatFH() { };

        virtual string getName() {return "LCE_SataFH";}

        void set_statService(StatServices* srv) {
            _statService = srv;

            // set the filenames
            string base = get_path() + get_service()->getBaseFileName();
            _statService->set_file_name_stats(base + "_stats.txt");
            _statService->set_file_name_mean(base + "_mean.txt");
            _statService->set_file_name_var(base + "_var.txt");
        }

        virtual void FHwrite();

        void set_save_choice(const int& i){_statService->set_save_choice(i);}

        void printStat_each_header(ostream& FH);
        void printStat_each    ( );
        void printStat_mean    ( );
        void printStat_variance( );
        void printStat_legend  ( );
};

// LCE_StatSH
//
/**StatHandler of the LCE_StatServiceNotifier class, adds a default StatRecorder to the recorders list (alive.rpl).*/
class LCE_StatSH : public StatHandler<LCE_StatSH> {

    public:
        LCE_StatSH () {}
        ~LCE_StatSH () {}

        virtual string getName() {return "LCE_StatSH";}
        virtual bool init ( ){
            StatHandler<LCE_StatSH>::init();
#ifdef _DEBUG
            message("  default (");
#endif

            add("Nb alive replicates","alive.rpl",GEN,ALL,0,&LCE_StatSH::get_isAlive,0,0,0);

#ifdef _DEBUG
            message(")\n");
#endif
            return true;
        }

        virtual bool setStatRecorders(const string& token) {return false;}

        double get_isAlive ();
};

// LCE_StatServiceNotifier
//
/**Initiates the StatServices' parameters (log time) when registering, calls StatServices::notify() when executing.
 * Registers the file handler used to save the stats to the '.txt' and '_mean.txt' files.*/
class LCE_StatServiceNotifier: public LCE
{
    private:
        StatServices* _service;

        unsigned int _occurrence;
        unsigned int _tot_occurrence;   // total number of values/generations to store
        unsigned int _save_choice;
        unsigned int* _index;

        string _arg;

        LCE_StatFH _fileHandler;
        LCE_StatSH _statHandler;

    public:

        LCE_StatServiceNotifier(int rank = my_NAN);

        virtual ~LCE_StatServiceNotifier ( ) {
            if(_index) delete[] _index;
        }

        virtual string getName() {return "LCE_StatServiceNotifier";}

        virtual void  execute ();

        virtual LCE_StatServiceNotifier* clone ( ) {return new LCE_StatServiceNotifier();}

        virtual bool init (Metapop* popPtr);
        virtual int get_save_choice(){return _save_choice;}

        //SimComponent overrides:
        virtual void loadFileServices ( FileServices* loader ) {loader->attach(&_fileHandler);}
        virtual void loadStatServices ( StatServices* loader );
        virtual age_t removeAgeClass ( ) {return 0;}
        virtual age_t addAgeClass ( ) {return 0;}
        virtual age_t requiredAgeClass () {return 0;}

        virtual void executeBeforeEachGeneration(const int& gen);

};

// LCE_Regulation
//
/** Regulates the population after dispersal.
 * Moves individuals from post-dispersal to adult containers according to Patch capacities.
 * Sets the age flag of the population to ADULTS.
 **/
class LCE_Regulation: public LCE{
    protected:
        void (LCE_Regulation::* regulation)();     // pointer to the correct regulation function

        string  _ageStr;
        age_idx _age;

#ifdef _DEBUG
        unsigned int _ex_cnt, _col_cnt, _ph_cnt;
#endif

    public:

        LCE_Regulation(age_idx age = ADLTx, int rank = my_NAN) : LCE(("regulation_"+ (string)(age==ADLTx ? "adults":"offspring")).c_str(),"",rank) {
            _age = age;
            _ageStr = age==ADLTx ? "adults":"offspring";
            this->add_parameter("regulation_model_"+_ageStr,INT2,false,true,0,1,"0");
        }

        virtual LCE_Regulation* clone () {return new LCE_Regulation();}

        virtual ~LCE_Regulation( ) { }

        virtual void execute ();
        virtual bool init(Metapop* popPtr);
        virtual void regulation_neutral();
        virtual void regulation_fitness_patch()  {_popPtr->regulate_selection_fitness_patch(_age);}
        virtual void regulation_fitness_metapop(){_popPtr->regulate_selection_fitness_metapop(_age);}
        virtual void regulation_fitness_hard()   {_popPtr->regulate_selection_fitness_hard(_age);}
        virtual void drawSuccessfullIndividuals(Patch* curPatch, const unsigned int& K, const sex_t& SEX);
        virtual void drawUnSuccessfullIndividuals(Patch* curPatch, const unsigned int& K, const sex_t& SEX);

        //SimComponent overrides:
        virtual void loadFileServices ( FileServices* loader ) {}
        virtual void loadStatServices ( StatServices* loader ) {}
        virtual age_t removeAgeClass ( ) {return 0;}
        virtual age_t addAgeClass ( ) {return 0;}
        virtual age_t requiredAgeClass () {return _age==ADLTx ? ADULTS : OFFSPRG;}
};

// LCE_Aging
//
/**Removes all adults from the patches adult containers.
 * This is the only LCE that actually removes the adults from the patches. It is thus
 * mandatory to use it once in the life cycle, preferably after breeding!**/
class LCE_Aging: public LCE
{
    public:

        LCE_Aging(int rank = my_NAN) : LCE("aging","",rank) {}

        virtual ~LCE_Aging( ) { }

        virtual void execute ();

        virtual LCE_Aging* clone ( ) {return new LCE_Aging();}

        //implementations:
        virtual void loadFileServices ( FileServices* loader ) {}
        virtual void loadStatServices ( StatServices* loader ) {}
        virtual age_t removeAgeClass ( ) {return ADULTS;}
        virtual age_t addAgeClass ( ) {return 0;}
        virtual age_t requiredAgeClass () {return ADULTS | OFFSPRG;}
};


// LCE_Patch_Extinction
//
/**Randomly extincts patches according to the extinction rate parameter.
 * Sets the patches extinction flag accordingly.**/
class LCE_Patch_Extinction: public LCE
{
    /**Patch extinction probability.*/
    double _Xtion_rate;

    public:

    LCE_Patch_Extinction(int rank = my_NAN) : LCE("extinction","",rank), _Xtion_rate(0)
    { this->add_parameter("extinction_rate",DBL,true,true,0,1,"0",true); }

    virtual ~LCE_Patch_Extinction( ) { }

    virtual bool init (Metapop* popPtr);

    virtual void execute ();

    virtual void executeBeforeEachGeneration(const int& gen);

    virtual LCE_Patch_Extinction* clone ( ) {return new LCE_Patch_Extinction();}

    //SimComponent overrides:
    virtual void loadFileServices ( FileServices* loader ) {}
    virtual void loadStatServices ( StatServices* loader ) {}
    virtual age_t removeAgeClass ( ) {return 0;}
    virtual age_t addAgeClass ( ) {return 0;}
    virtual age_t requiredAgeClass () {return ADULTS | OFFSPRG;}
};



#endif //LCEMISC_H
