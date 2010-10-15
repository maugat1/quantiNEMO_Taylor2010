/** @file filehandler.h
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

#ifndef filehandlerH
#define filehandlerH

#include "handler.h"
#include "fileservices.h"
#include "ttrait.h"
class Patch;
class Individual;

/** Interface to handle file input/output for any SimComponent.
 *  Stores the periodicity parameters and the file path and extension. The replicate file name
 *  is given by the FileServices. A file handler might be set to write output at specific generation
 *  of a specific replicate or at some periodic time during the simulation.
 */
class FileHandler : public Handler {

private:

  FileServices* _service;

  bool _isReplicatePeriodic;
  bool _isGenerationPeriodic;

  unsigned int _ReplicateOccurrence;
  unsigned int _GenerationOccurrence;
  unsigned int _TotGenerationOccurrence;
  unsigned int _current_replicate;
  unsigned int _current_generation;
  unsigned int _ExecRank;

  string  _path;
  string  _short_path;
  static string _base_path;
  string  _extension;
  string  _current_filename;
  int _fstat_choice;    // 0: no output; 1: FSTAT output; 2: FSTAT extended output
  string _script;       // a script name which should be executed after the write of a file


  void FHwrite (const age_idx& cur_age, const sex_t& cur_sex, ostream& FILE,
                Patch* current_patch, const int& patch_id, const int& nbPatchDigit, const int& position);
	void FHwriteIndividual2Stream(ofstream& FILE, Individual* ind, const int& ploidy, const int& position);

protected:
  vector<TTraitProto*> _trait;
  vector<int> _TTidx;
  int _nb_trait;

  age_t _age;           // which age should be outputted (ALL; ADULTS, OFFSPRG)?
  int   _sex;           // which sex should be outputet (0: both; 1: only female; 2: only male)?

public:

    FileHandler (const char* ext) : _service(0), _isReplicatePeriodic(0), _isGenerationPeriodic(0),
    _ReplicateOccurrence(0), _GenerationOccurrence(0), _current_replicate(0), _current_generation(0),
    _ExecRank(0), _extension(ext), _age(0), _sex(0) {
    }

  virtual ~FileHandler ( ) { }
  /**Called by notifier during simulation setup, performs file checking.
   *@return false if filename already exists on disk, true otherwise.
   */
  virtual bool  init ( );
  ///@name Accessors
  ///@{
  Metapop*      get_pop_ptr ( )      {return _service->get_pop_ptr();}

  FileServices* get_service () {return _service;}

	void set_service (FileServices* srv) {_service = srv;}

  string&       get_path () {return _path;}
  string&       get_short_path () {return _short_path;}

	void set_path (string& path) {
		_path = path;
		if(!_path.empty() && _path[_path.length()-1] != SEP) _path += SEP;
		_short_path = _path;
		if(!_service->getSimfolder().empty()){
			_path = _service->getSimfolder()+_path;
		}
		check_path(_path);
	}

  static void   check_path(string& path);

	static void   set_base_path(string s){
		if(s.empty())                       _base_path = "";  // nothing passed: to prevent problems
		else if(s[0] != '/' && s[1] != ':') _base_path = "";	// no absolute path passed, just the name
		else{                                                 // aboslute path passed: remove the file name from it
			assert(s.rfind(SEP) != string::npos);
			_base_path = s.substr(0, s.rfind(SEP));
			char lastChar = _base_path[_base_path.size()-1];    // add the directoy seperator if not present
			if(lastChar != '/' && lastChar != '\\') _base_path += SEP;
		}

		#ifdef _DEBUG
		 message("\n FileServices::_base_path: '%s'\n", _base_path.c_str());
		#endif
	}

	static string get_base_path()        {return _base_path;}

  string&       get_extension ( )  {return _extension;}

  void          set_extension (const char* ext) {_extension = ext;}

  string&       get_filename ();

  virtual bool get_isReplicatePeriodic ()         {return _isReplicatePeriodic;}
  virtual void set_isReplicatePeriodic (bool val) {_isReplicatePeriodic = val;}

  virtual int  get_ReplicateOccurrence  ()         {return _ReplicateOccurrence;}
  virtual void set_ReplicateOccurrence  (unsigned int val)  {_ReplicateOccurrence = val;}

  virtual bool get_isGenerationPeriodic ()         {return _isGenerationPeriodic;}
  virtual void set_isGenerationPeriodic (bool val) {_isGenerationPeriodic = val;}

  virtual int  get_GenerationOccurrence  ()         {return _GenerationOccurrence;}
  virtual void set_GenerationOccurrence  (unsigned int val)  {_GenerationOccurrence = val;}

  virtual int  get_TotGenerationOccurrence  ()         {return _TotGenerationOccurrence;}
  virtual void set_TotGenerationOccurrence  (unsigned int val)  {_TotGenerationOccurrence = val;}

  virtual int  get_ExecRank ()            {return _ExecRank;}
  virtual void set_ExecRank (int val)     {_ExecRank = val;}

  virtual int  get_fstat_choice()         {return _fstat_choice;}
  virtual void set_fstat_choice (int val) {_fstat_choice = val;}
	void write_individual_info_to_stream(ostream& FILE, Individual* ind, const age_idx& cur_age, const sex_t& cur_sex, char sep=' ');

  virtual string get_script ()            {return _script;}
  virtual void set_script (string val)    {_script = val;}
  virtual void set_sex(int i)             {_sex = i;}
	virtual void set_age(int i)             {
		switch(i){
			case 0: _age=ADULTS;    break;
			case 1: _age=OFFSPRG;   break;
			case 2: _age=ALL; 	    break;
		}
	}
  ///@}

  /**Sets the hanlder parameters.
    @param rpl_per replicate periodicity
    @param gen_per generation periodicity
    @param rpl_occ replicate occurence
    @param gen_occ generation occurence
    @param rank the rank in the life cycle, actualy unused...
    @param path the file path
    */
  virtual void set (bool rpl_per, bool gen_per, int rpl_occ, int gen_occ, int rank,
                    string path, string script, int sex, int age, TTraitProto* trait, const int& choice) {
	  set_isReplicatePeriodic(rpl_per);
    set_isGenerationPeriodic(gen_per);
    set_ReplicateOccurrence(rpl_occ);
	  set_GenerationOccurrence(gen_occ);
    set_ExecRank(rank);
    set_path(path);
    set_fstat_choice(choice);
    set_script(script);
    set_sex(sex);
    set_age(age);

    if(trait){
      _trait.push_back(trait);
      _TTidx.push_back(trait->get_absolute_index());
      _nb_trait = (int) _TTidx.size();
    }
  }

  virtual void set (TTraitProto* trait) {
    if(trait){
      _trait.push_back(trait);
      _TTidx.push_back(trait->get_absolute_index());
      _nb_trait = (int) _TTidx.size();
    }
  }

  /**Default behaviour of the class, called by Handler::update().**/
  virtual void  FHwrite ();

  virtual void update();

  virtual void execute_script(string script, string filename);
};
#endif //FILEHANDLER_H


