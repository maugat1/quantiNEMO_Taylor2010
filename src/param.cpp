/** @file param.cpp
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

#include <vector>
#include <sstream>
#include "param.h"
#include "output.h"
using namespace std;



// ----------------------------------------------------------------------------------------
// Param
// ----------------------------------------------------------------------------------------
Param::Param (string& Name,param_t Type,bool mandatory,bool bounded,double low_bnd,
        double up_bnd, string def, bool temp)
: _name(Name),_arg(def),_default_arg(def),_type(Type),_isSet(0),_isBounded(bounded),_isRequired(mandatory),_matrix(0),
    _matrixVar(0), _temporalParamAllowed(temp)
{
    _bounds[0] = low_bnd;
    _bounds[1] = up_bnd;
}

// ----------------------------------------------------------------------------------------
// update_arg
// ----------------------------------------------------------------------------------------
/** the argument is updated by the new temporal argument */
void Param::update_arg(const int& gen){
    assert(_temporalArgs.find(gen) != _temporalArgs.end());

    _arg = _temporalArgs.find(gen)->second;
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void Param::reset ()
{
    _arg = _default_arg;
    _isSet = false;
    _temporalArgs.clear();
}

// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
/** control if the value is within the bounds */
bool Param::check_arg (string& arg)
{
    //integer or decimal:
    if(_type == INT2 || _type == DBL || _type == INT_MAT || _type == DBL_MAT) {
        //!!a matrix may also be specified!! in that case, we don't check values here
        if(arg[0] != '{') {
            if(!STRING::is_number(arg)){
                fatal("The argument of parameter '%s' should be a number (arg: %s)!\n", arg.c_str());
                return false;
            }
            double val = STRING::str2int<double>(arg);

            // if the value is outside the bounds
            if(_isBounded && (val < _bounds[0] || val > _bounds[1]) ) {
                if(_type == INT2 || _type == INT_MAT){
                    fatal("The value of parameter '%s' is outside of the bounds (value: %i; bounds [%i - %i])!\n",
                            _name.c_str(), (int)val, (int)_bounds[0], (int)_bounds[1]);
                }
                else{
                    fatal("The value of parameter '%s' is outside of the bounds (value: %f; bounds [%f - %f])!\n",
                            _name.c_str(), val, _bounds[0], _bounds[1]);
                }
                return false;
            }
        }
    }
    else if(_type == MAT || _type == MAT_VAR || _type == INT_MAT || _type == DBL_MAT) {
        if(arg[0] != '{') {
            fatal("The argument of parameter '%s' should be a matrix!\n", _name.c_str());
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
/** control if the value is within the bounds */
bool Param::set (string& arg, ParamSet* parSet)
{
    _tot_arg = arg;
    _parSet  = parSet;

    // replace the sequences
    STRING::replace_seq(arg);

    // replace the repetition
    STRING::replace_rep(arg);

    // if it is a temporal argument
    if(arg[0] == '(') return set_temporal_args(arg);

    // if the string is enclosed with ""
    if(arg[0] == '\"'){
        assert(arg[arg.length()-1] == '\"');
        arg = arg.substr(1,arg.length()-2);
    }

    if(!check_arg(arg)) return false;

    _isSet = true;
    _arg = arg;
    return true;
}
// ----------------------------------------------------------------------------------------
// set
// ----------------------------------------------------------------------------------------
/** read a temporal parameter in the format int, string. The pairs are separated by , or ;
 * Comments have previously been removed
 */
bool Param::set_temporal_args (string& arg)
{
    if(!_temporalParamAllowed){
        fatal("Temporal parameter: Parameter %s can not change over time!\n", _name.c_str());
    }

    int gen, line, counter=0;
    string curArg;
    istringstream IN;
    char c;
    _temporalArgs.clear();

    IN.str(arg);
    IN.get(c);      // remove "("

    // for each temporal setting
    while(IN.good() && !IN.eof() && c != ')'){
        // read the generation time
        if(!STRING::removeCommentAndSpace(IN, line)) fatal("Reading temporal series!\n");
        IN >> gen;
        IN >> ws;

        // read the argument
        curArg = "";
        while((IN.get(c) && IN.good() && !IN.eof() && c != ';' && c != ',' && c != ')') || counter){
            if(c == '{')      ++counter;
            else if(c == '}') --counter;

            curArg += c;
        }

        if(gen < 1) gen=1;                // the first generation time has to be 1!!!

        // store the argument
        if(_temporalArgs.find(gen) != _temporalArgs.end()){
            warning("Temporal parameter: Generation %i of parameter %s already set!\n", gen, _name.c_str());
        }
        if(!check_arg(curArg)) return false;
        _temporalArgs[gen] = curArg;      // store the variables
    }

    // the first argument should be the one used to start the simulation (generations are sorted by the map)
    if(_temporalArgs.begin()->first != 1) fatal("Temporal parameter: The parameter %s is not defined for the first generation!\n", _name.c_str());
    _arg = _temporalArgs.begin()->second;

    // if the temporal parameter has only a single element it is indeed not a temporal parameter
    if(_temporalArgs.size() == 1){
        _temporalArgs.clear();
    }
    else{
        // add the generation times to the list in the paramSet
        for(map<int, string>::iterator pos = _temporalArgs.begin(); pos != _temporalArgs.end(); ++pos){
            _parSet->set_temporal_param(pos->first, this);
        }
    }

    _isSet = true;
    return true;
}

// ----------------------------------------------------------------------------------------
// get_value
// ----------------------------------------------------------------------------------------
double Param::get_value ()
{
    assert(!(is_matrix() || _type == STR));
    return STRING::str2int<double>(_arg.c_str());
}

// ----------------------------------------------------------------------------------------
// get_matrix
// ----------------------------------------------------------------------------------------
TMatrix* Param::get_matrix ()
{
    if( is_matrix() && _isSet){
        if(_matrix) delete _matrix;
        _matrix = parse_matrix();
        return _matrix;
    }
    else return NULL;
}

// ----------------------------------------------------------------------------------------
// get_matrix
// ----------------------------------------------------------------------------------------
TMatrixVar<double>* Param::get_matrixVar ()
{
    if( is_matrixVar() && _isSet ){
        if(_matrixVar) delete _matrixVar;
        _matrixVar = parse_matrixVarDbl();
        return _matrixVar;
    }
    else return NULL;
}

// ----------------------------------------------------------------------------------------
// parse_matrixVarDbl
// ----------------------------------------------------------------------------------------
/** creates and returns a TMatrixVar object
 * brackets are already "clean"
 */
TMatrixVar<double>* Param::parse_matrixVarDbl ()
{
    TMatrixVar<double>* mat = new TMatrixVar<double>();
    try{
        mat->read(_arg);
    }
    catch(const char* error){
        delete mat;
        fatal("Parameter '%s': %s\n", _name.c_str(), error);
    }

    return mat;
}

// ----------------------------------------------------------------------------------------
// parse_matrixVarStr
// ----------------------------------------------------------------------------------------
/** creates and returns a TMatrixVar object
 * brackets are already "clean"
 */
TMatrixVar<string>* Param::parse_matrixVarStr ()
{
    TMatrixVar<string>* mat = new TMatrixVar<string>();
    try{
        mat->read(_arg);
    }
    catch(const char* error){
        delete mat;
        fatal("Parameter '%s': %s\n", _name.c_str(), error);
    }

    return mat;
}

// ----------------------------------------------------------------------------------------
// parse_matrix
// ----------------------------------------------------------------------------------------
/*   */
TMatrix* Param::parse_matrix ()
{
    TMatrix* mat = new TMatrix();
    try{
        mat->read(_arg);
    }
    catch(const char* error){
        delete mat;
        fatal("Parameter '%s': %s\n", _name.c_str(), error);
    }

    return mat;
}

// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void Param::show_up ()
{
    message("\n%s\t%s\t%s\t%i", _name.c_str(), get_type_str().c_str(),
            _default_arg.empty() ? "\"\"":_default_arg.c_str(), _isRequired ? 1:0);
}

// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
string Param::get_type_str(){
    switch(_type){
        case DBL:       return "double";
        case INT2:      return "integer";
        case STR:       return "string";
        case MAT:       return "matrix";
        case DIST:      return "distribution";
        case MAT_VAR:   return "variable matrix";
        case INT_MAT:   return "integer/matrix";
        case DBL_MAT:   return "double/matrix";
        case STR_MAT:   return "string/matrix";
    }
    return "";
}

/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                             ******** ParamSet ********

// ----------------------------------------------------------------------------------------
// ParamSet
// ----------------------------------------------------------------------------------------
ParamSet::ParamSet(string name, bool isRequired) : _isSet(0){
_name=name;
_isRequired=isRequired;
#ifdef _DEBUG
message(" ParamSet::ParamSet(%s)\n", _name.c_str());
#endif
}

// ----------------------------------------------------------------------------------------
// ParamSet
// ----------------------------------------------------------------------------------------
ParamSet::ParamSet( ) : _isSet(0), _isRequired(0) {
#ifdef _DEBUG
message(" ParamSet::ParamSet(not set)\n", _name.c_str());
#endif
}

// ----------------------------------------------------------------------------------------
// ~ParamSet
// ----------------------------------------------------------------------------------------
ParamSet::~ParamSet ()
{
#ifdef _DEBUG
message(" ParamSet::~ParamSet(%s)\n", _name.c_str());
#endif
map<string, Param*>::iterator param = _params.begin();
for(; param != _params.end(); ++ param) {
if(param->second) delete param->second;
}
_params.clear();

multimap<int, map<string, Param*>* >::iterator pos = _temporalParams.begin();
for(; pos != _temporalParams.end(); ++pos){
if(pos->second) delete pos->second;
}
_temporalParams.clear();
}

// ----------------------------------------------------------------------------------------
// reset
// ----------------------------------------------------------------------------------------
void ParamSet::reset ()
{
_isSet = false;

map<string, Param*>::iterator param =  _params.begin();
for(; param != _params.end(); ++ param) {
param->second->reset();
}

multimap<int, map<string, Param*>* >::iterator pos = _temporalParams.begin();
for(; pos != _temporalParams.end(); ++pos){
if(pos->second) delete pos->second;
}
_temporalParams.clear();
}

// ----------------------------------------------------------------------------------------
// add_param
// ----------------------------------------------------------------------------------------
void ParamSet::add_param (Param* param)
{
assert(_params.find(param->get_name()) == _params.end());
_params[param->get_name()] = param;
}
// ----------------------------------------------------------------------------------------
// add_param
// ----------------------------------------------------------------------------------------
void ParamSet::add_param (string Name,param_t Type,bool mandatory,bool isBounded,
        double low_bnd,double up_bnd,string def,bool temp)
{
    // control that not another object is overridden if the names are identical
    map<string, Param*>::iterator cur_param = _params.find(Name);
    if(cur_param != _params.end()) {
        delete cur_param->second;
        fatal("The parameter '%s' was already set!", Name.c_str());
    }

    _params[Name] = new Param(Name, Type, mandatory, isBounded, low_bnd, up_bnd, def, temp);
}

// ----------------------------------------------------------------------------------------
// set_param
// ----------------------------------------------------------------------------------------
bool ParamSet::set_param (string& Name, string& Arg)
{
    map<string, Param*>::iterator param = _params.find(Name);

    if(param != _params.end())
        return param->second->set(Arg, this);
    else
        return false;
}

// ----------------------------------------------------------------------------------------
// set_param
// ----------------------------------------------------------------------------------------
/** returns true if the passed name is a prameter name, otherwise false */
bool ParamSet::check_param_name (const string& Name)
{
    return (_params.find(Name) != _params.end());
}

// ----------------------------------------------------------------------------------------
// set_param
// ----------------------------------------------------------------------------------------
void
ParamSet::set_temporal_param(const int& gen, Param* parm){
    multimap<int, map<string, Param*>* >::iterator pos = _temporalParams.find(gen);
    if(pos != _temporalParams.end()){
        pos->second->insert(pair<string, Param*>(parm->get_name(), parm));      // add the new parm
    }
    else{        // this generation has not yet been used
        map<string, Param*>* secMap = new map<string, Param*>;                  // make the second map with new
        secMap->insert(pair<string, Param*>(parm->get_name(), parm));           // set second map
        _temporalParams.insert(pair<int, map< string, Param*>* >(gen, secMap)); // set multimap
    }

}

// ----------------------------------------------------------------------------------------
// getTemporalParams
// ----------------------------------------------------------------------------------------
/** returns the map if the passed generation time is available, NULL if not */
map<string, Param*>*
ParamSet::getTemporalParams(const int& gen){
    multimap<int, map<string, Param*>* >::iterator pos;
    pos = _temporalParams.find(gen);

    // if present return the map, otherwise return NULL
    if(pos != _temporalParams.end()) return pos->second;
    return NULL;
}

// ----------------------------------------------------------------------------------------
// updateTemporalParams
// ----------------------------------------------------------------------------------------
/** updates the current value by the temporal value
 * returns the map with the params if a paramter was updated, NULL if no parameter was updated
 */
map<string, Param*>*
ParamSet::updateTemporalParams(const int& gen){
    multimap<int, map<string, Param*>* >::iterator pos;
    pos = _temporalParams.find(gen);

    // if not present go on
    if(pos == _temporalParams.end()) return NULL;

    // if present update the parameters and return true
    map<string, Param*>::iterator map_pos = pos->second->begin();
    for(; map_pos != pos->second->end(); ++map_pos){
        map_pos->second->update_arg(gen);
    }
    return pos->second;
}

// ----------------------------------------------------------------------------------------
// get_param
// ----------------------------------------------------------------------------------------
/** returns a pointer to the param if found and an error if not found */
Param* ParamSet::get_param (string Name)
{
    map<string, Param*>::iterator param = _params.find(Name);

    if(param == _params.end()) fatal("ParamSet::get_param could not find '%s'!", Name.c_str());
    return param->second;
}

// ----------------------------------------------------------------------------------------
// find_param
// ----------------------------------------------------------------------------------------
/** returns a pointer to the param if found and NULL if not found */
Param* ParamSet::find_param (string Name)
{
    map<string, Param*>::iterator param = _params.find(Name);

    if(param == _params.end()) return NULL;
    return param->second;
}
// ----------------------------------------------------------------------------------------
// check_consistency
// ----------------------------------------------------------------------------------------
bool ParamSet::check_consistency ()
{
    map<string, Param*>::iterator param = _params.begin();
    bool isOK = true;
    // bool touched = false;
    bool firstOK = true;

    for(int i=0; param != _params.end(); ++param, ++i) {
        //check if all required fields have been set properly
        if(!i)firstOK = param->second->isSet();
        if(param->second->isRequired()) {
            isOK &= param->second->isSet();
            if(!param->second->isSet() && firstOK){
                warning("Parameter set '%s' is disabled as parameter '%s' is required but not set!\n",
                        _params.begin()->second->get_name().c_str(), param->second->get_name().c_str());
            }
        }
        //else we don't care..
        //check if at least one param has been set
        // touched |= param->second->isSet();
    }
    _isSet = isOK;
    //return isOk or check if _isRequired in case no params are set (untouched params)
    return ( isOK | (!_isRequired));// & !touched) );
}
// ----------------------------------------------------------------------------------------
// show_up
// ----------------------------------------------------------------------------------------
void ParamSet::show_up ()
{
    message("\n%s",_name.c_str());
    map<string, Param*>::iterator param = _params.begin();
    while(param != _params.end()) {
        param->second->show_up();
        param++;
    }
}
// ----------------------------------------------------------------------------------------
// print
// ----------------------------------------------------------------------------------------
void ParamSet::print (ostream& FILE)
{
    map<string, Param*>::iterator param = _params.begin();
    while(param != _params.end()) {
        if(param->second->isSet()){
            FILE << param->second->get_name() << " " << param->second->get_tot_arg() << "\n";
        }
        param++;
    }
}

// ----------------------------------------------------------------------------------------
// print
// ----------------------------------------------------------------------------------------
void ParamSet::print_minimal (ostream& FILE)
{
    map<string, Param*>::iterator param = _params.begin();
    while(param != _params.end()) {
        if(param->second->isSet() && param->second->get_tot_arg() != param->second->get_default_arg()){
            FILE << param->second->get_name() << " " << param->second->get_tot_arg() << "\n";
        }
        param++;
    }
}

// ----------------------------------------------------------------------------------------
// print
// ----------------------------------------------------------------------------------------
void ParamSet::print_maximal (ostream& FILE)
{
    string name;

    // name of the param set
    FILE << "#### " << getName() << " ####\n";

    map<string, Param*>::iterator param = _params.begin();
    while(param != _params.end()) {
        // parameter + argument
        name = param->second->get_name() + "  ";
        if(param->second->isSet()) name += param->second->get_tot_arg();
        else                       name += param->second->get_default_arg();
        FILE.width(40);
        FILE.setf(ios::left,ios::adjustfield);
        FILE << name;

        // argument type
        FILE << "# type: ";
        name = param->second->get_type_str() + ";";
        FILE.width(16);
        FILE << name;

        // default value
        if(param->second->get_default_arg().empty()) name = "\"\";";
        else name = param->second->get_default_arg() + ";";
        FILE << " default: ";
        FILE.width(16);
        FILE << name;

        // temporal parameter
        FILE.width(11);
        if(param->second->isTemporalArgumentAllowed()) FILE << " temporal; ";
        else                                           FILE << "";

        // range
        if(param->second->isBounded()) FILE << " range: [" << param->second->get_bound(0)
            << "-"<< param->second->get_bound(1) << "];";

        FILE << "\n";

        param++;
    }
}
