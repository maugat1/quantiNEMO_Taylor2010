/** @file param.h
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

#ifndef paramH
#define paramH

#include <map>
#include "types.h"
//#include "output.h"
#include "tmatrix.h"

/**This structure stores one parameter, its definition and its string argument. **/
class ParamSet;
class Param
{

    private:
        string          _name;
        string          _arg;         // contains the current argument
        string          _tot_arg;     // contains as well the temporal parameters
        string          _default_arg; // contains the default argument
        param_t         _type;
        bool            _isSet;
        bool            _isBounded;
        bool            _isRequired;
        double          _bounds[2];
        TMatrix*        _matrix;
        TMatrixVar<double>*     _matrixVar;   // with variable size

        ParamSet*               _parSet;

        // parameters for temporal changes
        bool             _temporalParamAllowed;
        map<int, string> _temporalArgs;

    public:
        /**
         * @param Name the name of the parameter as read in the init file
         * @param Type the type of the parameter argument (see types.h), used to convert the argument string into a value
         * @param mandatory specifies whether this parameter is mandatory for the ParamSet owning it
         * @param bounded specifies whether the values this parameter can take are bounded
         * @param low_bnd the lower bound of the parameter value
         * @param up_bnd the upper bound
         * @param temp can tis parameter change over time?
         **/
        Param  (string& Name,param_t Type,bool mandatory = 0,bool bounded = 0,
                double low_bnd = 0,double up_bnd = 0, string def_arg = "",  bool temp = 0);

        ~Param () {
            if(_matrix)     delete _matrix;
            if(_matrixVar)  delete _matrixVar;
        }

        /**@brief clears the _isSet flag and the argument string**/
        void    reset   ();

        //accessors:
        /**@brief Sets the parameter's name**/
        void            set_name             (string value)          {_name = value;}

        /**@brief Sets the parameter's argument**/
        void            set_arg              (string value)          {_arg = value;}

        /**@brief Sets the parameter's type (see types.h)**/
        void            set_type             (param_t value)         {_type = value;}

        /**@brief Sets the _isSet flag**/
        void            set_isSet            (bool value)            {_isSet = value;}

        /**@brief Sets the _isBounded flag**/
        void            set_isBounded        (bool value)            {_isBounded = value;}

        /**@brief Sets the bounds**/
        void            set_bounds   (double low_bnd,double up_bnd)  {_bounds[0]=low_bnd; _bounds[1]=up_bnd;}

        /** the argument is updated by the new temporal argument */
        void            update_arg   (const int& gen);



        string&         get_name             ()                      {return _name;}
        string&         get_arg              ()                      {return _arg;}
        string&         get_default_arg      ()                      {return _default_arg;}
        string&         get_tot_arg          ()                      {return _tot_arg;}
        map<int,string>* get_temporal_args   ()                      {return &_temporalArgs;}
        param_t         get_type             ()                      {return _type;}
        bool            isSet                ()                      {return _isSet;}
        bool            isBounded            ()                      {return _isBounded;}
        bool            isRequired           ()                      {return _isRequired;}
        double          get_bound            (unsigned int i)        {return _bounds[i];}
        bool            isTemporalParam      ()                      {return !_temporalArgs.empty();}
        bool            isTemporalArgumentAllowed()                  {return _temporalParamAllowed;}

        /**@brief Sets the _isSet flag to true and _arg to arg if the arg is of the right type and whithin the bounds**/
        bool            set                  (string& arg, ParamSet* parSet);
        bool            set_temporal_args    (string& arg);   // set temporal arguments
        bool            check_arg            (string& arg);

        /**@brief Get the argument value according to its type
         *@return my_NAN if the parameter is not set or not of a the right type
         **/
        double          get_value            ();

        /**@brief Checks if the argument is of matrix type**/
        bool            is_matrix            () {return (_type == MAT || (_arg.size() ? _arg[0] == '{' : false));}

        /**@brief Gets the matrix from the argument string if the parameter is set and of matrix type
         *@return 0 if parameter not set or of the wrong type
         **/
        TMatrix*        get_matrix           ();

        /**@brief Parses the matrix from the argument string
         *@return a TMatrix ptr or abort on error
         **/
        TMatrix*        parse_matrix         ();

        /**@brief Checks if the argument is of matrix type**/
        bool            is_matrixVar            () {return (_type == MAT_VAR || (_arg.size() ? _arg[0] == '{' : false) );}

        /**@brief Gets the matrix from the argument string if the parameter is set and of matrix type
         *@return 0 if parameter not set or of the wrong type
         **/
        TMatrixVar<double>*     get_matrixVar           ();

        /**@brief Parses the matrix from the argument string
         *@return a TMatrix ptr or abort on error
         **/
        TMatrixVar<double>*     parse_matrixVarDbl         ();        // variable double matrix
        TMatrixVar<string>*     parse_matrixVarStr         ();        // variable string matrix

        /**@brief Print state to stdout**/
        void            show_up              ();
        string          get_type_str         ();

        ParamSet*       getParamSet          ()                      {return _parSet;}

};


/**Parameters container, implemented in each SimComponent. **/
class ParamSet
{
    private:
        string                  _name;
        bool                    _isSet;
        bool                    _isRequired;
        map<string, Param*>     _params;


        // temporal parameter
        multimap<int, map<string, Param*>* > _temporalParams;

    public:

        ParamSet                       ( );
        ParamSet                       (string name, bool isRequired);
        ~ParamSet                      ( );

        /**@brief Put the container in the unset state, reset each Param it contains**/
        void            reset                ( );

        /**@brief Sets the container's name**/
        void            set_name             (string value)          {
            _name = value;
        }

        /**@brief Sets the _isRequired flag meaning this container is mandatory and must be set in order to run a simulation**/
        void            set_isRequired       (bool value)            {_isRequired = value;}

        /**@brief Checks for the status of the required parameters
         *@return TRUE if all required parameters are or if nothing is set and the container is not required
         **/
        bool            check_consistency    ( );

        /**@brief print info to stdout**/
        void            show_up              ( );

        /**@brief print all set parameters to the outpout file stream**/
        void            print                (ostream& FILE);
        void            print_minimal        (ostream& FILE);
        void            print_maximal        (ostream& FILE);

        /**@brief Returns the number of parameters contained**/
        int             size                 ( )                     {return _params.size();}

        /**@brief Returns the complete list of parameters**/
        map<string, Param*>& getAllParams   ( )                     {return _params;}

        /** returns the map if the passed generation time is available, NULL if not */
        multimap<int, map<string, Param*>* >* getTemporalParams() {return &_temporalParams;}
        map<string, Param*>* getTemporalParams(const int& gen);
        map<string, Param*>* updateTemporalParams(const int& gen);

        ///@name Accessors to Param members.
        ///@{
        /**@brief Adds the param argument to the list**/
        void            add_param            (Param* param);

        /**@brief Adds a new param specified by arguments to the list.
         *@param Name the name of the parameter
         *@param Type the type of the parameter
         *@param mandatory specifies if this parameter is required and must be set for the container to gain the "set" status
         *@param isBounded specified whether this parameter is bounded
         *@param low_bnd the lower value the parameter can take, used if isBounded is true
         *@param up_bnd the upper value the parameter can take, used if isBounded is true
         *@param default_param the default value
         *@param temp can the parameter change over time?
         **/
        void            add_param            (string Name,param_t Type,bool mandatory,
                bool isBounded,double low_bnd,double up_bnd,
                string default_param, bool temp=false);

        /**@brief Look for a param named "Name" and try to set it with the "Arg" argument string.
         *@return TRUE if param Name has been found and set with Arg
         *@return FALSE otherwise
         *@param Name the name of the parameter to find in the list
         *@param Arg the argument string as found in the init params
         **/
        bool            set_param            (string& Name, string& Arg);
        bool            check_param_name     (const string& Name);
        void            set_temporal_param   (const int& gen, Param* parm);  // for temporal chainging parameters
        bool            set_general_parameter(map< string,string >* input);

        /**@brief Look for a param "name" in its parameters list.
         *@return NULL if no Param with _name = name exists
         **/
        Param*         get_param            (string name);     // returns an error if not found
        Param*         find_param           (string name);     // returns NULL if not found

        /**@brief Accessor to the status flag. **/
        bool            isSet                ()                      {return _isSet;}

        /**@brief Accessor to the parameters status flag. **/
        bool            isSet                (string name)           {return (get_param(name))->isSet();}

        /**@brief Accessor to the parameters status flag. **/
        void            set_isSet            (bool choice)           {_isSet = choice;}

        /**@brief Accessor to the mandatory flag. **/
        bool            isRequired           ()                      {return _isRequired;}

        bool            isTemporalParamSet   ()                      {return !_temporalParams.empty();}

        /**@brief Check if the parameter "name" is of matrix type. **/
        bool            isMatrix             (string name)           {return (get_param(name))->is_matrix();}

        /**@brief Name accessor.**/
        string&         getName              ()                      {return _name;}

        /**@brief Accessor to the parameters argument string.**/
        string&         getArg               (string name)           {return (get_param(name))->get_arg();}

        /**@brief Accessor the parameters value.**/
        double          getValue             (string name)           {return (get_param(name))->get_value();}

        /**@brief Accessor to the parameters matrix.**/
        TMatrix*        getMatrix            (string name)           {return (get_param(name))->get_matrix();}


        ///@}
};

#endif
