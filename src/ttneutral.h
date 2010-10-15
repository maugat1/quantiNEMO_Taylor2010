/** @file ttneutral.h
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



#ifndef ttneutralH
#define ttneutralH

#include "ttrait.h"
#include "types.h"
#include "filehandler.h"
#include "stathandler.h"
#include "metapop.h"

class TTNeutralFH;
class TTNeutralSH;
class TTNtrlPhenotyperFH;
class TTNeutralProto;

/**Microsatellites genome.*/
class TTNeutral : public TTrait
{

  TTNeutralProto* pProto;

private:

public:
  virtual void set_from_prototype(TTraitProto* T);

  TTNeutral () : pProto(0){}

  TTNeutral(const TTNeutral& T) : pProto(T.pProto){
    _copyTTraitParameters(T);  // copy the parameters of TTrait
  }

  virtual ~TTNeutral ();

  ///@}
  ///@name Implementations
  ///@{
  virtual TTNeutral& operator= (const TTrait& T);
  virtual bool    operator== (const TTrait& T);
  virtual bool    operator!= (const TTrait& T);
  virtual void    init                 ( );
  virtual void    ini_sequence         (Patch* patch);
  virtual void    reset                ( );
  virtual void*   set_trait            (void* value)           {return NULL;}
  inline virtual void**  get_sequence  ( )  const              {return (void**)sequence;}
  inline virtual void*   get_allele    (int loc, int all)  const;
  virtual void    set_sequence         (void** seq)            {reset();sequence = (unsigned char**)seq;}
  virtual void    set_value            ( )                     { }
  virtual void    set_value            (double value)          {return;}
  virtual double  get_value            ( )					           {return my_NAN;}
  virtual void    show_up              ( );
  virtual TTNeutral*  clone       ( )                      {return new TTNeutral(*this);}
  ///@}
};

/**Prototype class for the TTNeutral trait class.**/
class TTNeutralProto : public TTraitProto {
friend class TTNeutral; // we allow to access these parameters from TTNeutral directly
private:
  void ini_paramset();



public:

  static string _ini_fstat_file; // if the inital genotypes are given

  static TTNeutralFH* _writer;
  TTNeutralSH* _stats;


  TTNeutralProto ( );
  TTNeutralProto (int i);
  TTNeutralProto(const TTNeutralProto& T);

  virtual ~TTNeutralProto ( );

  //implementation of TTraitProto:
  virtual void                     init (Metapop* pMetapop);

  virtual TTNeutral*          hatch ();

  virtual TTNeutralProto*      clone () {return new TTNeutralProto(*this);}

  //implementation of SimComponent:
  virtual void loadFileServices ( FileServices* loader );

	virtual void loadStatServices ( StatServices* loader );

	string get_info();

  virtual void executeBeforeEachGeneration(const int& gen);
  virtual void executeBeforeEachReplicate(const int& rep){
    if(rep!= 1 && _random_genetic_map_matrix) setGeneticMapRandom(_random_genetic_map_matrix);
  }
};

////////////////////////////////////////////////////////////////////////////////
/**A file handler to save the neutral markers genotypes in the FSTAT format*/
class TTNeutralFH: public FileHandler {

public:

  TTNeutralFH (bool rpl_per = 0, bool gen_per = 0, int rpl_occ = 0, int gen_occ = 0, TTNeutralProto* TP = NULL)
    : FileHandler (".dat")
  {
    FileHandler::set(rpl_per,gen_per,rpl_occ,gen_occ, 0, "", "", 0,0,TP,1);
  }

  virtual ~TTNeutralFH ( ) { }

};

////////////////////////////////////////////////////////////////////////////////
/**The stat handler for neutral markers. */
class TTNeutralSH: public StatHandler<TTNeutralSH> {

public:

	TTNeutralSH (TTNeutralProto* TT){
    set(TT);
  }

  virtual ~TTNeutralSH ( ){}

  virtual bool init ( ) ;

  virtual bool setStatRecorders (const string& token);

};



#endif //TTNEUTRALGENES_H

