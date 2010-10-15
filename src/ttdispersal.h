/** @file ttdispersal.h
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

#ifndef ttdispersalH
#define ttdispersalH

#include "ttrait.h"
#include "types.h"
#include "stathandler.h"
#include "metapop.h"

class TTDispersalSH;
class TTDispersalProto;

/**Evolving dispersal trait, codes for female (_type = FDISP) or male (_type = MDISP) sex-specific dispersal rates.**/
class TTDispersal : public TTrait
{

  TTDispersalProto* pProto;

private:
  /** One diploid locus coding for a sex-specific dispersal rate.**/
  double** sequenceDisp;     // changed by Samuel 9-10-2006 /this makes life more easy as the functions remain the same as for the other traits

  double _phenotype;

  void (TTDispersal::* _inherit_func_ptr) (TTrait* mother, TTrait* father);

public:
  /** @param sex determines the type of this trait (FDISP for female dispersal, MDISP for male dispersal) **/
  TTDispersal (sex_t sex) : pProto(0){}

  TTDispersal (const TTDispersal& T) : pProto(T.pProto)
  {
    _copyTTraitParameters(T);  // copy the parameters of TTrait
  }

  virtual ~TTDispersal ();

  void  mutate();
  void  inherit(TTrait* mother, TTrait* father);


  virtual void set_from_prototype(TTraitProto* T);



  ///@}
  ///@name Implementations
  ///@{
  virtual void    init                 ( );
  virtual void    init_sequence        ( );
  virtual void    reset                ( )                 {/*init();*/}

  virtual void    set_value            ( )                 {_phenotype = (sequence[0][0] + sequence[0][1])/2.0;}
  virtual void    set_value            (double value)      {return;}
  /** @return the dispersal rate, mean of the 2 alleles beared at this dispersal locus. **/
  virtual void*   get_value            ( )                 {return &_phenotype;}

  virtual void**  get_sequence         ( )  const              {return (void**)sequence;}
//  virtual void*   get_allele           (int loc, int all) const {return ( !(all<2) ? NULL : (void*)&sequence[0][all] );}
  virtual void    set_sequence         (void** seq)        { }
  virtual void    show_up              ( );
  virtual void*   set_trait            (void* value)       {return value;}
  virtual TTDispersal*  clone          ( )                 {return new TTDispersal(*this);}
  virtual TTDispersal& operator= (const TTrait& TP);
  virtual bool operator== (const TTrait& TP);
  virtual bool operator!= (const TTrait& TP);
};
/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/*/

//                               ****** TTDispersalProto ******

/*_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\*/
/**Prototype of the evolving dispersal trait, defines the sex-specific trait type.*/
class TTDispersalProto : public TTraitProto {
friend class TTDispersal; // we allow to access these parameters from TTDispersal directly
public:

  double* _locus_position;
  void (TTDispersal::* _inherit_func_ptr)  (TTrait* mother, TTrait* father);

  TTDispersalProto(sex_t sex);
  TTDispersalProto(sex_t sex, int i);
  TTDispersalProto(const TTDispersalProto& TP);

  ~TTDispersalProto();

  //implements TTraitProto:
  virtual   void              init (Metapop* pMetapop);
  virtual   TTDispersal*     hatch ();
  virtual   TTDispersalProto* clone () {return new TTDispersalProto(*this);}

public:
  /**Mean mutation step. **/
  double _mut_mean;

  /**Initial allele for female dispersal.**/
  double _init_rate_fem;

  /**Initial allele for male dispersal.**/
  double _init_rate_mal;

  /**The gender of the trait, will determine its type.**/
  sex_t   _gender;

  /** The trait's StatHandler. **/
  TTDispersalSH* _stats;

  void init_paramset();

  /** Inheritance procedure, creates a new trait from mother's and father's traits
   * @param mother the mother's trait
   * @param father the father's trait
   **/
  void    (TTDispersalProto::* _inherit_disp_func_ptr)(TTrait* mother, TTrait* father, TTrait* child);
  void    inherit_low          (TTrait* mother, TTrait* father, TTrait* child);
  void    inherit_free         (TTrait* mother, TTrait* father, TTrait* child);

  ///@name Mutation
  ///@{
  void  mutate           (double** seq);
  ///@}

  //implements SimComponent:
  virtual void loadFileServices ( FileServices* loader ) {}
  virtual void loadStatServices ( StatServices* loader );

  virtual void executeAfterEachReplicate(const int& rep){
    if(_random_genetic_map_matrix) setGeneticMapRandom(_random_genetic_map_matrix);
  }
};


#endif //TTDISPERSALGENE_H

