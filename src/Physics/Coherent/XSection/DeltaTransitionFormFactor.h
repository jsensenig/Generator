//____________________________________________________________________________
/*!

\class    genie::DeltaTransitiontFormFactor

\brief    Transition form factor from Nucleon to Delta
          To be used for the evaluation of the Cross section of 
          COH NC Gamma production


\author  Marco Roda <mroda@liverpool.ac.uk>
         University of Liverpool

\created November 2019

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DELTA_TRANSITION_FORM_FACTOR_H_
#define _DELTA_TRANSITION_FORM_FACTOR_H_ 

#include "Framework/Algorithm/Algorithm.h"

namespace genie {

class DeltaTransitionFormFactor : public Algorithm {

public:
  DeltaTransitionFormFactor();
  DeltaTransitionFormFactor(string config);
  virtual ~DeltaTransitionFormFactor();

  double C3V( double Q2 ) const ;
  double C3VNC( double Q2 ) const ;
  double C5ANC( double Q2 ) const ;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);


private:

  void LoadConfig(void);

};

}       // genie namespace
#endif  // _DELTA_TRANSITION_FORM_FACTOR_H_ 
