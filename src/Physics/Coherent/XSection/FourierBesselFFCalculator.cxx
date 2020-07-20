//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda
         University of Liverpool

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <cmath>

#include "Physics/Coherent/XSection/FourierBesselFFCalculator.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

using namespace genie;

//____________________________________________________________________________
double DeVriesFormFactor::FormFactor( double Q ) const {

  double qr = Q * fRadius ;

  double aux_sum = 0.0, nu ;

  for (unsigned int i = 0 ; i < fFBCs.size() ; ++i ) {
     nu = i + 1. ;
     double pi_x_i = constants::kPi*nu ;
     aux_sum += pow( -1.0, i )*fFBCs[i]/( ( pi_x_i + qr )*( pi_x_i - qr ) ) ;
 }

 return 4.*constants::kPi*pow( fRadius/units::fm, 3)*aux_sum*(sin(qr)/(qr) ) ;

}
//____________________________________________________________________________
