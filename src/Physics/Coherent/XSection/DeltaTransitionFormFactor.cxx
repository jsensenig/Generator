//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location  (Elastic
   package)

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/Coherent/XSection/DeltaTransitionFormFactor.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

#include "Framework/Utils/StringUtils.h" 

using namespace genie;


//____________________________________________________________________________
DeltaTransitionFormFactor::DeltaTransitionFormFactor() :
Algorithm("genie::DeltaTransitionFormFactor")
{

}
//____________________________________________________________________________
DeltaTransitionFormFactor::DeltaTransitionFormFactor(string config) :
Algorithm("genie::DeltaTransitionFormFactor", config)
{

}
//____________________________________________________________________________
DeltaTransitionFormFactor::~DeltaTransitionFormFactor()
{

}
//____________________________________________________________________________
double DeltaTransitionFormFactor::FormFactor( double Q ) const {

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
void DeltaTransitionFormFactor::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeltaTransitionFormFactor::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DeltaTransitionFormFactor::LoadConfig(void)
{

  // load R33 parameters
  this -> GetParamVect( "DV-Coefficient", fFBCs ) ;

  GetParam( "DV-Radius", fRadius ) ;
  fRadius *= units::fm ;

  GetParam( "DV-Nucleus", fPDG ) ;
  
  LOG("DeltaTransitionFormFactor", pINFO) << "Loaded " << fFBCs.size() << " coeffictients for nucleus " << fPDG ; 
  
}
//____________________________________________________________________________
