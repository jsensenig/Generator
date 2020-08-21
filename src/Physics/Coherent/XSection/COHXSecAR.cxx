//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________
#include <chrono>

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Coherent/XSection/COHXSecAR.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

using namespace std::chrono ;
//____________________________________________________________________________
COHXSecAR::COHXSecAR() :
XSecIntegratorI("genie::COHXSecAR")
{

}
//____________________________________________________________________________
COHXSecAR::COHXSecAR(string config) :
XSecIntegratorI("genie::COHXSecAR", config)
{

}
//____________________________________________________________________________
COHXSecAR::~COHXSecAR()
{

}

//____________________________________________________________________________

double COHXSecAR::Integrate( const XSecAlgorithmI * model, const Interaction * in) const {

  if      ( fHasPion )   return IntegratePion( model, in ) ;
  else if ( fHasPhoton ) return IntegratePhoton( model, in ) ;
  
  return 0. ;

}
//____________________________________________________________________________
double COHXSecAR::IntegratePion(
      const XSecAlgorithmI * model, const Interaction * in) const
{
  const InitialState & init_state = in -> InitState();
  
  if(! model->ValidProcess(in) ) return 0.;
  
  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("COHXSecAR", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  
  Range1D_t y_lim = kps.Limits(kKVy);
  
  // Check this
  double Enu      = init_state.ProbeE(kRfLab);
  double Elep_min = (1.-y_lim.max) * Enu;
  double Elep_max = (1.-y_lim.min) * Enu;
  
  LOG("COHXSecAR", pINFO)
       << "Lepton energy integration range = [" << Elep_min << ", " << Elep_max << "]";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);
  
  double xsec = 0;
  if (fSplitIntegral) {
    
    utils::gsl::dXSec_dElep_AR_pion func(model, interaction, fGSLIntgType, fGSLRelTol, fGSLMaxEval);
    
    //~ ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kNONADAPTIVE;
    ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kADAPTIVE;
    
    double abstol = 1; // Pretty sure this parameter is unused by ROOT.
    int size = 1000;  // Max number of subintervals, won't reach nearly this.
    int rule = 2; // See https://www.gnu.org/software/gsl/manual/gsl-ref_17.html#SEC283
                  // Rule 2 is 21 points min
    ROOT::Math::Integrator ig( func,ig_type,abstol,fGSLRelTol,size,rule);
    
    xsec = ig.Integral(Elep_min, Elep_max) * (1E-38 * units::cm2);
  }
  else {
    double zero    = kASmallNum;
    double pi      = kPi-kASmallNum ;
    double twopi   = 2*kPi-kASmallNum ;
    
    //~ ROOT::Math::IBaseFunctionMultiDim * func = 
          //~ new utils::gsl::wrap::d5Xsec_dEldOmegaldOmegapi(model, interaction);
    //~ double kine_min[5] = { Elep_min, zero , zero    , zero, zero };
    //~ double kine_max[5] = { Elep_max, pi   , twopi   , pi  , twopi};
    
    ROOT::Math::IBaseFunctionMultiDim * func = 
          new utils::gsl::d4Xsec_dEldThetaldOmegapi(model, interaction);
    double kine_min[4] = { Elep_min, zero , zero    , zero    };
    double kine_max[4] = { Elep_max, pi   , pi      , twopi   };
    
    ROOT::Math::IntegrationMultiDim::Type ig_type = 
      utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
        
    double abstol = 1; //We mostly care about relative tolerance.
    ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
  
    xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
    delete func;
  }

  delete interaction;

  return xsec;
}
//____________________________________________________________________________
double COHXSecAR::IntegratePhoton( const XSecAlgorithmI * model, const Interaction * in) const {

  const InitialState & init_state = in -> InitState();

  if(! model->ValidProcess(in) ) { return 0.; }

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
    LOG("COHXSecAR", pDEBUG)  << "*** Below energy threshold";
    return 0;
  }
  
  // The time of this operation is monitored
  steady_clock::time_point start = steady_clock::now();

  // Check this
  double Enu      = init_state.ProbeE(kRfLab);
  double Egamma_min  = 0. ; 
  double Egamma_max = Enu;
  
  // LOG("COHXSecAR", pINFO)
  //      << "Lepton energy integration range = [" << Elep_min << ", " << Elep_max << "]";

  Interaction interaction(*in);
  interaction.SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);
  
  double xsec = 0;
  
  //for the time begin the option of splitting the integral is not there for photon

  // if (fSplitIntegral) {
  //   utils::gsl::dXSec_dElep_AR * func =
  //     new utils::gsl::dXSec_dElep_AR(model, interaction, fGSLIntgType, fGSLRelTol, fGSLMaxEval);
    
  //   //~ ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kNONADAPTIVE;
  //   ROOT::Math::IntegrationOneDim::Type ig_type = ROOT::Math::IntegrationOneDim::kADAPTIVE;
    
  //   double abstol = 1; // Pretty sure this parameter is unused by ROOT.
  //   int size = 1000;  // Max number of subintervals, won't reach nearly this.
  //   int rule = 2; // See https://www.gnu.org/software/gsl/manual/gsl-ref_17.html#SEC283
  //                 // Rule 2 is 21 points min
  //   ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,size,rule);
    
  //   xsec = ig.Integral(Elep_min, Elep_max) * (1E-38 * units::cm2);
  //   delete func;
  // }
  // else {

  double zero    = kASmallNum;
  double pi      = kPi-kASmallNum ;
  double twopi   = 2*kPi-kASmallNum ;
  
  //~ ROOT::Math::IBaseFunctionMultiDim * func = 
  //~ new utils::gsl::wrap::d5Xsec_dEldOmegaldOmegapi(model, interaction);
  //~ double kine_min[5] = { Elep_min, zero , zero    , zero, zero };
  //~ double kine_max[5] = { Elep_max, pi   , twopi   , pi  , twopi};
    
  ROOT::Math::IBaseFunctionMultiDim * func = nullptr ;
  if ( fOmegaIntegral ) func = new utils::gsl::d5Xsec_dEgdOmegaldOmegag(model, & interaction);
  else func = new utils::gsl::d4Xsec_dEgdThetaldThetagdPhig(model, & interaction);

  double min_theta = fOmegaIntegral ? -1. : zero ;
  double max_theta = fOmegaIntegral ?  1. : pi ;
  
  double kine_min[4] = { Egamma_min, min_theta , min_theta, zero    };
  double kine_max[4] = { Egamma_max, max_theta , max_theta, twopi   };
  
  ROOT::Math::IntegrationMultiDim::Type ig_type = 
    utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
  
  double abstol = 1; //We mostly care about relative tolerance.
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
  
  xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2) ;

  if ( fOmegaIntegral ) xsec *= 2 * constants::kPi ;

  delete func;
  //  }

  steady_clock::time_point end = steady_clock::now();

  duration<double> time_span = duration_cast<duration<double>>(end - start);
  LOG( "COHXSecAR", pINFO ) << "The integral was performed in " << time_span.count() << " sec" ;

  return xsec;
}
//____________________________________________________________________________
void COHXSecAR::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSecAR::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHXSecAR::LoadConfig(void)
{
  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("vegas") ) ;

  int max_eval;
  GetParamDef( "gsl-max-eval", max_eval, 4000 ) ;
  fGSLMaxEval    = (unsigned int) max_eval ;

  GetParamDef( "gsl-relative-tolerance", fGSLRelTol,  0.01) ;
  GetParamDef( "split-integral", fSplitIntegral, true ) ;

  GetParamDef( "IsCOHPion",  fHasPion,   false ) ;
  GetParamDef( "IsCOHGamma", fHasPhoton, false ) ;

  if ( fHasPhoton ) {
    GetParamDef( "OmegaPhaseSpace", fOmegaIntegral, true ) ;
  }

  bool error = false ;
  
  if ( !fHasPion && !fHasPhoton ) {
    LOG( "COHXSecAR", pERROR ) << "No pion nor gamma option has been requested" ;
    error = true ;
  }

  if ( fHasPion && fHasPhoton ) {
    LOG( "COHXSecAR", pERROR ) << "Pion and Gamma options have been requested at the same time" ;
    error = true ;
  }

  if ( error ) {
    LOG( "COHXSecAR", pFATAL ) << "Invalid configuration. Exiting" ;
    exit( 78 ) ;
  } 

}
//_____________________________________________________________________________


