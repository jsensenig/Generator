//____________________________________________________________________________
/*!

\program testXSec

\brief   test program used for testing/debugging differential xsec algorithms

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 20, 2004

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/DeepInelastic/XSection/DISStructureFunc.h"
#include "Physics/DeepInelastic/XSection/DISStructureFuncModelI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Interaction/InitialState.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"

#include <fstream>

using namespace genie;
using namespace genie::constants;
using namespace genie::units;
using namespace genie::utils;

void testRES();
void testRES2();
void testDeltaDecay() ;
//____________________________________________________________________________
int main(int argc, char** argv )
{
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);
  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gevgen", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();


  testDeltaDecay() ; 
  return 0;
}

//____________________________________________________________________________
void testDeltaDecay()  {

  Messenger * msg = Messenger::Instance();
  msg -> SetPriorityLevel( "ResonanceDecay", pDEBUG ) ;
  

  AlgFactory * algf = AlgFactory::Instance();
  const Algorithm * algo = algf -> GetAlgorithm("genie::BaryonResonanceDecayer", "BeforeHadronTransport") ;

}
//____________________________________________________________________________


void testRES()
{
  Messenger * msg = Messenger::Instance();
  msg->SetPriorityLevel("ReinSeghalRes", pDEBUG);
  msg->SetPriorityLevel("RSHAmpl", pDEBUG);
  msg->SetPriorityLevel("FKR", pDEBUG);

  double      Ev  = 1.8;
  double      Q2  = .1;
  double      W   = 1.5;

  const int nres=18;
  Resonance_t res[nres] = {
       kP33_1232,
       kS11_1535,
       kD13_1520,
       kS11_1650,
       kD13_1700,
       kD15_1675,
       kS31_1620,
       kD33_1700,
       kP11_1440,
       kP33_1600,
       kP13_1720,
       kF15_1680,
       kP31_1910,
       kP33_1920,
       kF35_1905,
       kF37_1950,
       kP11_1710,
       kF17_1970
  };
/*
  const int nres=1;
  Resonance_t res[nres] = {
       kP33_1232
  };
*/
  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * xsec_cc =
      dynamic_cast<const XSecAlgorithmI *> (
                   algf->GetAlgorithm("genie::ReinSeghalRESPXSec", "Default"));

  Interaction * interaction = Interaction::DISCC(kPdgTgtFreeN,kPdgNeutron,kPdgNuMu);
  interaction->InitStatePtr()->SetProbeE(Ev);
  Kinematics * kine = interaction->KinePtr();
  kine->SetW(W);
  kine->SetQ2(Q2);

  LOG("Main", pNOTICE) << *interaction;

  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  double sum=0;
  for(int i=0; i<nres; i++) {
     interaction->ExclTagPtr()->SetResonance(res[i]);
     double xsec = xsec_cc->XSec(interaction, kPSWQ2fE) / units::cm2;
     LOG("Main", pNOTICE) << "d2xsec/dWdQ2 [" << res::AsString(res[i]) << "] = " << xsec;
     sum+=xsec;
  }
  LOG("Main", pNOTICE) << "d2xsec/dWdQ2 [SUM] = " << sum;

//  LOG("Main", pNOTICE) << "integral{d2xsec/dWdQ2} = " << xsec_cc->Integral(interaction) / units::cm2;

  delete interaction;
}
//____________________________________________________________________________
void testRES2()
{
  Messenger * msg = Messenger::Instance();
  msg->SetPriorityLevel("ReinSeghalRes", pDEBUG);
  msg->SetPriorityLevel("RSHAmpl", pDEBUG);
  msg->SetPriorityLevel("FKR", pDEBUG);

  double Ev  = 1.8;

  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * xsec_cc =
      dynamic_cast<const XSecAlgorithmI *> (
                   algf->GetAlgorithm("genie::ReinSeghalRESPXSec", "Default"));

  Interaction * interaction = Interaction::RESCC(kPdgTgtFreeP,kPdgProton,kPdgNuMu);
  interaction->InitStatePtr()->SetProbeE(Ev);
  interaction->ExclTagPtr()->SetResonance(kP33_1232);

  int    nW = 500;
  double Wmin = 1.1;
  double Wmax = 1.7;
  int    nQ2 = 500;
  double Q2min = 0.001;
  double Q2max = 2.0;
//  double logQ2min = TMath::Log10(Q2min);
//  double logQ2max = TMath::Log10(Q2max);

  TH2D * h = new TH2D("h","",nW,Wmin,Wmax,nQ2,Q2min,Q2max);

  LOG("Main", pNOTICE) << *interaction;

  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  for(int i=0; i<nW; i++) {
   for(int j=0; j<nQ2; j++) {
     double W   = h->GetXaxis()->GetBinCenter(i+1);
     double Q2  = h->GetYaxis()->GetBinCenter(j+1);
     Kinematics * kine = interaction->KinePtr();
     kine->SetW(W);
     kine->SetQ2(Q2);
     double xsec = xsec_cc->XSec(interaction, kPSWQ2fE) / units::cm2;
//     LOG("Main", pNOTICE) << "d2xsec/dWdQ2 [" << res::AsString(res[i]) << "] = " << xsec;

     h->Fill(W,Q2,xsec);
   }
  }

  TFile f("res.out","recreate");
  h->Write();
  f.Close();

  delete interaction;
}
//____________________________________________________________________________
