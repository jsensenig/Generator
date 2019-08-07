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
#include "Physics/Hadronization/HadronizationModelI.h"
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
void testDIS();
void testQelCharm();
void testCOHGamma() ;
void testNCgamma();

//____________________________________________________________________________
int main(int argc, char** argv )
{
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);
  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gevgen", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  testNCgamma();
  return 0;
}


//____________________________________________________________________________
void testNCgamma()
{
  std::ofstream file;
  file.open("diff_XSec.dat");
	
  Messenger * msg = Messenger::Instance();

  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * xsec_model =
      dynamic_cast<const XSecAlgorithmI *> (
             algf->GetAlgorithm("genie::AlvarezRusoCOHPiPXSec", "Default"));


	  
//    ProcessInfo  proc_info(kScCoherent, kIntWeakNC);
// 
//    InitialState init_state(1000060120,14);
//    init_state.SetProbeE(1);
// 
	
  double Ev = 1. ;
  Interaction * interaction = Interaction::COHNC( 1000060120, 14, Ev ) ; 

  Kinematics * kine = interaction -> KinePtr() ;
   
   double gamma_th   = 10.0 * 3.14159/180.0;
   double gamma_phi  = 0.0 * 3.14159/180.0;
   double lepton_th  = 10.0 * 3.14159/180.0;
   double lepton_phi = 10.0 * 3.14159/180.0;
   
   double lepton_E_max   = Ev;
   double lepton_E_min = 0.0;
   double lepton_E = lepton_E_max;
   double npts = 50;
   
   std::cout<<std::endl<<"******************************************************"<<std::endl;
   for(int i = 0; i < npts; i++) {
	    
		
		double gamma_E = Ev - lepton_E ; 
		
		// convert in proper TLorentzVectors

		TLorentzVector lep_p( lepton_E * TMath::Sin( lepton_th ) * TMath::Cos( lepton_phi) , 
								lepton_E * TMath::Sin( lepton_th ) * TMath::Sin( lepton_phi) , 
								lepton_E * TMath::Cos( lepton_th ) , 
								lepton_E ) ;
								
		TLorentzVector gamma_p( gamma_E * TMath::Sin( gamma_th ) * TMath::Cos( gamma_phi) , 
								gamma_E * TMath::Sin( gamma_th ) * TMath::Sin( gamma_phi), 
								gamma_E * TMath::Cos( gamma_th ) , 
								gamma_E ) ;
		
		
		kine -> SetFSLeptonP4( lep_p ) ;
		kine -> SetHadSystP4( gamma_p ) ; 
		
		double val = xsec_model -> XSec( interaction, kPSElOlOpifE ) ;
		double units = 197.33e-16 * 197.33e-16;
		val *= units;
		
		std::cout<<"k0 = " << Ev <<"  ;    theta = "<< lepton_th <<"  ;    theta_g = "<< gamma_th <<"  ;    phi_g = "<< gamma_phi << "  ;    ";
		std::cout<<"E_gamma = " << gamma_E <<"  ;    XSec = "<<val<<std::endl;
		file<<gamma_E<<"  "<<val<<std::endl;
		
		lepton_E -= (lepton_E_max - lepton_E_min)/double(npts-1);
		
   }
  std::cout<<"******************************************************"<<std::endl<<std::endl;
  
  file.close();
  
  delete interaction;
}
//____________________________________________________________________________

//____________________________________________________________________________
void testQelCharm()
{
  Messenger * msg = Messenger::Instance();
  msg->SetPriorityLevel("QELCharmXSec", pDEBUG);

  AlgFactory * algf = AlgFactory::Instance();
  const XSecAlgorithmI * xsec_model =
      dynamic_cast<const XSecAlgorithmI *> (
             algf->GetAlgorithm("genie::KovalenkoQELCharmPXSec", "Default"));


   ProcessInfo  proc_info(kScQuasiElastic, kIntWeakCC);

   InitialState init_state(1000010010,14);
   init_state.SetProbeE(5);

   Interaction * interaction = new Interaction(init_state, proc_info);

   Target * target  = interaction->InitStatePtr()->TgtPtr();
   XclsTag * xcls = interaction->ExclTagPtr();

   target->SetHitNucPdg(kPdgProton);
   xcls->SetCharm(kPdgSigmaPPc);

   LOG("Main", pNOTICE)
     << "integral{dxsec/dQ2} = " << xsec_model->Integral(interaction) / (1E-40*units::cm2) << " 1E-40 cm2";;

  delete interaction;
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
void testDIS()
{
  Messenger * msg = Messenger::Instance();
//  msg->SetPriorityLevel("DISPXSec",  pDEBUG);
//  msg->SetPriorityLevel("DISXSec",  pDEBUG);
//  msg->SetPriorityLevel("DISSF",     pDEBUG);
//  msg->SetPriorityLevel("BodekYang", pDEBUG);
//  msg->SetPriorityLevel("PDFLIB",    pDEBUG);
/*
  double Ev = 10;
  double Q2 = 2;
  double x  = 0.7;
  double y  = Q2 / (2.*x*kNucleonMass*Ev);
*/
  double Ev = 80;
  double x  = 0.5;
  double y  = 0.5;
  double Q2 = 2*x*y*Ev*kNucleonMass;
  double W2 = kNucleonMass2 + 2*Ev*kNucleonMass*y*(1-x);
  double W  = TMath::Sqrt(TMath::Max(0., W2));

  const int Nq = 8;
  int  qrk[Nq] = {kPdgUQuark, kPdgUQuark, kPdgAntiUQuark, kPdgDQuark, kPdgDQuark, kPdgAntiDQuark, kPdgSQuark, kPdgAntiSQuark };
  bool sea[Nq] = {false,      true,       true,           false,      true,       true,           true,       true};

  // -- request the specified DIS SF model
  AlgFactory * algf = AlgFactory::Instance();

  const XSecAlgorithmI * xsec_cc =
      dynamic_cast<const XSecAlgorithmI *> (
                   algf->GetAlgorithm("genie::DISPartonModelPXSec", "CC-Default"));

  PDGLibrary * pdglib = PDGLibrary::Instance();

  double sum=0;

//  Interaction * interaction = Interaction::DISCC(kPdgTgtFe56,kPdgNeutron,kPdgNuMu);
//  Interaction * interaction = Interaction::DISCC(kPdgTgtFreeN,kPdgNeutron,kPdgNuMu);
//  Interaction * interaction = Interaction::DISNC(kPdgTgtFreeP,kPdgProton,kPdgNuMu);
//  Interaction * interaction = Interaction::DISNC(kPdgTgtFe56,kPdgProton,kPdgNuMu);
//  Interaction * interaction = Interaction::DISNC(kPdgTgtFreeP,kPdgProton,kPdgNuMu);
//  Interaction * interaction = Interaction::DISCC(kPdgTgtFreeN,kPdgNeutron,kPdgAntiNuMu);
//  Interaction * interaction = Interaction::DISNC(kPdgTgtFreeN,kPdgNeutron,kPdgNuMu);

//HERE
  Interaction * interaction = Interaction::DISNC(kPdgTgtFreeP,kPdgProton,kPdgAntiNuMu);

  interaction->InitStatePtr()->SetProbeE(Ev);

  Kinematics * kine = interaction->KinePtr();
  kine->Setx(x);
  kine->Sety(y);
  kine->SetQ2(Q2);
  kine->SetW(W);

  LOG("Main", pNOTICE) << *interaction;

  //interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  Interaction * interaction2 = new Interaction(*interaction);

/*
  LOG("Main", pNOTICE) << "d2xsec/dxdy [vN] = " << xsec_cc->XSec(interaction, kPSxyfE) / units::cm2;
  sum=0;
  for(int iq=0; iq<Nq; iq++) {
      LOG("Main", pNOTICE) << " **** USING: " << pdglib->Find(qrk[iq])->GetName() << ", sea:" << sea[iq] << "]";
      interaction2->InitStatePtr()->TgtPtr()->SetHitQrkPdg(qrk[iq]);
      interaction2->InitStatePtr()->TgtPtr()->SetHitSeaQrk(sea[iq]);
      //LOG("Main", pNOTICE) << *interaction;
      double xsec   = xsec_cc->XSec(interaction2, kPSxyfE) / units::cm2;
      LOG("Main", pNOTICE) << "d2xsec/dxdy [q:" << pdglib->Find(qrk[iq])->GetName()
                           << ", sea:" << sea[iq] << "]= " << xsec;
      sum+=xsec;
  }
  LOG("Main", pNOTICE) << "d2xsec/dxdy [SUM]= " << sum;
*/

  msg->SetPriorityLevel("DISXSec",   pNOTICE);
  msg->SetPriorityLevel("DISPXSec",  pNOTICE);
  msg->SetPriorityLevel("DISSF",     pNOTICE);
  msg->SetPriorityLevel("BodekYang", pNOTICE);
  msg->SetPriorityLevel("PDFLIB",    pNOTICE);

 interaction ->SetBit(kINoNuclearCorrection);
 interaction2->SetBit(kINoNuclearCorrection);

  LOG("Main", pNOTICE) << "Integral{d2xsec/dxdy [vN]} = " << xsec_cc->Integral(interaction) / units::cm2;
  sum=0;
  for(int iq=0; iq<Nq; iq++) {
      LOG("Main", pNOTICE) << " **** USING: " << pdglib->Find(qrk[iq])->GetName() << ", sea:" << sea[iq] << "]";
      interaction2->InitStatePtr()->TgtPtr()->SetHitQrkPdg(qrk[iq]);
      interaction2->InitStatePtr()->TgtPtr()->SetHitSeaQrk(sea[iq]);
      //LOG("Main", pNOTICE) << *interaction2;
      double xsec =  xsec_cc->Integral(interaction2) / units::cm2;
      LOG("Main", pNOTICE) << "integral{d2xsec/dxdy [q:" << pdglib->Find(qrk[iq])->GetName()
                           << ", sea:" << sea[iq] << "]}= " << xsec;
      sum+=xsec;
  }
  LOG("Main", pNOTICE) << "integral{d2xsec/dxdy} [SUM]= " << sum;

  delete interaction;
  delete interaction2;
}
//____________________________________________________________________________
