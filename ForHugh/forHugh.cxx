#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/EventGen/InteractionListGeneratorI.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/RunOpt.h"
#include "Physics/Resonance/XSection/MKSPPPXSec.h"
#include "Framework/Interaction/SppChannel.h"
#include "Framework/Conventions/KineVar.h"

#include <TLorentzVector.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>



namespace genie {
class XSecAlgorithmI;
class Interaction;

class d3XSecMK_dWQ2CosTheta_E: public ROOT::Math::IBaseFunctionMultiDim
{
public:
  d3XSecMK_dWQ2CosTheta_E(const XSecAlgorithmI * m, const Interaction * i, double wcut);
 ~d3XSecMK_dWQ2CosTheta_E();

  unsigned int                        NDim   (void)               const;
  double                              DoEval (const double * xin) const;
  ROOT::Math::IBaseFunctionMultiDim * Clone  (void)               const;

private:
  const XSecAlgorithmI * fModel;
  Interaction *    fInteraction;
//  bool isZero;
  KPhaseSpace * kps;
  double fW;
};

}

using namespace genie;

int main(int argc, char ** argv)
{
    
  double xsec, enu = 1., W = 1.3;
  
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);
  
  if ( ! RunOpt::Instance()->Tune() ) {
     std::cerr << " No TuneId in RunOption\n";
     exit(-1);
  }
  
  RunOpt::Instance()->BuildTune();
  AlgFactory * algf = AlgFactory::Instance();

  XSecAlgorithmI * alg = dynamic_cast<genie::XSecAlgorithmI*>(algf->AdoptAlgorithm("genie::MKSPPPXSec", "NoPauliBlock"));
  
  int tgtpdgc = kPdgTgtDeuterium;                         // ядерная мишень
  
  int nupdgc  = kPdgNuMu;                                // нейтрино
   
  ProcessInfo proc_info(kScResonant, kIntWeakCC);
  SppChannel_t spp_channel = kSpp_vp_cc_10100;
  //SppChannel_t spp_channel = kSpp_vn_cc_10010;
  //SppChannel_t spp_channel = kSpp_vn_cc_01100;
  
  //ProcessInfo proc_info(kScResonant, kIntWeakNC);
  //SppChannel_t spp_channel = kSpp_vp_nc_10010;
  //SppChannel_t spp_channel = kSpp_vp_nc_01100;
  //SppChannel_t spp_channel = kSpp_vn_nc_01010;
  //SppChannel_t spp_channel = kSpp_vn_nc_10001;
  
  
  //int nupdgc  = kPdgAntiNuMu;                          // антинейтрино
  
  //ProcessInfo proc_info(kScResonant, kIntWeakCC);
  //SppChannel_t spp_channel = kSpp_vbn_cc_01001;
  //SppChannel_t spp_channel = kSpp_vbp_cc_01010;
  //SppChannel_t spp_channel = kSpp_vbp_cc_10001;
  
  //ProcessInfo proc_info(kScResonant, kIntWeakNC);
  //SppChannel_t spp_channel = kSpp_vbp_nc_10010;
  //SppChannel_t spp_channel = kSpp_vbp_nc_01100;
  //SppChannel_t spp_channel = kSpp_vbn_nc_01010;
  //SppChannel_t spp_channel = kSpp_vbn_nc_10001;
  
  InitialState init_state = InitialState(tgtpdgc, nupdgc);
    
  Interaction * interaction = new Interaction(init_state, proc_info);
  
  int struck_nucleon = SppChannel::InitStateNucleon(spp_channel);
  
  Target * target = interaction->InitStatePtr()->TgtPtr();
  
  target->SetHitNucPdg(struck_nucleon);
  
  int nucpdg = SppChannel::FinStateNucleon(spp_channel);
  int pipdg  = SppChannel::FinStatePion(spp_channel);

  int nproton  = 0;
  int nneutron = 0;
  int npiplus  = 0;
  int npi0     = 0;
  int npiminus = 0;
  
  if       ( nucpdg == kPdgProton  ) nproton  = 1;
  else if  ( nucpdg == kPdgNeutron ) nneutron = 1;


  if       ( pipdg == kPdgPiP  ) npiplus  = 1;
  else if  ( pipdg == kPdgPi0  ) npi0     = 1;
  else if  ( pipdg == kPdgPiM  ) npiminus = 1;


  XclsTag exclusive_tag;

  exclusive_tag.SetNNucleons (nproton, nneutron);
  exclusive_tag.SetNPions    (npiplus, npi0, npiminus);

  interaction->SetExclTag(exclusive_tag);
  
  TLorentzVector p4(0,0,enu,enu);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);
  
  //const Interaction * in = dynamic_cast<const Interaction *>(interaction);
  double ml   = interaction->FSPrimLepton()->Mass();
  double ml2  = ml*ml;

  
  PDGLibrary * pdglib = PDGLibrary::Instance();
  
  double Mi   = pdglib->Find(SppChannel::InitStateNucleon(spp_channel))->Mass();
  double Mf   = pdglib->Find(SppChannel::FinStateNucleon(spp_channel))->Mass();
  double mpi = pdglib->Find(SppChannel::FinStatePion(spp_channel))->Mass();
  
//  double E_thr = (TMath::Power(Mf + ml + mpi, 2) - Mi*Mi)/2/Mi;
  KPhaseSpace * kps = interaction->PhaseSpacePtr();;
  if (enu < kps->Threshold_RSPP()) return -1;
  
  interaction->InitStatePtr()->SetProbeP4(p4);
  
  double s = Mi*(Mi + 2*enu);
  double sqrt_s = TMath::Sqrt(s);
 
  double Wlow = Mf + mpi;
  double Wup  = sqrt_s - ml;
  
  std::cout << "Wlow = " << Wlow << std::endl;
  std::cout << "Wup = " << Wup << std::endl;
  
  //double W = Wlow + (Wup - Wlow)*0.5;
  
  Range1D_t Wl  = kps->WLim_RSPP();
  if (W < Wl.min || W > Wl.max) return -1;
  
  interaction->InitStatePtr()->SetProbeP4(p4);
  
  d3XSecMK_dWQ2CosTheta_E func(alg, interaction, W ) ; 
  ROOT::Math::IntegrationMultiDim::Type ig_type = ROOT::Math::IntegrationMultiDim::kADAPTIVE;
  ROOT::Math::IntegratorMultiDim ig(ig_type,0,1e-9,1000000000);
  ig.SetFunction(func);
  double kine_min[2] = { 0., 0.};
  double kine_max[2] = { 1., 1.};
  xsec = ig.Integral(kine_min, kine_max);// * (1E-38 * units::cm2);
  
  std::cout << "xsec = " << xsec << std::endl;
}

d3XSecMK_dWQ2CosTheta_E::d3XSecMK_dWQ2CosTheta_E(
								    const XSecAlgorithmI * m, const Interaction * interaction, double  W) :
  ROOT::Math::IBaseFunctionMultiDim(),
  fModel(m),
  fW(W)
{

//  isZero = false;
  fInteraction = const_cast<Interaction*>(interaction);
  fInteraction->SetBit(kISkipProcessChk);
  fInteraction->SetBit(kISkipKinematicChk);
  
  kps = fInteraction->PhaseSpacePtr();
    
}
d3XSecMK_dWQ2CosTheta_E::~d3XSecMK_dWQ2CosTheta_E()
{

}
unsigned int d3XSecMK_dWQ2CosTheta_E::NDim(void) const
{
  return 2;
}
double d3XSecMK_dWQ2CosTheta_E::DoEval(const double * xin) const
{
  
  fInteraction->KinePtr()->SetW(fW);
   
  Range1D_t Q2l = kps->Q2Lim_W_RSPP(); 
   
  double Q2 = Q2l.min + (Q2l.max - Q2l.min)*xin[0];
  fInteraction->KinePtr()->SetQ2(Q2);
  
  fInteraction->KinePtr()->SetKV(kKVctp, -1. + 2.*xin[1]); //CosTheta
  
  
  double xsec = fModel->XSec(fInteraction, kPSWQ2ctpfE)*(Q2l.max-Q2l.min)*2;
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
d3XSecMK_dWQ2CosTheta_E::Clone() const
{
  return
    new d3XSecMK_dWQ2CosTheta_E(fModel,fInteraction,fW);
}
