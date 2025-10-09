#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>

using namespace ROOT;
using namespace ROOT::VecOps;

bool isAnalysisLepton(float Leading_Lepton_pdgId, float Leading_Lepton_pt, float Leading_Lepton_eta, float Leading_Lepton_phi){
  bool isAnaLepton = false;
  if(abs(Leading_Lepton_pdgId) == 11 && Leading_Lepton_pt > 35) isAnaLepton = true;
  if(abs(Leading_Lepton_pdgId) == 13 && Leading_Lepton_pt > 27) isAnaLepton = true;
  return isAnaLepton;
}

bool isVetoLepton(float nLepton, const RVec<Float_t>& Lepton_pt, const RVec<Int_t>& Lepton_isLoose){
  bool isVetoLepton = false;
  if (nLepton >1){
    for (int i=1; i<nLepton; i++){
      if (Lepton_pt[i] > 20 && Lepton_isLoose[i]) isVetoLepton = true;
    }
  }
  return isVetoLepton;
}

float deltaPhi(float phi1, float phi2)
{                                                        
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2)
{
  float deta = std::abs(eta1-eta2);
  float dphi = deltaPhi(phi1,phi2);
  return deta*deta + dphi*dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}


RVec<bool> isGoodFatjet(const RVec<Float_t>& FatJet_eta, const RVec<Float_t>& FatJet_phi,
                          const RVec<Float_t>& Lepton_eta, const RVec<Float_t>& Lepton_phi)
{
  RVec<bool> isGoodFatJet(FatJet_eta.size(), false);
  for(unsigned int iJet = 0; iJet < FatJet_eta.size(); iJet++){
    float dr = 999.;

    for (unsigned int iLep = 0; iLep < Lepton_eta.size(); iLep++){
      float tmp_dr  = deltaR(FatJet_eta.at(iJet), FatJet_phi.at(iJet),
	  Lepton_eta.at(iLep), Lepton_phi.at(iLep));
      if (tmp_dr < dr) dr = tmp_dr;
    }
    if (dr > 0.8) isGoodFatJet[iJet] = true;        
  }
  return isGoodFatJet;
}
