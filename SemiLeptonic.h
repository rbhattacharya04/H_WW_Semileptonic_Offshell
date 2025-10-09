#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include <string>

// Add these includes for mWW computation
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>

using namespace ROOT;
using namespace ROOT::VecOps;

// Optimized for speed - use inline and const
inline bool isAnalysisLepton(const float Leading_Lepton_pdgId, const float Leading_Lepton_pt, const float Leading_Lepton_eta, const float Leading_Lepton_phi){
  const int abs_pdgId = abs(Leading_Lepton_pdgId);
  return (abs_pdgId == 11 && Leading_Lepton_pt > 35.0f) || 
         (abs_pdgId == 13 && Leading_Lepton_pt > 27.0f);
}

inline bool isVetoLepton(const float nLepton, const RVec<Float_t>& Lepton_pt, const RVec<Int_t>& Lepton_isLoose){
  if (nLepton <= 1) return false;
  
  for (int i = 1; i < nLepton; i++){
    if (Lepton_pt[i] > 20.0f && Lepton_isLoose[i]) return true;
  }
  return false;
}

// Add mWW computation function
inline double computeMWW(const UInt_t& nLHEPart, 
                         const RVec<Float_t>& LHEPart_pt,
                         const RVec<Float_t>& LHEPart_eta, 
                         const RVec<Float_t>& LHEPart_phi,
                         const RVec<Float_t>& LHEPart_mass,
                         const RVec<Int_t>& LHEPart_pdgId,
                         const RVec<Int_t>& LHEPart_status) {
  
  ROOT::Math::PtEtaPhiMVector WW;
  
  for (UInt_t iPart = 0; iPart < nLHEPart; ++iPart) {
    if (LHEPart_status[iPart] != 1 || std::abs(LHEPart_pdgId[iPart]) != 24) continue;
    
    WW += ROOT::Math::PtEtaPhiMVector(LHEPart_pt[iPart], LHEPart_eta[iPart], 
                                      LHEPart_phi[iPart], LHEPart_mass[iPart]);
  }
  
  return WW.M();
}

// Smart mWW cut function - only applies cut for specific sample types
inline bool shouldApplyMWWCut(const std::string& sampleType) {
  return (sampleType.find("off") != std::string::npos && 
          (sampleType.find("sonly") != std::string::npos || sampleType.find("sand") != std::string::npos));
}

// Conditional mWW cut function
inline bool passMWWCutConditional(const UInt_t& nLHEPart, 
                                  const RVec<Float_t>& LHEPart_pt,
                                  const RVec<Float_t>& LHEPart_eta, 
                                  const RVec<Float_t>& LHEPart_phi,
                                  const RVec<Float_t>& LHEPart_mass,
                                  const RVec<Int_t>& LHEPart_pdgId,
                                  const RVec<Int_t>& LHEPart_status,
                                  const std::string& sampleType) {
  
  // If this sample type doesn't need mWW cut, always pass
  if (!shouldApplyMWWCut(sampleType)) return true;
  
  // Otherwise, apply the mWW > 160 cut
  ROOT::Math::PtEtaPhiMVector WW;
  
  for (UInt_t iPart = 0; iPart < nLHEPart; ++iPart) {
    if (LHEPart_status[iPart] != 1 || std::abs(LHEPart_pdgId[iPart]) != 24) continue;
    
    WW += ROOT::Math::PtEtaPhiMVector(LHEPart_pt[iPart], LHEPart_eta[iPart], 
                                      LHEPart_phi[iPart], LHEPart_mass[iPart]);
  }
  
  return WW.M() > 160.0f;
}
