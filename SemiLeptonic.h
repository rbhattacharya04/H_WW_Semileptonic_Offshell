#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include <string>

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"


using namespace ROOT;
using namespace ROOT::VecOps;

std::vector<std::array<float,7>> _values = {};
std::vector<std::array<float,9>> _wtagger_sfs = {};

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

float getLeptonIdSF(const float& Leading_Lepton_pdgId, const bool& Leading_Lepton_isTight, const RVec<Float_t>& Lepton_tightElectron_mvaFall17V2Iso_WP90_TotSF, const RVec<Float_t>& Lepton_tightMuon_cut_Tight_HWWW_TotSF){
  float weight = 1;
  if (Leading_Lepton_isTight){
    if (abs(Leading_Lepton_pdgId) == 11) weight =  Lepton_tightElectron_mvaFall17V2Iso_WP90_TotSF[0];
    else if (abs(Leading_Lepton_pdgId) == 13) weight = Lepton_tightMuon_cut_Tight_HWWW_TotSF[0];
  }
  return weight;
}

void initializeEleTriggerSF(){
  std::ifstream inputFile("Ele32_pt_eta_efficiency_withSys_Run2018.txt");
  std::string line;
  if (inputFile.is_open()){

    while(getline(inputFile, line)){
      std::stringstream ss(line);
      std::array<float,7> line_values{};
      int i = 0;
      float value;
      while (ss >> value)
      {
	line_values[i] = value;
	++i;
      }
      _values.push_back(line_values);    
    }
  }
}

void initializeWTaggerSF(std::string file){
  std::ifstream inputFile(file);
  std::string line;
  if (inputFile.is_open()){

    while(getline(inputFile, line)){
      std::stringstream ss(line);
      std::array<float,9> line_values{};
      int i = 0;
      float value;
      while (ss >> value)
      {
        line_values[i] = value;
        ++i;
      }
      _wtagger_sfs.push_back(line_values);
    }
  }
}


double getEleTriggerSF(float Lepton_pdgId, float Lepton_pt, float Lepton_eta){
  double weight = 1;
  if (abs(Lepton_pdgId) != 11) return weight;
  
  //handle overflow
  if (Lepton_eta < -2.5) Lepton_eta = -2.499;
  if (Lepton_eta > 2.5) Lepton_eta = 2.499;
  if (Lepton_pt > 100) Lepton_pt = 99.99;

  for (uint j = 0; j< _values.size(); j++){
    if (Lepton_eta >= _values[j][0] && Lepton_eta <= _values[j][1] && Lepton_pt >= _values[j][2] && Lepton_pt <= _values[j][3]){
      weight = _values[j][4];
      //output[1] = _values[j][4] + _values[j][5];
      //output[2] = _values[j][4] - _values[j][6] ;
      break;
    }
  }

  return weight;
}

double getWTaggerSF(float Jet_pt,float year = 2018, float wp = 0.5){
  double weight = 1;

  //handle overflow
  if (Jet_pt > 800) Jet_pt = 799.9;

  for (uint j = 0; j< _wtagger_sfs.size(); j++){
    if (_wtagger_sfs[j][1] != year) continue;
    if (_wtagger_sfs[j][3] != wp) continue;
    if (Jet_pt >= _wtagger_sfs[j][4] && Jet_pt <= _wtagger_sfs[j][5]){
      weight = _wtagger_sfs[j][6];
      //output[1] = _values[j][4] + _values[j][5];
      //output[2] = _values[j][4] - _values[j][6] ;
      break;
    }
  }

  return weight;
}

double getTriggerSF(float Lepton_pdgId, float Ele_Trigger_SF, float Mu_Trigger_SF){
  double weight  = 1;
  if (abs(Lepton_pdgId) == 11) weight = Ele_Trigger_SF;
  else if (abs(Lepton_pdgId) == 13) weight = Mu_Trigger_SF;
  return weight;
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

int isGoodFatjet_indx(const RVec<Float_t>& FatJet_eta, const RVec<Float_t>& FatJet_phi,
                          const RVec<Int_t>& FatJet_jetId,
                          const RVec<Float_t>& Lepton_eta, const RVec<Float_t>& Lepton_phi)
{
  for(unsigned int iJet = 0; iJet < FatJet_eta.size(); iJet++){
    float dr = 999.;
    if (FatJet_jetId.at(iJet) < 0) continue;
    if (abs(FatJet_eta.at(iJet)) > 2.4) continue;

    for (unsigned int iLep = 0; iLep < Lepton_eta.size(); iLep++){
      float tmp_dr  = deltaR(FatJet_eta.at(iJet), FatJet_phi.at(iJet),
	  Lepton_eta.at(iLep), Lepton_phi.at(iLep));
      if (tmp_dr < dr) dr = tmp_dr;
    }
    if (dr > 0.8) return iJet;         
  }
  return -1;
}




inline double computePUJetIdSF(const UInt_t& nJet,
                               const RVec<Int_t>& Jet_jetId,
                               const RVec<Int_t>& Jet_electronIdx1,
                               const RVec<Int_t>& Jet_muonIdx1,
                               const RVec<Float_t>& Jet_PUIDSF_loose,
                               const Int_t& Leading_Lepton_electronIdx,
                               const Int_t& Leading_Lepton_muonIdx) {

    double logSum = 0.0;
    for (UInt_t iJet = 0; iJet < nJet; ++iJet) {
        if (Jet_jetId[iJet] < 2) continue;
        if (Jet_electronIdx1[iJet] >= 0 && Jet_electronIdx1[iJet] == Leading_Lepton_electronIdx) continue;
        if (Jet_muonIdx1[iJet] >= 0 && Jet_muonIdx1[iJet] == Leading_Lepton_muonIdx) continue;
        if (Jet_PUIDSF_loose[iJet] > 0) logSum += TMath::Log(Jet_PUIDSF_loose[iJet]);
    }
    return TMath::Exp(logSum);
}




// Returns indices of CleanJets not overlapping with the selected FatJet (ΔR < 0.8)
inline RVec<int> getCleanJetNotOverlapping(
    float FatJet_eta,
    float FatJet_phi,
    const RVec<Float_t>& CleanJet_eta,
    const RVec<Float_t>& CleanJet_phi
) {
    RVec<int> nonOverlappingJets;
    for (size_t iJet = 0; iJet < CleanJet_eta.size(); ++iJet) {
        float dr = deltaR(FatJet_eta, FatJet_phi, CleanJet_eta[iJet], CleanJet_phi[iJet]);
        if (dr >= 0.8) nonOverlappingJets.push_back(iJet);
    }
    return nonOverlappingJets;
}


// bVeto_boo [DeepFlavB > 0.2783 medium WP]
inline bool bVeto_boo(
    const RVec<Float_t>& CleanJet_pt,
    const RVec<Float_t>& CleanJet_eta,
    const RVec<Int_t>& CleanJet_jetIdx,
    const RVec<Float_t>& Jet_btagDeepFlavB,
    const RVec<int>& CleanJet_notOverlapping,
    float bWP = 0.2783
) {
    for (auto idx : CleanJet_notOverlapping) {
        if (CleanJet_pt[idx] > 20 && std::abs(CleanJet_eta[idx]) < 2.5) {
            int jetIdx = CleanJet_jetIdx[idx];
            if (jetIdx >= 0 && jetIdx < Jet_btagDeepFlavB.size()) {
                if (Jet_btagDeepFlavB[jetIdx] > bWP) return false;
            }
        }
    }
    return true; // No b-tagged jets found
}

// Returns the b-veto SF for bVeto_boo
inline double bVeto_booSF(
    const RVec<Float_t>& CleanJet_pt,
    const RVec<Float_t>& CleanJet_eta,
    const RVec<Int_t>& CleanJet_jetIdx,
    const RVec<Float_t>& Jet_btagSF_deepjet_shape,
    const RVec<int>& CleanJet_notOverlapping,
    float ptCut = 20.0
) {
    double logSum = 0.0;
    for (auto idx : CleanJet_notOverlapping) {
        bool passJet = (CleanJet_pt[idx] > ptCut) && (std::abs(CleanJet_eta[idx]) < 2.5);
        int jetIdx = CleanJet_jetIdx[idx];
        if (jetIdx >= 0 && jetIdx < Jet_btagSF_deepjet_shape.size()) {
            if (passJet) {
                logSum += std::log(Jet_btagSF_deepjet_shape[jetIdx]);
            }
        }
    }
    return std::exp(logSum);
}


// for addition of bTagSF
inline bool bReq_boo(
    const RVec<Float_t>& CleanJet_pt,
    const RVec<Float_t>& CleanJet_eta,
    const RVec<Int_t>& CleanJet_jetIdx,
    const RVec<Float_t>& Jet_btagDeepFlavB,
    const RVec<int>& CleanJet_notOverlapping,
    float bWP = 0.2783
) {
    for (auto idx : CleanJet_notOverlapping) {
        if (CleanJet_pt[idx] > 30 && std::abs(CleanJet_eta[idx]) < 2.5) {
            int jetIdx = CleanJet_jetIdx[idx];
            if (jetIdx >= 0 && jetIdx < Jet_btagDeepFlavB.size()) {
                if (Jet_btagDeepFlavB[jetIdx] > bWP) return true;
            }
        }
    }
    return false;
}

inline double bReq_booSF(
    const RVec<Float_t>& CleanJet_pt,
    const RVec<Float_t>& CleanJet_eta,
    const RVec<Int_t>& CleanJet_jetIdx,
    const RVec<Float_t>& Jet_btagSF_deepjet_shape,
    const RVec<int>& CleanJet_notOverlapping,
    float ptCut = 30.0
) {
    double logSum = 0.0;
    for (auto idx : CleanJet_notOverlapping) {
        bool passJet = (CleanJet_pt[idx] > ptCut) && (std::abs(CleanJet_eta[idx]) < 2.5);
        int jetIdx = CleanJet_jetIdx[idx];
        if (jetIdx >= 0 && jetIdx < Jet_btagSF_deepjet_shape.size()) {
            if (passJet) {
                logSum += std::log(Jet_btagSF_deepjet_shape[jetIdx]);
            }
        }
    }
    return std::exp(logSum);
}

inline bool boosted_nocut_res(
    float PuppiMET_pt,
    int GoodFatJet_idx,
    const RVec<Float_t>& FatJet_pt,
    const RVec<Float_t>& FatJet_deepTag_WvsQCD,
    const RVec<Float_t>& FatJet_eta
) {
    if (
        PuppiMET_pt > 40 &&
        GoodFatJet_idx >= 0 &&
        FatJet_pt[GoodFatJet_idx] > 200 &&
        FatJet_deepTag_WvsQCD[GoodFatJet_idx] > 0.961 &&
        std::abs(FatJet_eta[GoodFatJet_idx]) < 2.4
    ) return true;
    return false;
}

// genjjMax -> Check carefully again if anything. is redundant
//  for addition of Samples Weight and cuts
// Computes the maximum mjj from all pairs of GenJets not overlapping with GenDressedLeptons
inline float genMjjmax(const UInt_t& nGenJet,
                       const RVec<Float_t>& GenJet_pt,
                       const RVec<Float_t>& GenJet_eta,
                       const RVec<Float_t>& GenJet_phi,
                       const RVec<Float_t>& GenJet_mass,
                       const UInt_t& nGenDressedLepton,
                       const RVec<Float_t>& GenDressedLepton_pt,
                       const RVec<Float_t>& GenDressedLepton_eta,
                       const RVec<Float_t>& GenDressedLepton_phi) {
    std::vector<int> cleanJetIdx;
    // Select GenJets with pt > 30 and |eta| < 4.7, not overlapping with GenDressedLeptons (ΔR > 0.4)
    for (UInt_t iJ = 0; iJ < nGenJet; ++iJ) {
        if (GenJet_pt[iJ] < 30. || std::abs(GenJet_eta[iJ]) > 4.7) continue;
        bool overlap = false;
        for (UInt_t iL = 0; iL < nGenDressedLepton; ++iL) {
            if (GenDressedLepton_pt[iL] < 10.) continue;
            float dr = deltaR(GenJet_eta[iJ], GenJet_phi[iJ], GenDressedLepton_eta[iL], GenDressedLepton_phi[iL]);
            if (dr < 0.4) {
                overlap = true;
                break;
            }
        }
        if (!overlap) cleanJetIdx.push_back(iJ);
    }
    float mjjmax = -999.;
    // Loop over all pairs of clean jets and compute mjj
    for (size_t i = 0; i < cleanJetIdx.size(); ++i) {
        TLorentzVector j1;
        j1.SetPtEtaPhiM(GenJet_pt[cleanJetIdx[i]], GenJet_eta[cleanJetIdx[i]], GenJet_phi[cleanJetIdx[i]], GenJet_mass[cleanJetIdx[i]]);
        for (size_t j = i+1; j < cleanJetIdx.size(); ++j) {
            TLorentzVector j2;
            j2.SetPtEtaPhiM(GenJet_pt[cleanJetIdx[j]], GenJet_eta[cleanJetIdx[j]], GenJet_phi[cleanJetIdx[j]], GenJet_mass[cleanJetIdx[j]]);
            float mjj = (j1 + j2).M();
            if (mjj > mjjmax) mjjmax = mjj;
        }
    }
    return mjjmax;
}

// Returns true if Gen_ZGstar_mass > 0 and < 4 (gstarLow)
inline bool gstarLow(float Gen_ZGstar_mass) {
    return (Gen_ZGstar_mass > 0. && Gen_ZGstar_mass < 4.);
}

// Returns true if Gen_ZGstar_mass < 0 or > 4 (gstarHigh)
inline bool gstarHigh(float Gen_ZGstar_mass) {
    return (Gen_ZGstar_mass < 0. || Gen_ZGstar_mass > 4.);
}

// DY photon filter: returns true if event passes the DY photon veto
inline bool DYPhotonFilter(const UInt_t& nPhotonGen,
                           const RVec<Float_t>& PhotonGen_pt,
                           const RVec<Float_t>& PhotonGen_eta,
                           const RVec<Int_t>& PhotonGen_isPrompt,
                           const UInt_t& nLeptonGen,
                           const RVec<Float_t>& LeptonGen_pt,
                           const RVec<Int_t>& LeptonGen_isPrompt) {
    int nPromptPhoton = 0;
    int nPromptLepton = 0;
    for (UInt_t i = 0; i < nPhotonGen; ++i) {
        if (PhotonGen_isPrompt[i] == 1 && PhotonGen_pt[i] > 15 && std::abs(PhotonGen_eta[i]) < 2.6)
            nPromptPhoton++;
    }
    for (UInt_t i = 0; i < nLeptonGen; ++i) {
        if (LeptonGen_isPrompt[i] == 1 && LeptonGen_pt[i] > 15)
            nPromptLepton++;
    }
    // Passes filter if NOT (at least one prompt photon and at least two prompt leptons)
    return !(nPromptPhoton > 0 && nPromptLepton >= 2);
}


// Wjets photon filter: returns true if there are NO prompt photons with pt > 10 and |eta| < 2.5
inline bool WjetsPhotonFilter(const UInt_t& nPhotonGen,
                              const RVec<Float_t>& PhotonGen_pt,
                              const RVec<Float_t>& PhotonGen_eta,
                              const RVec<Int_t>& PhotonGen_isPrompt) {
    for (UInt_t i = 0; i < nPhotonGen; ++i) {
        if (PhotonGen_isPrompt[i] == 1 && PhotonGen_pt[i] > 10 && std::abs(PhotonGen_eta[i]) < 2.5)
            return false; // Event fails filter if any such photon exists
    }
    return true; // Event passes filter if no such photon exists
}

inline float Top_pTrw(const ROOT::VecOps::RVec<int>& GenPart_pdgId,
                      const ROOT::VecOps::RVec<unsigned int>& GenPart_statusFlags,
                      const ROOT::VecOps::RVec<float>& GenPart_pt) {
    float topGenPtOTF = 0.0;
    float antitopGenPtOTF = 0.0;
    for (size_t i = 0; i < GenPart_pdgId.size(); ++i) {
        if (GenPart_pdgId[i] == 6 && ((GenPart_statusFlags[i] >> 13) & 1)) {
            topGenPtOTF += GenPart_pt[i];
        }
        if (GenPart_pdgId[i] == -6 && ((GenPart_statusFlags[i] >> 13) & 1)) {
            antitopGenPtOTF += GenPart_pt[i];
        }
    }
    if (topGenPtOTF * antitopGenPtOTF > 0.) {
        float w1 = 0.103 * std::exp(-0.0118 * topGenPtOTF) - 0.000134 * topGenPtOTF + 0.973;
        float w2 = 0.103 * std::exp(-0.0118 * antitopGenPtOTF) - 0.000134 * antitopGenPtOTF + 0.973;
        return std::sqrt(w1 * w2);
    } else {
        return 1.0;
    }
}
