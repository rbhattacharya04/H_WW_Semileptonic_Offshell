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

inline RVec<bool> getCleanJetNotOverlapping(
    float FatJet_eta,
    float FatJet_phi,
    const RVec<Float_t>& CleanJet_eta,
    const RVec<Float_t>& CleanJet_phi
) {
    RVec<bool> mask(CleanJet_eta.size(), false);
    for (size_t iJet = 0; iJet < CleanJet_eta.size(); ++iJet) {
        float dr = deltaR(FatJet_eta, FatJet_phi, CleanJet_eta[iJet], CleanJet_phi[iJet]);
        mask[iJet] = (dr >= 0.8);
    }
    return mask;
}

inline double getBTagSF(const RVec<Float_t>& CleanJet_btagSF_notOverlap){
  double sf = 1.0;
  for (int i=0; i<CleanJet_btagSF_notOverlap.size(); i++){
    sf *= CleanJet_btagSF_notOverlap.at(i);
  }
  return sf;
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

inline double getHiggsCandidate (const Float_t& Lepton_pt, const Float_t& Lepton_eta, const Float_t& Lepton_phi,
				const Float_t& Jet_pt, const Float_t& Jet_eta, const Float_t& Jet_phi, const Float_t& Jet_mass, int var){
  ROOT::Math::PtEtaPhiMVector Lepton = ROOT::Math::PtEtaPhiMVector(Lepton_pt, Lepton_eta,
                                      Lepton_phi, 0);
  ROOT::Math::PtEtaPhiMVector Jet = ROOT::Math::PtEtaPhiMVector(Jet_pt, Jet_eta, Jet_phi, Jet_mass);
  ROOT::Math::PtEtaPhiMVector H_vis = Lepton + Jet;
  if(var == 0 ) return H_vis.M();
  if(var == 1) return H_vis.Pt();
  if(var == 2) return H_vis.Eta();
  if(var == 3) return H_vis.Phi();
}
