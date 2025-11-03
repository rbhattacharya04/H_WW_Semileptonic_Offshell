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

std::vector<std::array<double,7>> _values = {};
std::vector<std::array<double,9>> _wtagger_sfs = {};

bool isHole_ex(const double cand_eta, const double cand_phi){
  bool Hole_ex = false;
  if ((cand_eta < -1.3 && cand_eta > -2.5) && (cand_phi > -1.57 && cand_phi < -0.87)) Hole_ex = true;
  return Hole_ex;

}

bool isAnalysisLepton(double Leading_Lepton_pdgId, double Leading_Lepton_pt, double Leading_Lepton_eta, double Leading_Lepton_phi){
  bool isAnaLepton = false;
  if(abs(Leading_Lepton_pdgId) == 11 && Leading_Lepton_pt > 35) isAnaLepton = true;
  if(abs(Leading_Lepton_pdgId) == 13 && Leading_Lepton_pt > 27) isAnaLepton = true;
  return isAnaLepton;
}

inline bool isVetoLepton(double nLepton, const RVec<Double_t>& Lepton_pt, const RVec<Int_t>& Lepton_isLoose){
  bool isVetoLepton = false;
  if (nLepton >1){
    for (int i=1; i<nLepton; i++){
      if (Lepton_pt[i] > 20 && Lepton_isLoose[i]) isVetoLepton = true;
    }
  }
  return isVetoLepton;
}

inline double computeMWW(const UInt_t& nLHEPart, 
                         const RVec<Double_t>& LHEPart_pt,
                         const RVec<Double_t>& LHEPart_eta, 
                         const RVec<Double_t>& LHEPart_phi,
                         const RVec<Double_t>& LHEPart_mass,
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
/*bool shouldApplyMWWCut(const std::string& sampleType) {
  return (sampleType.find("off") != std::string::npos && 
          (sampleType.find("sonly") != std::string::npos || sampleType.find("sand") != std::string::npos));
}*/

inline double getLeptonIdSF(const double& Leading_Lepton_pdgId, const bool& Leading_Lepton_isTight, const RVec<Double_t>& Lepton_tightElectron_mvaFall17V2Iso_WP90_TotSF, const RVec<Double_t>& Lepton_tightMuon_cut_Tight_HWWW_TotSF){
  double weight = 1;
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
      std::array<double,7> line_values{};
      int i = 0;
      double value;
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
      std::array<double,9> line_values{};
      int i = 0;
      double value;
      while (ss >> value)
      {
        line_values[i] = value;
        ++i;
      }
      _wtagger_sfs.push_back(line_values);
    }
  }
}


double getEleTriggerSF(double Lepton_pdgId, double Lepton_pt, double Lepton_eta){
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

double getWTaggerSF(double Jet_pt,double year = 2018, double wp = 0.5){
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

double getTriggerSF(double Lepton_pdgId, double Ele_Trigger_SF, double Mu_Trigger_SF){
  double weight  = 1;
  if (abs(Lepton_pdgId) == 11) weight = Ele_Trigger_SF;
  else if (abs(Lepton_pdgId) == 13) weight = Mu_Trigger_SF;
  return weight;
}

double deltaPhi(double phi1, double phi2)
{                                                        
  double result = phi1 - phi2;
  while (result >= double(M_PI)) result -= double(2*M_PI);
  while (result < -double(M_PI)) result += double(2*M_PI);
  return result;
}

double deltaR2(double eta1, double phi1, double eta2, double phi2)
{
  double deta = std::abs(eta1-eta2);
  double dphi = deltaPhi(phi1,phi2);
  return deta*deta + dphi*dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

inline int isGoodFatjet_indx(const RVec<Double_t>& FatJet_eta, const RVec<Double_t>& FatJet_phi,
                          const RVec<Int_t>& FatJet_jetId,
                          const RVec<Double_t>& Lepton_eta, const RVec<Double_t>& Lepton_phi)
{
  int matchedJet = -1;
  for(unsigned int iJet = 0; iJet < FatJet_eta.size(); iJet++){
    double dr = 999.;
    if (FatJet_jetId.at(iJet) <= 0) continue;
    if (abs(FatJet_eta.at(iJet)) >= 2.4) continue;
    if (isHole_ex(FatJet_eta.at(iJet),FatJet_phi.at(iJet))) continue;
    //matchedJet = iJet;

    for (unsigned int iLep = 0; iLep < Lepton_eta.size(); iLep++){
      double tmp_dr  = deltaR(FatJet_eta.at(iJet), FatJet_phi.at(iJet),
    	  Lepton_eta.at(iLep), Lepton_phi.at(iLep));
      if ( dr > tmp_dr ) dr = tmp_dr;
    }
    if (dr >= 0.8) return iJet;       
  }
  return matchedJet;
}

inline RVec<bool> getCleanJetNotOverlapping(
    double FatJet_eta,
    double FatJet_phi,
    const RVec<Double_t>& CleanJet_eta,
    const RVec<Double_t>& CleanJet_phi,
    const RVec<Double_t>& CleanJet_pt
    ) {
  RVec<bool> mask(CleanJet_eta.size(), false);
  for (size_t iJet = 0; iJet < CleanJet_eta.size(); ++iJet) {
    if (CleanJet_pt.at(iJet) < 20) continue;
    double dr = deltaR(FatJet_eta, FatJet_phi, CleanJet_eta[iJet], CleanJet_phi[iJet]);
    mask[iJet] = (dr >= 0.8);
  }
  return mask;
}

inline double getBTagSF(const RVec<Double_t>& CleanJet_btagSF_notOverlap){
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
                               const RVec<Double_t>& Jet_PUIDSF_loose,
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

double getHiggsCandidate (const Double_t& Lepton_pt, const Double_t& Lepton_eta, const Double_t& Lepton_phi,
				const Double_t& Jet_pt, const Double_t& Jet_eta, const Double_t& Jet_phi, const Double_t& Jet_mass, int var){
  ROOT::Math::PtEtaPhiMVector Lepton = ROOT::Math::PtEtaPhiMVector(Lepton_pt, Lepton_eta,
                                      Lepton_phi, 0);
  ROOT::Math::PtEtaPhiMVector Jet = ROOT::Math::PtEtaPhiMVector(Jet_pt, Jet_eta, Jet_phi, Jet_mass);
  ROOT::Math::PtEtaPhiMVector H_vis = Lepton + Jet;
  if(var == 0 ) return H_vis.M();
  if(var == 1) return H_vis.Pt();
  if(var == 2) return H_vis.Eta();
  if(var == 3) return H_vis.Phi();
}

// genjjMax -> Check carefully again if anything. is redundant
//  for addition of Samples Weight and cuts
// Computes the maximum mjj from all pairs of GenJets not overlapping with GenDressedLeptons
inline double genMjjmax(const UInt_t& nGenJet,
                       const RVec<Double_t>& GenJet_pt,
                       const RVec<Double_t>& GenJet_eta,
                       const RVec<Double_t>& GenJet_phi,
                       const RVec<Double_t>& GenJet_mass,
                       const UInt_t& nGenDressedLepton,
                       const RVec<Double_t>& GenDressedLepton_pt,
                       const RVec<Double_t>& GenDressedLepton_eta,
                       const RVec<Double_t>& GenDressedLepton_phi) {
    std::vector<int> cleanJetIdx;
    // Select GenJets with pt > 30 and |eta| < 4.7, not overlapping with GenDressedLeptons (ΔR > 0.4)
    for (UInt_t iJ = 0; iJ < nGenJet; ++iJ) {
        if (GenJet_pt[iJ] < 30. || std::abs(GenJet_eta[iJ]) > 4.7) continue;
        bool overlap = false;
        for (UInt_t iL = 0; iL < nGenDressedLepton; ++iL) {
            if (GenDressedLepton_pt[iL] < 10.) continue;
            double dr = deltaR(GenJet_eta[iJ], GenJet_phi[iJ], GenDressedLepton_eta[iL], GenDressedLepton_phi[iL]);
            if (dr < 0.4) {
                overlap = true;
                break;
            }
        }
        if (!overlap) cleanJetIdx.push_back(iJ);
    }
    double mjjmax = -999.;
    // Loop over all pairs of clean jets and compute mjj
    for (size_t i = 0; i < cleanJetIdx.size(); ++i) {
        TLorentzVector j1;
        j1.SetPtEtaPhiM(GenJet_pt[cleanJetIdx[i]], GenJet_eta[cleanJetIdx[i]], GenJet_phi[cleanJetIdx[i]], GenJet_mass[cleanJetIdx[i]]);
        for (size_t j = i+1; j < cleanJetIdx.size(); ++j) {
            TLorentzVector j2;
            j2.SetPtEtaPhiM(GenJet_pt[cleanJetIdx[j]], GenJet_eta[cleanJetIdx[j]], GenJet_phi[cleanJetIdx[j]], GenJet_mass[cleanJetIdx[j]]);
            double mjj = (j1 + j2).M();
            if (mjj > mjjmax) mjjmax = mjj;
        }
    }
    return mjjmax;
}

// Returns true if Gen_ZGstar_mass > 0 and < 4 (gstarLow)
bool gstarLow(double Gen_ZGstar_mass) {
    return (Gen_ZGstar_mass > 0. && Gen_ZGstar_mass < 4.);
}

// Returns true if Gen_ZGstar_mass < 0 or > 4 (gstarHigh)
bool gstarHigh(double Gen_ZGstar_mass) {
    return (Gen_ZGstar_mass < 0. || Gen_ZGstar_mass > 4.);
}

// Top pT reweighting function (Top PAG)
double Top_pTrw(double topGenPtOTF, double antitopGenPtOTF) {
    if (topGenPtOTF * antitopGenPtOTF > 0.) {
        double w1 = 0.103 * std::exp(-0.0118 * topGenPtOTF) - 0.000134 * topGenPtOTF + 0.973;
        double w2 = 0.103 * std::exp(-0.0118 * antitopGenPtOTF) - 0.000134 * antitopGenPtOTF + 0.973;
        return std::sqrt(w1 * w2);
    } else {
        return 1.0;
    }
}

// DY photon filter: returns true if event passes the DY photon veto
inline bool DYPhotonFilter(const UInt_t& nPhotonGen,
                           const RVec<Double_t>& PhotonGen_pt,
                           const RVec<Double_t>& PhotonGen_eta,
                           const RVec<Int_t>& PhotonGen_isPrompt,
                           const UInt_t& nLeptonGen,
                           const RVec<Double_t>& LeptonGen_pt,
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
                              const RVec<Double_t>& PhotonGen_pt,
                              const RVec<Double_t>& PhotonGen_eta,
                              const RVec<Int_t>& PhotonGen_isPrompt) {
    for (UInt_t i = 0; i < nPhotonGen; ++i) {
        if (PhotonGen_isPrompt[i] == 1 && PhotonGen_pt[i] > 10 && std::abs(PhotonGen_eta[i]) < 2.5)
            return false; // Event fails filter if any such photon exists
    }
    return true; // Event passes filter if no such photon exists
}

bool isHoleLepton(const double cand_eta, const double cand_phi, const double pdgId){
  if(abs(pdgId) == 13) return false;
  if(abs(pdgId) == 11) return isHole_ex(cand_eta,cand_phi);
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

float Top_pTrw(const ROOT::VecOps::RVec<int>& GenPart_pdgId,
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

int GenLHE(const ROOT::VecOps::RVec<int>& LHEPart_pdgId) {
    int sum = 0;
    for (auto pdgId : LHEPart_pdgId) {
        if (pdgId == 21) sum++;
    }
    return sum == 0;
}
