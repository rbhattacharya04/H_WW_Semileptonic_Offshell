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

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"


using namespace ROOT;
using namespace ROOT::VecOps;

std::vector<std::array<float,7>> _values = {};

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

double getTriggerSF(float Lepton_pdgId, float Ele_Trigger_SF, float Mu_Trigger_SF){
  double weight  = 1;
  if (abs(Lepton_pdgId) == 11) weight = Ele_Trigger_SF;
  else if (abs(Lepton_pdgId) == 13) weight = Mu_Trigger_SF;
  return weight;
}
