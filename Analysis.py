from Dataset import dataset
import pickle

import ROOT
ROOT.gInterpreter.ProcessLine(".O3")
ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "SemiLeptonic.h"')


def makeRDF(dataset_name):
    results = {}

    # Get files and isMC from dataset
    files = dataset[dataset_name]["files"]
    isMC = dataset[dataset_name]["isMC"]
    df = ROOT.RDataFrame("Events", files)

    #df = df.Range(1000)
    ROOT.RDF.Experimental.AddProgressBar(df)

    # Define mWW first 
    if isMC:
        df = df.Define("mWW", "computeMWW(nLHEPart, LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId, LHEPart_status)")
    
   # Use weight from Dataset.py (this includes mWW cuts automatically!)
    df = df.Define("weight", dataset[dataset_name]["weight"])
    
    
    #comment out the following two lines as I have defined weight as above
    # df = df.Redefine("weight","weight*Generator_weight")
    
    df = df.Define("cutflow_stage","0")
    results["Cutflow1"] = df.Histo1D(("h_cutflow_1","Cutflow 1",1,-0.5,0.5),"cutflow_stage","weight")
    
    # Using direct HLT filter
    df = df.Filter("HLT_IsoMu24 || HLT_Ele32_WPTight_Gsf","HLT Cut")
    results["Cutflow2"] = df.Histo1D(("h_cutflow_2","Cutflow 2",1,-0.5,0.5),"cutflow_stage","weight")
    #if isMC: df = df.Redefine("weight","weight*puWeight*EMTFbug_veto")
    #df = df.Define("weight","1")
    
    # Original lepton tight definitions
    ele_tight = "(abs(Lepton_pdgId) == 11 && Lepton_isTightElectron_mvaFall17V2Iso_WP90)"
    mu_tight = "(abs(Lepton_pdgId) == 13 && Lepton_isTightMuon_cut_Tight_HWWW)"
    lepton_tight = ele_tight + " || " + mu_tight
    
    # Original lepton tight definitions (simpler approach)
    # The try-except block above will catch and handle missing columns
    ele_tight = "(abs(Lepton_pdgId) == 11 && Lepton_isTightElectron_mvaFall17V2Iso_WP90)"
    mu_tight = "(abs(Lepton_pdgId) == 13 && Lepton_isTightMuon_cut_Tight_HWWW)"
    lepton_tight = ele_tight + " || " + mu_tight
    df = df.Define("Leading_Lepton_pt","Lepton_pt[0]")
    df = df.Define("Leading_Lepton_isLoose","Lepton_isLoose[0]")
    df = df.Define("Lepton_isTight",lepton_tight)
    df = df.Define("Leading_Lepton_isTight","Lepton_isTight[0]")
    df = df.Define("Leading_Lepton_eta","Lepton_eta[0]")
    df = df.Define("Leading_Lepton_phi","Lepton_phi[0]")
    df = df.Define("Leading_Lepton_pdgId","Lepton_pdgId[0]")
    df = df.Define("Leading_Lepton_electronIdx","Lepton_electronIdx[0]")
    df = df.Define("Leading_Lepton_muonIdx","Lepton_muonIdx[0]")
    df = df.Define("isAnalysisLepton","isAnalysisLepton(Leading_Lepton_pdgId,Leading_Lepton_pt,Leading_Lepton_eta,Leading_Lepton_phi)")
    df = df.Filter("isAnalysisLepton", "Analysis Lepton Selection")
    results["Cutflow3"] = df.Histo1D(("h_cutflow_3","Cutflow 3",1,-0.5,0.5),"cutflow_stage","weight")
    df = df.Define("isVetoLepton","isVetoLepton(nLepton,Lepton_pt,Lepton_isLoose)")
    df = df.Filter("!isVetoLepton", "Veto Lepton Cut")
   
    results["Cutflow4"] = df.Histo1D(("h_cutflow_4","Cutflow 4",1,-0.5,0.5),"cutflow_stage","weight")
      
    # Generate report
    report = df.Report()
    report.Print()
 
    return results

histograms = {}


#for keys in dataset:
#histograms["ggH_sonly_off"] = makeRDF(dataset["ggH_sonly_off"],True)
histograms["ggH_sonly_off"] = makeRDF("ggH_sonly_off")
#print(histograms)

#file_path = "my_histograms_ggH_sonly_off_Step_3.pkl"

#print(f"Saving dictionary of histograms to {file_path}")
#with open(file_path, "wb") as f:
    # Use pickle.dump() to save the dictionary
#    pickle.dump(histograms, f)

output_file = ROOT.TFile("output.root", "RECREATE")
histograms["ggH_sonly_off"]["Cutflow1"].Write()
histograms["ggH_sonly_off"]["Cutflow2"].Write()
histograms["ggH_sonly_off"]["Cutflow3"].Write()
histograms["ggH_sonly_off"]["Cutflow4"].Write()
output_file.Close()

