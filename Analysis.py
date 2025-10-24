from Dataset import dataset
import pickle
import argparse

import ROOT
ROOT.gInterpreter.ProcessLine(".O3")
ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "SemiLeptonic.h"')

def makeRDF(dataset_name, wtagger="Nominal"):
    results = {}
    # Get files and isMC from dataset
    #sample_dict = dataset[dataset_name]
    files = dataset[dataset_name]["files"]
    isMC = dataset[dataset_name]["isMC"]
    isSignal = dataset[dataset_name]["isSignal"]
    isOffshell = dataset[dataset_name]["isOffshell"]
    
    df = ROOT.RDataFrame("Events", files)
    #df = df.Range(1000)
    ROOT.RDF.Experimental.AddProgressBar(df)
    if isSignal:
        #XS_Weight_cal = (1000*0.4357)/genEventSumw
        XS_Weight_cal = "XSWeight"
    df = df.Define("weight","1")
    # Define only needed variables for weights/filters
    

 
    if isMC:
        df = df.Define("passDYPhotonFilter", "DYPhotonFilter(nPhotonGen, PhotonGen_pt, PhotonGen_eta, PhotonGen_isPrompt, nLeptonGen, LeptonGen_pt, LeptonGen_isPrompt)" )
        df = df.Define("passWjetsPhotonFilter", "WjetsPhotonFilter(nPhotonGen, PhotonGen_pt, PhotonGen_eta, PhotonGen_isPrompt)")
        df = df.Define("mjjGenmax","genMjjmax(nGenJet, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, nGenDressedLepton, GenDressedLepton_pt, GenDressedLepton_eta, GenDressedLepton_phi)")
        df = df.Define("Top_pTrw", "Top_pTrw(GenPart_pdgId, GenPart_statusFlags, GenPart_pt)")
        df = df.Define("gstarLowWeight", "0.94 * float(gstarLow(Gen_ZGstar_mass))")
        df = df.Define("gstarHighWeight", "1.14 * float(gstarHigh(Gen_ZGstar_mass))")
        df = df.Define("VgWeight", "gstarLowWeight + gstarHighWeight")
        df = df.Define("GenLHE", "GenLHE(LHEPart_pdgId)")
        if isSignal:
            df = df.Define("Lhe_mWW", "computeMWW(nLHEPart, LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId, LHEPart_status)")
            if isOffshell:
                df = df.Filter("Lhe_mWW > 160")
            else:
                df = df.Filter("Lhe_mWW < 160")
        if isSignal:
           #df = df.Redefine("weight",f"weight*genWeight*{XS_Weight_cal}*METFilter_MC*puWeight*EMTFbug_veto") #XSWeight is genweight*baseW https://github.com/sv3048/LatinoAnalysis/blob/SemilepOFFSHEL    L/NanoGardener/python/data/formulasToAdd_MCnoSF_Full2018v9.py#L29-L31
           df = df.Redefine("weight",f"weight*genWeight*{XS_Weight_cal}*METFilter_MC*puWeight*EMTFbug_veto") #XSWeight is genweight*baseW https://github.com/sv3048/LatinoAnalysis/blob/SemilepOFFSHELL    /NanoGardener/python/data/formulasToAdd_MCnoSF_Full2018v9.py#L29-L31
        else:
            df = df.Redefine("weight","weight*XSWeight*METFilter_MC*puWeight*EMTFbug_veto")

        # # MC weights and corrections
        # df = df.Redefine("weight","weight*XSWeight*METFilter_MC*puWeight*EMTFbug_veto")
    
    # Apply sample-specific weight/filter
    if dataset[dataset_name]["sample_weights"]:
        df = df.Redefine("weight", f"weight*({dataset[dataset_name]['sample_weights']})")
    if dataset[dataset_name]["sample_filters"]:
        df = df.Filter(dataset[dataset_name]["sample_filters"])

    # Cutflow 1: Initialization
    df = df.Define("cutflow_stage","0")
    results["Cutflow1"] = df.Histo1D(("h_cutflow_1","Cutflow 1",1,-0.5,0.5),"cutflow_stage","weight")

    # Using direct HLT filter
    df = df.Filter("HLT_IsoMu24 || HLT_Ele32_WPTight_Gsf","HLT Cut")
    results["Cutflow2"] = df.Histo1D(("h_cutflow_2","Cutflow 2",1,-0.5,0.5),"cutflow_stage","weight")

     
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
    
    if isMC:
        df = df.Define("Leading_Lepton_promptgenmatched","Lepton_promptgenmatched[0]")
        df = df.Filter("Leading_Lepton_promptgenmatched", "Gen Matching of the leading Lepton")    
        df = df.Define("Lepton_ID_SF","getLeptonIdSF(Leading_Lepton_pdgId,Leading_Lepton_isTight,Lepton_tightElectron_mvaFall17V2Iso_WP90_TotSF,Lepton_tightMuon_cut_Tight_HWWW_TotSF)")
        df = df.Redefine("weight","weight*Lepton_ID_SF")
        
        if hasattr(ROOT, "initializeEleTriggerSF"):
            ROOT.initializeEleTriggerSF()
        
        df = df.Define("EleTriggerSF","getEleTriggerSF(Leading_Lepton_pdgId,Leading_Lepton_pt,Leading_Lepton_eta)")
        df = df.Define("LepTriggerSF","getTriggerSF(Leading_Lepton_pdgId,EleTriggerSF,TriggerEffWeight_1l)")
        df = df.Redefine("weight", "weight*LepTriggerSF")
       

    df = df.Define("isAnalysisLepton","isAnalysisLepton(Leading_Lepton_pdgId,Leading_Lepton_pt,Leading_Lepton_eta,Leading_Lepton_phi)")
    df = df.Filter("isAnalysisLepton", "Analysis Lepton Selection")


    results["Cutflow3"] = df.Histo1D(("h_cutflow_3","Cutflow 3",1,-0.5,0.5),"cutflow_stage","weight")
    
    df = df.Define("isVetoLepton","isVetoLepton(nLepton,Lepton_pt,Lepton_isLoose)")
    df = df.Filter("!isVetoLepton", "Veto Lepton Cut")

    df = df.Filter("PuppiMET_pt > 30", "PuppiMET_pt > 30 GeV cut")
   
    results["Cutflow4"] = df.Histo1D(("h_cutflow_4","Cutflow 4",1,-0.5,0.5),"cutflow_stage","weight")

    # Jet Selection
    df = df.Filter("nFatJet>=1","At Least 1 Fat Jet")
    if isMC:
        df = df.Define("JetPUIDSF","computePUJetIdSF(nJet,Jet_jetId,Jet_electronIdx1,Jet_muonIdx1,Jet_PUIDSF_loose,Leading_Lepton_electronIdx,Leading_Lepton_muonIdx)")
        df = df.Redefine("weight","weight*JetPUIDSF")
    df = df.Define("GoodFatJet_idx","isGoodFatjet_indx(FatJet_eta,FatJet_phi,FatJet_jetId,Lepton_eta,Lepton_phi)")
    df = df.Filter("GoodFatJet_idx != -1", "Good Fat Jet cut")
    df = df.Define("AnaFatJet_pt","FatJet_pt[GoodFatJet_idx]")
    df = df.Define("AnaFatJet_eta","FatJet_eta[GoodFatJet_idx]")
    df = df.Define("AnaFatJet_phi", "FatJet_phi[GoodFatJet_idx]")
    df = df.Define("AnaFatJet_jetId","FatJet_jetId[GoodFatJet_idx]")
    df = df.Define("CleanJet_notOverlapping", "getCleanJetNotOverlapping(FatJet_eta[GoodFatJet_idx], FatJet_phi[GoodFatJet_idx], CleanJet_eta, CleanJet_phi)")
    df = df.Define("bVeto_boo", "bVeto_boo(CleanJet_pt, CleanJet_eta, CleanJet_jetIdx, Jet_btagDeepFlavB, CleanJet_notOverlapping)")
    df = df.Filter("bVeto_boo", "bjet veto")
    df = df.Define("bReq_boo", "bReq_boo(CleanJet_pt, CleanJet_eta, CleanJet_jetIdx, Jet_btagDeepFlavB, CleanJet_notOverlapping)")
    df = df.Define("bReq_booSF", "bReq_booSF(CleanJet_pt, CleanJet_eta, CleanJet_jetIdx, Jet_btagSF_deepjet_shape, CleanJet_notOverlapping, 30.0)")
    df = df.Define("boosted_nocut_res", "boosted_nocut_res(PuppiMET_pt, GoodFatJet_idx, FatJet_pt, FatJet_deepTag_WvsQCD, FatJet_eta)"
)
    #df = df.Define("boosted_nocut_res", "PuppiMET_pt > 40 && GoodFatJet_idx >= 0 && FatJet_pt[GoodFatJet_idx] > 200 && FatJet_deepTag_WvsQCD[GoodFatJet_idx] > 0.961 && abs(FatJet_eta[GoodFatJet_idx]) < 2.4")

    if isMC:
        df = df.Define("bVeto_booSF","bVeto_booSF(CleanJet_pt, CleanJet_eta, CleanJet_jetIdx, Jet_btagSF_deepjet_shape, CleanJet_notOverlapping, 20.0)")
        df = df.Define("btagSF", "boosted_nocut_res * bVeto_booSF * bVeto_boo + bReq_booSF * bReq_boo")
        df = df.Redefine("weight","weight*bVeto_booSF")
        df = df.Redefine("weight","weight*btagSF")
   
    if wtagger == "Nominal":
        df = df.Define("AnaFatJet_nom_wtag","FatJet_particleNet_WvsQCD[GoodFatJet_idx]")
        if isMC:
            if hasattr(ROOT, "initializeWTaggerSF"):
                ROOT.initializeWTaggerSF("W_Nominal_Run2_SF.csv")
    elif wtagger == "MD":
        df = df.Define("AnaFatJet_md_wtag","(FatJet_particleNetMD_Xcc[GoodFatJet_idx] + FatJet_particleNetMD_Xqq[GoodFatJet_idx])/(FatJet_particleNetMD_Xcc[GoodFatJet_idx] + FatJet_particleNetMD_Xqq[GoodFatJet_idx] + FatJet_particleNetMD_QCD[GoodFatJet_idx])")
        if isMC:
            if hasattr(ROOT, "initializeWTaggerSF"):
                ROOT.initializeWTaggerSF("W_MD_Run2_SF.csv")
    #df = df.Filter("AnaFatJet_jetId == 2", "Tight Jet Id cut")
    df = df.Filter("AnaFatJet_pt>200","Jet pT cut")
    if wtagger == "Nominal":
        df = df.Filter("AnaFatJet_nom_wtag > 0.94","Wtagger Nominal 0p5 cut")
        if isMC:
            df = df.Define("WTagger_SF","getWTaggerSF(AnaFatJet_pt)")
            df = df.Redefine("weight","weight*WTagger_SF")
    elif wtagger == "MD":
        df = df.Filter("AnaFatJet_md_wtag > 0.90"," Wtagger MD 0p5 cut")
        if isMC:
            df = df.Define("WTagger_SF","getWTaggerSF(AnaFatJet_pt)")
            df = df.Redefine("weight","weight*WTagger_SF")
    

    results["Cutflow5"] = df.Histo1D(("h_cutflow_5","Cutflow 5",1,-0.5,0.5),"cutflow_stage","weight")
    
    report = df.Report()
    report.Print()
    

 
    return results

parser = argparse.ArgumentParser()
parser.add_argument("-w","--wtag", help="WTagger option Nominal, MD", type=str,default = "Nominal")
args = parser.parse_args()


histograms = {}
#for keys in dataset:
#histograms["ggH_sonly_off"] = makeRDF(dataset["ggH_sonly_off"],True)
#histograms["DY_else"] = makeRDF("DY_else",args.wtag)
#histograms["DY_else"] = makeRDF("DY_else",args.wtag)
histograms["WW"] = makeRDF("WW",args.wtag)


#

# output_file = ROOT.TFile("output.root", "RECREATE")
# histograms["ggH_sonly_off"]["Cutflow1"].Write()
# histograms["ggH_sonly_off"]["Cutflow2"].Write()
# histograms["ggH_sonly_off"]["Cutflow3"].Write()
# histograms["ggH_sonly_off"]["Cutflow4"].Write()
# histograms["ggH_sonly_off"]["Cutflow5"].Write()
# output_file.Close()

# Example snippet after makeRDF
# ...existing code...

# results = makeRDF("ggH_sonly_off", args.wtag)  # Call makeRDF and get results

# output_file = ROOT.TFile("output.root", "RECREATE")
# for h in results.values():
#     h.Write()
# output_file.Close()

# Example main loop

output_file = ROOT.TFile("output.root", "RECREATE")
histograms["WW"]["Cutflow1"].Write()
histograms["WW"]["Cutflow2"].Write()
histograms["WW"]["Cutflow3"].Write()
histograms["WW"]["Cutflow4"].Write()
histograms["WW"]["Cutflow5"].Write()
output_file.Close()



