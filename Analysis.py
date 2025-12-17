from Dataset import dataset
import pickle
import argparse

import ROOT
ROOT.gInterpreter.ProcessLine(".O3")
ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "SemiLeptonic.h"')

def makeRDF(dataset_name, wtagger="Nominal"):
    print(dataset_name)
    results = {}
    # Get files and isMC from dataset
    #sample_dict = dataset[dataset_name]
    files = dataset[dataset_name]["files"]
    isMC = dataset[dataset_name]["isMC"]
    isSignal = dataset[dataset_name]["isSignal"]
    isOffshell = dataset[dataset_name]["isOffshell"]
    #files = ["/afs/cern.ch/work/r/rbhattac/public/DY_v9_NanoAOD/13D0AD97-6B32-CB4C-BA87-5E37BA4CF20E.root"]
    if isMC:
        df_run = ROOT.RDataFrame("Runs", files)
        ROOT.RDF.Experimental.AddProgressBar(df_run)
        sum_result = df_run.Sum("genEventSumw")
        genEventSumw = sum_result.GetValue()
        #genEventSumw = 109864.115
        print(f"genEventSumw = {genEventSumw}")
        h1 = ROOT.TH1F("genEventWeight", "Example Histogram;X-axis Label;Y-axis Label", 5, -0.5, 4.5)
        h1.SetBinContent(1, genEventSumw)
 
    
    df = ROOT.RDataFrame("Events", files)
    #df = df.Range(1000)
    ROOT.RDF.Experimental.AddProgressBar(df)
    #XS_Weight_cal = (59.7*1000*0.775)/genEventSumw
    if "xsec" in dataset[dataset_name]: 
        XS_Weight_cal = (1000*dataset[dataset_name]["xsec"])/genEventSumw
    
    df = df.Define("weight","1")
    
    #count = df.Count()
    #print(f"No. of events : {count.GetValue()}")
    if isMC:
        sum_genWeight = df.Sum("genWeight")

        if not isSignal:
            df = df.Define("DYPhotonFilter", "DYPhotonFilter(nPhotonGen, PhotonGen_pt, PhotonGen_eta, PhotonGen_isPrompt, nLeptonGen, LeptonGen_pt, LeptonGen_isPrompt)")
            df = df.Define("WjetsPhotonFilter", "WjetsPhotonFilter(nPhotonGen, PhotonGen_pt, PhotonGen_eta, PhotonGen_isPrompt)")
            df = df.Define("mjjGenmax","genMjjmax(nGenJet, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, nGenDressedLepton, GenDressedLepton_pt, GenDressedLepton_eta, GenDressedLepton_phi)")
            df = df.Define("Top_pTrw", "Top_pTrw(GenPart_pdgId, GenPart_statusFlags, GenPart_pt)")
            df = df.Define("gstarLowWeight", "0.94 * float(gstarLow(Gen_ZGstar_mass))")
            df = df.Define("gstarHighWeight", "1.14 * float(gstarHigh(Gen_ZGstar_mass))")
            df = df.Define("VgWeight", "gstarLowWeight + gstarHighWeight")
            if dataset_name not in ["ZZ", "WZ"]:
                df = df.Define("GenLHE", "GenLHE(LHEPart_pdgId)")
            if dataset_name in ["ggH_bonly_on", "ggH_bonly_off", "ggH_sand_off", "ggH_sand_on"]:
                df = df.Define("Lhe_mWW", "computeMWW(nLHEPart, LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId, LHEPart_status)")
                if dataset_name in ["ggH_bonly_off","ggH_sand_off"]:
                    df = df.Filter("Lhe_mWW > 160")
                else:
                    df = df.Filter("Lhe_mWW < 160")
        
        if isSignal:
            df = df.Define("Lhe_mWW", "computeMWW(nLHEPart, LHEPart_pt, LHEPart_eta, LHEPart_phi, LHEPart_mass, LHEPart_pdgId, LHEPart_status)")
            if isOffshell:
                df = df.Filter("Lhe_mWW > 160")
            else:
                df = df.Filter("Lhe_mWW < 160")
        #comment out the following two lines as I have defined weight as above
    
        if "xsec" in dataset[dataset_name]:
            #df = df.Redefine("weight",f"weight*genWeight*{XS_Weight_cal}*METFilter_MC*puWeight*EMTFbug_veto") #XSWeight is genweight*baseW https://github.com/sv3048/LatinoAnalysis/blob/SemilepOFFSHELL/NanoGardener/python/data/formulasToAdd_MCnoSF_Full2018v9.py#L29-L31
            print("xsec = ",dataset[dataset_name]["xsec"])
            print("XS_Weight_cal = ", XS_Weight_cal)
            df = df.Redefine("weight",f"weight*genWeight*{XS_Weight_cal}*METFilter_MC*puWeight*EMTFbug_veto") #XSWeight is genweight*baseW https://github.com/sv3048/LatinoAnalysis/blob/SemilepOFFSHELL/NanoGardener/python/data/formulasToAdd_MCnoSF_Full2018v9.py#L29-L31
        else:
            df = df.Redefine("weight","weight*XSWeight*METFilter_MC*puWeight*EMTFbug_veto")

   
    ## CHECK ?? 
    # Apply sample-specific weight/filter
    #if isMC:
    #    if dataset[dataset_name]["sample_weights"]:
    #        df = df.Redefine("weight", f"weight*({dataset[dataset_name]['sample_weights']})")
    #    if dataset[dataset_name]["sample_filters"]:
    #        df = df.Filter(dataset[dataset_name]["sample_filters"])

 
    if dataset_name in ["ggH_sonly_off","ggH_bonly_off","ggH_sand_off"]:
        results["Lhe_mWW"] = df.Histo1D(("h_Lhe_mWW","Lhe_mWW",50,0,1000),"Lhe_mWW","weight")
        df = df.Define("Gen_mWW","computeMWWGen(nGenPart, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)")
        results["Gen_mWW"] = df.Histo1D(("h_Gen_mWW","Gen_mWW",50,0, 1000),"Gen_mWW", "weight") 
    df = df.Define("cutflow_stage","0")
    results["Cutflow_Start"] = df.Histo1D(("h_cutflow_Start","Cutflow Start",1,-0.5,0.5),"cutflow_stage","weight")

    # Using direct HLT filter
    df = df.Filter("HLT_IsoMu24 || HLT_Ele32_WPTight_Gsf","HLT Cut")
    results["Cutflow_Trigger"] = df.Histo1D(("h_cutflow_Trigger","Cutflow Trigger",1,-0.5,0.5),"cutflow_stage","weight")

     
    results["All_Lepton_pt"] = df.Histo1D(("h_All_Lepton_pt","All Lepton pt",50,0,1000),"Lepton_pt","weight")
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
    results["Lepton_pt"] = df.Histo1D(("h_Lepton_pt","Lepton pt",50,0,1000),"Leading_Lepton_pt","weight")
    
    if isMC:
        df = df.Define("Leading_Lepton_promptgenmatched","Lepton_promptgenmatched[0]")
        df = df.Filter("Leading_Lepton_promptgenmatched", "Gen Matching of the leading Lepton")
        results["Cutflow_LeptonGenMatching"] = df.Histo1D(("h_cutflow_LeptonGenMatching","Cutflow LeptonGenMatching",1,-0.5,0.5),"cutflow_stage","weight")
        df = df.Define("Lepton_ID_SF","getLeptonIdSF(Leading_Lepton_pdgId,Leading_Lepton_isTight,Lepton_tightElectron_mvaFall17V2Iso_WP90_TotSF,Lepton_tightMuon_cut_Tight_HWWW_TotSF)")
        df = df.Redefine("weight","weight*Lepton_ID_SF")
        
        if hasattr(ROOT, "initializeEleTriggerSF"):
            ROOT.initializeEleTriggerSF()
        
        df = df.Define("EleTriggerSF","getEleTriggerSF(Leading_Lepton_pdgId,Leading_Lepton_pt,Leading_Lepton_eta)")
        df = df.Define("LepTriggerSF","getTriggerSF(Leading_Lepton_pdgId,EleTriggerSF,TriggerEffWeight_1l)")
        df = df.Redefine("weight", "weight*LepTriggerSF")
       

    df = df.Define("isAnalysisLepton","isAnalysisLepton(Leading_Lepton_pdgId,Leading_Lepton_pt,Leading_Lepton_eta,Leading_Lepton_phi)")
    df = df.Filter("isAnalysisLepton", "Analysis Lepton Selection")


    results["Cutflow_AnaLepton"] = df.Histo1D(("h_cutflow_AnaLepton","Cutflow AnaLepton",1,-0.5,0.5),"cutflow_stage","weight")
    
    df = df.Define("isLeptonHole_ex","isHoleLepton(Leading_Lepton_eta,Leading_Lepton_phi,Leading_Lepton_pdgId)")
    df = df.Filter("!isLeptonHole_ex","Is Hole Lepton")
    results["Cutflow_notHoleLepton"] = df.Histo1D(("h_cutflow_notHoleLepton","Cutflow notHoleLepton",1,-0.5,0.5),"cutflow_stage","weight")

    
    df = df.Define("isVetoLepton","isVetoLepton(nLepton,Lepton_pt,Lepton_isLoose)")
    df = df.Filter("!isVetoLepton", "Veto Lepton Cut")
    results["Cutflow_Veto_Lepton"] = df.Histo1D(("h_cutflow_Veto_Lepton","Cutflow Veto Lepton",1,-0.5,0.5),"cutflow_stage","weight")

    df = df.Filter("PuppiMET_pt > 30", "PuppiMET_pt > 30 GeV cut")
    results["Cutflow_MET"] = df.Histo1D(("h_cutflow_MET","Cutflow MET",1,-0.5,0.5),"cutflow_stage","weight")


    results["nCleanedFatJet"] = df.Histo1D(("nCleanedFatJet","nCleanedFatJet",6,-0.5,5.5),"nCleanFatJet")
    
    results["nGenJetAK8"] = df.Histo1D(("h_nGenJetAK8","nGenJetAK8",6,-0.5,5.5),"nGenJetAK8","weight")
    results["nFatJet"] = df.Histo1D(("h_nFatJet","nGenJetAK8",6,-0.5,5.5),"nFatJet","weight")

    #df = df.Filter("nGenJetAK8 >=1","GenJet cut")
    #df = df.Define("Lead_GenJetAK8_mass","GenJetAK8_mass.at(0)")
    #results["h_Lead_GenJetAK8_mass"] = df.Histo1D(("h__GenJetAK8_mass","Lead_GenJetAK8_mass",20, 0,100),"Lead_GenJetAK8_mass", "weight")
    results["All_Jet_pt_1"] = df.Histo1D(("h_All_Jet_pt_1","All Jet pt",50,0,1000),"FatJet_pt","weight") 

    # Jet Selection
    df = df.Filter("nFatJet>=1","At Least 1 Fat Jet")
    if isMC:
        df = df.Define("JetPUIDSF","computePUJetIdSF(nJet,Jet_jetId,Jet_electronIdx1,Jet_muonIdx1,Jet_PUIDSF_loose,Leading_Lepton_electronIdx,Leading_Lepton_muonIdx)")
        df = df.Redefine("weight","weight*JetPUIDSF")
    
    results["Cutflow_nFatJet"] = df.Histo1D(("h_cutflow_nFatJet","Cutflow nFatjet",1,-0.5,0.5),"cutflow_stage","weight")
    
    df = df.Define("Lead_jet_pt","FatJet_pt[0]")
    results["Jet_pt_1"] = df.Histo1D(("h_Jet_pt_1","Jet pt",50,0,1000),"Lead_jet_pt","weight") 

 
    df = df.Define("GoodFatJet_idx","isGoodFatjet_indx(FatJet_eta,FatJet_phi,FatJet_jetId,Lepton_eta,Lepton_phi)")
    df = df.Filter("!(GoodFatJet_idx == -1)", "Good Fat Jet cut")
    results["Cutflow_JetCleaning"] = df.Histo1D(("h_cutflow_JetCleaning","Cutflow JetCleaning",1,-0.5,0.5),"cutflow_stage","weight")
    df = df.Define("AnaFatJet_pt","FatJet_pt[GoodFatJet_idx]")
    df = df.Define("AnaFatJet_eta","FatJet_eta[GoodFatJet_idx]")
    df = df.Define("AnaFatJet_phi", "FatJet_phi[GoodFatJet_idx]")
    df = df.Define("AnaFatJet_mass", "FatJet_mass[GoodFatJet_idx]")
    df = df.Define("AnaFatJet_msoftdrop", "FatJet_msoftdrop[GoodFatJet_idx]")
    df = df.Define("AnaFatJet_jetId","FatJet_jetId[GoodFatJet_idx]")
    
    results["Jet_pt"] = df.Histo1D(("h_Jet_pt","Jet pt",50,0,1000),"AnaFatJet_pt","weight") 

    
   
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
    
    #df = df.Define("isJetHole_ex","isHole_ex(AnaFatJet_eta,AnaFatJet_phi)")
    #df = df.Filter("!isJetHole_ex","Is Hole Jet")
    results["Cutflow_notHoleJet"] = df.Histo1D(("h_cutflow_notHoleJet","Cutflow notHoleJet",1,-0.5,0.5),"cutflow_stage","weight")
    
    df = df.Filter("AnaFatJet_pt>200","Jet pT cut")
    results["Cutflow_Jet_Pt"] = df.Histo1D(("h_cutflow_Jet_Pt","Cutflow Jet Pt",1,-0.5,0.5),"cutflow_stage","weight")
    
    df = df.Filter("(AnaFatJet_msoftdrop > 65 && AnaFatJet_msoftdrop < 105)","Mass cut")
    results["Cutflow_Jet_mass"] = df.Histo1D(("h_cutflow_Mass_cut","Cutflow Jet Mass",1,-0.5,0.5),"cutflow_stage","weight")
    
    
    df = df.Define("CleanJet_btag", "Take(Jet_btagDeepFlavB, CleanJet_jetIdx)")
    if isMC:
        df = df.Define("CleanJet_btag_SF","Take(Jet_btagSF_deepjet_shape,CleanJet_jetIdx)")
    
    df = df.Define("CleanJet_notOverlapping", "getCleanJetNotOverlapping(AnaFatJet_eta, AnaFatJet_phi, CleanJet_eta, CleanJet_phi, CleanJet_pt)")

    df = df.Define("CleanJet_btag_notOverlap","CleanJet_btag[CleanJet_notOverlapping]")
    
    if isMC:
        df = df.Define("CleanJet_btagSF_notOverlap","CleanJet_btag_SF[CleanJet_notOverlapping]")
        df = df.Define("bTagSF","getBTagSF(CleanJet_btagSF_notOverlap)")
        df = df.Redefine("weight","weight*bTagSF")

    btagWP = 0.2783
    df = df.Define("bTagged_jets",f"CleanJet_btag_notOverlap[CleanJet_btag_notOverlap > {btagWP}]")
    df = df.Filter("bTagged_jets.size() == 0", "bVeto cut")
    results["Cutflow_bVeto"] = df.Histo1D(("h_cutflow_bVeto","Cutflow bVeto",1,-0.5,0.5),"cutflow_stage","weight")

    if wtagger == "Nominal":
        results["Jet_Nominal_WTagger"] = df.Histo1D(("h_Nominal_WTagger", "WTagger Nominal", 10, 0, 1), "AnaFatJet_nom_wtag","weight")
    elif wtagger == "MD":
        results["Jet_MD_WTagger"] = df.Histo1D(("h_MD_WTagger", "WTagger MD", 10, 0, 1), "AnaFatJet_md_wtag","weight") 
    
    if wtagger == "Nominal":
        df = df.Filter("AnaFatJet_nom_wtag > 0.98","Wtagger Nominal 0p5 cut")
        if isMC:
            df = df.Define("WTagger_SF","getWTaggerSF(AnaFatJet_pt)")
            df = df.Redefine("weight","weight*WTagger_SF")
    elif wtagger == "MD":
        df = df.Filter("AnaFatJet_md_wtag > 0.90"," Wtagger MD 0p5 cut")
        if isMC:
            df = df.Define("WTagger_SF","getWTaggerSF(AnaFatJet_pt)")
            df = df.Redefine("weight","weight*WTagger_SF")
    results["Cutflow_WTagger"] = df.Histo1D(("h_cutflow_WTagger","Cutflow WTagger",1,-0.5,0.5),"cutflow_stage","weight")
    results["Lepton_pt_selection"] = df.Histo1D(("h_Lepton_pt_selection","Lepton pt",50,0,1000),"Leading_Lepton_pt","weight")
    results["Lepton_eta_selection"] = df.Histo1D(("h_Lepton_eta_selection", "Lepton eta", 25, -2.5, 2.5), "Leading_Lepton_eta", "weight")
    results["Lepton_phi_selection"] = df.Histo1D(("h_Lepton_phi_selection", "Lepton phi", 16, -3.2, 3.2), "Leading_Lepton_phi", "weight")    

    results["Jet_pt_selection"] = df.Histo1D(("h_Jet_pt_selection","Jet pt",50,0,1000),"AnaFatJet_pt","weight")
    results["Jet_eta_selection"] = df.Histo1D(("h_Jet_eta_selection", "Jet eta", 25, -2.5, 2.5), "AnaFatJet_eta", "weight")
    results["Jet_phi_selection"] = df.Histo1D(("h_Jet_phi_selection", "Jet phi", 16, -3.2, 3.2), "AnaFatJet_phi", "weight")    
    results["Jet_mass_selection"] = df.Histo1D(("h_Jet_mass_selection","Jet mass",25,0,250),"AnaFatJet_msoftdrop","weight")
    if wtagger == "Nominal":
        results["Jet_Nominal_WTagger_selection"] = df.Histo1D(("h_Nominal_WTagger_selection", "WTagger Nominal", 10, 0, 1), "AnaFatJet_nom_wtag","weight")
    elif wtagger == "MD":
        results["Jet_MD_WTagger_selection"] = df.Histo1D(("h_MD_WTagger_selection", "WTagger MD", 10, 0, 1), "AnaFatJet_md_wtag","weight")
    
    df = df.Define("H_vis_m","getHiggsCandidate(Leading_Lepton_pt,Leading_Lepton_eta,Leading_Lepton_phi,AnaFatJet_pt,AnaFatJet_eta,AnaFatJet_phi,AnaFatJet_msoftdrop,0)") 
    df = df.Define("H_vis_pt","getHiggsCandidate(Leading_Lepton_pt,Leading_Lepton_eta,Leading_Lepton_phi,AnaFatJet_pt,AnaFatJet_eta,AnaFatJet_phi,AnaFatJet_msoftdrop,1)") 
    df = df.Define("H_vis_eta","getHiggsCandidate(Leading_Lepton_pt,Leading_Lepton_eta,Leading_Lepton_phi,AnaFatJet_pt,AnaFatJet_eta,AnaFatJet_phi,AnaFatJet_msoftdrop,2)") 
    df = df.Define("H_vis_phi","getHiggsCandidate(Leading_Lepton_pt,Leading_Lepton_eta,Leading_Lepton_phi,AnaFatJet_pt,AnaFatJet_eta,AnaFatJet_phi,AnaFatJet_msoftdrop,3)") 
    results["Higgs_pt"] = df.Histo1D(("h_Higgs_pt","Higgs pt",50,0,1000),"H_vis_pt","weight")
    results["Higgs_eta"] = df.Histo1D(("h_Higgs_eta", "Higgs eta", 25, -2.5, 2.5), "H_vis_eta", "weight")
    results["Higgs_phi"] = df.Histo1D(("h_Higgs_phi", "Higgs phi", 16, -3.2, 3.2), "H_vis_phi", "weight")
    results["Higgs_mass"] = df.Histo1D(("h_Higgs_mass", "Higgs mass", 75, 0, 1500), "H_vis_m", "weight")    
   
    results["Met_pt"] = df.Histo1D(("h_MET_pt", "MET pt", 50,0, 1000), "PuppiMET_pt", "weight")
    results["Met_phi"] = df.Histo1D(("h_MET_phi", "MET phi", 16, -3.2, 3.2), "PuppiMET_phi", "weight")

    if isMC:
        sum_genWeight_value = sum_genWeight.GetValue()
        print(f"sum_genWeight = {sum_genWeight_value}")
        h1.SetBinContent(2,sum_genWeight_value)
        results["genWeightSum"] = h1
        sum_Weight = df.Sum("weight")
        print(f"sum_Weight = {sum_Weight.GetValue()}")
      
    report = df.Report()
    report.Print() 
 
    return results

parser = argparse.ArgumentParser()
parser.add_argument("-w","--wtag", help="WTagger option Nominal, MD", type=str,default = "Nominal")
parser.add_argument("-r","--run", help="run option", choices=["sig", "sig+sbi", "MC", "bkg", "data"], type=str,default = "sig")
args = parser.parse_args()


histograms = {}

if args.run == "MC":
    del dataset["data"]
elif args.run == "bkg":
    del dataset["data"]
    del dataset["ggH_sand_off"]
    del dataset["ggH_sonly_off"]

print(f"{args.run}")
if args.run == "data":
    print("Wrong1")
#    histograms["data"] = makeRDF("data",args.wtag)
elif args.run == "sig+sbi":
    print("Wrong2")
#    histograms["ggH_sonly_off"] = makeRDF("ggH_sonly_off",args.wtag)
#    histograms["ggH_sand_off"] = makeRDF("ggH_sand_off",args.wtag)
elif args.run == "sig":
    print("Wrong3")
    #histograms["ggH_sonly_off"] = makeRDF("ggH_sonly_off",args.wtag)
    #histograms["ggH_bonly_off"] = makeRDF("ggH_bonly_off",args.wtag)
    #histograms["ggH_sand_off"] = makeRDF("ggH_sand_off",args.wtag)
    #histograms["WW"] = makeRDF("WW",args.wtag)
    histograms["DYJetsToLL_M-50"] = makeRDF("DYJetsToLL_M-50",args.wtag)
    histograms["DY_else"] = makeRDF("DY_else",args.wtag)
else:
    print("Right")
    for keys in dataset:
        print(f"Running for {keys}")
        histograms[keys] = makeRDF(keys,args.wtag)
#print(histograms)

#file_path = "my_histograms_ggH_sonly_off_Step_3.pkl"

#print(f"Saving dictionary of histograms to {file_path}")
#with open(file_path, "wb") as f:
    # Use pickle.dump() to save the dictionary
#    pickle.dump(histograms, f)

output_file = ROOT.TFile("output_8_Dec.root", "RECREATE")
for key1 in histograms:
    new_directory = output_file.mkdir(key1)
    new_directory.cd()
    for key2 in histograms[key1]:
        histograms[key1][key2].Write()
    output_file.cd()
output_file.Close()



