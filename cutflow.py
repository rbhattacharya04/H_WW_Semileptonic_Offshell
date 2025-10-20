import ROOT
import pandas as pd

cutflows = ["h_cutflow_Start","h_cutflow_Trigger","h_cutflow_LeptonGenMatching","h_cutflow_AnaLepton","h_cutflow_Veto_Lepton","h_cutflow_MET","h_cutflow_nFatJet","h_cutflow_JetCleaning","h_cutflow_bVeto","h_cutflow_Jet_Pt","h_cutflow_WTagger"]

bkg_samples = ["W + jets","Top","DY","WW","Vg","VgS","VZ","VVV","WWewk","VBF_V","ggH_sonly_on","ggH_sand_on"]
signal = ["ggH_sonly_off"]
SBI = ["ggH_sand_off"]
cutflow_bkg = {}
cutflow_bkg_weighted = {}
cutflow_signal = {}
cutflow_signal_weighted = {}
cutflow_sbi = {}
cutflow_sbi_weighted = {}

root_file = ROOT.TFile.Open("output.root","READ")

pd.set_option('display.float_format', '{:.2f}'.format)

for cut in cutflows:
    #bkg_integral = 0
    #bkg_nentries = 0
    #for bkg in bkg_samples:
    #    histogram = root_file.Get(f"{bkg}/{cut}")
    #    bkg_nentries += histogram.GetEntries()
    #    bkg_integral += histogram.Integral()
    #cutflow_bkg[f"{cut}"] = {"value" : bkg_nentries}
    #cutflow_bkg_weighted[f"{cut}"] = {"value" : bkg_integral}
    histo_signal = root_file.Get(f"ggH_sonly_off/{cut}")
    cutflow_signal[f"{cut}"] = {"value" : histo_signal.GetEntries()}
    cutflow_signal_weighted[f"{cut}"] = {"value" : histo_signal.Integral()}
    #histo_sbi = root_file.Get(f"ggH_sand_off/{cut}")
    #cutflow_sbi[f"{cut}"] = {"value" : histo_sbi.GetEntries()}
    #cutflow_sbi_weighted[f"{cut}"] = {"value" : histo_sbi.Integral()}


df_signal = pd.DataFrame.from_dict(cutflow_signal,orient='index')
df_signal['cum_eff'] = df_signal['value']/df_signal['value'].shift(1)
start = df_signal['value'].iloc[0]
df_signal['abs_eff'] = df_signal['value']/start
print("Signal_CutFlow")
print(df_signal)

df_signal_weighted = pd.DataFrame.from_dict(cutflow_signal_weighted,orient='index')
df_signal_weighted['cum_eff'] = df_signal_weighted['value']/df_signal_weighted['value'].shift(1)
start = df_signal_weighted['value'].iloc[0]
df_signal_weighted['abs_eff'] = df_signal_weighted['value']/start
print("Signal_CutFlow_Weighted")
print(df_signal_weighted)

#df_bkg = pd.DataFrame.from_dict(cutflow_bkg,orient='index')
#df_bkg['cum_eff'] = df_bkg['value']/df_bkg['value'].shift(1)
#start = df_bkg['value'].iloc[0]
#df_bkg['abs_eff'] = df_bkg['value']/start
#print("Bkg_CutFlow")
#print(df_bkg)

#df_bkg_weighted = pd.DataFrame.from_dict(cutflow_bkg_weighted,orient='index')
#df_bkg_weighted['cum_eff'] = df_bkg_weighted['value']/df_bkg_weighted['value'].shift(1)
#start = df_bkg_weighted['value'].iloc[0]
#df_bkg_weighted['abs_eff'] = df_bkg_weighted['value']/start
#print("Bkg_CutFlow_weighted")
#print(df_bkg_weighted)

#df_sbi = pd.DataFrame.from_dict(cutflow_sbi,orient='index')
#df_sbi['cum_eff'] = df_sbi['value']/df_sbi['value'].shift(1)
#start = df_sbi['value'].iloc[0]
#df_sbi['abs_eff'] = df_sbi['value']/start
#print("SBI_CutFlow")
#print(df_sbi)

#df_sbi_weighted = pd.DataFrame.from_dict(cutflow_sbi_weighted,orient='index')
#df_sbi_weighted['cum_eff'] = df_sbi_weighted['value']/df_sbi_weighted['value'].shift(1)
#start = df_sbi_weighted['value'].iloc[0]
#df_sbi_weighted['abs_eff'] = df_sbi_weighted['value']/start
#print("SBI_CutFlow_Weighted")
#print(df_sbi_weighted)
