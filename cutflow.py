import ROOT
import pandas as pd
#from Dataset import dataset

#dataset_list = list(dataset.keys())
#bkg_samples = dataset_list
#bkg_samples.remove("data")
#bkg_samples.remove("ggH_sonly_off")
#bkg_samples.remove("ggH_sand_off")
cutflows = ["h_cutflow_Start","h_cutflow_Trigger","h_cutflow_LeptonGenMatching","h_cutflow_AnaLepton","h_cutflow_notHoleLepton", "h_cutflow_Veto_Lepton","h_cutflow_MET","h_cutflow_nFatJet","h_cutflow_JetCleaning","h_cutflow_notHoleJet", "h_cutflow_Jet_Pt","h_cutflow_Mass_cut","h_cutflow_bVeto","h_cutflow_WTagger"]

#bkg_samples = ["W + jets","Top","DY","WW","Vg","VgS","VZ","VVV","WWewk","VBF_V","ggH_sonly_on","ggH_sand_on"]
signal = ["ggH_sonly_off"]
SBI = ["ggH_sand_off"]
cutflow_bkg = {}
cutflow_bkg_weighted = {}
cutflow_signal = {}
cutflow_signal_weighted = {}
cutflow_sbi = {}
cutflow_sbi_weighted = {}


lumi = 59.7
#root_file = ROOT.TFile.Open("output.root","READ")
root_file = ROOT.TFile.Open("output_19_Nov_correct_WTagger_full_bkg_v2.root","READ")

#pd.set_option('display.float_format', '{:.2f}'.format)
bkg_samples = ["TTToSemiLeptonic"]#,"TTTo2L2Nu","TTWJetsToLNu","ST_s-channel","ST_t-channel_antitop","ST_t-channel_top","ST_tW_antitop","ST_tW_top","DYJetsToLL_M-50","DY_else","ggH_bonly_off","qqWWqq","WJetsToLNu-LO","WGToLNuG","ZGToLLG","ZZ","WZ","VVV"]
#bkg_samples = []
#bkg_samples = ["ggH_bonly_off"]#, "WJetsToLNu_HT70To100", "WJetsToLNu_HT100To200", "WJetsToLNu_HT200To400", "WJetsToLNu_HT400To600" , "WJetsToLNu_HT600To800" , "WJetsToLNu_HT800To1200", "WJetsToLNu_HT1200To2500", "WJetsToLNu_HT2500ToInf"]
cutflow_bkg_ind = {}
default_value = {}
cutflow_bkg_ind_weighted = dict.fromkeys(bkg_samples,default_value)
for cut in cutflows:
    print(f"{cut}")
    bkg_integral = 0
    bkg_nentries = 0
    for bkg in bkg_samples:
        histogram = root_file.Get(f"{bkg}/{cut}")
        bkg_nentries += histogram.GetEntries()
        bkg_integral += histogram.Integral()
        #cutflow_bkg_ind[bkg][f"{cut}"] = {"value": histogram.GetEntries()}
        #cutflow_bkg_ind_weighted[f"{bkg}"][f"{cut}"] = {"value": lumi*histogram.Integral()}
        #print(cutflow_bkg_ind_weighted[f"{bkg}"][f"{cut}"]["value"])
    cutflow_bkg[f"{cut}"] = {"value" : bkg_nentries}
    cutflow_bkg_weighted[f"{cut}"] = {"value" : lumi*bkg_integral}
    histo_signal = root_file.Get(f"ggH_sonly_off/{cut}")
    cutflow_signal[f"{cut}"] = {"value" : histo_signal.GetEntries()}
    cutflow_signal_weighted[f"{cut}"] = {"value" : lumi*histo_signal.Integral()}
    #histo_signal = root_file.Get(f"VVV/{cut}")
    #cutflow_signal[f"{cut}"] = {"value" : histo_signal.GetEntries()}
    #cutflow_signal_weighted[f"{cut}"] = {"value" : lumi*histo_signal.Integral()}
    #histo_sbi = root_file.Get(f"ggH_sand_off/{cut}")
    #cutflow_sbi[f"{cut}"] = {"value" : histo_sbi.GetEntries()}
    #cutflow_sbi_weighted[f"{cut}"] = {"value" : histo_sbi.Integral()}


#df_signal = pd.DataFrame.from_dict(cutflow_signal,orient='index')
#df_signal['rel_eff'] = df_signal['value']/df_signal['value'].shift(1)
#start = df_signal['value'].iloc[0]
#df_signal['abs_eff'] = df_signal['value']/start
#print("Signal_CutFlow")
#print(df_signal)

df_signal_weighted = pd.DataFrame.from_dict(cutflow_signal_weighted,orient='index')
df_signal_weighted = df_signal_weighted.rename(index={"h_cutflow_Start": "h_cutflow_TightLepton_METFilter_EMTFBugVeto"})
start_value = {"h_start": {"value" : 26011.29}}
new_row_df = pd.DataFrame.from_dict(start_value,orient='index')
df_signal_weighted = pd.concat([new_row_df,df_signal_weighted])
df_signal_weighted['rel_eff'] = df_signal_weighted['value']/df_signal_weighted['value'].shift(1)
start = df_signal_weighted['value'].iloc[0]
df_signal_weighted['abs_eff'] = df_signal_weighted['value']/start
print("Signal_CutFlow_Weighted")
print(df_signal_weighted)

#df_signal = pd.DataFrame.from_dict(cutflow_signal,orient='index')
#df_signal['cum_eff'] = df_signal['value']/df_signal['value'].shift(1)
#start = df_signal['value'].iloc[0]
#df_signal['abs_eff'] = df_signal['value']/start
#print("DY_CutFlow")
#print(df_signal)

#df_signal_weighted = pd.DataFrame.from_dict(cutflow_signal_weighted,orient='index')
#df_signal_weighted['cum_eff'] = df_signal_weighted['value']/df_signal_weighted['value'].shift(1)
#start = df_signal_weighted['value'].iloc[0]
#df_signal_weighted['abs_eff'] = df_signal_weighted['value']/start
#print("DY_CutFlow_Weighted")
#print(df_signal_weighted)

#df_bkg = pd.DataFrame.from_dict(cutflow_bkg,orient='index')
#df_bkg['cum_eff'] = df_bkg['value']/df_bkg['value'].shift(1)
#start = df_bkg['value'].iloc[0]
#df_bkg['abs_eff'] = df_bkg['value']/start
#print("Bkg_CutFlow")
#print(df_bkg)

df_bkg_weighted = pd.DataFrame.from_dict(cutflow_bkg_weighted,orient='index')
df_bkg_weighted['rel_eff'] = df_bkg_weighted['value']/df_bkg_weighted['value'].shift(1)
start = df_bkg_weighted['value'].iloc[0]
df_bkg_weighted['abs_eff'] = df_bkg_weighted['value']/start
print("Bkg_CutFlow_weighted")
print(df_bkg_weighted)

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
#print(cutflow_bkg_ind_weighted)

#for keys in cutflow_bkg_ind_weighted:
#    cutflow = cutflow_bkg_ind_weighted[keys]
#    df = pd.DataFrame.from_dict(cutflow,orient='index')
#    df['rel_eff'] = df['value']/df['value'].shift(1)
#    start = df['value'].iloc[0]
#    df['abs_eff'] = df['value']/start
#    print(f"Cutflow for {keys}")
#    print(df)
