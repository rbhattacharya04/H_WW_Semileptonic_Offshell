# ROOT imports 
import os, ROOT
import cmsstyle as CMS
#import narf

# python imports 
import matplotlib.pyplot as plt  # matplotlib library
import mplhep as hep  # HEP (CMS) extensions/styling on top of mpl

# For constructing examples
import hist  # histogramming library
import numpy as np 
import uproot
import math
import shutil

from matplotlib.colors import ListedColormap
#petroff10 = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]
petroff10 = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70"]
#plots = ["h_Jet_pt","h_Nominal_WTagger","h_Lepton_pt_selection","h_Lepton_eta_selection","h_Lepton_phi_selection","h_Jet_pt_selection","h_Jet_eta_selection","h_Jet_phi_selection","h_Jet_mass_selection","h_Nominal_WTagger_selection","h_Higgs_pt","h_Higgs_eta","h_Higgs_phi","h_Higgs_mass", "h_MET_pt","h_MET_phi"]

plots = {
    #"h_Jet_pt_selection" : {
    #    "x_axis" : r"$p_{T}^{Jet}$", 
    #    "y_axis": "Events/20 GeV",
    #} ,
    #"h_Jet_eta_selection" : {
    #    "x_axis" : r"$\eta^{Jet}$",
    #    "y_axis" : "Events/0.2",
    #},
    #"h_Jet_phi_selection" : {
    #    "x_axis" : r"$\phi^{Jet}$",
    #    "y_axis" : "Events/0.4",
    #},
    "h_Jet_mass_selection" : {
        "x_axis" : r"$m^{Jet}$",
        "y_axis" : "Events/10 GeV",
    },
    "h_Nominal_WTagger_selection" : {
        "x_axis" : "WTagger_Score",
        "y_axis" : "Events/0.1",
    },
    "h_Nominal_TopTagger_selection": {
        "x_axis" : "TopTagger_Score",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_tau1_selection": {
        "x_axis" : r"$\tau_{1}$",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_tau2_selection": {
        "x_axis" : r"$\tau_{2}$",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_tau3_selection": {
        "x_axis" : r"$\tau_{3}$",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_tau4_selection": {
        "x_axis" : r"$\tau_{4}$",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_tau2_tau1_selection": {
        "x_axis" : r"$\frac{\tau_{2}}{\tau_{1}}$",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_tau3_tau2_selection": {
        "x_axis" : r"$\frac{\tau_{3}}{\tau_{2}}$",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_subJetIdx1_selection": {
        "x_axis" : "Jet_subJetIdx1",
        "y_axis" : "Events/1",
    },
    "h_Jet_n2b1_selection": {
        "x_axis" : "n2b1",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_n3b1_selection": {
        "x_axis" : "n3b1",
        "y_axis" : "Events/0.1",
    },
    "h_Jet_rawFactor_selection": {
        "x_axis" : "Jet_rawFactor",
        "y_axis" : "Events/10",
    },
    "h_Jet_electronIdx3SJ": {
        "x_axis" : "Jet_electronIdx3SJ",
        "y_axis" : "Events/1",
    },
    "h_Jet_muonIdx3SJ": {
        "x_axis" : "Jet_muonIdx3SJ",
        "y_axis" : "Events/1",
    },
    "h_Jet_nConstituents": {
        "x_axis" : "Jet_nConstituents",
        "y_axis" : "Events/1",
    },
    "h_Jet_genJetAK8Idx": {
        "x_axis" : "Jet_genJetAK8Idx",
        "y_axis" : "Events/1",
    },
    "h_Jet_hadronFlavour": {
        "x_axis" : "Jet_hadronFlavour",
        "y_axis" : "Events/1",
    },
    "h_Jet_nBHadrons": {
        "x_axis" : "Jet_nBHadrons",
        "y_axis" : "Events/1",
    },
    "h_Jet_nChadrons": {
        "x_axis" : "Jet_nCHadrons",
        "y_axis" : "Events/1",
    },
    
}

#bkg_samples = ["W + jets","Top","DY","Multiboson","ggH_sonly_on","ggH_sand_on"]
#bkg_samples.reverse()
#Multi_bosons = ["WW","Vg","VgS","VZ","VVV","WWewk","VBF_V"]
#signal = ["ggH_sonly_off"]
#SBI = ["ggH_sand_off"]
#bkg_samples = []
#SBI = []i

#bkg_samples = {"Top" : ["TTToSemiLeptonic", "TTTo2L2Nu","ST_tW_top","ST_tW_antitop","ST_t-channel_top","ST_t-channel_antitop", "TTWJetsToLNu", "ST_s-channel"]}
bkg_samples = {"Top" : ["TTToSemiLeptonic"]}
bkg_samples["Top"].reverse()
uproot_file = uproot.open("output_15Nov_v3.root")
 

output_path = "/eos/user/r/rbhattac/www/Marc_TNP_Plots/H_WW/Debug_Plots/test_plot_18_Nov_2025_High_WTagger_TTToSemiLeptonic"
os.makedirs(output_path,exist_ok=True)
shutil.copy("index.php", output_path)

for plot in plots:
    #plot = "h_Jet_phi_selection"
    print(plot)
    #Styling
    hep.style.use("CMS")
    fig, ax = plt.subplots()
    hep.cms.label("Preliminary", com = 13, lumi = 59.7, data = True, loc=0, ax=ax);
    ax.set_ylabel(plots[plot]["y_axis"])
    #ax.set_yscale('log')
    #plt.yscale('log')
    axis_label = plot.replace("h_","")
    axis_label = axis_label.replace("_selection","")
    ax.set_xlabel(plots[plot]["x_axis"])
    histo_list_plot = []
    legend_list = []

    for key in bkg_samples:
        #histolist = []
        for bkg in bkg_samples[key]:
            print(bkg)
            h = uproot_file[f"{bkg}/{plot}"].to_boost()
            h = h * 59.7
            integral = h.view(flow=True).sum().value
            #histolist.append(histogram)
            histo_list_plot.append(h)
            legend_list.append(f"{bkg}")
        #h_sum = histolist[0].Clone()
        #for i in range(len(histolist)):
        #    if i==0:
        #        continue
        #    h_sum.Add(histolist[i])
        #integral = h_sum.Integral()
        #h = narf.root_to_hist(h_sum)
        #histo_list_plot.append(h)
        #legend_list.append(f"{key}") 

    hep.histplot(histo_list_plot, ax=ax, stack=True, histtype='fill', label=legend_list)#,color=petroff10)
    h_signal = uproot_file[f"ggH_sonly_off/{plot}"].to_boost()
    #hist_signal.Rebin(2)
    h_signal = h_signal * 59.7
    #integral_signal = hist_signal.Integral(0,hist_signal.GetNbinsX()+1)
    integral_signal = h_signal.view(flow=True).sum().value
    hep.histplot(h_signal, ax=ax, histtype='step', linewidth=3, label=f"Signal")#, color="#717581")
    #hist_sbi = root_file.Get(f"ggH_sand_off/{plot}")
    #integral_sbi = hist_sbi.Integral()
    #h_sbi = narf.root_to_hist(hist_sbi)
    #hep.histplot(h_sbi, ax=ax, histtype='step', linewidth=3,color='blue', label=f"SBI: {integral_sbi:.2f}")
    #ax.set_yscale("log")
    ax.legend()
    plt.savefig(f'{output_path}/{plot}.png')  
    
        
        

