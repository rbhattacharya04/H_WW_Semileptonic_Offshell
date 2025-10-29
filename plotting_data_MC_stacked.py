# ROOT imports 
import os, ROOT
import cmsstyle as CMS
import narf

# python imports 
import matplotlib.pyplot as plt  # matplotlib library
import mplhep as hep  # HEP (CMS) extensions/styling on top of mpl

# For constructing examples
import hist  # histogramming library
import numpy as np 
import uproot
import math
import shutil


plots = ["h_Jet_pt","h_Nominal_WTagger","h_Lepton_pt_selection","h_Lepton_eta_selection","h_Lepton_phi_selection","h_Jet_pt_selection","h_Jet_eta_selection","h_Jet_phi_selection","h_Nominal_WTagger_selection","h_Higgs_pt","h_Higgs_eta","h_Higgs_phi","h_Higgs_mass"]

bkg_samples = ["W + jets","Top","DY","Multiboson","ggH_sonly_on","ggH_sand_on"]
bkg_samples.reverse()
Multi_bosons = ["WW","Vg","VgS","VZ","VVV","WWewk","VBF_V"]
signal = ["ggH_sonly_off"]
SBI = ["ggH_sand_off"]

#root_file = ROOT.TFile.Open("output_saved.root","READ")
root_file = ROOT.TFile.Open("output.root","READ")



output_path = "/eos/user/s/sverma/www/Marc_TNP_Plots/H_WW/28Oct"
os.makedirs(output_path,exist_ok=True)
shutil.copy("index.php", output_path)
for plot in plots:
    print(plot)
    #Styling
    hep.style.use("CMS")
    fig, ax = plt.subplots()
    hep.cms.label("Preliminary", com = 13, lumi = 59.7, data = True, loc=0, ax=ax);
    ax.set_ylabel("Events")
    ax.set_yscale('log')
    ax.set_xlabel(f"{plot}")
    histo_list = []
    legend_list = []
    for bkg in bkg_samples:
        print(bkg)
        if bkg == "Multiboson":
                h_ww = root_file.Get(f"WW/{plot}")
                histogram = h_ww.Clone()
                for sample in Multi_bosons:
                    if sample == "WW": 
                        continue
                    h = root_file.Get(f"{sample}/{plot}")
                    histogram.Add(h)
        else:
            histogram = root_file.Get(f"{bkg}/{plot}")
        integral = histogram.Integral()
        h= narf.root_to_hist(histogram)
        histo_list.append(h)
        legend_list.append(f"{bkg}: {integral:.2f}") 
    hep.histplot(histo_list, ax=ax, stack=True, histtype='fill', label=legend_list)
    hist_signal = root_file.Get(f"ggH_sonly_off/{plot}")
    integral_signal = hist_signal.Integral()
    h_signal = narf.root_to_hist(hist_signal)
    hep.histplot(h_signal, ax=ax, histtype='step', linewidth=3,color='green', label=f"Signal: {integral_signal:.2f}")
    hist_sbi = root_file.Get(f"ggH_sand_off/{plot}")
    integral_sbi = hist_sbi.Integral()
    h_sbi = narf.root_to_hist(hist_sbi)
    hep.histplot(h_sbi, ax=ax, histtype='step', linewidth=3,color='blue', label=f"SBI: {integral_sbi:.2f}")
    ax.legend()
    plt.savefig(f'{output_path}/{plot}.png')  
    
        
        

