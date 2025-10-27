import os

def list_files_in_folder(folder_path, search_string=None):
    found_files = []
    try:
        for filename in os.listdir(folder_path):
            if filename.endswith(".swp"):
                continue
            if filename.endswith(".swx"):
                continue
            file_path = os.path.join(folder_path, filename)
            if os.path.isfile(file_path):
                if search_string is None or search_string in filename:
                    found_files.append(file_path)
    except FileNotFoundError:
        print(f"Error: Folder not found at '{folder_path}'")
    return found_files
#Step_3
mc_path = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer20UL18_106x_nAODv9_Full2018v9/MCl1loose2018v9__MCCorr2018v9NoJERInHorn__MCCombJJLNu2018"
#Step_2
#mc_path = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer20UL18_106x_nAODv9_Full2018v9/MCl1loose2018v9__MCCorr2018v9NoJERInHorn"

#Step_1
#mc_path = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer20UL18_106x_nAODv9_Full2018v9/MCl1loose2018v9"

#Signal NanoAOD
sig_path = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Sig"

dataset = {
    "data": {
        "files": list_files_in_folder("/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Run2018_UL2018_nAODv9_Full2018v9/DATAl1loose2018v9__DATACombJJLNu2018/"),
        "isMC": False,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None,
        "sample_filters": None
    },

    "DYJetsToLL_M-50": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "(LHE_HT < 70)",
        "sample_filters": "DYPhotonFilter"
    },

    "DYJetsToLL_M-50_HT-70to100": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-70to100"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.0",
        "sample_filters": "DYPhotonFilter"
    },
    "DYJetsToLL_M-50_HT-100to200": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-100to200"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.000",
        "sample_filters": "DYPhotonFilter"
    },
    "DYJetsToLL_M-50_HT-200to400": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-200to400"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "0.999",
        "sample_filters": "DYPhotonFilter"
    },
    "DYJetsToLL_M-50_HT-400to600": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-400to600"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "0.990",
        "sample_filters": "DYPhotonFilter"
    },
    "DYJetsToLL_M-50_HT-600to800": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-600to800"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "0.975",
        "sample_filters": "DYPhotonFilter"
    },
    "DYJetsToLL_M-50_HT-800to1200": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-800to1200"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "0.907",
        "sample_filters": "DYPhotonFilter"
    },
    "DYJetsToLL_M-50_HT-1200to2500": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-1200to2500"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "0.833",
        "sample_filters": "DYPhotonFilter"
    },
    "DYJetsToLL_M-50_HT-2500toInf": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-2500toInf"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.015",
        "sample_filters": "DYPhotonFilter"
    },
    "DYJetsToLL_M-10to50-LO": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-10to50-LO"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.0",
        "sample_filters": "DYPhotonFilter"
    },
    "DY_else": {
        "files": list_files_in_folder(mc_path, "DYJetsToLL_M-10to50-LO"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.0",
        "sample_filters": "DYPhotonFilter"
    },

    "TTToSemiLeptonic": {
        "files": list_files_in_folder(mc_path, "TTToSemiLeptonic"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "Top_pTrw",
        "sample_filters": None
    },
    "TTTo2L2Nu": {
        "files": list_files_in_folder(mc_path, "TTTo2L2Nu"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "Top_pTrw",
        "sample_filters": None
    },
    "TTWJetsToLNu": {
        "files": list_files_in_folder(mc_path, "TTWJetsToLNu"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None,
        "sample_filters": None
    },
    "ST_s-channel": {
        "files": list_files_in_folder(mc_path, "ST_s-channel"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None,
        "sample_filters": None
    },
    "ST_t-channel_antitop": {
        "files": list_files_in_folder(mc_path, "ST_t-channel_antitop"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "100. / 32.4",
        "sample_filters": None
    },
    "ST_t-channel_top": {
        "files": list_files_in_folder(mc_path, "ST_t-channel_top"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "100. / 32.4",
        "sample_filters": None
    },
    "ST_tW_antitop": {
        "files": list_files_in_folder(mc_path, "ST_tW_antitop"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None,
        "sample_filters": None
    },
    "ST_tW_top": {
        "files": list_files_in_folder(mc_path, "ST_tW_top"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None,
        "sample_filters": None
    },

    "WW": {
        "files": (
            list_files_in_folder(mc_path, "WmToLNu_WmTo2J_QCD") +
            list_files_in_folder(mc_path, "WpToLNu_WpTo2J_QCD")
        ),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "(!((mjjGenmax > 150) * GenLHE))", #fix it later
        "sample_filters": None
    },
    "ggH_bonly_on": {
        "files": list_files_in_folder(mc_path, "GluGluToWWToQQ_Cont_private"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "Lhe_mWW<160", #fix it later
        "sample_filters": None
    },
    "ggH_bonly_off": {
        "files": list_files_in_folder(mc_path, "GluGluToWWToQQ_Cont_private"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": True,
        "sample_weights": "Lhe_mWW > 160", #fix it later
        "sample_filters": None
    },
    "qqWWqq": {
        "files": (
            list_files_in_folder(mc_path, "WpTo2J_WmToLNu_QCD") +
            list_files_in_folder(mc_path, "WpToLNu_WmTo2J_QCD")
        ),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "(mjjGenmax > 150) * (GenLHE)", #added sample weight
        "sample_filters": None
    },
    "WJetsToLNu-LO": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu-LO"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "(LHE_HT < 70)", # Apply HT cut to W+jets LO sample as weight ? # check with Rajarshi
        "sample_filters": "WjetsPhotonFilter" # Apply WjetsPhotonFilter to W+jets LO sample
    },
    "WJetsToLNu_HT70To100": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu_HT70To100"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.21 * 0.95148",
        "sample_filters": "WjetsPhotonFilter"
    },
    "WJetsToLNu_HT100To200": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu_HT100To200"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "0.9471",
        "sample_filters": "WjetsPhotonFilter"
    },
    "WJetsToLNu_HT200To400": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu_HT200To400"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "0.9515",
        "sample_filters": "WjetsPhotonFilter"
    },
    "WJetsToLNu_HT400To600": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu_HT400To600"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "0.9581",
        "sample_filters": "WjetsPhotonFilter"
    },
    "WJetsToLNu_HT600To800": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu_HT600To800"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.0582",
        "sample_filters": "WjetsPhotonFilter"
    },
    "WJetsToLNu_HT800To1200": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu_HT800To1200"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.1285",
        "sample_filters": "WjetsPhotonFilter"
    },
    "WJetsToLNu_HT1200To2500": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu_HT1200To2500"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "1.3268",
        "sample_filters": "WjetsPhotonFilter"
    },
    "WJetsToLNu_HT2500ToInf": {
        "files": list_files_in_folder(mc_path, "WJetsToLNu_HT2500ToInf"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "2.7948",
        "sample_filters": "WjetsPhotonFilter"
    },
    "Vg": {
        "files": (
            list_files_in_folder(mc_path, "WGToLNuG") +
            list_files_in_folder(mc_path, "ZGToLLG")
        ),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "(!(Gen_ZGstar_mass > 0))",
        "sample_filters": None
    },
    "WGToLNuG": {
        "files": list_files_in_folder(mc_path, "WGToLNuG"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "VgWeight * (Gen_ZGstar_mass > 0 && Gen_ZGstar_mass < 0.1)",
        "sample_filters": None
    },
    "ZGToLLG": {
        "files": list_files_in_folder(mc_path, "ZGToLLG"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "VgWeight * (Gen_ZGstar_mass > 0)",  # *0.448 XS correction can be added if needed
        "sample_filters": None
    },
    "WZTo3LNu_mllmin0p1": {
        "files": list_files_in_folder(mc_path, "WZTo3LNu_mllmin0p1"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "VgWeight * (Gen_ZGstar_mass > 0.1)",
        "sample_filters": None
    },

    "ZZ": {
        "files": list_files_in_folder(mc_path, "ZZ"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None,
        "sample_filters": None
    },
    "WZ": {
        "files": list_files_in_folder(mc_path, "WZ"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "(mjjGenmax < 150)",  # Uses mjjGenmax defined in Analysis.py
        "sample_filters": None
    },
    "WmToLNu_ZTo2J_QCD": {
        "files": list_files_in_folder(mc_path, "WmToLNu_ZTo2J_QCD"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "(mjjGenmax > 150)",  # Uses mjjGenmax defined in Analysis.py
        "sample_filters": None
    },
    "WpToLNu_ZTo2J_QCD": {
        "files": list_files_in_folder(mc_path, "WpToLNu_ZTo2J_QCD"),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": "(mjjGenmax > 150)",  # Uses mjjGenmax defined in Analysis.py
        "sample_filters": None
    },
    "VVV": {
        "files": (
            list_files_in_folder(mc_path, "ZZZ") +
            list_files_in_folder(mc_path, "WZZ") +
            list_files_in_folder(mc_path, "WWZ") +
            list_files_in_folder(mc_path, "WWW")
        ),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None, #no sample weight for VVV
        "sample_filters": None
    },
    "WWewk": {
        "files": (
            list_files_in_folder(mc_path, "WpTo2J_WmToLNu") +
            list_files_in_folder(mc_path, "WpToLNu_WmTo2J") +
            list_files_in_folder(mc_path, "WpToLNu_WpTo2J") +
            list_files_in_folder(mc_path, "WmToLNu_WmTo2J")
        ),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None, #no sample weight for WWewk
        "sample_filters": None
    },
    "VBF_V": {
        "files": (
            list_files_in_folder(mc_path, "Wm_LNuJJ_EWK") +
            list_files_in_folder(mc_path, "Wp_LNuJJ_EWK")
        ),
        "isMC": True,
        "isSignal": False,
        "isOffshell": False,
        "sample_weights": None, #no sample weight for VBF_V
        "sample_filters": None
    },
    "ggH_sonly_on": {
        "files": list_files_in_folder(mc_path, "GluGluToWWToQQ_Sig_private"),
        "isMC": True,
        "isSignal": True,
        "isOffshell": False,
        "sample_weights": None,
        "sample_filters": "Lhe_mWW < 160"
    },
    "ggH_sonly_off": {
        "files": list_files_in_folder(mc_path, "GluGluToWWToQQ_Sig_private"),
        "isMC": True,
        "isSignal": True,
        "isOffshell": True,
        "sample_weights": None,
        "sample_filters": "Lhe_mWW > 160"
    },
    "ggH_sand_off": {
        "files": list_files_in_folder(mc_path, "GluGluToWWToQQ_SBI_private"),
        "isMC": True,
        "isSignal": True,
        "isOffshell": True,
        "sample_weights": None,
        "sample_filters": "Lhe_mWW > 160"
    },
    "ggH_sand_on": {
        "files": list_files_in_folder(mc_path, "GluGluToWWToQQ_SBI_private"),
        "isMC": True,
        "isSignal": True,
        "isOffshell": False,
        "sample_weights": None,
        "sample_filters": "Lhe_mWW < 160"
    }
}































