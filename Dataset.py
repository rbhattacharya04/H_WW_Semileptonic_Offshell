import os 

def list_files_in_folder(folder_path, search_string=None):
    """
    Lists files in a given folder, optionally filtering by a string in the filename.

    Args:
        folder_path (str): The path to the folder to search.
        search_string (str, optional): A string to search for in the filenames. 
                                       If None, all files are returned. Defaults to None.

    Returns:
        list: A list of filenames (including their full paths) that match the criteria.
    """
    found_files = []
    try:
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            if os.path.isfile(file_path):  # Ensure it's a file, not a directory
                if search_string is None or search_string in filename:
                    found_files.append(file_path)
    except FileNotFoundError:
        print(f"Error: Folder not found at '{folder_path}'")
    return found_files


#mc_path = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer20UL18_106x_nAODv9_Full2018v9/MCl1loose2018v9"
mc_path = "/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer20UL18_106x_nAODv9_Full2018v9/MCl1loose2018v9__MCCorr2018v9NoJERInHorn__MCCombJJLNu2018"


dataset = { 
"data" : list_files_in_folder("/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Run2018_UL2018_nAODv9_Full2018v9/DATAl1loose2018v9__DATACombJJLNu2018/"),
#"DY" : list_files_in_folder(mc_path, "DYJetsToLL_M-50,"),
"DY" : list_files_in_folder(mc_path, "DYJetsToLL_M-50") + list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-70to100") +list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-100to200") + list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-200to400") + list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-400to600") + list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-600to800") + list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-800to1200") + list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-1200to2500") + list_files_in_folder(mc_path, "DYJetsToLL_M-50_HT-2500toInf"),
"DY_else" : list_files_in_folder(mc_path, "DYJetsToLL_M-10to50-LO"),
"Top" : list_files_in_folder(mc_path, "TTToSemiLeptonic") + list_files_in_folder(mc_path, 'TTTo2L2Nu') + list_files_in_folder(mc_path,'TTWJetsToLNu') + list_files_in_folder(mc_path,'TTWjets') + list_files_in_folder(mc_path,'ST_s-channel') + list_files_in_folder(mc_path,'ST_t-channel_antitop') + list_files_in_folder(mc_path,'ST_t-channel_top') + list_files_in_folder(mc_path, 'ST_tW_antitop') + list_files_in_folder(mc_path,'ST_tW_top'),
"WW" : list_files_in_folder(mc_path,'WmToLNu_WmTo2J_QCD') + list_files_in_folder(mc_path,'WpToLNu_WpTo2J_QCD'),
"ggH_bonly_on": list_files_in_folder(mc_path,'GluGluToWWToQQ_Cont_private'),
"ggH_bonly_off": list_files_in_folder(mc_path,'GluGluToWWToQQ_Cont_private'),
"qqWWqq": list_files_in_folder(mc_path,'WpTo2J_WmToLNu_QCD') + list_files_in_folder(mc_path,'WpToLNu_WmTo2J_QCD'),
"W + jets":list_files_in_folder(mc_path,'WJetsToLNu-LO') + list_files_in_folder(mc_path,'WJetsToLNu_HT70To100') + list_files_in_folder(mc_path,'WJetsToLNu_HT100To200') + list_files_in_folder(mc_path,'WJetsToLNu_HT200To400') + list_files_in_folder(mc_path,'WJetsToLNu_HT400To600') + list_files_in_folder(mc_path,'WJetsToLNu_HT600To800') + list_files_in_folder(mc_path,'WJetsToLNu_HT800To1200') + list_files_in_folder(mc_path,'WJetsToLNu_HT1200To2500') + list_files_in_folder(mc_path,'WJetsToLNu_HT2500ToInf'),
"Vg": list_files_in_folder(mc_path,'WGToLNuG') + list_files_in_folder(mc_path,'ZGToLLG'),
"VgS": list_files_in_folder(mc_path,'WGToLNuG') + list_files_in_folder(mc_path,'ZGToLLG') + list_files_in_folder(mc_path,'WZTo3LNu_mllmin0p1'),
"VZ": list_files_in_folder(mc_path, 'ZZ') + list_files_in_folder(mc_path,'WZ') + list_files_in_folder(mc_path,'WmToLNu_ZTo2J_QCD') + list_files_in_folder(mc_path,'WpToLNu_ZTo2J_QCD'),
"VVV": list_files_in_folder(mc_path,'ZZZ') + list_files_in_folder(mc_path,'WZZ') + list_files_in_folder(mc_path,'WWZ') + list_files_in_folder(mc_path,'WWW'),
"WWewk": list_files_in_folder(mc_path,'WpTo2J_WmToLNu') + list_files_in_folder(mc_path,'WpToLNu_WmTo2J') + list_files_in_folder(mc_path,'WpToLNu_WpTo2J') + list_files_in_folder(mc_path,'WmToLNu_WmTo2J'),
"VBF_V": list_files_in_folder(mc_path,'Wm_LNuJJ_EWK') + list_files_in_folder(mc_path,'Wp_LNuJJ_EWK'),
"ggH_sonly_on": list_files_in_folder(mc_path,'GluGluToWWToQQ_Sig_private'),
"ggH_sonly_off" : list_files_in_folder(mc_path,'GluGluToWWToQQ_Sig_private'),
"ggH_sand_off" : list_files_in_folder(mc_path,'GluGluToWWToQQ_SBI_private'),
"ggH_sand_on" : list_files_in_folder(mc_path,'GluGluToWWToQQ_SBI_private'),
}


