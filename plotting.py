import numpy as np
import hist, uproot
from array import array
import cmsstyle
import pickle
import matplotlib.pyplot as plt
import mplhep as hep
import os
import shutil

with open("my_histograms_Top_Step_1.pkl", "rb") as file:
    # Unpickle the object from the file
    histograms = pickle.load(file)

print(histograms)
top_histo = histograms["Top"]
histo_nlepton = top_histo["nLepton"]
integral = histo_nlepton.sum
print("Integral =",integral)
print(histo_nlepton)

eos_path = "/eos/user/r/rbhattac/www/Marc_TNP_Plots"
output_dir = "/H_WW/Debug_Plots/Top_Step_1"
output_path = eos_path + output_dir
os.makedirs(output_path,exist_ok=True)
shutil.copy("index.php", output_path)


hep.style.use("CMS")
fig, ax = plt.subplots()
hep.cms.label("Preliminary", loc=0, ax=ax);
hep.histplot(histo_nlepton, ax=ax, histtype='fill', label=f"Events = {integral}")
ax.set_ylabel("Events")
ax.set_xlabel("nLepton")
#ax.legend()


# Style
#plt.legend()
hep.cms.label();

plot_filename = os.path.join(output_path, "nLepton.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "nLepton.pdf")
plt.savefig(plot_filename)
plt.close(fig)

histo_LeptonPt_0_3D = top_histo["Lepton_pt_0"]
print(histo_LeptonPt_0_3D[:,hist.sum,hist.sum])
histo_LeptonPt_1_3D = top_histo["Lepton_pt_1"]
print(histo_LeptonPt_1_3D[:,hist.sum,hist.sum])

fig, ax1 = plt.subplots()
hep.cms.label("Preliminary", data = True, loc=0, ax=ax1);
ax1.set_ylabel("Events")
ax1.set_xlabel("Lepton_pt_0")
histoLeptonPt_0 = histo_LeptonPt_0_3D[:,hist.sum,hist.sum]
integral = histoLeptonPt_0.sum

hep.histplot(histoLeptonPt_0, ax=ax1, histtype='fill', label = f"Events = {integral}")
#ax1.legend()

plot_filename = os.path.join(output_path, "Lepton_0_Pt.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "Lepton_0_Pt.pdf")
plt.savefig(plot_filename)
plt.close(fig)

fig, ax2 = plt.subplots()
hep.cms.label("Preliminary", data = True, loc=0, ax=ax2);
ax2.set_ylabel("Events")
ax2.set_xlabel("Lepton_pt_1")

histoLeptonPt_1 = histo_LeptonPt_1_3D[:,hist.sum,hist.sum]
integral = histoLeptonPt_1.sum
hep.histplot(histoLeptonPt_1, ax=ax2, histtype='fill', label = f"Events = {integral}")
#ax2.legend()

plot_filename = os.path.join(output_path, "Lepton_1_Pt.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "Lepton_1_Pt.pdf")
plt.savefig(plot_filename)
plt.close(fig)


hist_LeptonPt_0_isLoose = histo_LeptonPt_0_3D[:, :, hist.sum]

fig, ax3 = plt.subplots()
hist_LeptonPt_0_isLoose.plot(ax=ax3, cbarextend=True, flow='none');
fig.get_axes()[-1].set_ylabel('Events', fontsize=20)
hep.cms.label(ax=ax3)

plot_filename = os.path.join(output_path, "hist_LeptonPt_0_isLoose.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "hist_LeptonPt_0_isLoose.pdf")
plt.savefig(plot_filename)
plt.close(fig)

hist_LeptonPt_0_isTight = histo_LeptonPt_0_3D[:, hist.sum, :]
fig, ax4 = plt.subplots()
hist_LeptonPt_0_isTight.plot(ax=ax4, cbarextend=True, flow='none');
fig.get_axes()[-1].set_ylabel('Events', fontsize=20)
hep.cms.label(ax=ax4)

plot_filename = os.path.join(output_path, "hist_LeptonPt_0_isTight.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "hist_LeptonPt_0_isTight.pdf")
plt.savefig(plot_filename)
plt.close(fig)


hist_LeptonPt_1_isLoose = histo_LeptonPt_1_3D[:, :, hist.sum]

fig, ax5 = plt.subplots()
hist_LeptonPt_1_isLoose.plot(ax=ax5, cbarextend=True, flow='none');
fig.get_axes()[-1].set_ylabel('Events', fontsize=20)
hep.cms.label(ax=ax5)

plot_filename = os.path.join(output_path, "hist_LeptonPt_1_isLoose.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "hist_LeptonPt_1_isLoose.pdf")
plt.savefig(plot_filename)
plt.close(fig)

hist_LeptonPt_1_isTight = histo_LeptonPt_1_3D[:, hist.sum, :]
fig, ax6 = plt.subplots()
hist_LeptonPt_1_isTight.plot(ax=ax6, cbarextend=True, flow='none');
fig.get_axes()[-1].set_ylabel('Events', fontsize=20)
hep.cms.label(ax=ax6)
plot_filename = os.path.join(output_path, "hist_LeptonPt_1_isTight.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "hist_LeptonPt_1_isTight.pdf")
plt.savefig(plot_filename)
plt.close(fig)

hist_LeptonPt_1_isLoose_Tight = histo_LeptonPt_1_3D[:, :, 1j]
fig, ax7 = plt.subplots()
hist_LeptonPt_1_isLoose_Tight.plot(ax=ax7, cbarextend=True, flow='none');
fig.get_axes()[-1].set_ylabel('Events', fontsize=20)
hep.cms.label(ax=ax7)
plot_filename = os.path.join(output_path, "hist_LeptonPt_1_isLoose_Tight.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "hist_LeptonPt_1_isLoose_Tight.pdf")
plt.savefig(plot_filename)
plt.close(fig)


hist_LeptonPt_1_isLoose_NoTight = histo_LeptonPt_1_3D[:, :, 0j]
fig, ax8 = plt.subplots()
hist_LeptonPt_1_isLoose_NoTight.plot(ax=ax8, cbarextend=True, flow='none');
fig.get_axes()[-1].set_ylabel('Events', fontsize=20)
hep.cms.label(ax=ax8)
plot_filename = os.path.join(output_path, "hist_LeptonPt_1_isLoose_NoTight.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "hist_LeptonPt_1_isLoose_NoTight.pdf")
plt.savefig(plot_filename)
plt.close(fig)


hist_LeptonPt_0_isLoose_NoTight = histo_LeptonPt_0_3D[:, :, 0j]
fig, ax9 = plt.subplots()
hist_LeptonPt_0_isLoose_NoTight.plot(ax=ax9, cbarextend=True, flow='none');
fig.get_axes()[-1].set_ylabel('Events', fontsize=20)
hep.cms.label(ax=ax9)
plot_filename = os.path.join(output_path, "hist_LeptonPt_0_isLoose_NoTight.png")
plt.savefig(plot_filename)
plot_filename = os.path.join(output_path, "hist_LeptonPt_0_isLoose_NoTight.pdf")
plt.savefig(plot_filename)
plt.close(fig)
