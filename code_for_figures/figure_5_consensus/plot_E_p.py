'''

% Title: Figure_5_Plot_E_P_final
% Description: Plots the final precision error for different configuraiton of the
% DMP, the data is shown in the Figure 5 of the paper
% Author: Mohsen Raoufi
% Contact: mohsenraoufi@icloud.com
% Affiliation: Research Cluster of Excellence, "Science of Intelligence"
% Date Created: July, 2024
% Version: 1.0
% Usage: Run this script with MATLAB R2023a or later. Run it in the directory that contains your data.
% License: Distributed under the MIT License. See LICENSE.txt for more information.
% Citation: please cite our work if you use this code
% "Messengers: Breaking Echo Chambers in Collective Opinion Dynamics with Homophily"
% NOTE: We recommend you to visualize the data using the MATLAB code provided, otherwise, the visualization might not be as expected.

'''

import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import cm

# Load the .mat files
folder_str = 'data/'
name_str = 'MC_1_Arena__50k_tf_cone_BasicMarkov_initENumMsngr_LargeParReg'

mat_data = scipy.io.loadmat(folder_str + name_str + '.mat')
mat_data_ss = scipy.io.loadmat(folder_str + name_str + '_SS.mat')

# Access the required data from the .mat files
p2expltArr = mat_data['p2expltArr'].flatten()
p2msngrArr = mat_data['p2msngrArr'].flatten()
E_P_p_SS_MC = mat_data_ss['E_P_p_SS_MC']

# Save figure in the end?
saveBool = False

# Initialize figure
fig = plt.figure(figsize=(12, 10))

# use colormap
plt.set_cmap('pink')

p2e = np.log10(p2expltArr[:-1])
p2m = np.log10(p2msngrArr)

P2M, P2E = np.meshgrid(p2e, p2m)

var2Plot = E_P_p_SS_MC / E_P_p_SS_MC[0, -1]  # Normalize the final precision error w.r.t the initial
var2Plot = var2Plot[:-1, :]

# Plot surface
surf = plt.contourf(P2M, P2E, var2Plot.T, cmap=cm.hot, alpha=1.0)
plt.colorbar(surf, label='E_P^O')

# Smooth the data (ONLY FOR CONTOURS!!)
from scipy.ndimage import gaussian_filter
# var2Plot_smoothed = var2Plot #  if the following gaussian filter does not work properly
var2Plot_smoothed = gaussian_filter(var2Plot, sigma=1)

# Add contour lines
lw_1 = 2
lw_2 = 3
plt.contour(P2M, P2E, var2Plot_smoothed.T + 0.1, levels=[0.1 + 0.275], colors='w', linewidths=lw_1)
plt.contour(P2M, P2E, var2Plot_smoothed.T + 0.1, levels=[0.1 + 0.25], colors='g', linewidths=lw_1)
plt.contour(P2M, P2E, var2Plot_smoothed.T + 0.1, levels=[0.1 + 0.12], colors='y', linewidths=lw_1)
plt.contour(P2M, P2E, var2Plot_smoothed.T + 0.2, levels=[0.2 + 0.65], colors='b', linewidths=lw_2)
plt.contour(P2M, P2E, var2Plot_smoothed.T + 0.2, levels=[0.2 + 0.65], colors='r', linewidths=lw_2)

# Update colors of the contours (R1, R2, R3 regions)
delta = 0.1
plt.plot(p2e, p2e, '--', linewidth=lw_1, color="#64E3FA")

# Add text labels for regions
txt_fontsz = 40
txt_fontsz_2 = 30
plt.text(-1.8, -7, "R1", color='b', fontsize=txt_fontsz, fontweight='bold')
plt.text(-3.5, -2.8, "R2", color='w', fontsize=txt_fontsz, fontweight='bold')
plt.text(-7.8, -1.2, "R3", color='r', fontsize=txt_fontsz, fontweight='bold')
plt.text(-1.66, -2.26, "R2-b", color='g', fontsize=txt_fontsz_2, fontweight='bold')
plt.text(-7.81, -8.1, "R2-a", color=[0.6]*3, fontsize=txt_fontsz_2, fontweight='bold')
plt.text(-4.52, -1.06, "R2-c", color='y', fontsize=txt_fontsz_2, fontweight='bold')

plt.xlabel("$\log_{10} p_{E}$", fontsize=30, fontweight='bold')
plt.ylabel("$\log_{10} p_{M}$", fontsize=30, fontweight='bold')
plt.xticks(np.arange(-8, 0))
plt.yticks(np.arange(-8, 0))
plt.gca().set_aspect('equal', adjustable='box')
plt.gca().tick_params(labelsize=22)
plt.grid(True)


plt.title("Final Precision Error Visualization", fontsize=22, fontweight='bold')

if saveBool:
    plt.savefig("E_p_o_contours.png", dpi=300)
plt.show()