% Title: Figure_7_Plot_E_s_final
% Description: Plots the final precision error (of opinion) for different configuraiton of the
% DMP, the data is shown in the Figure 7 of the paper
% Author: Mohsen Raoufi
% Contact: mohsenraoufi@icloud.com
% Affiliation: Research Cluster of Excellence, "Science of Intelligence"
% Date Created: July, 2024
% Version: 1.0
% Usage: Run this script with MATLAB R2023a or later. Run it in the directory that contains your data.
% License: Distributed under the MIT License. See LICENSE.txt for more information.
% Citation: please cite our work if you use this code
% "Messengers: Breaking Echo Chambers in Collective Opinion Dynamics with Homophily"

clc;
clear all;
close all;

folder_str = 'data/';
% folder_str = '/Users/mohsenraoufi/Project/echo_chambers/echo-chambers/MesExp_2D/MATLAB/for_paper/figure_7_init_M/data/';

name_str = 'MC_1_Arena__50k_tf_cone_BasicMarkov_initENumMsngr_LargeParReg';
name_str = 'MCs_13_Jan_2022_PC_1.0x1.0Arena__50k_tf__lowRand4Explt_cone_BasicMarkov_initAllExploiters_LargeParReg';

load([folder_str, name_str, '.mat']);
load([folder_str, name_str, '_SS.mat']);

% save figure in the end?
saveBool = false;

% not imp: used only for naming the saved figures
init_expected = false; true;

close all

fig = figure(1);
fig.Position = [719 176 866 704];

p2e = log10(p2expltArr(1:end-1));
p2m = log10(p2msngrArr);

[P2M, P2E] = meshgrid(p2e,p2m);

var2Plot = E_P_s_SS_MC./E_P_s_SS_MC(1,end); % normalize the final precision error w.r.t the baseline
var2Plot = var2Plot(1:end-1,:);

SH = surf(P2M, P2E, var2Plot','EdgeAlpha',0.,'FaceAlpha',1.0);
hold on;

hold on;
axis([min(p2e) max(p2e) min(p2m) max(p2m)])
view(0,90);
colormap('hot');
axis square

set(gca,'FontSize',22)
%
cbar = colorbar;

cbar.Label.set('String',"E_P^S") % ,'Interpreter','latex','FontSize',30,'FontWeight','bold')
cbar.Label.set("FontSize",30)
cbar.Label.Rotation = 0;

xlabel("$\log_{10} p_{E}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')
ylabel("$\log_{10} p_{M}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')

set(gcf,'color','w')
xticks(-8:-1)
set(gca,'YTick',get(gca,'XTick'))

cmap = colormap('Pink');
colormap(cmap);

f = gcf;

file_strng = "E_p_o_contours";

if(saveBool)
    savefig(file_strng+'.fig')
    exportgraphics(f,file_strng+'_v_sq.pdf','ContentType','vector')
    exportgraphics(f,file_strng+'_sq.pdf')
    exportgraphics(f,file_strng+"_sq.png","Resolution",300)
end
