% Title: Plotting Analytical tau_c from DMP Result
% Description: This script processes and visualizes data of the DMP from
% checkStateMachine 
% Author: Mohsen Raoufi
% Contact: mohsenraoufi@icloud.com
% Affiliation: Research Cluster of Excellence, "Science of Intelligence"
% Date Created: May 27, 2024
% Version: 2.0
% Usage: Run this script with MATLAB R2023a or later. Run it in the directory that contains your data.
% License: Distributed under the MIT License. See LICENSE.txt for more information.
% Citation: please cite our work if you use this code 
% "Messengers: Breaking Echo Chambers in Collective Opinion Dynamics with Homophily"

clc; clear all; close all

load Res_CheckStateMachine_expInit.mat

bool_save_fig = true;

logp2e = log10(p2expltArr);
logp2m = log10(p2msngrArr);

[P2E, P2M] = meshgrid(p2expltArr,p2msngrArr);

fig888 = figure(888);

hold off
plot(nan,nan)
fig888.Color = 'w';
fig888.Position = [341 200 691 661];

% [C,H] = contour(LogP2E,LogP2M, log10(1./(1./(10.^LogP2E + 10.^LogP2M))),'LineWidth',3);
% [C,H] = contour(LogP2E,LogP2M, log10((1./(10.^LogP2E + 10.^LogP2M))),'LineWidth',3);


% Tau_s : the sojourn time
% [C,H] = contour(log10(P2E),log10(P2M), log10(0.5*(P2E + P2M)./(P2E.*P2M)),'LineWidth',3);

% [C,H] = contour(log10(P2E),log10(P2M), log10(2*0.5*(P2M+P2E)./(P2E.*P2M)),'LineWidth',3);

% Tau_c : the relaxation time
[C,H] = contour(log10(P2E),log10(P2M), log10((1./(P2E + P2M))),'LineWidth',3);

% colormap hot
colormap copper
% colormap abyss
% colormap autumn
% colormap(flipud(copper))

set(gca,'FontSize',24)

xlabel("$\log_{10} p_{E}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')
ylabel("$\log_{10} p_{M}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')

axis square

clabel(C,H,'fontSize',18,'labelspacing', 300); % 'Color','k');
% title("$\log_{10}\tau_{c}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')

title("Relaxation Time $\log_{10}\tau_{c}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')

colorbar
% print('waiting_time_P2E_P2M_abyss','-dpng')
% print('waiting_time_P2E_P2M_copper','-dpng')

% exportgraphics(fig888, "DMP_relaxation_time_P2E_P2M.png","Resolution",300)

