% Title: Post-processing Plotting Analytical and simulation tau_s (sojourn time) from DMP Result
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

clc;
% clear all;
close all;

% load Res_CheckStateMachine_expInit.mat
% load('Res_CheckStateMachine_expInit_2024.mat')

bool_save_fig = true;

logp2e = log10(p2expltArr(1:end-1)); % to make it square
logp2m = log10(p2msngrArr);

[LogP2E, LogP2M] = meshgrid(logp2e,logp2m);
[P20Arr, P21Arr] = meshgrid(p2expltArr,p2msngrArr);
[P2E, P2M] = meshgrid(p2expltArr,p2msngrArr);

mExpltTime = mean(expltDur(:,:,1:end-1,:),2);
mmExpltTime = mean(mExpltTime,1); 
mmExpltTime = reshape(mmExpltTime,size(mmExpltTime,3),size(mmExpltTime,4));

mExpltNum = mean(numExpl(:,:,1:end-1,:),2);
mmExpltNum = mean(mExpltNum,1);
mmExpltNum = reshape(mmExpltNum,size(mmExpltNum,3),size(mmExpltNum,4));

mNumEdge = mean(numEdges(:,:,1:end-1,:),2);
mmNumEdge = mean(mNumEdge,1);
mmNumEdge = reshape(mmNumEdge,size(mmNumEdge,3),size(mmNumEdge,4));

mean_explt_waiting_time = (tf*(expltDur) + tf*(1-expltDur))./(numEdges+1);
mean_explt_waiting_time = mean_explt_waiting_time(:,:,1:end-1,:);
mMean_explt_wt = mean(mean_explt_waiting_time,2);
mmMean_explt_wt = mean(mMean_explt_wt,1);
mmMean_explt_wt = reshape(mmMean_explt_wt,size(mmMean_explt_wt,3),size(mmMean_explt_wt,4));

fig = figure(888);

hold off
plot(nan,nan)
fig.Color = 'w';
fig.Position = [336 188 782 661];

% color analytical
col_analyt = 0*[255, 51, 51]/255;
col_analyt = [0, 150, 0]/255;

% Analytical Results
% contour(log10(P2E),log10(P2M), log10(0.5*(P2M+P2E)./(P2E.*P2M)),'-','color',col_analyt,'LineWidth',1.5);
contour(log10(P2E),log10(P2M), log10(0.5*(P2M+P2E)./(P2E.*P2M)),'-','LineWidth',1.5);

hold on

% surf(LogP2E, LogP2M, log10(mmMean_explt_wt)','EdgeAlpha',0.1);
% [C,H] = contour(LogP2E, LogP2M, log10(mmMean_explt_wt)','b');
[C,H] = contour(LogP2E,LogP2M, log10(mmMean_explt_wt)','--','LineWidth',6,'EdgeAlpha',0.6);

hold on

% [C,H] = contour(LogP2E,LogP2M, log10(1./(1./(10.^LogP2E + 10.^LogP2M))),'LineWidth',3);
% [C,H] = contour(LogP2E,LogP2M, log10((1./(10.^LogP2E + 10.^LogP2M))),'LineWidth',3);
% [C,H] = contour(LogP2E,LogP2M, log10((1./(10.^LogP2E + 10.^LogP2M))),':b','LineWidth',0.2);


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
% title("$\log_{10}\tau_{S}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')
title("Sojourn Time",'Interpreter','latex','FontSize',22,'FontWeight','bold')
cbar = colorbar;
cbar.Label.String = "$\log_{10}\tau_{S}$";
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 30;

set(gca,'YTick',get(gca,'XTick'))

file_strng = 'DMP_sojourn_time_analytic_numeric_same_cmap_sq';
if(bool_save_fig)
    savefig(fig, [file_strng,'.fig'])
    % exportgraphics(fig,[file_strng,'.pdf'])
    exportgraphics(fig,[file_strng,'.png'],"Resolution",300)
end

% 
% clabel(C,H,'fontSize',12,'labelspacing', 200); 
% % surf(1./(P20Arr), 1./(P21Arr), mmExpltTime./(mmNumEdge+1),'EdgeAlpha',0.1);
% xlabel('log P_{2E}');
% ylabel('log P_{2M}');
% view(0,90);
% axis square;
% colorbar
% axis([min(logp2e) max(logp2e) min(logp2m) max(logp2m)])
% % colormap jet
% title('mean exploit duration');