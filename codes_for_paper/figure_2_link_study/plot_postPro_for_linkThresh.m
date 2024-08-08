% Title: Plotting Figure 2 of the paper (number of clusters and Precision
% error vs communication range)
% Description: This script processes and visualizes data of the simulations
% for no Messenger (baseline) setup
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
close all;
% clear all

% %% LOAD the following file that contains the relevant data
% load data/MCs_27_May_2024__NPop_100_Arena_1__tf_20k__lowRand4Explt__cone__BasicMarkov__NoMessenger_sensRang_Study_ONLY_INITIAL__nMC_8.mat
% load data/MCs_29_May_2024__NPop_100_Arena_1__tf_50k__cone__BasicMarkov__NoMessenger_sensRang_Study__nMC_40.mat
load data/MCs_31_May_2024__NPop_100_Arena_1__tf_50k__cone__BasicMarkov__NoMessenger_sensRang_Study__nMC_40.mat
save_fig_bool = false;

plot_shade_bool = false;

fig = figure(1);
fig.Position = [1063 495 800 400];
hold on
axis([0,1,0,1])

if(plot_shade_bool)
    start_shade_x = graph_data.commRangeArr(1:end-1);
    finish_shade_x = graph_data.commRangeArr(2:end);
    giant_comp_arr = graph_data.giantComp;
    shade(start_shade_x, finish_shade_x, mean(graph_data.giantComp(1:end-1,:),2)'/100)
    % shade(start_shade_x, finish_shade_x, min(graph_data.giantComp(1:end-1,:),[],2)'/100)
end

% %% %%% NUMBER OF CLUSTERS

% %% Initial clusters
% plot_color = [0, 98, 255]/255;
% plot_shaded_with_std(nClustMat_init/100, linkThreshArr, 2, plot_color, true)
% p1 = plot(linkThreshArr,mean(nClustMat_init,2)/100,'--','color',plot_color,'LineWidth',2);
% 
% scatter(linkThreshArr,reshape(nClustMat_init,length(linkThreshArr),[])/100,'MarkerFaceColor',plot_color,'MarkerEdgeColor',plot_color,'SizeData',0.2)
% 
% % %% Final clusters
% plot_color = [0, 70, 180]/255;
% plot_shaded_with_std(nClustMat_final/100, linkThreshArr, 2, plot_color, true)
% p2 = plot(linkThreshArr,mean(nClustMat_final,2)/100,'-','color',plot_color,'LineWidth',2);
% 
% scatter(linkThreshArr,reshape(nClustMat_final/100,length(linkThreshArr),[]),'MarkerFaceColor',plot_color,'MarkerEdgeColor',plot_color,'SizeData',0.2)
% 
% yyaxis left
% set(gca,'color','w')
% set(gcf,'color','w')
% xlabel('Communication Range')
% ylabel('Number of Clusters/N')
% set(gca,'YColor',plot_color)
% ylim([0, 1])


yyaxis right
E_p_p_init = std(zpArr,0,2);
E_p_p = std(zpArr,1,2);
avg_ep_f = mean(E_p_p_init(:,1,end,:),4);
plot_color = [194, 39, 0]/255;
plot_shaded_with_std(E_p_p_init(:,1,end,:), linkThreshArr, 4, plot_color, true)
plot(linkThreshArr,avg_ep_f,'--','color',plot_color,'LineWidth',2);
scatter(linkThreshArr,reshape(E_p_p_init(:,1,1,:),length(linkThreshArr),[]),'MarkerFaceColor',plot_color,'MarkerEdgeColor',plot_color,'SizeData',0.2)


plot_color = [140, 30, 0]/255;
plot_shaded_with_std(E_p_p(:,1,1,:), linkThreshArr, 4, plot_color, true)
plot(linkThreshArr,mean(E_p_p(:,1,1,:),4),'-','color',plot_color,'LineWidth',2);
scatter(linkThreshArr,reshape(E_p_p(:,1,1,:),length(linkThreshArr),[]),'MarkerFaceColor',plot_color,'MarkerEdgeColor',plot_color,'SizeData',0.2)

ylabel('Precision Error')
set(gca,'YColor',plot_color)

plot_color = [130, 130, 130]/255;
p1 = plot(nan,nan,'--', 'color',plot_color,'LineWidth',2);
plot_color = [30, 30, 30]/255;
p2 = plot(nan,nan,'-', 'color',plot_color,'LineWidth',2);
lgnd = legend([p1;p2],{'initial';'final'},'FontSize',14,'Location','east','Box','on');


set(gca,'fontSize', 14)
ax = gca;
strng = "not_to_push_data/";




if(save_fig_bool)
    exportgraphics(ax,[strng+'__nClust_vs_com_Range__init_final.png'],'Resolution',300)
    exportgraphics(ax,[strng+'__nClust_vs_com_Range__init_final.pdf'])
    savefig([strng+'__nClust_vs_com_Range__init_final'])
end

% exportgraphics(ax,[strng+'__clustOrd_vs_com_Range__init_final.png'],'Resolution',300)
% exportgraphics(ax,[strng+'__clustOrd_vs_com_Range__init_final.pdf'])
% savefig([strng+'__clustOrd_vs_com_Range__init_final'])


% print([strng+'__nClust_vs_com_Range__init_final'],'-dpng','-r300')
% print([strng+'__nClust_vs_com_Range__init_final'],'-dpdf')
% savefig([strng+'__nClust_vs_com_Range__init_final'])
% print([strng+'__nClust_vs_com_Range__init_final'],'-dmeta')


%
%
% nPlots = size(nClustMat,2);
% for iTime=1:nPlots
%     %     plot(linkThreshArr,mean(nClustMat(:,iTime,:),3),'color',[0.1,0.1,(iTime/nPlots)^0.2],'LineWidth',2)
%     hold on
%         plot(iTime/nPlots*tf,reshape(mean(nClustMat(:,iTime,:)./nClustMat(:,end,:),3),1,[]),'.b','markerSize',0.1)%,'SizeData',1)
% %     scatter(iTime/nPlots*tf,reshape(nClustMat(11,iTime,:)./nClustMat(11,end,:),1,[]),'b','SizeData',1)
% %     plot(iTime/nPlots*tf,reshape(nClustMat(11,iTime,:),1,[]),'.b')%,'SizeData',1)
%     drawnow
%     %     pause(0.1)
%     %     hold off
% end
% set(gca,'color','w')
% set(gcf,'color','w')
% xlabel('Time')
% ylabel('# Clusters/Final # Clusters')


function plot_shaded_with_std(data, x_data, nth_dim_for_std, color, constraint_to_zero)
avg_data = mean(data,nth_dim_for_std);
std_data = std(data,[],nth_dim_for_std);
low_bound_data = avg_data - std_data;
hi_bound_data = avg_data + std_data;

if(constraint_to_zero)
    low_bound_data = max(low_bound_data, 0);
end

% make sure it is vertical!
low_bound_data = reshape(low_bound_data, [], 1);
hi_bound_data = reshape(hi_bound_data,[], 1);

fill([x_data, fliplr(x_data)], [low_bound_data', fliplr(hi_bound_data')], color, 'FaceAlpha',0.2, 'EdgeColor','none')
end


function []=shade(Start,Finish,color_variable)
curax=axis;
y=[curax(3) curax(4) curax(4) curax(3)];
for i=1:length(Start)
    x=[Start(i) Start(i) Finish(i) Finish(i)];
    % color_rgb = (1-color_variable(i))*[1,1,1];
    color_rgb = [1,1,(color_variable(i))];
    h=fill(x,y,color_rgb,'EdgeAlpha',0);
    set(h,'facealpha',.5)
end
end

