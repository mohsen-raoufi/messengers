% Title: Post-processing Plotting Analytical and simulation Messenger_ratio from DMP Result
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
close all;
% clear all;

% load Res_CheckStateMachine_expInit.mat
% load('Res_CheckStateMachine_expInit_2024.mat')
% load Res_CheckStateMachine_expInit_7_bug_fixed.mat
load Res_CheckStateMachine_expInit_2024_new.mat

bool_save_fig = false; true;

fig = figure(333);

hold off
plot(nan,nan)
fig.Color = 'w';
fig.Position = [336 188 782 661];

logp2e = log10(p2expltArr(1:end-1));
logp2m = log10(p2msngrArr);

[LogP2E, LogP2M] = meshgrid(logp2e,logp2m);
[P20Arr, P21Arr] = meshgrid(p2expltArr,p2msngrArr);
[P2E, P2M] = meshgrid(p2expltArr,p2msngrArr);

% mExpltTime = mean(expltDur(:,:,1:end-1,:),2);
% mmExpltTime = mean(mExpltTime,1); % abs(mExpltTime(1,:,:,:)-mExpltTime(2,:,:,:));%
% mmExpltTime = reshape(mmExpltTime,size(mmExpltTime,3),size(mmExpltTime,4));
% 
% mExpltTime = mean(expltDur(:,:,1:end-1,:),2);

mExpltNum = mean(numExpl(:,:,1:end-1,:),2);
mmExpltNum = mean(mExpltNum,1);
mmExpltNum = reshape(mmExpltNum,size(mmExpltNum,3),size(mmExpltNum,4));

hold on;

% surf(LogP2E, LogP2M, 1-mmExpltNum', 'EdgeAlpha',0.0)

line_width = 4.0;
edge_alpha = 0.5;

my_CMap = [];

analytic_m = (P2M./(P2E+P2M));

% col_analyt = [0, 0, 0]/255;
col_analyt = [0, 150, 0]/255;
lw_analyt = 1.5;

n_lines = 10;

v_e = linspace(0,0.48,n_lines);
for ii=1:length(v_e)
    % line_color = [1-v_e(ii), 0, v_e(ii)];
    line_color = [v_e(ii), 0, 1-v_e(ii)];
    my_CMap = [my_CMap; line_color];

    % line_color = 'k';

    % contour3(LogP2E, LogP2M, 1-(mmExpltNum'),[v_e(ii), v_e(ii)],':','EdgeColor',line_color,'LineWidth',line_width)

    % contour3(log10(P2E),log10(P2M), analytic_m , [v_e(ii), v_e(ii)], '-','EdgeColor',col_analyt,'LineWidth',lw_analyt);
    % using the same colormap
    contour3(log10(P2E),log10(P2M), analytic_m , [v_e(ii), v_e(ii)], '-','EdgeColor',line_color,'LineWidth',lw_analyt);

    % contour3(LogP2E, LogP2M, 1-(mmExpltNum'),[v_e(ii), v_e(ii)],'--','EdgeColor',line_color,'LineWidth',line_width,'EdgeAlpha',edge_alpha)
    

    mean_avg_num_expl = mean(avgNumExpl(:,1:end-1,:),1);
    mean_avg_num_expl = reshape(mean_avg_num_expl,size(mean_avg_num_expl,2),[]);
    contour3(LogP2E, LogP2M, 1-(mean_avg_num_expl'),[v_e(ii), v_e(ii)],'--','EdgeColor',line_color,'LineWidth',line_width,'EdgeAlpha',edge_alpha)


    % contour3(LogP2E, LogP2M, 1-(mmExpltNum'),[v_e(ii), v_e(ii)],'--','EdgeColor',line_color,'LineWidth',line_width,'EdgeAlpha',edge_alpha)
    
    tmp = reshape(mean(avgNumExpl(:,1:end-1,:),1),size(avgNumExpl,3),[]);
    contour3(LogP2E, LogP2M, 1-(tmp'),[v_e(ii), v_e(ii)],'--','EdgeColor',line_color,'LineWidth',line_width,'EdgeAlpha',edge_alpha)
    
end
% [C,H] = contour3(LogP2E, LogP2M, (mmExpltNum'),v_e,'b','LineWidth',1.5);




% % % % % % % [C,H] = contour3(LogP2E, LogP2M, (mmExpltNum'),10.^([-2:0.01:log10(0.49)]),'b');
% % % % % % % [C,H] = contour3(LogP2E, LogP2M, (mmExpltNum'),10.^(linspace(-2, log10(0.49), 15)),'b');
% % % % % % % [C,H] = contour3(LogP2E, LogP2M, (mmExpltNum'),10.^(linspace(-2, 0, 25)),'b');
% % % % % % % [C,H] = contour3(LogP2E, LogP2M, (mmExpltNum'),logspace(0,0.49,10),'b');


% hold on;
% H.ZData = H.ZData + 0.1;
% clabel(C,H,'Color','b');
% contour3(LogP2E, LogP2M, (mmExpltNum')+eps,[0.5 0.5],'--k','LineWidth',2);

line_color = 0.5*[1,0,1];
my_CMap = [my_CMap; line_color];

% line_color = 'w';
% contour3(LogP2E, LogP2M, 1-(mmExpltNum')+eps,[0.5 0.5],'--','color', line_color,'LineWidth',line_width);


v_m = linspace(0.52,1,n_lines);
for ii=1:length(v_m)
    % line_color = [1-v_m(ii), 0, v_m(ii)];
    line_color = [v_m(ii), 0, 1-v_m(ii)];
    my_CMap = [my_CMap; line_color];

    % line_color = 'k';

    % contour3(LogP2E, LogP2M, 1-(mmExpltNum'),[v_m(ii), v_m(ii)],':','EdgeColor',line_color,'LineWidth',line_width)


    % contour3(log10(P2E),log10(P2M), analytic_m , [v_m(ii), v_m(ii)], '-','EdgeColor',col_analyt,'LineWidth',lw_analyt);

    % using the same colormap
    contour3(log10(P2E),log10(P2M), analytic_m , [v_m(ii), v_m(ii)], '-','EdgeColor',line_color,'LineWidth',lw_analyt);


    contour3(LogP2E, LogP2M, 1-(mmExpltNum'),[v_m(ii), v_m(ii)],'--','EdgeColor',line_color,'LineWidth',line_width, 'EdgeAlpha',edge_alpha)

end
% [C,H] = contour3(LogP2E, LogP2M, (mmExpltNum'),v_m,'EdgeColor',line_colors','LineWidth',1.5);


% clabel(C,H,v_m([4, 14]),'Color','r');
% x0 = log(P20Arr); x1 = log(P21Arr);
% contour(x0, x1, exp(x0)./(exp(x0)+exp(x1)),'k')
% surf(log(P20Arr), log(P20Arr./(P20Arr+P21Arr)), mmExpltNum','EdgeAlpha',0.1);

% p_mid = plot(nan,nan,'--','color',0.5*[1,0,1],'LineWidth',line_width);
p_mid = plot(nan,nan,'--','color','w','LineWidth',line_width);
p_expl = plot(nan,nan,'-b','LineWidth',2);
p_mess = plot(nan,nan,'-r','LineWidth',2);


colormap(my_CMap)
% colormap(fliplr(my_CMap))
colorbar()
caxis([0 1]);
cb = colorbar('YTick',[0:0.2:1], 'FontSize', 24);
% ylabel(cb,'Messenger Population Ratio')
% ylabel(cb,'m')

% legend([p_expl, p_mid, p_mess], {'Exploiters in Majority', 'Equal Distribution', 'Messengers in Majority'}, ...
    % 'Location','best', 'FontSize',18)
% p_ = plot(nan,nan,'--k','LineWidth',2);

set(gca,'FontSize',24)

xlabel("$\log_{10} p_{E}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')
ylabel("$\log_{10} p_{M}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')

% legend([p_mid], {'Equal Distribution'}, 'Location','southeast', ...
%     'FontSize',22)

lgnd = legend([p_mid], {'\color{white} Equal Distribution'}, 'Location','southeast', ...
    'FontSize',22,'color','b','box','off');
% lgnd.Title.Color = 'w';

axis square

set(gca,'YTick',get(gca,'XTick'))
title('Expected Messenger Population Ratio (m)','Interpreter','latex','FontWeight','normal', 'fontsize',20) %,'Interpreter','latex')
cbar = colorbar;
cbar.Label.String = "$m$";
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 30;

file_strng = 'DMP_messenger_ratio_analytic_numeric_same_cmap_sq';
if(bool_save_fig)
    savefig(fig, [file_strng,'.fig'])
    % exportgraphics(fig,[file_strng,'.pdf'])
    exportgraphics(fig,[file_strng,'.png'],"Resolution",300)
end


% 
% ii = 1;
% set(findobjc(gca,'Type','patch','UserData',v_m(ii)),'EdgeColor',[0,0,0])

% print('DMP_messenger_ratio_contours_analytic_numeric_transparent_less_lines','-dpng','-r300')
% print('DMP_messenger_ratio_surf','-dpng')
