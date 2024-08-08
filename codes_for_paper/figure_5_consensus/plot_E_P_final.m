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

clc;
clear all;
close all;

folder_str = 'data/';

name_str = 'MC_1_Arena__50k_tf_cone_BasicMarkov_initENumMsngr_LargeParReg';

load([folder_str, name_str, '.mat']);
load([folder_str, name_str, '_SS.mat']);

% save figure in the end?
saveBool = false;

% not imp: used only for naming the saved figures
init_expected = true;

close all

fig = figure(1);
fig.Position = [719 176 866 704];

p2e = log10(p2expltArr(1:end-1));
p2m = log10(p2msngrArr);

[P2M, P2E] = meshgrid(p2e,p2m);

var2Plot = E_P_p_SS_MC./E_P_p_SS_MC(1,end); % normalize the final precision error w.r.t the baseline
var2Plot = var2Plot(1:end-1,:);

SH = surf(P2M, P2E, var2Plot','EdgeAlpha',0.,'FaceAlpha',1.0);
hold on;

lw_1 = 2; % line width config
lw_2 = 3; % line width config

% Smooth the data (ONLY FOR CONTOURS!!)
var2Plot = smoothdata(var2Plot,'gaussian',4);
[Qc1_k,Hc1_k] = contour3(P2M, P2E, var2Plot'+0.1,0.1+0.275*ones(1,2) ,'w','LineWidth',lw_1);
[Qc1_g,Hc1_g] = contour3(P2M, P2E, var2Plot'+0.1,0.1+0.25*ones(1,2) ,'g','LineWidth',lw_1);

[Qc2_y,Hc2_y] = contour3(P2M, P2E, var2Plot'+0.1,0.1+0.12*ones(1,2) ,'y','LineWidth',lw_1);
gray_col = 0.6*[1,1,1];

[Qc3_b,Hc3_b] = contour3(P2M, P2E, var2Plot'+0.2,0.2+0.65*ones(1,2) ,'b','LineWidth',lw_2);
[Qc3_r,Hc3_r] = contour3(P2M, P2E, var2Plot'+0.2,0.2+0.65*ones(1,2) ,'r','LineWidth',lw_2);

% % %% update colours of the contours

% %% R1
delta = 0.1;
Hc3_b.XData(Hc3_b.XData<Hc3_b.YData + delta) = nan;
Hc3_b.YData(Hc3_b.XData<Hc3_b.YData + delta) = nan;

% %% R3
delta = 0.2;
Hc3_r.XData(Hc3_r.XData>Hc3_r.YData + delta) = nan;
Hc3_r.YData(Hc3_r.XData>Hc3_r.YData + delta) = nan;
Hc3_r.XData(Hc3_r.XData>-1) = nan;
Hc3_r.YData(Hc3_r.XData>-1) = nan;

% %% R2

% %% %% R2-a
delta = -0.55;
Hc1_g.XData(Hc1_g.XData<Hc1_g.YData + delta) = nan;
Hc1_g.YData(Hc1_g.XData<Hc1_g.YData + delta) = nan;
Hc1_g.XData(Hc1_g.XData<-2.3) = nan;
Hc1_g.YData(Hc1_g.XData<-2.3) = nan;

% %% %% R2-c
Hc2_y.XData(Hc2_y.YData<-2) = nan;
Hc2_y.YData(Hc2_y.YData<-2) = nan;
Hc2_y.XData(Hc2_y.XData>Hc2_y.YData) = nan;
Hc2_y.YData(Hc2_y.XData>Hc2_y.YData) = nan;

% draw the off-diagonal line
off_diag_x = p2e;
off_diag_y = p2e;
plot3(off_diag_x, off_diag_y, 10*ones(size(off_diag_y)),'--','LineWidth',lw_1,'Color',"#64E3FA")


txt_fontsz = 40;
txt_fontsz_2 = 30;
% eps = 0.0;
txt_R1 = text(-1.8,-7,10,"R1",'Color','b','FontSize',txt_fontsz,'FontWeight','bold');
txt_R2 = text(-3.5,-2.8,10,"R2",'Color','w','FontSize',txt_fontsz,'FontWeight','bold');
txt_R3 = text(-7.8,-1.2,10,"R3",'Color','r','FontSize',txt_fontsz,'FontWeight','bold');

txt_R2_b = text(-1.66,-2.26,10,"R2-b",'Color','g','FontSize',txt_fontsz_2,'FontWeight','bold');
txt_R2_a = text(-7.81,-8.1,10,"R2-a",'Color',gray_col*1.5,'FontSize',txt_fontsz_2,'FontWeight','bold');
txt_R2_c = text(-4.52,-1.06,10,"R2-c",'Color','y','FontSize',txt_fontsz_2,'FontWeight','bold');

hold on;
axis([min(p2e) max(p2e) min(p2m) max(p2m)])
view(0,90);
colormap('hot');
axis square

set(gca,'FontSize',22)
%
cbar = colorbar;

cbar.Label.set('String',"E_P^O") % ,'Interpreter','latex','FontSize',30,'FontWeight','bold')
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
