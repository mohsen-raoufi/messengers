clc;
% clear all;
close all;
%
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_1__i_p2m_49.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.1__i_p2e_8__i_p2m_12.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_100k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_43.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_100k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_43_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_100k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_100k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_100k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_46.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_100k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_46_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat_final_pos.png
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_43.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_43_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_3__i_p2m_4.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_3__i_p2m_45.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_3__i_p2m_45_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_3__i_p2m_4_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_41__i_p2m_48.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_41__i_p2m_48.mat_final_pos.png
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_41__i_p2m_48_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_45__i_p2m_3.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_45__i_p2m_3_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_49__i_p2m_1.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_49__i_p2m_1.mat_final_pos.png
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_49__i_p2m_1_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_4__i_p2m_3.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_4__i_p2m_3_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_6__i_p2m_7.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_6__i_p2m_7_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_8__i_p2m_12.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_8__i_p2m_12_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_9__i_p2m_10.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_9__i_p2m_10_Params.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.1__i_p2e_8__i_p2m_12.mat
% % % single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.1__i_p2e_8__i_p2m_12_Params.mat
% % single_run_01_Nov_2023__NPop_200_Arena_1__tf_10k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_41__i_p2m_48.mat


% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_15__i_p2m_2.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_21__i_p2m_30.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_49__i_p2m_1.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_1__i_p2m_49.mat

% load not_to_push_data/single_run_01_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_6__i_p2m_7.mat
% load not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_7__i_p2m_12.mat
% % single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_8__i_p2m_12.mat
% % single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_9__i_p2m_10.mat
% % single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.1__i_p2e_8__i_p2m_12.mat
% % single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_1__i_p2m_49.mat

% % single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_2__i_p2m_16.mat



% strng = 'not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat';
% do_the_job(strng)
% strng = 'not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_6__i_p2m_7.mat';
% do_the_job(strng)
% strng = 'not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_7__i_p2m_12.mat';
% do_the_job(strng)
% strng = 'not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_21__i_p2m_30.mat';
% do_the_job(strng)
strng = 'not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat';
do_the_job(strng)
% strng = 'not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_15__i_p2m_2.mat';
% do_the_job(strng)
% strng = 'not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_49__i_p2m_1.mat';
% do_the_job(strng)


function do_the_job(strng)
close all

path_string = '../../codes_for_paper/';

load([path_string,strng]);

height = 0.3;
height_narrow = 0.1;
y_lowest = 0.11;
gap = 0.05;

axPos3 = [0.1300 y_lowest                       0.7750 height_narrow];
axPos2 = [0.1300 axPos3(2)+height_narrow+gap    0.7750 height];
axPos1 = [0.1300 axPos2(2)+height+gap           0.7750 height];

cbarPos = [0.92, y_lowest, 0.02, axPos1(2)+axPos1(4)-y_lowest];

% custom_cmap_map = load('Berlin_RdBu_denser.mat');
% custom_cmap_map = custom_cmap_map.RdBu;
% cMAP = colormap(custom_cmap_map);


tmp = 0.5*linspace(0,1.0,100)';
ones_tmp = ones(100,1);
cMAP = [flipud(tmp)+0.5, 0*tmp, tmp; flipud(tmp), 0*tmp, tmp+0.5];
colormap(cMAP)


N = NPop;

maxZ = max(zz(:));
minZ = min(zz(:));
lenZ = maxZ - minZ;

fig = figure(1);
fig.Position = [449 300 758 320];

hold on

nSkip = 2;

tf = 10000;

disp(tf)
timeArr = linspace(1,tf,nTVars);
timeArr = timeArr(1:nSkip:end);

ax1 = subplot('Position',axPos1);%  3,1,1);
plot(nan,nan);
hold on;

title_str = sprintf("log_{10} p_e: %1.2f, p_m: %1.2f", log10(p2expltArr(i_p2e)),log10(p2msngrArr(i_p2m)));
title(title_str)

var2Plot = zsArr;
var2Plot = reshape(var2Plot, [], nTVars);
var2Plot = var2Plot(:,1:nSkip:end);

timeArr_Mat = repmat(timeArr, 100, 1);
scatter(timeArr_Mat,var2Plot,1)
% plot(timeArr_Mat,var2Plot)
for i = 1:length(timeArr)/1
    t_ind = timeArr(i);

    z = var2Plot(:,i);
    indCol = round((z-minZ)/(lenZ)*(length(cMAP)-1));
    indCol = max(1,indCol);
    col = cMAP(indCol,:);

    scatter(t_ind*ones(N,1), z, 2, col,'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0)
end
ylabel("Positions")

ylim([0 1])

ax1.XTickLabel = {};
set(gca,'fontsize',12);




ax2 = subplot('Position',axPos2); % 3,1,2);
plot(nan,nan);
hold on;

var2Plot = zpArr;
var2Plot = reshape(var2Plot, [], nTVars);
var2Plot = var2Plot(:,1:nSkip:end);

timeArr_Mat = repmat(timeArr, 100, 1);
% scatter(timeArr_Mat,var2Plot,1)
% plot(timeArr_Mat,var2Plot)
ylabel("Opinions")
%
for i = 1:length(timeArr)/1
    t_ind = timeArr(i);

    z = var2Plot(:,i);
    indCol = round((z-minZ)/(lenZ)*(length(cMAP)-1));
    indCol = max(1,indCol);
    col = cMAP(indCol,:);

    scatter(t_ind*ones(N,1), z, 2, col,'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0)
end

ylim([0 1])
ax2.XTickLabel = {};
set(gca,'fontsize',12);

ax3 = subplot('Position',axPos3); %3,1,3);
plot(nan,nan);
hold on;
% Opinion Plots with fixed with

statesArr = stateArr;
statesArr = reshape(statesArr, [], nTVars);
statesArr = statesArr(:,1:nSkip:end);

for i = 1:length(timeArr)/1
    t_ind = timeArr(i);

    z = var2Plot(:,i);
    
    z = sort(z);

    indCol = round((z-minZ)/(lenZ)*(length(cMAP)-1));
    indCol = max(1,indCol);
    col = cMAP(indCol,:);

    scatter(t_ind*ones(N,1), [1:100]', 2, col,'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',0)
    
    % if(rem(t_ind,2)==1)
        states = statesArr(:,i);

        messengers = states==0;
        messengers_id = find(messengers);
        scatter(t_ind*ones(1,sum(messengers)), messengers_id, 1.5,'filled','MarkerFaceColor','w','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0);
    % end
end



set(gca,'fontsize',12);
set(gcf,'color','w');

xlabel('Time Step')
ylabel({'Collective';'Opinion'})

cbar = colorbar();
cbar.Position = cbarPos; % [0.9200 0.10 0.02 0.800];

% subplot(2,1,1)
% title_str = sprintf("log_{10} p_e: %1.2f, p_m: %1.2f", log10(p2expltArr(i_p2e)),log10(p2msngrArr(i_p2m)));
% title(title_str)


file_str = sprintf("data/opinions_i_p2e__%i__i_p2m__%i_customCMap",i_p2e,i_p2m);
exportgraphics(fig,file_str+".png",'resolution',300)

end

% var2Plot = zpArr;
% scatter(timeArr_Mat,var2Plot,1,'k')

