clc;
clear all;
close all;

addpath("../functions/")
% strng = 'data/single_run_28_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat';
% strng = 'data/single_run_28_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_21__i_p2m_30.mat';
% strng = 'data/single_run_28_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_7__i_p2m_12.mat';
% strng = 'data/single_run_28_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat';
% strng = 'data/single_run_28_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat';

strng = 'not_to_push_data/single_run_24_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_6__i_p2m_7.mat';
strng = 'not_to_push_data/single_run_24_May_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_5__i_p2m_6.mat';

% strng = '../../codes_for_paper/not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat';
load(strng)

nSkip = 1;

NRobots = NPop;

fig = figure(1);
fig.Position = [76 345 1743 153];
% fig.Position = [1 64 1512 773]; % wide

set(gca,'FontSize',14)
set(gcf,'color','w');

plot(nan,nan);
hold on;
axis([0 tf, 0 100])


% custom_cmap_map = load('Berlin_RdBu_denser.mat');
% custom_cmap_map = custom_cmap_map.RdBu;
% cMAP = colormap(custom_cmap_map);

% cMAP = colormap('parula'); % turbo
% colormap(cMAP)

tmp = 0.5*linspace(0,1.0,100)';
ones_tmp = ones(100,1);
cMAP = [flipud(tmp)+0.5, 0*tmp, tmp; flipud(tmp), 0*tmp, tmp+0.5];
colormap(cMAP)





maxZ = max(zz(:));
minZ = min(zz(:));
lenZ = maxZ - minZ;

% pos = posArr(:,:,1);
tmp = posArr(1,1,:,:,1,1);
pos = reshape(tmp,[2,size(tmp,4)]);

zp_0 = zpArr(1,1,:,1);
[~, sorting_index] = sort(zp_0);

marker_size = 1.5;
marker_size = 15;
% marker_size = 50;

for time=1:nSkip:nTVars
    
    tt = (time-1)/nSkip + 1;
    
    timeStep = time*nSkipSave;
    
    opinions = reshape(zpArr(1,1,:,time),1,[]);

    states = reshape(stateArr(1,1,:,time),1,[]);

    % based on the initial position
    % sorted_opinions = opinions(sorting_index);

    % based on the current position
    % sorted_opinions = sort(opinions);
    

    % based on clusters and opinion
    positions = reshape(posArr(1,1,:,:,time),2,[]);
    groupIndxArr = numCluster_rad(positions(1,:)', positions(2,:)', linkThresh);
    % combi_val = [groupIndxArr, opinions'];
    combi_val = [opinions', groupIndxArr];
    sorted_combi_val = sortrows(combi_val);
    sorted_opinions = sorted_combi_val(:,1);

    % scatter(timeStep*ones(1,100), [1:100], marker_size, sorted_opinions, 'filled');
    scatter(timeStep*ones(1,100), [1:100], marker_size, sorted_opinions, 'filled','MarkerFaceAlpha',0.25);
    
    
    % 
    % messengers = states==0;
    % messengers_id = find(messengers);
    % 
    % scatter(timeStep*ones(1,sum(messengers)), messengers_id, marker_size/2,'MarkerFaceColor','w','MarkerEdgeAlpha',0); % ,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
    messengers = states==0;
    messengers_id = find(messengers);


    % scatter(timeStep*ones(1,sum(messengers)), messengers_id, 2,'filled','MarkerFaceColor','w','MarkerEdgeAlpha',0);

    
    % if(rem)
    % drawnow

end

set(gca,'FontSize',14)

xlabel('Time Step');
ylabel('ID');
cbar = colorbar;
cbar.Label.String = 'Opinion';

save_mode = true; false; % 
if(save_mode)
    % file_strng = [strng, '_opinions_only_with_dots_BerlinCMap'];
    % file_strng = [strng, '_opinions_only_BerlinCMap'];
    file_strng = [strng, '_opinions_only_RdBuCMap'];
    % file_strng = [strng, '_opinions_only_wide'];
    
    exportgraphics(fig,[file_strng,'.png'],"Resolution",300)
    % exportgraphics(fig,[file_strng,'.pdf'])
    % savefig([file_strng,'.fig'])
end
