function fig = opinion_plot_func(strng)
clc;
% clear all;


load(strng,"-mat");

tf = 50000; % bug in loading .mat file!
% disp(tf)

close all;

nSkip = 5;

NRobots = NPop;

fig = figure();
fig.Position = [76 345 1743 153];

set(gca,'FontSize',14)
set(gcf,'color','w');

plot(nan,nan);
hold on;
axis([0 tf, 0 100]);


custom_cmap_map = load('Berlin_RdBu_denser.mat');
custom_cmap_map = custom_cmap_map.RdBu;
cMAP = colormap(custom_cmap_map);

% cMAP = colormap('parula'); % turbo
% colormap(cMAP)

% tmp = 0.5*linspace(0,1.0,100)';
% ones_tmp = ones(100,1);
% cMAP = [flipud(tmp)+0.5, 0*tmp, tmp; flipud(tmp), 0*tmp, tmp+0.5];
% colormap(cMAP)

cbar = colorbar;
cbar.Limits = [0, 1];
caxis([0, 1]);



maxZ = max(zz(:));
minZ = min(zz(:));
lenZ = maxZ - minZ;

% pos = posArr(:,:,1);
tmp = posArr(1,1,:,:,1,1);
pos = reshape(tmp,[2,size(tmp,4)]);

zp_0 = zpArr(1,1,:,1);
[~, sorting_index] = sort(zp_0);

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

    scatter(timeStep*ones(1,100), [1:100], 5, sorted_opinions, 'filled','MarkerEdgeAlpha',0);
    
    messengers = states==0;
    messengers_id = find(messengers);


    scatter(timeStep*ones(1,sum(messengers)), messengers_id, 2,'filled','MarkerFaceColor','w','MarkerEdgeAlpha',0);

    
    % if(rem)
    % drawnow

end

drawnow();

set(gca,'FontSize',14)

xlabel('Time Step');
ylabel('ID');
cbar = colorbar;
cbar.Label.String = 'Opinion';

save_mode = false; % true; 
if(save_mode)
    file_strng = [strng, '_opinions'];
    % savefig([file_strng,'.fig'])
    % exportgraphics(fig,[file_strng,'.pdf'])
    exportgraphics(fig,[file_strng,'.png'],"Resolution",300)
end
