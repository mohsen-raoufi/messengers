%
clc;
clear all; close all;

addpath('../functions/');

folder_path = 'data/';
files = dir(fullfile(folder_path, '*.mat'));


custom_cmap_map = load('Berlin_RdBu_denser.mat');
custom_cmap_map = custom_cmap_map.RdBu;
cMAP = colormap(custom_cmap_map);


for k = 1:length(files)

    % for k = 3:3 % length(files):length(files) % 3:4
    filename = fullfile(folder_path, files(k).name);
    fprintf('Processing %s\n', filename);


    load(filename);

    close all;


    fig = figure(1);
    fig.Position = [248 396 600 600];

    set(gcf,'color','w');
    hold on;
    box on;
    axis(1.0*[min(x) max(x), min(y) max(y)]);
    axis square

    % xlabel('x'); ylabel('y');
    ax = gca;
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';

    colormap gray


    [~, cntr1] = contour(xx,yy,zz,'LineStyle','--','EdgeAlpha',0.2);
    [~, cntr2] = contour(xx,yy,zz,'LineStyle','--','EdgeAlpha',0.2);



    maxZ = max(zz(:));
    minZ = min(zz(:));
    lenZ = maxZ - minZ;
    

    time = size(posArr,5);
    [sc, sc2] = plot_pos_net(posArr, zpArr, stateArr, time, cMAP, minZ, lenZ, NPop, linkThresh, true, 40);

    [~, cntr1] = contour(xx,yy,zz,[zEnv zEnv],'-.k','LineWidth',1.5,'EdgeAlpha',0.5);
    

    file_strng = [strng,'_final_pos'];
    savefig(fig, [file_strng,'.fig'])
    % exportgraphics(fig,[file_strng,'.pdf'])
    exportgraphics(fig,[file_strng,'.png'],"Resolution",300)

    beep;
    pause(2);
    % beep;
end


beep;
% beep;
pause(2);
beep;
pause(1);
beep
