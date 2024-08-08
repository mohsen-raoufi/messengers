% Title: Animation of Analytical Parameter space for Messenger_ratio and Tau_s from DMP Result
% Description: ***
% Author: Mohsen Raoufi
% Contact: mohsenraoufi@icloud.com
% Affiliation: Research Cluster of Excellence, "Science of Intelligence"
% Date Created: July 20, 2024
% Version: 1.0
% Usage: Run this script with MATLAB R2023a or later. Run it in the directory that contains your data.
% License: Distributed under the MIT License. See LICENSE.txt for more information.
% Citation: please cite our work if you use this code
% "Messengers: Breaking Echo Chambers in Collective Opinion Dynamics with Homophily"

clc;
close all;
clear all;

saveMode = true; % false; % 
if(saveMode)
    video_name = 'DMP_switching_dynamic_2';
    if(isunix)
        if(true)%~ismac)
            v = VideoWriter(video_name ,'Motion JPEG AVI'); % 'MPEG-4'); %     % For Saving Animation
        else
            v = VideoWriter(video_name ,'MPEG-4'); % 'MPEG-4'); %     % For Saving Animation
        end
    end
    open(v);
end

nSkipPar = 1;

p2expltArr = exp(-[2:nSkipPar:50]./5);
p2msngrArr = p2expltArr;

fig = figure(333);

hold off
plot(nan,nan)
fig.Color = 'w';
fig.Position = [336 188 782 661];


logp2e = log10(p2expltArr(1:end));
logp2m = log10(p2msngrArr(1:end));

[LogP2E, LogP2M] = meshgrid(logp2e,logp2m);
[P20Arr, P21Arr] = meshgrid(p2expltArr,p2msngrArr);
[P2E, P2M] = meshgrid(p2expltArr,p2msngrArr);

analytic_m = (P2M./(P2E+P2M));

analytic_tau = floor(0.5*(P2E + P2M)./(P2E.*P2M));

tf_sim = 10^3;

plot(nan,nan)
tmp = 0.5*(linspace(0,1.0,200)');
ones_tmp = ones(100,1);
cMAP = [flipud(tmp)+0.5, 0*tmp, tmp; flipud(tmp), 0*tmp, tmp+0.5];
% cMAP
cMAP = flipud(cMAP);
cMAP = [cMAP; 1 1 1];
% cMAP = (cMAP).^.2;
% cMAP = cMAP./sum(cMAP,2);
% cMAP = normr(cMAP);
colormap(cMAP)
% colormap([0.2 0.2 0.2; 1 1 1]);

plot(nan,nan);

sf = surf(LogP2E, LogP2M, (analytic_m));
hold on;

view(0,90)

swch_Mat = zeros(size(analytic_m));
counter_Mat = swch_Mat;

colorbar()
caxis([0 1.01]);
cb = colorbar('YTick',[0:0.2:1], 'FontSize', 24);
set(gca,'FontSize',24)
cb.Label.String = 'Messenger Ratio';

xlabel("$\log_{10} p_{E}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')
ylabel("$\log_{10} p_{M}$",'Interpreter','latex','FontSize',30,'FontWeight','bold')

axis square
axis([min(logp2e), max(logp2e), min(logp2m), max(logp2m)])

for t=1:tf_sim

    delete(sf)

    title(sprintf('Time %3i',t))

    counter_Mat = counter_Mat + ones(size(analytic_m));

    swch_Mat = counter_Mat > analytic_tau;

    % reset counter_mat
    counter_Mat = counter_Mat.*not(swch_Mat);

    % disp(counter_Mat)
    % surf(LogP2E, LogP2M, counter_Mat)
    empty_mat = (swch_Mat);
    % empty_mat = double(empty_mat);
    % empty_mat(empty_mat == 1) = 1;

    display_mat = (analytic_m);
    display_mat(empty_mat == 1) = 1.1;

    % sf = surf(LogP2E, LogP2M, display_mat,'EdgeAlpha',0);%,'FaceColor','flat');
    sf = imagesc(logp2e, logp2m, display_mat);

    drawnow;

    pause(0.001)

    if(saveMode)
        frame = getframe(gcf);    % For Saving Animation
        writeVideo(v,frame);
    end

end

if(saveMode)
    close(v)
end
