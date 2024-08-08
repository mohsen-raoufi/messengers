
%%
clc; clearvars -except strng;
close all;

% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_10k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1.mat
% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_15k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1.mat
% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_15k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1.mat
% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_20k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1

% %% For presentation
% load data/single_run_22_Jul_2024__NPop_100_Arena_1__tf_1k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat
load data/single_run_22_Jul_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat
% load data/single_run_22_Jul_2024__NPop_100_Arena_1__tf_1k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_1__i_p2m_49.mat
% load data/single_run_23_Jul_2024__NPop_100_Arena_1__tf_2k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat
nSkip = 2;
nTrace = 10;
nSkipTrace = 1; 0.2*nSkip;


load(strng)

% the ones I made for the animation
% load
% not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_1__i_p2m_49.mat
load(strng);

% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_20k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1
debugMode  = 1;     % 1: all 3 plots, 2: just error, 3: just pos

NRobots = NPop;
% strng = 'Test_01';
% % load(strng);% finalTestForVideo.mat

saveMode = true; %  false; %
if(saveMode)
    video_name = ['video_sample_echo_chambers__i_p2e_',num2str(i_p2e), '__i_p2m_',num2str(i_p2m),'_newCMap_2_withCBar_and_opinions'];
    if(isunix)
        if(true) % ~ismac)
            v = VideoWriter(video_name ,'Motion JPEG AVI'); % 'MPEG-4'); %     % For Saving Animation
        else
            v = VideoWriter(video_name ,'MPEG-4'); % 'MPEG-4'); %     % For Saving Animation
        end
    end
    open(v);
end

fig = figure(1);
width = 700;
height = 600;
if(ismac)
    fig.Position = [636 125 width height];
else
    fig.Position = [2407 292 width height];
end


pos_ax1 = [0.1 0.1 0.7 0.7*width/height];
pos_ax2 = [0.9 0.1 0.07 0.7*width/height];

% load from external file
custom_cmap_map = load('Berlin_RdBu_denser.mat');
custom_cmap_map = custom_cmap_map.RdBu;
cMAP = colormap(custom_cmap_map);


maxZ = max(zz(:));
minZ = min(zz(:));
lenZ = maxZ - minZ;

% pos = posArr(:,:,1);
tmp = posArr(1,1,:,:,1,1);
pos = reshape(tmp,[2,size(tmp,4)]);

Adjc = zeros(NRobots); % Adjacancy matrix, previously called Links!
for jj=1:NRobots-1
    for ii=jj+1:NRobots
        if(norm(pos(:,jj)-pos(:,ii))<linkThresh)
            Adjc(ii,jj) = 1;
        end
    end
end
Adjc = Adjc + Adjc';

ax1 = subplot('Position',pos_ax1,'Units','pixels');
ax2 = subplot('Position',pos_ax2);

subplot(ax1);

sc = plot(nan,nan);% scatter(pos(1,:),pos(2,:),'filled','ro');
sc2 = sc;
hold on;
hold on;
ax1.Units = 'pixels';


xlabel('x'); ylabel('y');
scSW = plot(nan,nan);
colormap gray
%         axis(1.0*[-1 1 -1 1])
axis(1.0*[min(x) max(x), min(y) max(y)]);
axis square


% box('on');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 20;
fontsize(fig, 20, "pixels")

set(gca,'Units','pixels')
set(gca,'Color','w')
set(gcf,'Color','w')


% plt(1) = plot(
% lgnd = legend

envCent = 0.5*[min(x)+max(x); min(y)+max(y)];
envR = envCent(1) - min(x);

objFArrOld = 0;
% ss = imagesc(x,y,zz);
rng('shuffle');


swch = zeros(1,NRobots);


robCol = 'r';
tSwch = 0;%tSw*tf;
swchBool = 0;
% nSkipTime = 10*20;

%
% time = nTVars;
% close all
% plot_pos_net(posArr, zpArr, stateArr, time, cMAP, minZ, lenZ, NRobots, linkThresh, true, 40);
% hold on
% plot_pos_net(posArr, zpArr, stateArr, time, cMAP, minZ, lenZ, NRobots, linkThresh, true, 40);
% axis square
% axis([-1 1 -1 1])

% beep
% beep
% pause
%
% hold off

% title_str = sprintf('P_E: %1.2f, P_M: %1.2f',log10(p2expltArr), log10(p2msngrArr));
% title(title_str)


% for time=tSwch:nSkipTime:tSwch+1000%


% set(gca,'fontsize',12);
% set(gca,'fontsize',40); % MAC requires larger font-size
% fontsize(18, "inches")
%
% set(gcf,'color','w');

[~, cntr1] = contour(xx,yy,zz,'LineStyle','--','EdgeAlpha',0.2);
[~, cntr2] = contour(xx,yy,zz,'LineStyle','--','EdgeAlpha',0.2);

scs_trace = repmat(sc, [1,nTrace]);
scs_trace2 = scs_trace;

hold on;

subplot(ax2);
hold on;

sw_plt = swarmchart(ones(NPop,1),zeros(NPop,1),15,zeros(NPop,1),'filled');

ax2.Units = 'pixels';
% ax2.Position = [0.85 0.1100 0.07 0.8150];
% ax2.Position = [720 50 40 450];
set(ax2,'Box','off');
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';
ax2.YLim = [0 1];

colormap(cMAP);
cbar = colorbar;
cbar.Position
cbar.Units = 'pixels';
cbar.Position = [600 ax2.Position(2) 20 ax2.Position(4)];
cbar.FontSize = 14;
cbar.Limits = [0, 1];
set(cbar.Title, 'String', "Opinions")

for time=1:nSkip:nTVars% nSkipTime:tf
    % for time=nTVars:nTVars% nSkipTime:tf
    %     if(debugMode), sgtitle(sprintf('Time: %4i',time)); end
    subplot(ax1);

    delete(sc);
    delete(sc2);
    delete(cntr1);
    delete(cntr2);
    %
    %
    delete(scs_trace);
    delete(scs_trace2);

    delete(findobj(gca, 'type', 'line'))

    % delete(scs_trace);

    % hold off;
    % plot(nan,nan);
    % hold on;
    % axis square
    % axis([-1 1 -1 1]*ArenaScale)

    [~, cntr1] = contour(xx,yy,zz,[zEnv zEnv],'-.k','LineWidth',1.5,'EdgeAlpha',0.5);
    tt = (time-1)/nSkip + 1;

    timeStep = time*nSkipSave;
    title(['Time Step: ',num2str(timeStep),' / ', num2str(tf)],'FontSize',15);


    marker_size_big = 60;
    for i_trace=1:nTrace
        % delete(scs_trace(i_trace));
        % delete(scs_trace2(i_trace));
        time_trace = max(1, time-i_trace*nSkipTrace);
        marker_size = (1 - (i_trace-1)/nTrace) * marker_size_big/4;
        [scs_trace(i_trace), scs_trace2(i_trace)] = plot_pos_net(posArr, zpArr, stateArr, time_trace, cMAP, minZ, lenZ, NRobots, linkThresh, false, marker_size);
        % scs_trace(time_trace) = tmp_sc;
    end

    [sc, sc2] = plot_pos_net(posArr, zpArr, stateArr, time, cMAP, minZ, lenZ, NRobots, linkThresh, true, 40*2);

    %
    % axis off
    % box on
    % % box off
    % xlm = xlim;
    % ylm = ylim;
    % line(xlm,[ylm(2), ylm(2)],'Color','k')
    % line(xlm,[ylm(1), ylm(1)],'Color','k')
    % line([xlm(2), xlm(2)],ylm,'Color','k')
    % line([xlm(1), xlm(1)],ylm,'Color','k')

    [~, cntr2] = contour(xx,yy,zz,'LineStyle','--','EdgeAlpha',0.2);

    subplot(ax2);
    delete(sw_plt);
    
    opinions = zpArr(1,1,:,time);
    opinions = reshape(opinions,1,[]);
    sw_plt = swarmchart(ones(size(opinions)),opinions,15,opinions,'filled');

    clim([0,1])

    % colormap(cMAP)



    % cbar = colorbar;
    %     % Get the position of the colorbar
    % cbarPos = cbar.Position;
    %
    % tmp = zpArr(1,1,:,time,1);
    % z = reshape(tmp,[1,size(tmp,3)]);
    %
    % % Create x and y coordinates for scatter plot on the colorbar
    % xScatter = ones(size(z)) * (cbarPos(1) + cbarPos(3) / 2); % x coordinates (center of the colorbar)
    % yScatter = cbarPos(2) + z * cbarPos(4); % y coordinates (mapped along the colorbar length)
    % scatter(xScatter, yScatter, 30, z, 'filled'); % Create scatter plot
    %



    % title(title_str)

    % colorbar(cMAP)
    drawnow

    %     pos = posArr(:,:,time);

    if(saveMode)
        frame = getframe(gcf);    % For Saving Animation
        writeVideo(v,frame);
    end

end

th = 0:0.02:2*pi;
r = zEnv;
% plot(r*cos(th),r*sin(th),'--g','LineWidth',2)
drawnow
if(saveMode), close(v); end
% print([strng,'_final_pos.png'],'-dpng','-r300','-fillpage')

ax = gca;
exportgraphics(ax,[strng,'_final_pos.png'],'Resolution',300);

beep

% figure();
% sc = scatter(pos(1,:),pos(2,:),10,'filled','MarkerFaceColor',0*[1,1,1]);
% axis(ArenaScale*[-1 1 -1 1])
% axis square
