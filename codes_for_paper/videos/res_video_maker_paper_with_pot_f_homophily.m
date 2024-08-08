
%%
clc; clearvars -except strng;
close all;

% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_10k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1.mat
% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_15k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1.mat
% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_15k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1.mat
% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_20k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1


% load data/single_run_22_Jul_2024__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_36__i_p2m_45.mat
load data/single_run_22_Jul_2024__NPop_100_Arena_1__tf_1k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_1__i_p2m_49.mat
% load data/single_run_23_Jul_2024__NPop_100_Arena_1__tf_2k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.15__i_p2e_1__i_p2m_49.mat


nSkip = 2;
nTrace = 10;
nSkipTrace = 4; 0.2*nSkip;

% the ones I made for the animation
% load
% not_to_push_data/single_run_02_Nov_2023__NPop_100_Arena_1__tf_50k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.4__i_p2e_1__i_p2m_49.mat
load(strng);

% load not_to_push_data/single_run_04_Oct_2023__NPop_100_Arena_1__tf_20k__lowRand4Explt__cone__BasicMarkov__initENumMsngr_sensRang_0.2__nMC_1__i_p2e_49__i_p2m_1
debugMode  = 1;     % 1: all 3 plots, 2: just error, 3: just pos

NRobots = NPop;

saveMode = true; %  false; % 
if(saveMode)
    video_name = ['video_sample_echo_chambers__i_p2e_',num2str(i_p2e), '__i_p2m_',num2str(i_p2m),'_newCMap_2_withCBar_and_pot_field_link_4'];
    if(isunix)
        if(~ismac)
            v = VideoWriter(video_name ,'Motion JPEG AVI'); % 'MPEG-4'); %     % For Saving Animation
        else
            v = VideoWriter(video_name ,'MPEG-4'); % 'MPEG-4'); %     % For Saving Animation
        end
    end
    open(v);
end

switch debugMode
    case 1
        f = figure(1);
        f.Position = [248 396 641 523];%[209 394 1493 383];%[209 394 1113 383];%[300 180 416 383]; % PC

        %         f.Position = [212 69 1100 1100];
    case 2
        f = figure(1);
        f.Position = [277.6667 138.3333 732.6667 465.3333]; % Laptop
        % f.Position = [2433 34.3333 732.6667 465.3333]; % Laptop + Monitor
        %         f.Position = [708 388 722 485];%[209 394 1113 383]; % PC SCIoI %[300 180 416 383];
        legend('off')
    case 3
        f = figure(1);
        f.Position = [1121 401 546 514];%[0 0 546 514];%[1121 401 546 514]; %
end

colormap('hot');
colormap('winter')

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

switch debugMode
    case 1
        sc = plot(nan,nan);% scatter(pos(1,:),pos(2,:),'filled','ro');
        sc2 = sc;
        hold on;

        xlabel('x'); ylabel('y');
        scSW = plot(nan,nan);
        %         axis(1.0*[-1 1 -1 1])
        axis(1.0*[min(x) max(x), min(y) max(y)]);
        axis square
end

% box('on');
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.FontSize = 20;
fontsize(f, 20, "pixels")

set(gcf,'color','w');

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

% get the dissonance homophily potential field
% xx, yy, zz, hh
time = 1;

% reshape xx,yy,zz for higher resolution
landBounds = ArenaScale*[-1 1; -1 1]; %
initBounds = landBounds;

nDiscreteSpace = 200;
x = linspace(landBounds(1,1),landBounds(1,2),nDiscreteSpace);
y = linspace(landBounds(2,1),landBounds(2,2),nDiscreteSpace);
[xx,yy] = meshgrid(x,y);
zz = landFunc(xx,yy);
% landFunc = @(x,y)landFunc(x,y);
landFunc = @(x,y) (landFunc(x,y)-min(zz(:)))/(max(zz(:))-min(zz(:)));
zz = landFunc(xx,yy);

hh = get_dissonance(time,xx,yy,zz,NPop, posArr, linkThresh, zpArr);

sf = surf(xx,yy,hh);

hold on;

cbar = colorbar;
view(0,90);
cbar.Label.String = 'Dissonance';

clim([0, 0.3])
for time=1:nSkip:nTVars
    delete(sc);
    delete(sf);

    tt = (time-1)/nSkip + 1;

    timeStep = time*nSkipSave;
    title(['Time Step: ',num2str(timeStep),' / ', num2str(tf)],'FontSize',15);
    
    hh = get_dissonance(time,xx,yy,zz,NPop, posArr, linkThresh, zpArr);
    
    sf = surf(xx,yy,hh,'FaceColor','flat','EdgeAlpha',0.0);
    
    % limits = [nanmin(hh(:)),11];%nanmax(hh(:))];
    % zlim(limits);

    % clim(limits);
    
    pos = posArr(1,1,:,:,time);
    pos = reshape(pos,[2,NPop]);

    sc = scatter3(pos(1,:),pos(2,:),10*ones(length(pos)),20,'MarkerFaceColor','k','LineWidth',1.5,'MarkerEdgeColor','w','MarkerEdgeAlpha',1.0);

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

function hh = get_dissonance(time, xx, yy, zz, NPop, posArr, linkThresh, zpArr)
hh = nan(size(zz));
for ii=1:size(xx,1)
    for jj=1:size(yy,1)
        pnt = [xx(ii,jj); yy(ii,jj)];
        input_pnt = zz(ii,jj);

        neighbors = zeros(NPop,1);
        pos = posArr(1,1,:,:,time);
        pos = reshape(pos,[2,NPop]);
        for iii=1:NPop
            if(norm(pnt-pos(:,iii))<linkThresh)
                neighbors(iii) = 1;
            end
        end
        cur_zp = zpArr(1,1,:,time);
        cur_zp = reshape(cur_zp, [],1);
        num_neighb = sum(neighbors);
        if(num_neighb>0)
            soc_sig = cur_zp'*neighbors/num_neighb;
        else
            soc_sig = nan;
        end

        disonance = 0.5*abs(input_pnt-soc_sig)^0.5;
        hh(ii,jj) = disonance;

    end
end

end