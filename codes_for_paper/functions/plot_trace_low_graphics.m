close all;
clc;

% load single_run_28_Sep_2023__NPop_200_Arena_1__tf_4k__lowRand4Explt__sine2__BasicMarkov__initENumMsngr_sensRang_0.4__nMC_1__i_p2e_49__i_p2m_1.mat

fig = figure(2);
fig.Position = [680 458 560 420];
hold on;
axis([-1 1 -1 1])
% axis([-0 2 -0 2])
axis square



for time=1:nTVars
    tmp = posArr(1,1,:,:,time,1);
    pos = reshape(tmp,[2,size(tmp,4)]);

    if(rem(time,2)==0)
        % delete(sct)
        clr = (time/nTVars).^1*[1 0 0];
        sct = scatter(pos(1,:), pos(2,:), 10, clr, 'filled', 'MarkerFaceAlpha', 0.01, 'MarkerEdgeAlpha', 0);
%     scatter(posArr(1,:,t), posArr(2,:,t), 10, clr, 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0);
        % drawnow();
    end
    
    
end


set(gca,'Color','w')
set(gcf,'Color','w')
axis off

% print(strng+"_trace.png","-dpng",'-r300')

